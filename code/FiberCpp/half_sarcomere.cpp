/**
 * @file    half_sarcomere.cpp
 * @brief   Source file for the half_sarcomere class
 * @author  Ken Campbell
 */

#include <cstdio>

#include "half_sarcomere.h"
#include "thick_filament.h"
#include "kinetic_scheme.h"
#include "transition.h"
#include "m_state.h"
#include "thin_filament.h"
#include "muscle.h"
#include "FiberSim_model.h"
#include "FiberSim_options.h"
#include "FiberSim_protocol.h"
#include "FiberSim_data.h"
#include "JSON_functions.h"

#include "gsl_math.h"
#include "gsl_vector.h"
#include "gsl_matrix.h"
#include "gsl_spmatrix.h"
#include "gsl_splinalg.h"
#include "gsl_linalg.h"
#include "gsl_rng.h"
#include "gsl_roots.h"

// Structure used for root finding for force control mode
struct force_control_params
{
    double target_force;
    half_sarcomere* p_hs;
};

// Structure used for handling lattice binding / unbinding events
struct lattice_event
{
    char mol_type;          // m (myosin) or c (mybpc)
    int m_f;                // thick filament index
    int m_n;                // thick filament myosin
    int a_f;                // thin filament index
    int a_n;                // thin filament binding site
    transition * p_trans;   // transition
};

// Constructor
half_sarcomere::half_sarcomere(
    FiberSim_model* set_p_fs_model,
    FiberSim_options* set_p_fs_options,
    FiberSim_protocol* set_p_fs_protocol,
    muscle* set_p_parent_m,
    int set_hs_id)
{
    // Initialise

    // Set the pointers
    p_fs_model = set_p_fs_model;
    p_fs_options = set_p_fs_options;
    p_fs_protocol = set_p_fs_protocol;
    p_parent_m = set_p_parent_m;
    
    // Set the id
    hs_id = set_hs_id;

    // Set the class pointers to the kinetic scheme for myosin
    for (int i = 0; i < p_fs_model->m_no_of_isotypes; i++) {
        p_m_scheme[i] = p_fs_model->p_m_scheme[i];
    }

    // and MyBPC
    for (int i = 0; i < p_fs_model->c_no_of_isotypes; i++)
    {
        p_c_scheme[i] = p_fs_model->p_c_scheme[i];
    }

    // Set up the attachment span
    adjacent_bs = p_fs_options->adjacent_bs;
    m_attachment_span = 1 + (2 * adjacent_bs);

    // Create an array of lattice events
    for (int i = 0; i < (MAX_NO_OF_ADJACENT_BS * MAX_NO_OF_TRANSITIONS); i++)
    {
        p_event[i] = new lattice_event{-1,-1,-1,-1};
    }

    // Initialize macroscopic state variables
    time_s = 0.0;
    hs_length = p_fs_model->initial_hs_length;
    hs_force = 0.0;
    pCa = 0.0;
    hs_titin_force = 0.0;
    hs_extracellular_force = 0.0;

    // Zero the step_counter
    step_counter = 0;

    // Initialise the random number generator
    // This needs to be done before the thick filaments are allocated to allow for
    // lambda jitter
    const gsl_rng_type* T;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    rand_generator = gsl_rng_alloc(T);
    gsl_rng_set(rand_generator,
        unsigned long(100 * (p_parent_m->muscle_id + 1) + (hs_id + 1)));

    // Log
    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file, "In half_sarcomere constructor for hs[%i]\n", hs_id);
    }

    // Set some parameters from the FiberSim_model object

    // Start with the hs_length
    hs_length = p_fs_model->initial_hs_length;

    // Now the number of thick filaments
    m_n = p_fs_model->m_n;

    // There must be twice as many thin as thick filaments
    a_n = 2 * m_n;

    // Make new thick filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        p_mf[m_counter] = new thick_filament(
            p_fs_model,
            p_fs_options,
            this,
            m_counter);
    }

    // Make new thin filaments
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        p_af[a_counter] = new thin_filament(
            p_fs_model,
            p_fs_options,
            this,
            a_counter);
    }

    // Allocate space for the nearest_m_matrix
    nearest_actin_matrix = gsl_matrix_short_alloc(m_n, 6);

    // Calculate the y and z coordinates
    initialise_filament_y_z_coordinates(m_n);

    // Initialise the nearest_a_matrix
    initialise_nearest_actin_matrix();

    // Set some values
    m_nodes_per_thick_filament = p_fs_model->m_crowns_per_filament;
    m_cbs_per_thick_filament = p_mf[0]->m_no_of_cbs;
    m_cbs_per_node = p_fs_model->m_hubs_per_crown * p_fs_model->m_myosins_per_hub;
    m_inter_crown_rest_length = p_fs_model->m_inter_crown_rest_length;

    a_nodes_per_thin_filament = p_fs_model->a_bs_per_unit *
        p_fs_model->a_regulatory_units_per_strand;
    a_bs_per_thin_filament = p_af[0]->a_no_of_bs;
    a_bs_per_node = p_fs_model->a_strands_per_filament;
    a_inter_bs_rest_length = p_fs_model->a_inter_bs_rest_length;

    a_regulatory_units_per_strand = p_fs_model->a_regulatory_units_per_strand;
    a_bs_per_unit = p_fs_model->a_bs_per_unit;
    a_strands_per_filament = p_fs_model->a_strands_per_filament;

    t_attach_a_node = p_fs_model->t_attach_a_node;
    t_attach_m_node = p_fs_model->t_attach_m_node;

    // Set up for position calculations
    hs_total_nodes = (m_n * m_nodes_per_thick_filament) +
        (a_n * a_nodes_per_thin_filament);

    // Copy parameters from p_fs_model
    a_no_of_bs_states = p_fs_model->a_no_of_bs_states;
    a_k_stiff = p_fs_model->a_k_stiff;
    a_k_on = p_fs_model->a_k_on;
    a_k_off = p_fs_model->a_k_off;
    a_k_coop = p_fs_model->a_k_coop;

    m_no_of_cb_states = p_fs_model->p_m_scheme[0]->no_of_states;
    m_k_stiff = p_fs_model->m_k_stiff;

    c_no_of_pc_states = p_fs_model->p_c_scheme[0]->no_of_states;
    c_no_of_pcs = p_mf[0]->c_no_of_pcs;

    t_k_stiff = p_fs_model->t_k_stiff;
    t_slack_length = p_fs_model->t_slack_length;
    sprintf_s(t_passive_mode, _MAX_PATH, p_fs_model->t_passive_mode);

    // Extracellular parameters
    sprintf_s(e_passive_mode, _MAX_PATH, p_fs_model->e_passive_mode);
    e_sigma = p_fs_model->e_sigma;
    e_L = p_fs_model->e_L;
    e_slack_length = p_fs_model->e_slack_length;

    m_k_cb = p_fs_model->m_k_cb;

    c_k_stiff = p_fs_model->c_k_stiff;

    // Allocate space for tri-diagonal solve used by iterative technique
    tri_d_vector = gsl_vector_alloc(hs_total_nodes);
    tri_e_vector = gsl_vector_alloc((size_t)hs_total_nodes-1);
    tri_f_vector = gsl_vector_alloc((size_t)hs_total_nodes-1);

    // Intialise them
    initialize_tridiagonal_vectors();

    // Allocate space for f0_vector and initialize it
    f0_vector = gsl_vector_alloc(hs_total_nodes);
    initialise_f0_vector();

    // Allocate space for g_vector and df_vector
    g_vector = gsl_vector_alloc(hs_total_nodes);
    df_vector = gsl_vector_alloc(hs_total_nodes);

    // Allocate space for the x vector
    x_vector = gsl_vector_alloc(hs_total_nodes);

    // Calculate the x positions
    calculate_x_positions();

    // Unpack x_vector to set p_af[]->bs_x and p_mf[]->cb_x
    unpack_x_vector();

    // Map the filaments
    set_cb_nearest_a_f();
    set_cb_nearest_a_n();

    set_pc_nearest_a_f();
    set_pc_nearest_a_n();

    // Initialise and zero the state population vectors
    a_pops = gsl_vector_alloc(a_no_of_bs_states);
    m_pops = gsl_vector_alloc(m_no_of_cb_states);
    c_pops = gsl_vector_alloc(c_no_of_pc_states);

    gsl_vector_set_zero(a_pops);
    gsl_vector_set_zero(m_pops);
    gsl_vector_set_zero(c_pops);
}

// Destructor
half_sarcomere::~half_sarcomere()
{
    // Tidy up
    if (p_fs_options->log_mode > 0)
    {
        fprintf_s(p_fs_options->log_file, "In half_sarcomere destructor for hs[%i]\n", hs_id);
    }

    // Deallocate the nearest_a_matrix
    gsl_matrix_short_free(nearest_actin_matrix);


    // Delete the lattice events
    for (int i = 0; i < (MAX_NO_OF_ADJACENT_BS * MAX_NO_OF_TRANSITIONS); i++)
    {
        delete p_event[i];
    }

    // Delete thick filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        delete p_mf[m_counter];
    }

    // Delete thin filaments
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        delete p_af[a_counter];
    }

    // Delete gsl_vectors
    gsl_vector_free(x_vector);

    // Delete gsl_vectors used for tridiagonal solve
    gsl_vector_free(f0_vector);
    gsl_vector_free(tri_d_vector);
    gsl_vector_free(tri_e_vector);
    gsl_vector_free(tri_f_vector);
    gsl_vector_free(g_vector);
    gsl_vector_free(df_vector);

    // Delete state vectors
    gsl_vector_free(a_pops);
    gsl_vector_free(m_pops);
    gsl_vector_free(c_pops);

    // Delete the random number generator
    gsl_rng_free(rand_generator);
}

// Functions

size_t half_sarcomere::implement_time_step(double time_step,
    double delta_hsl, double sim_mode, double set_pCa)
{
    //! Code runs the half-sarcomere through a single time-step

    // Variables
    size_t x_solve_iterations;

    // Code

    // Update the time and pCa
    time_s = time_s + time_step;
    pCa = set_pCa;

    // Update the step_counter
    step_counter = step_counter + 1;

    // Map the filaments and run kinetics
    set_cb_nearest_a_n();
    set_pc_nearest_a_n();

    thin_filament_kinetics(time_step, pow(10, -pCa));
    thick_filament_kinetics(time_step);

    // Branch depending on sim_mode
    if (sim_mode > 0.0)
    {
        // Replace the delta_hsl with the delta_hsl found by
        // an iterative search
        delta_hsl = calculate_delta_hsl_for_force(sim_mode);
    }

    // Apply length change
    if (fabs(delta_hsl) > 0.0)
    {
        hs_length = hs_length + delta_hsl;
        update_f0_vector(delta_hsl);
    }

    // Calculate positions and deduce force
    x_solve_iterations = calculate_x_positions();

    unpack_x_vector();
    hs_titin_force = calculate_titin_force();
    hs_extracellular_force = calculate_extracellular_force();
    hs_force = calculate_force();


    // Calculate mean filament lengths
    calculate_mean_filament_lengths();

    // Hold state variables
    calculate_a_pops();
    calculate_m_pops();
    calculate_c_pops();

    // Return
    return x_solve_iterations;
}

double half_sarcomere::calculate_delta_hsl_for_force(double target_force)
{
    //! Returns the delta_hsl required for the half-sarcomere to generate the target force

    // Variables
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type* T;
    gsl_root_fsolver* s;
    double r = 0.0;
    double x_lo = -150, x_hi = 150.0;
    struct force_control_params params = { target_force, this };

    gsl_function F;
    F.function = &test_force_wrapper;
    F.params = &params;

    // Test

    // Code
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.01);

//        printf("%5d [%f, %f] %f %f %f\n",
  //          iter, x_lo, x_hi, r, r - 0, x_hi - x_lo);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    return r;
}

double half_sarcomere::test_force_wrapper(double delta_hsl, void* params)
{
    struct force_control_params* p =
        (struct force_control_params*) params;
    half_sarcomere* p_hs = p->p_hs;

    return p_hs->test_force_for_delta_hsl(delta_hsl, params);
}

double half_sarcomere::test_force_for_delta_hsl(double delta_hsl, void *params)
{
    //! Returns the difference between the calculated force and the target force
    //  for a given delta_hsl

    // Variables
    struct force_control_params* p =
        (struct force_control_params*) params;
    double target_force = p->target_force;
    double test_force;

    // Current state - need to return to this at the end
    double original_hs_length = hs_length;
    gsl_vector* original_x_vector = gsl_vector_alloc(hs_total_nodes);

    // Code

    // Save the x_vector which we need later
    gsl_vector_memcpy(original_x_vector, x_vector);

    // Apply the length length
    hs_length = hs_length + delta_hsl;
    update_f0_vector(delta_hsl);

    // Solve the x positions and calculate force
    calculate_x_positions();
    test_force = calculate_force();
    printf("test_force: %g\n", test_force);

    // Restore the half_sarcomere
    hs_length = original_hs_length;
    gsl_vector_memcpy(x_vector, original_x_vector);
    update_f0_vector(-delta_hsl);

    // Recover memory
    gsl_vector_free(original_x_vector);

    // Return the difference
    return (test_force - target_force);
}

void half_sarcomere::initialize_tridiagonal_vectors(void)
{
    //! Code sets values for tridiagonal vectors used for iterative solve

    // Variables

    // Code
    int row_index = 0;

    // Zero vectors
    gsl_vector_set_zero(tri_d_vector);
    gsl_vector_set_zero(tri_e_vector);
    gsl_vector_set_zero(tri_f_vector);

    // Loop through actin filaments
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        for (int node_counter = 0; node_counter < a_nodes_per_thin_filament; node_counter++)
        {
            if (node_counter == 0)
            {
                gsl_vector_set(tri_d_vector, row_index, 2.0 * a_k_stiff);
                gsl_vector_set(tri_e_vector, row_index, -1.0 * a_k_stiff);
            }

            if ((node_counter > 0) && (node_counter < (a_nodes_per_thin_filament - 1)))
            {
                gsl_vector_set(tri_f_vector, ((size_t)row_index - 1), -1.0 * a_k_stiff);
                gsl_vector_set(tri_d_vector, row_index, 2.0 * a_k_stiff);
                gsl_vector_set(tri_e_vector, row_index, -1.0 * a_k_stiff);
            }

            if (node_counter == (a_nodes_per_thin_filament - 1))
            {
                gsl_vector_set(tri_f_vector, ((size_t)row_index - 1), -1.0 * a_k_stiff);
                gsl_vector_set(tri_d_vector, row_index, 1.0 * a_k_stiff);
            }

            row_index = row_index + 1;
        }
    }

    // Loop through myosin
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int node_counter = 0; node_counter < m_nodes_per_thick_filament; node_counter++)
        {
            if (node_counter == 0)
            {
                gsl_vector_set(tri_d_vector, row_index, 2.0 * m_k_stiff);
                gsl_vector_set(tri_e_vector, row_index, -1.0 * m_k_stiff);
            }

            if ((node_counter > 0) && (node_counter < (m_nodes_per_thick_filament - 1)))
            {
                gsl_vector_set(tri_f_vector, ((size_t)row_index - 1), -1.0 * m_k_stiff);
                gsl_vector_set(tri_d_vector, row_index, 2.0 * m_k_stiff);
                gsl_vector_set(tri_e_vector, row_index, -1.0 * m_k_stiff);
            }

            if (node_counter == (m_nodes_per_thick_filament - 1))
            {
                gsl_vector_set(tri_f_vector, ((size_t)row_index - 1), -1.0 * m_k_stiff);
                gsl_vector_set(tri_d_vector, row_index, 1.0 * m_k_stiff);
            }

            row_index = row_index + 1;
        }
    }
}

void half_sarcomere::calculate_g_vector(gsl_vector* x_trial)
{
    //! Sets g_vector which is the part of f in kx = f that is
    //  due to titin, cross-bridges, or mybpc and depends on x

    // Code

    // Zero g_vector
    gsl_vector_set_zero(g_vector);

    // Add in titin
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        int thick_node_index = (a_n * a_nodes_per_thin_filament) +
            (m_counter * m_nodes_per_thick_filament) +
            t_attach_m_node - 1;

        // Loop through surrounding actins
        for (int a_counter = 0; a_counter < 6; a_counter++)
        {
            int thin_node_index =
                (gsl_matrix_short_get(nearest_actin_matrix, m_counter, a_counter) *
                    a_nodes_per_thin_filament) +
                t_attach_a_node - 1;

            double g_adjustment = 0.0;

            if (!strcmp(t_passive_mode, "linear"))
            {
                double x_a = gsl_vector_get(x_trial, thin_node_index);
                double x_m = gsl_vector_get(x_trial, thick_node_index);

                g_adjustment = t_k_stiff * (x_a - x_m);
            }

            gsl_vector_set(g_vector, thin_node_index,
                gsl_vector_get(g_vector, thin_node_index) + g_adjustment);

            gsl_vector_set(g_vector, thick_node_index,
                gsl_vector_get(g_vector, thick_node_index) - g_adjustment);
        }
    }

    // Add in myosins
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int cb_counter = 0; cb_counter < m_cbs_per_thick_filament; cb_counter++)
        {
            // Check for a link
            if (gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_f, cb_counter) >= 0)
            {
                int thick_node_index = node_index('m', m_counter, cb_counter);
                int thin_node_index = node_index('a',
                    gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_f, cb_counter),
                    gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_n, cb_counter));

                double g_adjustment =
                    m_k_cb * (gsl_vector_get(x_trial, thin_node_index) - gsl_vector_get(x_trial, thick_node_index));

                gsl_vector_set(g_vector, thin_node_index,
                    gsl_vector_get(g_vector, thin_node_index) + g_adjustment);

                gsl_vector_set(g_vector, thick_node_index,
                    gsl_vector_get(g_vector, thick_node_index) - g_adjustment);
            }
        }
    }

    // Add in mybpc
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int pc_counter = 0; pc_counter < p_mf[m_counter]->c_no_of_pcs; pc_counter++)
        {
            // Check for a link
            if (gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_f, pc_counter) >= 0)
            {
                int thick_node_index = (a_n * a_nodes_per_thin_filament) +
                    (m_counter * m_nodes_per_thick_filament) +
                    gsl_vector_short_get(p_mf[m_counter]->pc_node_index, pc_counter);

                int thin_node_index = node_index('a',
                    gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_f, pc_counter),
                    gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_n, pc_counter));

                double g_adjustment =
                    c_k_stiff * (gsl_vector_get(x_trial, thin_node_index) - gsl_vector_get(x_trial, thick_node_index));

                gsl_vector_set(g_vector, thin_node_index,
                    gsl_vector_get(g_vector, thin_node_index) + g_adjustment);

                gsl_vector_set(g_vector, thick_node_index,
                    gsl_vector_get(g_vector, thick_node_index) - g_adjustment);
            }
        }
    }
}

void half_sarcomere::calculate_df_vector(gsl_vector* x_trial)
{
    //! Sets df_vector which is the part of f in kx = f that is
    //  due to titin or cross-bridges and does not depend on x

    // Code

    // Zero df_vector
    gsl_vector_set_zero(df_vector);

    // Add in titin
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        int thick_node_index = (a_n * a_nodes_per_thin_filament) +
            (m_counter * m_nodes_per_thick_filament) +
            t_attach_m_node - 1;

        // Loop through surrounding actins
        for (int a_counter = 0; a_counter < 6; a_counter++)
        {
            int thin_node_index =
                (gsl_matrix_short_get(nearest_actin_matrix, m_counter, a_counter) *
                    a_nodes_per_thin_filament) +
                t_attach_a_node - 1;

            double df_adjustment =
                t_k_stiff * t_slack_length;

            gsl_vector_set(df_vector, thin_node_index,
                gsl_vector_get(df_vector, thin_node_index) - df_adjustment);

            gsl_vector_set(df_vector, thick_node_index,
                gsl_vector_get(df_vector, thick_node_index) + df_adjustment);
        }
    }

    // Add in myosins
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int cb_counter = 0; cb_counter < m_cbs_per_thick_filament; cb_counter++)
        {
            // Check for a link
            if (gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_f, cb_counter) >= 0)
            {
                // Check whether there is an extension
                int cb_state = gsl_vector_short_get(p_mf[m_counter]->cb_state, cb_counter);
                int cb_iso = gsl_vector_short_get(p_mf[m_counter]->cb_iso, cb_counter);
                double ext = p_m_scheme[cb_iso-1]->p_m_states[cb_state-1]->extension;

                if (fabs(ext) > 0.0)
                {
                    int thick_node_index = node_index('m', m_counter, cb_counter);
                    int thin_node_index = node_index('a',
                        gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_f, cb_counter),
                        gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_n, cb_counter));

                    double df_adjustment = m_k_cb * ext;

                    gsl_vector_set(df_vector, thin_node_index,
                        gsl_vector_get(df_vector, thin_node_index) + df_adjustment);

                    gsl_vector_set(df_vector, thick_node_index,
                        gsl_vector_get(df_vector, thick_node_index) - df_adjustment);
                }
            }
        }
    }
}

size_t half_sarcomere::calculate_x_positions()
{
    //! Calculates the x positions by solving kx = f

    // Variables

    size_t no_of_iterations;
    int keep_going;

    double min_deviation;
    double max_deviation;

    // Vectors used in the iterative calculation
    gsl_vector* x_last = gsl_vector_alloc(hs_total_nodes);
    gsl_vector* x_worker = gsl_vector_alloc(hs_total_nodes);
    gsl_vector* x_temp = gsl_vector_alloc(hs_total_nodes);
    gsl_vector* h_vector = gsl_vector_alloc(hs_total_nodes);
    gsl_vector* f_temp = gsl_vector_alloc(hs_total_nodes);

    // Code
    
    // Copy x_vector to x_last for comparison
    gsl_vector_memcpy(x_last, x_vector);

    // Calculate g_vector
    calculate_g_vector(x_vector);

    // Calculate df_vector
    calculate_df_vector(x_vector);

    // Loop until convergence
    no_of_iterations = 0;
    keep_going = 1;
    while (keep_going)
    {
        // Calculate f_temp which is the part of f that does not depend on x
        // f_temp = f0_vector + df_vector
        gsl_vector_memcpy(f_temp, f0_vector);
        gsl_vector_add(f_temp, df_vector);

        // Calculate h = k0\f_temp
        gsl_linalg_solve_tridiag(tri_d_vector, tri_e_vector, tri_f_vector,
            f_temp, h_vector);

        // Now calculate x_worker which is -k0\g + h
        // Store k0/g as x_temp
        gsl_linalg_solve_tridiag(tri_d_vector, tri_e_vector, tri_f_vector,
            g_vector, x_temp);

        // Copy h to x_worker, and subtract x_temp
        gsl_vector_memcpy(x_worker, h_vector);
        gsl_vector_sub(x_worker, x_temp);

        // Find the largest difference between x_worker and x_last
        // Start by copy x_worker to x_temp, then subtract x_last
        gsl_vector_memcpy(x_temp, x_worker);
        gsl_vector_sub(x_temp, x_last);

        // Get the smallest and largest values
        gsl_vector_minmax(x_temp, &min_deviation, &max_deviation);

        // Cope with signs
        max_deviation = GSL_MAX(fabs(min_deviation), fabs(max_deviation));

        // Save x_worker as x_last
        gsl_vector_memcpy(x_last, x_worker);

        // Update iterations
        no_of_iterations = no_of_iterations + 1;

        if (max_deviation < p_fs_options->x_pos_rel_tol)
        {
            keep_going = 0;
        }
        else
        {
            // Update vectors for next loop
            calculate_g_vector(x_worker);
            calculate_df_vector(x_worker);
        }
    }

    // Copy x_worker to x_vector
    gsl_vector_memcpy(x_vector, x_worker);

    // Tidy up
    gsl_vector_free(x_last);
    gsl_vector_free(x_temp);
    gsl_vector_free(x_worker);
    gsl_vector_free(h_vector);
    gsl_vector_free(f_temp);

    return no_of_iterations;
}

void half_sarcomere::initialise_filament_y_z_coordinates(int m_n)
{
    //! Sets the y and z coordinates for the filaments

    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file, "In hs[%i]::initialise_filament_y_z_coordinates\n", hs_id);
    }

    // Check m_n is a square
    if (!gsl_fcmp(sqrt((double)m_n), 1.0, 0.001))
    {
        printf("m_n is not a perfect square. Filament array cannot be constructed\n");
        exit(1);
    }

    // Check for valid m_n
    if (m_n != 9)
    {
        printf("m_n is not == 3. Intialize filament_y_z_coordinates has not been checked\n");
        //exit(1);
    }

    // Variables
    int n_rows = (int)round(sqrt(m_n));         /**< no_of_rows of myosins, y_direction */
    int n_cols = n_rows;                        /**< no of cols of myosins, z direction */

    double y_distance = 0.5;                    /**< distance between actin and myosin
                                                     planes in y direction */
    double z_distance = (sqrt(3.0) / 6.0);      /**< distance between actin and myosin
                                                     planes in z direction */

    // Variables used in loops
    int counter;                                // Counter
    double y;                                   // y coord
    double z;                                   // z coord

    // Code

    // Start with myosins

    counter = 0;
    y = 0.0;                                    // coordinate of first myosin
    z = 0.0;                                    // coordinate of first myosin

    for (int r = 0 ; r < n_rows ; r++)
    {
        for (int c = 0; c < n_cols; c++)
        {
            p_mf[counter]->m_y = y;
            p_mf[counter]->m_z = z;

            y = y + (2.0 * y_distance);
            counter = counter + 1;
        }

        z = z + (3.0 * z_distance);

        if (GSL_IS_ODD(r))
            y = 0.0;
        else
            y = y_distance;
    }

    // Now actins
    counter = 0;
    y = 0.0;
    z = (-2.0 * z_distance);
    int c;
    
    for (int r = 0; r < (2 * n_rows); r++)
    {
        for (c = 0; c < n_cols; c++)
        {
            p_af[counter]->a_y = y;
            p_af[counter]->a_z = z;

            y = y + (2.0 * y_distance);
            counter = counter + 1;
        }

        if (GSL_IS_ODD(r))
            z = z + (2.0 * z_distance);
        else
            z = z + z_distance;

        if (((r % 4) == 0) || ((r % 4) == 1))
            y = y_distance;
        else
            y = 0;
    }
}

void half_sarcomere::initialise_nearest_actin_matrix(void)
{
    //! Sets the nearest actin matrix

    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file,
            "In hs[%i]::initialise_nearest_actin_matrix\n", hs_id);
    }

    // Variables
    int m_rows = (int)round(sqrt(m_n));         /**< no_of_rows of myosins, y_direction */
    int m_cols = m_rows;                        /**< no of cols of myosins, z direction */

    short int** a_m_matrix;                     /**< short int matrix
                                                     positive entries show positions of thick filaments
                                                     negative entries show positions of thin filaments */
    int n_rows = (2 * m_rows) + 1;              /**< size of a_m_matrix */
    int n_cols = (2 * m_rows) + 2;

    int r, c;
    int row_index;
    int col_start;
    int col_index;
    short counter;
    int temp;
    
    // Code

    // Initialize the a_m_matrix
    a_m_matrix = new short int* [n_rows];
    for (int i = 0 ; i < n_rows ; i++)
        a_m_matrix[i] = new short int[n_cols];

    // Zero the array
    for (r = 0 ; r < n_rows ; r++)
    {
        for (c = 0; c < n_cols; c++)
            a_m_matrix[r][c] = 0;
    }

    // Loop through m_rows filling in m_cols of myosins on alternate rows
    counter = 1;
    for (r = 0 ; r < m_rows ; r++)
    {
        if (GSL_IS_EVEN(r))
        {
            if (GSL_IS_EVEN(m_cols))
                col_start = 2;
            else
                col_start = 1;
        }
        else
        {
            if (GSL_IS_EVEN(m_cols))
                col_start = 1;
            else
                col_start = 2;
        }

        // Fill in the entries
        row_index = (2 * r) + 1;

        for (c = 0; c < m_cols; c++)
        {
            col_index = col_start + (c * 2);
            a_m_matrix[row_index][col_index] = counter;
            counter = counter + 1;
        }
    }

    // Now loop through the actins
    counter = 1;
    for (r = 0; r < m_rows; r++)
    {
        row_index = (2 * r);

        for (c = 0 ; c < m_rows ; c++)
        {
            if (GSL_IS_EVEN(r))
                col_start = 1;
            else
                col_start = 2;
            
            col_index = (c * 2) + col_start;

            a_m_matrix[row_index][col_index] = -counter;
            
            counter = counter + 1;
        }
        for (c = 0; c < m_rows; c++)
        {
            if (GSL_IS_EVEN(r))
                col_start = 2;
            else
                col_start = 1;
            
            col_index = (c * 2) + col_start;

            a_m_matrix[row_index][col_index] = -counter;

            counter = counter + 1;
        }
    }
    
    // Now add in mirrors
    for (r = 1; r < (n_rows-1); r++)
    {
        // If there is a myosin in col[1], copy actins from right-hand side
        if (a_m_matrix[r][1] > 0)
        {
            a_m_matrix[r + 1][0] = a_m_matrix[r + 1][n_cols - 2];
            a_m_matrix[r - 1][0] = a_m_matrix[r - 1][n_cols - 2];
        }

        // If there is a myosin in col[end-1], copy actins from left-hand side
        if (a_m_matrix[r][n_cols - 2] > 0)
        {
            a_m_matrix[r + 1][n_cols - 1] = a_m_matrix[r + 1][1];
            a_m_matrix[r - 1][n_cols - 1] = a_m_matrix[r - 1][1];
        }
    }
    
    // Now copy bottom row to top with an offset
    for (c = 0; c < (n_cols - 3); c++)
    {
        a_m_matrix[n_rows - 1][c+1] = a_m_matrix[0][c];
    }
    for (c = 1 ; c < (n_cols - 1); c++)
    {
        a_m_matrix[n_rows - 1][c - 1] = a_m_matrix[0][c];
    }

    if (GSL_IS_EVEN(m_n))
    {
        for (c = 0; c < m_rows; c++) {
            col_index = 2 * c;
            temp = a_m_matrix[n_rows - 1][col_index];
            a_m_matrix[n_rows - 1][col_index] = a_m_matrix[n_rows - 1][col_index + 1];
            a_m_matrix[n_rows - 1][col_index + 1] = temp;
        }
    }
    
    a_m_matrix[n_rows - 1][n_cols - 2] = a_m_matrix[n_rows - 1][0];
    
        
    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file, "hs[%i]: a_m_matrix\n", hs_id);
        for (r = (n_rows-1) ; r >= 0 ; r--)
        {
            for (c = 0; c < n_cols; c++)
            {
                fprintf_s(p_fs_options->log_file, "%5i\t", a_m_matrix[r][c]);
                if (c == (n_cols - 1))
                    fprintf_s(p_fs_options->log_file, "\n");
            }
        }
    }

    // Now we have the a_m_matrix, fill in the nearest_actin_matrix
    for (int thick_counter = 0; thick_counter < m_n; thick_counter++)
    {
        for (r = 0; r < n_rows; r++)
        {
            for (c = 0; c < n_cols; c++)
            {
                if (a_m_matrix[r][c] == (thick_counter + 1))
                {
                    gsl_matrix_short_set(nearest_actin_matrix, thick_counter, 0,
                        -a_m_matrix[r + 1][c] - 1);
                    gsl_matrix_short_set(nearest_actin_matrix, thick_counter, 1,
                        -a_m_matrix[r + 1][c + 1] - 1);
                    gsl_matrix_short_set(nearest_actin_matrix, thick_counter, 2,
                        -a_m_matrix[r - 1][c + 1] - 1);
                    gsl_matrix_short_set(nearest_actin_matrix, thick_counter, 3,
                        -a_m_matrix[r - 1][c] - 1);
                    gsl_matrix_short_set(nearest_actin_matrix, thick_counter, 4,
                        -a_m_matrix[r - 1][c - 1] - 1);
                    gsl_matrix_short_set(nearest_actin_matrix, thick_counter, 5,
                        -a_m_matrix[r + 1][c - 1] - 1);
                }
            }
        }
    }

    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file, "hs[%i]: nearest_actin_matrix\n", hs_id);
        for (r = 0 ; r < m_n ; r++)
        {
            fprintf_s(p_fs_options->log_file, "m_%i: ", r);
            for (c = 0; c < 6; c++)
            {
                fprintf_s(p_fs_options->log_file, "%2i\t",
                    gsl_matrix_short_get(nearest_actin_matrix, r, c));
                if (c == 5)
                    fprintf_s(p_fs_options->log_file, "\n");
            }
        }
    }

    // Recover space
    for (int i = 0; i < n_rows; i++)
        delete[] a_m_matrix[i];
    delete[] a_m_matrix;
}

void half_sarcomere::initialise_f0_vector(void)
{
    //! Initialises bare f0_vector which holds right-hand side of kx = f
    //! without cross-bridges or titin

    // Variables
    int ind;                                    // index

    // First zero f_vector
    gsl_vector_set_zero(f0_vector);

    // Loop through actins
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        ind = node_index('a', a_counter, a_bs_per_thin_filament - 1);
        gsl_vector_set(f0_vector, ind, a_k_stiff * a_inter_bs_rest_length);
    }

    // Loop through myosins
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        ind = node_index('m', m_counter, 0);
        gsl_vector_set(f0_vector, ind,
            m_k_stiff * (hs_length - p_mf[m_counter]->m_lambda));

        ind = node_index('m', m_counter, m_cbs_per_thick_filament - 1);
        gsl_vector_set(f0_vector, ind, (-m_k_stiff * m_inter_crown_rest_length));
    }
}

void half_sarcomere::unpack_x_vector(void)
{
    //! Unpacks x_vector into p_af and p_mf structures

    // Variables
    int row_counter = 0;
    int ind;

    // Loop through actins
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        for (int node_counter = 0; node_counter < a_nodes_per_thin_filament; node_counter++)
        {
            for (int bs_counter = 0; bs_counter < a_bs_per_node; bs_counter++)
            {
                ind = (node_counter * a_bs_per_node) + bs_counter;
                gsl_vector_set(p_af[a_counter]->bs_x, ind, gsl_vector_get(x_vector, row_counter));
            }
            row_counter = row_counter + 1;
        }
    }

    // Loop through myosins
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int node_counter = 0; node_counter < m_nodes_per_thick_filament; node_counter++)
        {
            for (int cb_counter = 0; cb_counter < m_cbs_per_node; cb_counter++)
            {
                ind = (node_counter * m_cbs_per_node) + cb_counter;
                gsl_vector_set(p_mf[m_counter]->cb_x, ind, gsl_vector_get(x_vector, row_counter));
            }
            row_counter = row_counter + 1;
        }
    }
}


void half_sarcomere::calculate_mean_filament_lengths(void)
{
    //! Updates the values of a_mean_fil_length and m_mean_fil_length

    // Variables
    double holder;

    // Code
    
    // First the a filaments
    holder = 0.0;
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        holder = holder + gsl_vector_get(p_af[a_counter]->bs_x,
                                (size_t)a_bs_per_thin_filament - 1);
    }
    a_mean_fil_length = holder / (double)a_n;

    // Now the m filaments
    holder = 0.0;
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        holder = holder +
            (hs_length - gsl_vector_get(p_mf[m_counter]->cb_x,
                            (size_t)m_cbs_per_thick_filament - 1));
    }
    m_mean_fil_length = holder / (double)m_n;
}

void half_sarcomere::calculate_a_pops(void)
{
    //! Updates the vector holding the proportion of 
    // binding sites in each state

    // Variables
    int bs_ind;

    // Code

    // Zero the vector 
    gsl_vector_set_zero(a_pops);
    
    // Loop through thin filaments
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        for (int bs_counter = 0; bs_counter < a_bs_per_thin_filament; bs_counter++)
        {
            bs_ind = gsl_vector_short_get(p_af[a_counter]->bs_state, bs_counter) - 1;
            gsl_vector_set(a_pops, bs_ind,
                gsl_vector_get(a_pops, bs_ind) + 1.0);
        }
    }

    // Turn into proportions
    gsl_vector_scale(a_pops, 1.0 / ((double)a_n * (double)a_bs_per_thin_filament));
}

void half_sarcomere::calculate_m_pops(void)
{
    //! Updates the vector holding the proportion of 
    //! cross-bridges sites in each state

    // Variables
    int cb_ind;

    // Code

    // Zero the vector 
    gsl_vector_set_zero(m_pops);

    // Loop through thick filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int cb_counter = 0; cb_counter < m_cbs_per_thick_filament ; cb_counter++)
        {
            cb_ind = gsl_vector_short_get(p_mf[m_counter]->cb_state, cb_counter) - 1;
            gsl_vector_set(m_pops, cb_ind,
                gsl_vector_get(m_pops, cb_ind) + 1.0);
        }
    }

    // Turn into proportions
    gsl_vector_scale(m_pops, 1.0 / ((double)m_n * (double)m_cbs_per_thick_filament));
}

void half_sarcomere::calculate_c_pops(void)
{
    //! Updates the vector holding the proportion of MyBPC in each state

    // Variables
    int pc_ind;

    // Code

    // Zero the vector
    gsl_vector_set_zero(c_pops);

    // Loop through thick filaments

    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int pc_counter = 0; pc_counter < p_mf[m_counter]->c_no_of_pcs; pc_counter++)
        {
            pc_ind = gsl_vector_short_get(p_mf[m_counter]->pc_state, pc_counter) - 1;
            gsl_vector_set(c_pops, pc_ind,
                gsl_vector_get(c_pops, pc_ind) + 1.0);
        }
    }

    // Turn into proportions
    gsl_vector_scale(c_pops, 1.0 / ((double)m_n * (double)(p_mf[0]->c_no_of_pcs)));
}

double half_sarcomere::calculate_force(void)
{
    //! Calculate the force from the average strain in the last myosin node

    double holder = 0.0;

    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        holder = holder + (m_k_stiff *
            (hs_length - p_mf[m_counter]->m_lambda - m_inter_crown_rest_length -
                gsl_vector_get(x_vector, node_index('m', m_counter, 0))));
    }

    // Adjust for nm scale of filaments, proportion of non-fibrosis muscle,
    // proportion of myofilaments and density of thick filaments
    // Normalize to the number of thick filaments in the calculation
    // Return force in N m^-2
    holder = (holder * (1.0 - p_fs_model->prop_fibrosis) *
                    p_fs_model->prop_myofilaments * p_fs_model->m_filament_density *
                    1e-9 / (double)m_n) +
                hs_extracellular_force;

    // Return
    return holder;
}

double half_sarcomere::calculate_titin_force(void)
{
    //! Calculate the titin contribution to total force 

    double holder = 0.0;

    // Loop through thick filaments

    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        int thick_node_index = (a_n * a_nodes_per_thin_filament) +
            (m_counter * m_nodes_per_thick_filament) +
            t_attach_m_node - 1;

        double x_m = gsl_vector_get(x_vector, thick_node_index);

        // Loop through surrounding actins

        for (int a_counter = 0; a_counter < 6; a_counter++)
        {
            int thin_node_index =
                (gsl_matrix_short_get(nearest_actin_matrix, m_counter, a_counter) *
                    a_nodes_per_thin_filament) +
                t_attach_a_node - 1;

            double x_a = gsl_vector_get(x_vector, thin_node_index);

            if (!strcmp(t_passive_mode, "linear"))
            {
                holder = holder + t_k_stiff * (x_m - x_a - t_slack_length);
            }
        }
    }

    // Adjust for nm scale of filaments and density of thick filaments
    // Normalize to the number of thick filaments in the calculation
    // Return force in N m^-2
    holder = holder * (1.0 - p_fs_model->prop_fibrosis) *
                p_fs_model->prop_myofilaments * p_fs_model->m_filament_density *
                1e-9 / (double)m_n;

    // Return
    return holder;
}

double half_sarcomere::calculate_extracellular_force(void)
{
    //! Calculate the extracellular contribution to total force
    
    // Variables
    double pas_force = 0.0;

    // Code

    if (!strcmp(e_passive_mode, "exponential"))
    {
        if (hs_length >= e_slack_length)
            pas_force = p_fs_model->prop_fibrosis *
                e_sigma * (exp((hs_length - e_slack_length) / e_L) - 1.0);
        else
            pas_force = p_fs_model->prop_fibrosis *
                e_sigma * (exp(-(hs_length - e_slack_length) / e_L) - 1.0);
    }

    return pas_force;
}

void half_sarcomere::update_f0_vector(double delta_hsl)
{
    //! Updates the f0_vector, the part of the right-hand side that holds information
    // on half-sarcomere length and titin does not depend on cross-bridges

    // Code
    if (fabs(delta_hsl) > 0.0)
    {
        // Implement half-sarcomere length change
        for (int m_counter = 0; m_counter < m_n; m_counter++)
        {
            int ind = node_index('m', m_counter, 0);
            gsl_vector_set(f0_vector, ind,
                gsl_vector_get(f0_vector, ind) + (delta_hsl * m_k_stiff));
        }
    }
}

int half_sarcomere::node_index(char molecule_type, int filament_index, int n_index)
{
    /** Return the node index for use in calculating positions */

    int node_index = 0;

    // Code
    switch (molecule_type)
    {
        case 'a':
        {
            // It's an actin binding site
            node_index = (filament_index * a_nodes_per_thin_filament) +
                            (n_index / a_bs_per_node);
            break;
        }
        case 'm':
        {
            // It's a myosin
            node_index = (a_n * a_nodes_per_thin_filament) +
                            (filament_index * m_nodes_per_thick_filament) +
                            (n_index / m_cbs_per_node);
            break;
        }
        case 'c':
        {
            // It's a mybpc
            node_index = (a_n * a_nodes_per_thin_filament) +
                            (filament_index * m_nodes_per_thick_filament) +
                                gsl_vector_short_get(
                                    p_mf[filament_index]->pc_node_index, n_index);
            break;
        }
        default:
        {
            printf("Error, half_sarcomere::node_index called with %c", molecule_type);
            exit(1);
        }
    }

    return node_index;
}

void half_sarcomere::set_cb_nearest_a_f(void)
{
    //! Code sets the nearest_f

    // Variables

    gsl_vector* thin_angles = gsl_vector_alloc(6);
    gsl_vector* angle_differences = gsl_vector_alloc(6);

    size_t nearest_fil_index;

    double temp_diff;

    // Code

    // Set the thin_angles - [0, 60, 120, 180, 240, 300] - these are the angles at which the
    // thin filaments surround each thick filament
    for (int i = 0; i < 6; i++)
        gsl_vector_set(thin_angles, i, (double)i * 60.0);

    // Loop through thick filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int cb_counter = 0; cb_counter < p_mf[m_counter]->m_no_of_cbs; cb_counter++)
        {
            for (int i = 0; i < 6; i++)
            {
                // Calculate the angular difference, where difference between angles a and b is calculated as
                // norm_deg = mod(a-b,360)
                // diff_deg = min(norm_deg, 360-norm_deg)
                temp_diff = fabs(fmod(gsl_vector_get(p_mf[m_counter]->cb_angle, cb_counter) - 
                                    gsl_vector_get(thin_angles, i), 360.0));
                temp_diff = GSL_MIN(temp_diff, 360.0 - temp_diff);
                gsl_vector_set(angle_differences, i, temp_diff);
            }

            // Find the index of the nearest thin filament
            nearest_fil_index = gsl_vector_min_index(angle_differences);

            // Set that
            gsl_vector_short_set(p_mf[m_counter]->cb_nearest_a_f, cb_counter,
                 gsl_matrix_short_get(nearest_actin_matrix, m_counter, nearest_fil_index));
        }
    }

    // Tidy up
    gsl_vector_free(thin_angles);
    gsl_vector_free(angle_differences);
}

void half_sarcomere::set_cb_nearest_a_n(void)
{
    // Code sets the cb_nearest_a_n for each cb_x in p_mf
    // This is a matrix with the first dimension showing cb_indices
    // and the second dimension holding the m_atachment span

    // Variables
    int cb_nearest_a_f;
    int nearest_thin_node;
    int bs_ind;

    double bs_angle;
    double temp_diff;

    // Variables
    gsl_vector* cb_to_thin_node_x = gsl_vector_alloc(a_nodes_per_thin_filament);
    gsl_vector* angle_differences = gsl_vector_alloc(a_bs_per_node);

    // Code

    // Loop through thick filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        // Reset the matrix for the thick filament
        gsl_matrix_short_set_all(p_mf[m_counter]->cb_nearest_a_n, -1);

        for (int cb_counter = 0; cb_counter < p_mf[m_counter]->m_no_of_cbs; cb_counter++)
        {
            cb_nearest_a_f = gsl_vector_short_get(
                p_mf[m_counter]->cb_nearest_a_f, cb_counter);

            // Fill vector with distance from cb to nearest_node
            for (int node_counter = 0; node_counter < a_nodes_per_thin_filament; node_counter++)
            {
                double x1 = gsl_vector_get(p_mf[m_counter]->cb_x, cb_counter);
                double x2 = gsl_vector_get(x_vector,
                    (size_t)((cb_nearest_a_f * a_nodes_per_thin_filament) + node_counter));

                gsl_vector_set(cb_to_thin_node_x, node_counter, fabs(x1 - x2));
            }
            nearest_thin_node = (int)gsl_vector_min_index(cb_to_thin_node_x);

            // Now loop throug the attachment span
            for (int span_i = -adjacent_bs; span_i <= adjacent_bs; span_i++)
            {
                int thin_node_ind = nearest_thin_node + span_i;

                if ((thin_node_ind < 0) ||
                    (thin_node_ind >= a_nodes_per_thin_filament))
                {
                    // Node is not on thin filament
                    continue;
                }

                // There are a_bs_per_node binding sites at this node. Find the one pointing to the cb
                for (int bs_counter = 0; bs_counter < a_bs_per_node; bs_counter++)
                {
                    bs_ind = (thin_node_ind * a_bs_per_node) + bs_counter;
                    bs_angle = gsl_vector_get(p_af[cb_nearest_a_f]->bs_angle, bs_ind);

                    temp_diff = fabs(fmod(gsl_vector_get(p_mf[m_counter]->cb_angle, cb_counter) - bs_angle, 360.0));
                    temp_diff = GSL_MIN(temp_diff, 360.0 - temp_diff);

                    gsl_vector_set(angle_differences, bs_counter, temp_diff);
                }

                // Set the index as the bs that has the biggest angular difference from the cb
                gsl_matrix_short_set(p_mf[m_counter]->cb_nearest_a_n,
                    cb_counter, (span_i + adjacent_bs),
                        (short)((thin_node_ind * a_bs_per_node) +
                            gsl_vector_max_index(angle_differences)));

                // Note the angular difference
                gsl_matrix_set(p_mf[m_counter]->cb_nearest_bs_angle_diff,
                    cb_counter, (size_t)(span_i + adjacent_bs),
                    gsl_vector_max(angle_differences));
            }
        }
    }

    // Tidy up
    gsl_vector_free(cb_to_thin_node_x);
    gsl_vector_free(angle_differences);
}

void half_sarcomere::set_pc_nearest_a_f(void)
{
    //! Code sets the nearest_f

    // Variables

    gsl_vector* thin_angles = gsl_vector_alloc(6);
    gsl_vector* angle_differences = gsl_vector_alloc(6);

    size_t nearest_fil_index;

    double temp_diff;

    // Code

    // Set the thin_angles - [0, 60, 120, 180, 240, 300] - these are the angles at which the
    // thin filaments surround each thick filament
    for (int i = 0; i < 6; i++)
        gsl_vector_set(thin_angles, i, (double)i * 60.0);

    // Loop through thick filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int pc_counter = 0; pc_counter < p_mf[m_counter]->c_no_of_pcs; pc_counter++)
        {
            for (int i = 0; i < 6; i++)
            {
                // Calculate the angular difference, where difference between angles a and b is calculated as
                // norm_deg = mod(a-b,360)
                // diff_deg = min(norm_deg, 360-norm_deg)
                temp_diff = fabs(fmod(gsl_vector_get(p_mf[m_counter]->pc_angle, pc_counter) -
                    gsl_vector_get(thin_angles, i), 360.0));
                temp_diff = GSL_MIN(temp_diff, 360.0 - temp_diff);
                gsl_vector_set(angle_differences, i, temp_diff);
            }

            // Find the index of the nearest thin filament
            nearest_fil_index = gsl_vector_min_index(angle_differences);

            // Set that
            gsl_vector_short_set(p_mf[m_counter]->pc_nearest_a_f, pc_counter,
               gsl_matrix_short_get(nearest_actin_matrix, m_counter, nearest_fil_index));
        }
    }

    // Tidy up
    gsl_vector_free(thin_angles);
    gsl_vector_free(angle_differences);
}

void half_sarcomere::set_pc_nearest_a_n(void)
{
    //! Code sets the pc_nearest_a_n for each pc_x in p_mf

    // Variables
    int pc_nearest_a_f;
    int nearest_thin_node;
    int bs_ind;

    double bs_angle;
    double temp_diff;

    // Variables
    gsl_vector* pc_to_thin_node_x = gsl_vector_alloc(a_nodes_per_thin_filament);
    gsl_vector* angle_differences = gsl_vector_alloc(a_bs_per_node);

    // Loop through thick filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        // Reset the matrix for the thick filament
        gsl_matrix_short_set_all(p_mf[m_counter]->pc_nearest_a_n, -1);

        for (int pc_counter = 0; pc_counter < p_mf[m_counter]->c_no_of_pcs; pc_counter++)
        {
            pc_nearest_a_f = gsl_vector_short_get(
                p_mf[m_counter]->pc_nearest_a_f, pc_counter);

            // Fill vector with distance from cb to nearest_node
            for (int node_counter = 0; node_counter < a_nodes_per_thin_filament; node_counter++)
            {
                int pc_index = (a_n * a_nodes_per_thin_filament) +
                    (m_counter * m_nodes_per_thick_filament) +
                    gsl_vector_short_get(p_mf[m_counter]->pc_node_index, pc_counter);
                double x1 = gsl_vector_get(x_vector, pc_index);
                double x2 = gsl_vector_get(x_vector,
                    (size_t)((pc_nearest_a_f * a_nodes_per_thin_filament) + node_counter));

                gsl_vector_set(pc_to_thin_node_x, node_counter, fabs(x1 - x2));
            }
            nearest_thin_node = (int)gsl_vector_min_index(pc_to_thin_node_x);

            // Now loop through the attachment span
            for (int span_i = -adjacent_bs; span_i <= adjacent_bs; span_i++)
            {
                int thin_node_ind = nearest_thin_node + span_i;

                if ((thin_node_ind < 0) ||
                    (thin_node_ind >= a_nodes_per_thin_filament))
                {
                    // Node is not on filament
                    continue;
                }

                // There are a_bs_per_node_binding sites at this node. Find the one
                // pointing to the pc
                for (int bs_counter = 0; bs_counter < a_bs_per_node; bs_counter++)
                {
                    bs_ind = (thin_node_ind * a_bs_per_node) + bs_counter;
                    bs_angle = gsl_vector_get(p_af[pc_nearest_a_f]->bs_angle, bs_ind);

                    temp_diff = fabs(fmod(gsl_vector_get(p_mf[m_counter]->pc_angle, pc_counter) - bs_angle, 360.0));
                    temp_diff = GSL_MIN(temp_diff, 360.0 - temp_diff);

                    gsl_vector_set(angle_differences, bs_counter, temp_diff);
                }

                // Set the index as the bs that has the biggest angular difference the cb
                gsl_matrix_short_set(p_mf[m_counter]->pc_nearest_a_n,
                    pc_counter, size_t(span_i + adjacent_bs),
                    (thin_node_ind * a_bs_per_node) +
                        (short)gsl_vector_max_index(angle_differences));
                
                // Set the angular difference
                gsl_matrix_set(p_mf[m_counter]->pc_nearest_bs_angle_diff,
                    pc_counter, (size_t)(span_i + adjacent_bs),
                        gsl_vector_max(angle_differences));
            }
        }
    }

    // Tidy up
    gsl_vector_free(pc_to_thin_node_x);
    gsl_vector_free(angle_differences);
}


void half_sarcomere::thick_filament_kinetics(double time_step)
{
    //! Code implements thick filament kinetics

    // Update node forces
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        p_mf[m_counter]->calculate_node_forces();
    }

    myosin_kinetics(time_step);


//    mybpc_kinetics(time_step);
  
}

void half_sarcomere::myosin_kinetics(double time_step)
{
    //! Code implements myosin kienetics

    // Variables

    int cb_state;                       // cb_state (>=1)
    int cb_isotype;                     // cb_isotoype(>=0)

    int new_state;                      // New cb state

    char old_type;                      // Existing state type
    char new_type;                      // New state type

    int transition_index;               // index to an m_transition

    m_state* p_m_state;                 // pointer to a myosin state

    // Code

    // Cycle through filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        // Now cycle through dimers
        if (p_mf[m_counter]->m_myosins_per_hub != 2)
        {
            printf("half_sarcomere::myosin_kinetics() does not work without myosin dimers");
            exit(1);
        }

        for (int cb_counter = 0; cb_counter < m_cbs_per_thick_filament;
            cb_counter = cb_counter + 2)
        {
            transition_index = return_m_transition(time_step, m_counter, cb_counter);

            if (transition_index >= 0)
            {
                // Transition occurred
                cb_state = gsl_vector_short_get(p_mf[m_counter]->cb_state, cb_counter);
                cb_isotype = gsl_vector_short_get(p_mf[m_counter]->cb_iso, cb_counter);
                p_m_state = p_m_scheme[cb_isotype - 1]->p_m_states[cb_state - 1];

                old_type = p_m_state->state_type;

                new_state = p_event[transition_index]->p_trans->new_state;
                new_type = p_m_scheme[cb_isotype - 1]->p_m_states[new_state - 1]->state_type;

                // Implement transition
                handle_lattice_event(p_event[transition_index]);

                // If the head is transitioning into or out of a S state, do the same
                // for the partner head
                if ((old_type == 'S') || (new_type == 'S'))
                {
                    p_event[transition_index]->m_n = p_event[transition_index]->m_n + 1;
                    handle_lattice_event(p_event[transition_index]);
                }
                else
                {
                    if (false)
                    {
                        // Check whether the other head will undergo a transition
                        transition_index = return_m_transition(time_step, m_counter, cb_counter + 1);

                        if (transition_index >= 0)
                        {
                            // Get the potential transition
                            cb_state = gsl_vector_short_get(p_mf[m_counter]->cb_state, (size_t)cb_counter + 1);
                            cb_isotype = gsl_vector_short_get(p_mf[m_counter]->cb_iso, (size_t)cb_counter + 1);
                            p_m_state = p_m_scheme[cb_isotype - 1]->p_m_states[cb_state - 1];

                            old_type = p_m_state->state_type;

                            new_state = p_event[transition_index]->p_trans->new_state;
                            new_type = p_m_scheme[cb_isotype - 1]->p_m_states[new_state - 1]->state_type;

                            // Exclude transitions to or from 'S'
                            if ((old_type != 'S') && (new_type != 'S'))
                            {
                                handle_lattice_event(p_event[transition_index]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void half_sarcomere::mybpc_kinetics(double time_step)
{
    //! Code implements mybpc kinetics

    int transition_index;           // index to a transition
                                    // -1, if nothing happens

    // Code

    // Cycle through filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int pc_counter = 0; pc_counter < p_mf[m_counter]->c_no_of_pcs; pc_counter++)
        {
            // Check for event
            transition_index = return_c_transition(time_step, m_counter, pc_counter);

            if (transition_index >= 0)
            {
                handle_lattice_event(p_event[transition_index]);
            }
        }
    }
}

int half_sarcomere::return_m_transition(double time_step, int m_counter, int cb_counter)
{
    // Code returns the transition a cb undergoes

    // Variables

    int cb_state;                       // cb_state (>=1)
    int cb_isotype;                     // cb_isotoype(>=1)

    int crown_index;                    // index of crown with cb

    int a_f;                            // relevant actin filament
    int a_n;                            // relevant binding site

    int bs_ind;                         // index of binding site

    int a_f_partner;                    // actin filament for the partner head
                                        // will be -1 if unattached

    int mybpc_state;                    // state number for MyBPC controlling cb
    int mybpc_iso;                      // isotype number for MyBPC controlling cb

    m_state* p_m_state;                 // pointer to a myosin state
    transition* p_trans;                // pointer to a transition

    int new_state;                      // new cb_state after transition

    int max_transitions;

    double x;                           // distance between cn and relevant bs
    double node_f;                      // node_force

    gsl_vector* transition_probs;

    double alignment_factor;            // double from 0 to 1 that adjusts
                                        // attachment probability based on angle
                                        // between cb and bs

    double prob;                        // doubles to do with transition probabilities

    double angle;                       // aiignment angle between head and bs

    int prob_index;                     // integer noting the transition index

    int event_index;                    // the event which occurred

    // Code

    // Set values
    cb_state = gsl_vector_short_get(p_mf[m_counter]->cb_state, cb_counter);
    cb_isotype = gsl_vector_short_get(p_mf[m_counter]->cb_iso, cb_counter);
    p_m_state = p_m_scheme[cb_isotype - 1]->p_m_states[cb_state - 1];

    max_transitions = p_m_scheme[cb_isotype - 1]->max_no_of_transitions;

    // Allocate and zero the transition vector
    // This vector is max_transitions * m_attachment_span
    // Binding events to different actin nodes are in different elements
    // Non-binding transitions are single elements
    transition_probs = gsl_vector_alloc((size_t)max_transitions * (size_t)m_attachment_span);
    gsl_vector_set_zero(transition_probs);

    // Get the a_f and the a_n for the myosin head
    if (p_m_state->state_type == 'A')
    {
        a_f = gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_f, cb_counter);
        a_n = gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_n, cb_counter);
    }
    else
    {
        a_f = gsl_vector_short_get(p_mf[m_counter]->cb_nearest_a_f, cb_counter);
    }

    // Get the a_f for the partner dimer
    if (cb_counter < (m_cbs_per_thick_filament - 2))
    {
        a_f_partner = gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_f,
            (size_t)cb_counter + 1);
    }
    else
    {
        a_f_partner = -1;
    }

    // Deduce node force
    crown_index = cb_counter / (p_fs_model->m_hubs_per_crown * p_fs_model->m_myosins_per_hub);
    node_f = gsl_vector_get(p_mf[m_counter]->node_forces, crown_index);

    // Deduce state and isotype of controlling MyBPC
    if (gsl_vector_short_get(p_mf[m_counter]->cb_controlling_pc_index, cb_counter) == -1)
    {
        // Set to 0, as no MyBPC control
        mybpc_state = 0;
        mybpc_iso = 0;
    }
    else
    {
        // Pull the state and isotype
        mybpc_state = gsl_vector_short_get(p_mf[m_counter]->pc_state,
            gsl_vector_short_get(p_mf[m_counter]->cb_controlling_pc_index, cb_counter));
        mybpc_iso = gsl_vector_short_get(p_mf[m_counter]->pc_iso,
            gsl_vector_short_get(p_mf[m_counter]->cb_controlling_pc_index, cb_counter));
    }

    // Prepare for calculating rates
    gsl_vector_set_zero(transition_probs);

    // Cycle through transitions
    for (int t_counter = 0; t_counter < max_transitions; t_counter++)
    {
        p_trans = p_m_state->p_transitions[t_counter];
        new_state = p_trans->new_state;

        if (new_state > 0)
        {
            // It's a possible transition
            if (p_trans->transition_type == 'a')
            {
                // Need to cycle through adjacent_bs, searching for binding sites
                for (int bs_counter = 0; bs_counter < m_attachment_span; bs_counter++)
                {
                    bs_ind = gsl_matrix_short_get(
                        p_mf[m_counter]->cb_nearest_a_n, cb_counter, bs_counter);

                    if ((bs_ind < 0) || (bs_ind >= a_bs_per_thin_filament))
                    {
                        continue;           // binding site is not on filament
                    }

                    if (gsl_vector_short_get(p_af[a_f]->bound_to_m_f, bs_ind) >= 0)
                    {
                        continue;           // binding site is already occupied
                    }

                    if (gsl_vector_short_get(p_af[a_f]->bs_state, bs_ind) == 1)
                    {
                        continue;           // binding site is off
                    }

                    if ((p_m_scheme[cb_isotype - 1]->p_m_states[new_state - 1]->state_type == 'S') &&
                        (a_f_partner >= 0))
                    {
                        continue;           // transition into S state is prevented by partner head
                                            // being attached
                    }

                    // Transition is possible
                    x = gsl_vector_get(p_mf[m_counter]->cb_x, cb_counter) -
                        gsl_vector_get(p_af[a_f]->bs_x, bs_ind);

                    angle = gsl_matrix_get(p_mf[m_counter]->cb_nearest_bs_angle_diff, cb_counter, bs_ind);
                    alignment_factor = -cos(angle * M_PI / 180.0);

                    prob = (1.0 - exp(-time_step * alignment_factor *
                        p_trans->calculate_rate(x, node_f, mybpc_state, mybpc_iso)));

                    // Update the probability vector
                    prob_index = (t_counter * m_attachment_span) + bs_counter;
                    gsl_vector_set(transition_probs, prob_index, prob);

                    // Note the transition
                    p_event[prob_index]->mol_type = 'm';
                    p_event[prob_index]->m_f = m_counter;
                    p_event[prob_index]->m_n = cb_counter;
                    p_event[prob_index]->a_f = a_f;
                    p_event[prob_index]->a_n = bs_ind;
                    p_event[prob_index]->p_trans = p_trans;
                }
            }
            else
            {
                // It's a simpler event, entry goes into a single row

                // If head is attached, calculate a_n, and x
                if (a_f >= 0)
                {
                    a_n = gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_n, cb_counter);
                    x = gsl_vector_get(p_mf[m_counter]->cb_x, cb_counter) -
                        gsl_vector_get(p_af[a_f]->bs_x, a_n);
                }
                else
                {
                    a_n = -1;
                    x = 0.0;
                }

                prob = (1.0 - exp(-time_step *
                    p_trans->calculate_rate(x, node_f, mybpc_state, mybpc_iso)));

                // Update the probability vector
                prob_index = (t_counter * m_attachment_span);
                gsl_vector_set(transition_probs, prob_index, prob);

                // Note the transition
                p_event[prob_index]->mol_type = 'm';
                p_event[prob_index]->m_f = m_counter;
                p_event[prob_index]->m_n = cb_counter;
                p_event[prob_index]->a_f = a_f;
                p_event[prob_index]->a_n = a_n;
                p_event[prob_index]->p_trans = p_trans;
            }
        }
    }

    // Use random number to determine which event (if any) occurred
    event_index = return_event_index(transition_probs);

    // Tidy up
    gsl_vector_free(transition_probs);

    // Return
    return event_index;
}

int half_sarcomere::return_c_transition(double time_step, int m_counter, int pc_counter)
{
    // Code returns the transition a mybpc undergoes

    // Variables
    int c_state;                    // state of mybpc
    int c_isotype;                  // isotype of mybpc

    int a_f;                        // index of actin filament
    int a_n;                        // index of actin binding site

    int bs_ind;                     // index of binding site

    int new_state;                  // state after transition

    int max_transitions;            // max number of transitions

    int prob_index;                 // index in the probability vector

    int event_index;                // -1 if no event occurs
                                    // otherwise, indices to a p_event

    double node_force;              // force at the node the mypbc is
                                    // associated with

    double pc_x;                    // position of mybpc

    double x;                       // distance between mybpc and actin bs

    double angle;                   // angle between mybpc and actin bs

    double alignment_factor;        // allows for angle between mybpc and bs

    double prob;                    // probability of transition

    m_state* p_c_state;             // pointer to a state
    transition* p_trans;            // pointer to a transition

    gsl_vector* transition_probs;   // vector of transition probabilities

    // Code

    // Set values
    c_state = gsl_vector_short_get(p_mf[m_counter]->pc_state, pc_counter);
    c_isotype = gsl_vector_short_get(p_mf[m_counter]->pc_iso, pc_counter);
    
    p_c_state = p_c_scheme[c_isotype - 1]->p_m_states[c_state - 1];

    max_transitions = p_c_scheme[c_isotype - 1]->max_no_of_transitions;

    // Allocate and zero the transition vector
    // This vector is max_transitions * m_attachment_span
    // Binding events to different actin nodes are in different elements
    // Non-binding transitions are single elements
    transition_probs = gsl_vector_alloc((size_t)max_transitions * (size_t)m_attachment_span);
    gsl_vector_set_zero(transition_probs);

    // Allocate and zero the transition vector
    // This vector is max_transitions * m_attachment_span
    // Binding events to different actin nodes are in different elements
    // Non-binding transitions are single elements
    transition_probs = gsl_vector_alloc((size_t)max_transitions * (size_t)m_attachment_span);
    gsl_vector_set_zero(transition_probs);

    // Get the a_f and the a_n for the mybpc
    if (p_c_state->state_type == 'A')
    {
        a_f = gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_f, pc_counter);
        a_n = gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_n, pc_counter);
    }
    else
    {
        a_f = gsl_vector_short_get(p_mf[m_counter]->pc_nearest_a_f, pc_counter);
    }

    // Deduce the node force
    node_force = gsl_vector_get(p_mf[m_counter]->node_forces,
        gsl_vector_short_get(p_mf[m_counter]->pc_node_index, pc_counter));

    // Cycle through transitions
    for (int t_counter = 0; t_counter < max_transitions; t_counter++)
    {
        p_trans = p_c_state->p_transitions[t_counter];
        new_state = p_trans->new_state;

        if (new_state > 0)
        {
            // It's a possible transition
            
            if (p_trans->transition_type == 'a')
            {
                // Need to cycle through adjacent bs, searching for possibilities
                for (int bs_counter = 0; bs_counter < m_attachment_span; bs_counter++)
                {
                    bs_ind = gsl_matrix_short_get(
                        p_mf[m_counter]->pc_nearest_a_n, pc_counter, bs_counter);
                    
                    if ((bs_ind < 0) || (bs_ind >= a_bs_per_thin_filament))
                    {
                        continue;       // binding site is not on filament
                    }

                    if (gsl_vector_short_get(p_af[a_f]->bound_to_m_f, bs_ind) >= 0)
                    {
                        continue;       // binding site is already occupied
                    }

                    if (gsl_vector_short_get(p_af[a_f]->bs_state, bs_ind) == 1)
                    {
                        continue;       // binding site is off
                    }

                    // Transition is possible
                    // Deduce the x and the angle
                    pc_x = gsl_vector_get(x_vector,
                            node_index('c', m_counter, pc_counter));
                    x = pc_x - gsl_vector_get(p_af[a_f]->bs_x, bs_ind);

                    angle = gsl_matrix_get(p_mf[m_counter]->pc_nearest_bs_angle_diff,
                        pc_counter, bs_counter);
                    alignment_factor = -cos(angle * M_PI / 180.0);

                    prob = (1.0 - exp(-time_step * alignment_factor *
                            p_trans->calculate_rate(x, node_force, c_state, c_isotype)));

                    // Update the probability vector
                    prob_index = (t_counter * m_attachment_span) + bs_counter;
                    gsl_vector_set(transition_probs, prob_index, prob);

                    // Note the transition
                    p_event[prob_index]->mol_type = 'c';
                    p_event[prob_index]->m_f = m_counter;
                    p_event[prob_index]->m_n = pc_counter;
                    p_event[prob_index]->a_f = a_f;
                    p_event[prob_index]->a_n = bs_ind;
                    p_event[prob_index]->p_trans = p_trans;
                }
            }
            else
            {
                // It's a simpler event, entry goes into a single row

                // If head is attached, calculate a_n and x
                if (a_f >= 0)
                {
                    a_n = gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_n,
                            pc_counter);
                    pc_x = gsl_vector_get(x_vector,
                            node_index('c', m_counter, pc_counter));
                    x = pc_x - gsl_vector_get(p_af[a_f]->bs_x, a_n);
                }
                else
                {
                    a_n = -1;
                    x = 0.0;
                }

                prob = (1.0 - exp(-time_step *
                    p_trans->calculate_rate(x, node_force, c_state, c_isotype)));

                // Update the probability vector
                prob_index = (t_counter * m_attachment_span);
                gsl_vector_set(transition_probs, prob_index, prob);

                // Note the transition
                p_event[prob_index]->mol_type = 'c';
                p_event[prob_index]->m_f = m_counter;
                p_event[prob_index]->m_n = pc_counter;
                p_event[prob_index]->a_f = a_f;
                p_event[prob_index]->a_n = a_n;
                p_event[prob_index]->p_trans = p_trans;
            }
        }
    }

    // Use random number to determine which event (if any) occurred
    event_index = return_event_index(transition_probs);

    // Tidy up
    gsl_vector_free(transition_probs);

    // Return
    return event_index;
}

void half_sarcomere::handle_lattice_event(lattice_event* p_event)
{
    //! Handles lattice event

    // Variables
    int current_state;
    int new_state;

    // Code

    if (p_event->mol_type == 'm')
    {
        // It's a myosin transition
        
        // Pull the states
        current_state = gsl_vector_short_get(p_mf[p_event->m_f]->cb_state, p_event->m_n);
        new_state = p_event->p_trans->new_state;

        // Set the new one
        gsl_vector_short_set(p_mf[p_event->m_f]->cb_state, p_event->m_n, new_state);

        // Handle lattice interactions
        switch (p_event->p_trans->transition_type)
        {
            case 'a':
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_type, p_event->a_n, 1);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_f, p_event->a_n, p_event->m_f);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_n, p_event->a_n, p_event->m_n);

                gsl_vector_short_set(p_mf[p_event->m_f]->cb_bound_to_a_f, p_event->m_n, p_event->a_f);
                gsl_vector_short_set(p_mf[p_event->m_f]->cb_bound_to_a_n, p_event->m_n, p_event->a_n);

                break;

            case 'd':
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_type, p_event->a_n, 0);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_f, p_event->a_n, -1);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_n, p_event->a_n, -1);

                gsl_vector_short_set(p_mf[p_event->m_f]->cb_bound_to_a_f, p_event->m_n, -1);
                gsl_vector_short_set(p_mf[p_event->m_f]->cb_bound_to_a_n, p_event->m_n, -1);

                break;
        }
    }

    if (p_event->mol_type == 'c')
    {
        // It's a mybpc transition

        // Pull the states
        current_state = gsl_vector_short_get(p_mf[p_event->m_f]->pc_state, p_event->m_n);
        new_state = p_event->p_trans->new_state;

        // Set the new one
        gsl_vector_short_set(p_mf[p_event->m_f]->pc_state, p_event->m_n, new_state);

        // Handle lattice interactions
        switch (p_event->p_trans->transition_type)
        {
            case 'a':
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_type, p_event->a_n, 2);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_f, p_event->a_n, p_event->m_f);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_n, p_event->a_n, p_event->m_n);

                gsl_vector_short_set(p_mf[p_event->m_f]->pc_bound_to_a_f, p_event->m_n, p_event->a_f);
                gsl_vector_short_set(p_mf[p_event->m_f]->pc_bound_to_a_n, p_event->m_n, p_event->a_n);
    
                break;

            case 'd':
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_type, p_event->a_n, 0);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_f, p_event->a_n, -1);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_n, p_event->a_n, -1);

                gsl_vector_short_set(p_mf[p_event->m_f]->pc_bound_to_a_f, p_event->m_n, -1);
                gsl_vector_short_set(p_mf[p_event->m_f]->pc_bound_to_a_n, p_event->m_n, -1);
        }
    }
}

void half_sarcomere::thin_filament_kinetics(double time_step, double Ca_conc)
{
    //! Code implements thin filament kinetics

    // Variables
    int unit_occupied;

    int down_neighbor_status;
    int up_neighbor_status;

    double rate;
    double rand_double;

    double coop_boost;

    gsl_vector_short * bs_indices;

    // Code
    bs_indices = gsl_vector_short_alloc(a_bs_per_unit);

    // Loop through thin filaments
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        // Loop through strands
        for (int str_counter = 0; str_counter < a_strands_per_filament; str_counter++)
        {
            // Loop through regulatory units
            for (int unit = 0; unit < a_regulatory_units_per_strand; unit++)
            {
                int unit_counter = str_counter + (unit * a_strands_per_filament);

                // Deduce the status of the neighbors
                if (unit_counter >= a_strands_per_filament)
                    down_neighbor_status =
                        gsl_vector_short_get(p_af[a_counter]->unit_status,
                                ((size_t)unit_counter - (size_t)a_strands_per_filament));
                else
                    down_neighbor_status = -1;

                if (unit_counter <= (a_strands_per_filament * a_regulatory_units_per_strand - a_strands_per_filament - 1))
                    up_neighbor_status =
                        gsl_vector_short_get(p_af[a_counter]->unit_status,
                            ((size_t)unit_counter + (size_t)a_strands_per_filament));
                else
                    up_neighbor_status = -1;

                // Set the indices for the unit
                p_af[a_counter]->set_regulatory_unit_indices(unit_counter, bs_indices);
                
                if (gsl_vector_short_get(p_af[a_counter]->unit_status, unit_counter) == 1)
                {
                    // Site is off and can turn on
                    coop_boost = 0.0;
                    if (down_neighbor_status == 2)
                        coop_boost = coop_boost + a_k_coop;
                    if (up_neighbor_status == 2)
                        coop_boost = coop_boost + a_k_coop;

                    rate = a_k_on * Ca_conc * (1.0 + coop_boost);

                    // Test event with a random number
                    rand_double = gsl_rng_uniform(rand_generator);

                    if (rand_double > exp(-rate * time_step))
                    {
                        // Unit activates
                        for (int i = 0; i < a_bs_per_unit; i++)
                            gsl_vector_short_set(p_af[a_counter]->bs_state,
                                gsl_vector_short_get(bs_indices, i), 2);
                    }
                }
                else
                {
                    // Site might turn off if it is empty
                    unit_occupied = 0;
                    for (int i = 0; i < a_bs_per_unit; i++)
                        if (gsl_vector_short_get(p_af[a_counter]->bound_to_m_f,
                                gsl_vector_short_get(bs_indices, i)) != -1)
                        {
                            unit_occupied = 1;
                            break;
                        }

                    if (unit_occupied == 0)
                    {
                        coop_boost = 0.0;
                        if (down_neighbor_status == 1)
                            coop_boost = coop_boost + a_k_coop;
                        if (up_neighbor_status == 1)
                            coop_boost = coop_boost + a_k_coop;

                        rate = a_k_off * (1.0 + coop_boost);

                        // Test event with a random number
                        rand_double = gsl_rng_uniform(rand_generator);

                        if (rand_double > exp(-rate * time_step))
                        {
                            // Unit deactivates
                            for (int i = 0; i < a_bs_per_unit; i++)
                                gsl_vector_short_set(p_af[a_counter]->bs_state,
                                    gsl_vector_short_get(bs_indices, i), 1);
                        }
                    }
                }
            }
        }

    }

    // Update unit status for next round
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        p_af[a_counter]->set_unit_status();
    }

    // Tidy up
    gsl_vector_short_free(bs_indices);
}

int half_sarcomere::return_event_index(gsl_vector* prob)
{
    // Returns -1 if no event occurs or the index of an event
    // prob holds individual probabilities

    // Variables
    int n = (int)prob->size;        // the size of the probability array
    int event_index;                // the index of the event

    double holder;                  // used for running total
    double rand_number;             // uniformly distributed between 0 and 1

    gsl_vector* cum_prob;           // cumulative probabiity

    // Code

    // Allocate the cum prob vector and fill the values
    cum_prob = gsl_vector_alloc(n);

    holder = 0.0;
    for (int i = 0; i < n; i++)
    {
        holder = holder + gsl_vector_get(prob, i);
        gsl_vector_set(cum_prob, i, holder);
    }

    // Scale if requried (where the probabilities are very high)
    if (holder > 1.0)
        gsl_vector_scale(cum_prob, 1.0 / holder);

    // Set the event index to -1 (nothing happened)
    event_index = -1;

    // Get a random number, keep looping until it is less than cum_prob
    rand_number = gsl_rng_uniform(rand_generator);
    for (int i = 0; i < n; i++)
    {
        if (rand_number < gsl_vector_get(cum_prob, i))
        {
            event_index = i;
            break;
        }
    }

    // Tidy up
    gsl_vector_free(cum_prob);

    // Return
    return event_index;
}

void half_sarcomere::write_gsl_spmatrix_to_file(gsl_spmatrix* p_sparse_matrix, char output_file_string[])
{
    //! Writes gsl_sparse_triplet_matrix to file

    FILE* output_file = NULL;

    // Check file can be opened, abort if not
    errno_t err = fopen_s(&output_file, output_file_string, "w");
    if (err != 0)
    {
        printf("half_sarcomere::write_gsl_vector_to_file: %s\ncould not be opened\n",
            output_file_string);
        exit(1);
    }

    // Output
    for (size_t r = 0; r < p_sparse_matrix->size1; r++)
    {
        for (size_t c = 0; c < p_sparse_matrix->size2; c++)
        {
            fprintf_s(output_file, "%g",
                gsl_spmatrix_get(p_sparse_matrix, r, c));
            if (c < (p_sparse_matrix->size2 - 1))
                fprintf_s(output_file, "\t");
            else
                fprintf_s(output_file, "\n");
        }
    }

    // Tidy up
    if (output_file != NULL)
        fclose(output_file);
}

void half_sarcomere::write_gsl_vector_to_file(gsl_vector* p_vector, char output_file_string[])
{
    //! Writes gsl_vector to file

    FILE* output_file;

    // Check file can be opened, abort if not
    errno_t err = fopen_s(&output_file, output_file_string, "w");
    if (err != 0)
    {
        printf("half_sarcomere::write_gsl_vector_to_file: %s\ncould not be opened\n",
            output_file_string);
        exit(1);
    }

    gsl_vector_fprintf(output_file, p_vector, "%g");

    // Tidy up
    fclose(output_file);
}


void half_sarcomere::write_hs_status_to_file(char output_file_string[])
{
    /** Writes hs status to file */

    // Variables
    FILE* output_file;
    char temp_string[_MAX_PATH];
    int i;

    // Code

    // Check file can be opened, abort if not
    errno_t err = fopen_s(&output_file, output_file_string, "w");
    if (err != 0)
    {
        printf("Options log file file: %s\ncould not be opened\n",
            output_file_string);
        exit(1);
    }

    fprintf_s(output_file, "{\n\"hs_data\": {\n");
    fprintf_s(output_file, "\t\"hs_id\": %i,\n", hs_id);
    fprintf_s(output_file, "\t\"time\": %g,\n", time_s);
    fprintf_s(output_file, "\t\"hs_length\": %.*g,\n", p_fs_options->dump_precision, hs_length);
    fprintf_s(output_file, "\t\"hs_force\": %.*g,\n", p_fs_options->dump_precision, hs_force);
    fprintf_s(output_file, "\t\"pCa\": %.*g,\n", p_fs_options->dump_precision, pCa);
    fprintf_s(output_file, "\t\"m_nodes_per_thick_filament\": %i,\n",
        m_nodes_per_thick_filament);
    fprintf_s(output_file, "\t\"a_nodes_per_thin_filament\": %i,\n",
        a_nodes_per_thin_filament);

    // CB extension parameters

    fprintf_s(output_file, "\t\"cb_extensions\": [");

    for (int j = 0; j < p_fs_model->m_no_of_isotypes; j++)
    {
        for (i = 0; i < p_m_scheme[j]->no_of_states; i++)
        {
            fprintf_s(output_file, "%g", p_m_scheme[j]->p_m_states[i]->extension);
            if ((i == (p_m_scheme[j]->no_of_states - 1)) &&
                (j == (p_fs_model->m_no_of_isotypes - 1)))
            {
                fprintf_s(output_file, "],\n");
            }
            else
            {
                fprintf_s(output_file, ", ");
            }
        }
    }

    // Titin parameters

    fprintf(output_file, "\"titin\": {\n");
    fprintf(output_file, "\t\"t_k_stiff\": %.*F,\n", p_fs_options->dump_precision, t_k_stiff);
    fprintf(output_file, "\t\"t_slack_length\": %.*F,\n", p_fs_options->dump_precision, t_slack_length);
    fprintf(output_file, "\t\"t_attach_a_node\": %i,\n", t_attach_a_node);
    fprintf(output_file, "\t\"t_attach_m_node\": %i", t_attach_m_node);
    fprintf_s(output_file, "},\n");
    
    fprintf_s(output_file, "\"thick\": [\n");

    for (int thick_counter = 0; thick_counter < m_n; thick_counter++)
    {

        fprintf_s(output_file, "{\n\t\"thick_id\": %i,\n", p_mf[thick_counter]->thick_id);
        fprintf_s(output_file, "\t\"m_y\": %.*F,\n", p_fs_options->dump_precision,
            p_mf[thick_counter]->m_y);
        fprintf_s(output_file, "\t\"m_z\": %.*F,\n", p_fs_options->dump_precision,
            p_mf[thick_counter]->m_z);
        fprintf_s(output_file, "\t\"m_no_of_cbs\": %i,\n", p_mf[thick_counter]->m_no_of_cbs);
        fprintf_s(output_file, "\t\"m_k_stiff\": %.*F,\n", p_fs_options->dump_precision,
            m_k_stiff);
        fprintf_s(output_file, "\t\"m_k_cb\": %.*F, \n", p_fs_options->dump_precision, m_k_cb);
        fprintf_s(output_file, "\t\"c_k_stiff\": %.*F,\n", p_fs_options->dump_precision, c_k_stiff);
        fprintf_s(output_file, "\t\"m_inter_crown_rest_length\": %.*F,\n",
            p_fs_options->dump_precision, m_inter_crown_rest_length);
        fprintf_s(output_file, "\t\"m_cbs_per_node\": %i,\n", m_cbs_per_node);
        fprintf_s(output_file, "\t\"m_lambda\": %.*F,\n", p_fs_options->dump_precision,
            p_mf[thick_counter]->m_lambda);
        fprintf_s(output_file, "\t\"c_no_of_pcs\": %i,\n", c_no_of_pcs);

        sprintf_s(temp_string, _MAX_PATH, "%s", "nearest_actin_filaments");
        JSON_functions::write_gsl_matrix_short_as_JSON_array(
            nearest_actin_matrix,
            output_file, temp_string,
            false, 1);

        sprintf_s(temp_string, _MAX_PATH, "cb_x");
        JSON_functions::write_gsl_vector_as_JSON_array(
            p_mf[thick_counter]->cb_x, output_file,
            temp_string, false, p_fs_options->dump_precision);

        sprintf_s(temp_string, _MAX_PATH, "cb_angle");
        JSON_functions::write_gsl_vector_as_JSON_array(
            p_mf[thick_counter]->cb_angle,
            output_file,
            temp_string, false, p_fs_options->dump_precision);

        sprintf_s(temp_string, _MAX_PATH, "cb_state");

        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->cb_state,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "cb_iso");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->cb_iso,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "cb_bound_to_a_f");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->cb_bound_to_a_f,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "cb_bound_to_a_n");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->cb_bound_to_a_n,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "cb_nearest_a_f");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->cb_nearest_a_f,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "cb_nearest_a_n");
        JSON_functions::write_gsl_matrix_short_as_JSON_array(
            p_mf[thick_counter]->cb_nearest_a_n,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "cb_nearest_bs_angle_diff");
        JSON_functions::write_gsl_matrix_as_JSON_array(
            p_mf[thick_counter]->cb_nearest_bs_angle_diff, output_file,
            temp_string, false, p_fs_options->dump_precision);

        sprintf_s(temp_string, _MAX_PATH, "cb_controlling_pc_index");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->cb_controlling_pc_index,
            output_file,
            temp_string, false);


        sprintf_s(temp_string, _MAX_PATH, "node_forces");
        JSON_functions::write_gsl_vector_as_JSON_array(
            p_mf[thick_counter]->node_forces, output_file,
            temp_string, false, p_fs_options->dump_precision);

        // MyBPC
        sprintf_s(temp_string, _MAX_PATH, "pc_node_index");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->pc_node_index,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "pc_angle");
        JSON_functions::write_gsl_vector_as_JSON_array(
            p_mf[thick_counter]->pc_angle, output_file,
            temp_string, false, p_fs_options->dump_precision);

        sprintf_s(temp_string, _MAX_PATH, "pc_state");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->pc_state,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "pc_iso");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->pc_iso,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "pc_bound_to_a_f");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->pc_bound_to_a_f,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "pc_bound_to_a_n");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->pc_bound_to_a_n,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "pc_nearest_a_f");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_mf[thick_counter]->pc_nearest_a_f,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "pc_nearest_a_n");
        JSON_functions::write_gsl_matrix_short_as_JSON_array(
            p_mf[thick_counter]->pc_nearest_a_n,
            output_file,
            temp_string, true);

        if (thick_counter == (m_n - 1))
        {
            // Last filament
            fprintf_s(output_file, "}\n");
        }
        else
            fprintf_s(output_file, "},\n");
    }
    // Closes thick array
    fprintf_s(output_file, "],\n");

    // Opens thin structure
    fprintf_s(output_file, "\"thin\": [\n");

    for (int thin_counter = 0; thin_counter < a_n; thin_counter++)
    {

        fprintf_s(output_file, "{\n\t\"thin_id\": %i,\n", p_af[thin_counter]->thin_id);
        fprintf_s(output_file, "\t\"a_y\": %.*F,\n", p_fs_options->dump_precision,
            p_af[thin_counter]->a_y);
        fprintf_s(output_file, "\t\"a_z\": %.*F,\n", p_fs_options->dump_precision,
            p_af[thin_counter]->a_z);
        fprintf_s(output_file, "\t\"a_no_of_bs\": %i,\n", p_af[thin_counter]->a_no_of_bs);
        fprintf_s(output_file, "\t\"a_bs_per_node\": %i,\n", a_bs_per_node);
        fprintf_s(output_file, "\t\"a_k_stiff\": %.*F,\n", p_fs_options->dump_precision,
            a_k_stiff);
        fprintf_s(output_file, "\t\"a_inter_bs_rest_length\": %.*F,\n",
            p_fs_options->dump_precision, a_inter_bs_rest_length);

        sprintf_s(temp_string, _MAX_PATH, "bs_x");
        JSON_functions::write_gsl_vector_as_JSON_array(
            p_af[thin_counter]->bs_x, output_file,
            temp_string, false, p_fs_options->dump_precision);

        sprintf_s(temp_string, _MAX_PATH, "bs_angle");
        JSON_functions::write_gsl_vector_as_JSON_array(
            p_af[thin_counter]->bs_angle, output_file,
            temp_string, false, p_fs_options->dump_precision);

        sprintf_s(temp_string, _MAX_PATH, "bs_unit");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_af[thin_counter]->bs_unit,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "bs_state");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_af[thin_counter]->bs_state,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "bs_isoform");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_af[thin_counter]->bs_isoform,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "unit_status");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_af[thin_counter]->unit_status,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "bound_to_m_f");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_af[thin_counter]->bound_to_m_f,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "bound_to_m_n");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_af[thin_counter]->bound_to_m_n,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "bound_to_m_type");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_af[thin_counter]->bound_to_m_type,
            output_file,
            temp_string, false);

        if (thin_counter == (a_n - 1))
        {
            // Last filament
            fprintf_s(output_file, "}\n");
        }
        else
            fprintf_s(output_file, "},\n");
    }
    // Closes thin array
    fprintf_s(output_file, "]\n");
         
    // closes hs_status
    fprintf_s(output_file, "}\n");

    // closes JSON file
    fprintf_s(output_file, "}\n");


    // Tidy up
    fclose(output_file);
}
