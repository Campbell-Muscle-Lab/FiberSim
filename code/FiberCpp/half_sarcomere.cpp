/**
 * @file    half_sarcomere.cpp
 * @brief   Source file for the half_sarcomere class
 * @author  Ken Campbell
 */

#include <cstdio>
#include <chrono>
#include <iostream>
#include <regex>

#include "half_sarcomere.h"
#include "thick_filament.h"
#include "kinetic_scheme.h"
#include "transition.h"
#include "m_state.h"
#include "thin_filament.h"
#include "iso_scheme.h"
#include "iso_type.h"
#include "iso_transition.h"
#include "muscle.h"
#include "FiberSim_model.h"
#include "model_hs_variation.h"
#include "FiberSim_options.h"
#include "FiberSim_protocol.h"
#include "FiberSim_data.h"
#include "JSON_functions.h"

#include "gsl_math.h"
#include "gsl_vector.h"
#include "gsl_matrix.h"
#include "gsl_spmatrix.h"
#include "gsl_spblas.h"
#include "gsl_splinalg.h"
#include "gsl_linalg.h"
#include "gsl_rng.h"
#include "gsl_randist.h"
#include "gsl_roots.h"

using namespace::std;

// Structure used for root finding for force control mode
struct force_control_params
{
    double target_force;
    double time_step;
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

    // Loop through potential variations, changing the half-sarcomere's model
    for (int i = 0; i < p_fs_model->no_of_model_hs_variations; i++)
    {
        handle_hs_variation(p_fs_model->p_model_hs_variation[i]);
    }

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
        p_event[i] = new lattice_event{'x',1,-1,-1,-1,NULL};
    }

    // Initialize macroscopic state variables
    time_s = 0.0;
    hs_length = p_fs_model->initial_hs_length;
    hs_force = 0.0;
    pCa = 0.0;
    hs_titin_force = 0.0;
    viscosity = p_fs_model->viscosity;
    hs_viscous_force = 0.0;
    hs_extracellular_force = 0.0;

    hs_inter_hs_titin_force_effect = 0.0;

    // Zero the step_counter
    step_counter = 0;

    // Set the command length
    hs_command_length = hs_length;

    // Set the slack length to default
    hs_slack_length = GSL_NAN;

    // Initialise the random number generator
    // This needs to be done before the thick filaments are allocated to allow for
    // lambda jitter

    // Set time_seed to 0
    long time_seed = 0;

    // If rand_seed is not empty, set the time_seed based on the string, or to an unpredictable
    // value based on the system clock
    if (abs(strcmp(p_fs_options->rand_seed, "")))
    {
        if (!strcmp(p_fs_options->rand_seed, "random"))
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            auto t2 = t1 - floor<std::chrono::seconds>(t1);
            time_seed = (long)std::chrono::duration_cast<std::chrono::microseconds>(t2).count();
        }
        else
        {
            time_seed = atol(p_fs_options->rand_seed);
        }
    }

    const gsl_rng_type* T;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    rand_generator = gsl_rng_alloc(T);
    gsl_rng_set(rand_generator,
        unsigned long(100 * (p_parent_m->muscle_id + 1) + (hs_id + 1) + time_seed));

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

    // Now make vectors to hold the transition probabilities and the
    // cumulative probabilities. We need them before initialising
    // the thick filaments
    // Work out the maximum number of transitions

    max_m_transitions = 0;
    for (int i = 0; i < p_fs_model->m_no_of_isotypes; i++)
    {
        max_m_transitions = GSL_MAX(max_m_transitions, p_fs_model->p_m_scheme[i]->max_no_of_transitions);
    }

    max_c_transitions = 0;
    for (int i = 0; i < p_fs_model->c_no_of_isotypes; i++)
    {
        max_c_transitions = GSL_MAX(max_c_transitions, p_fs_model->p_c_scheme[i]->max_no_of_transitions);
    }

    max_transitions = GSL_MAX(max_m_transitions, max_c_transitions);

    // Now make the vectors
    transition_probs = gsl_vector_alloc(max_transitions * m_attachment_span);

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

    // Now set the nearest thick filaments for each actin
    for (short int a_counter = 0; a_counter < a_n; a_counter++)
    {
        int surrounding_thick_index = 0;

        // Now scan through the nearest_actin_matrix
        for (int thick_counter = 0; thick_counter < m_n; thick_counter++)
        {
            for (int thin_pos = 0; thin_pos < 6; thin_pos++)
            {
                if (gsl_matrix_short_get(nearest_actin_matrix, thick_counter, thin_pos) == a_counter)
                {
                    gsl_vector_short_set(p_af[a_counter]->nearest_thick_filaments,
                        surrounding_thick_index, thick_counter);

                    // Increment the counter
                    surrounding_thick_index = surrounding_thick_index + 1;

                    break;
                }
            }

            // If we have all the surrounding thick filaments, break out
            if (surrounding_thick_index == 3)
                break;
        }
    }

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
    a_gamma_coop = p_fs_model->a_gamma_coop;

    m_no_of_cb_states = p_fs_model->p_m_scheme[0]->no_of_states;
    m_k_stiff = p_fs_model->m_k_stiff;

    c_no_of_pc_states = p_fs_model->p_c_scheme[0]->no_of_states;
    c_no_of_pcs = p_mf[0]->c_no_of_pcs;

    sprintf_s(t_passive_mode, _MAX_PATH, p_fs_model->t_passive_mode);
    t_k_stiff = p_fs_model->t_k_stiff;
    t_offset = p_fs_model->t_offset;
    if (!strcmp(t_passive_mode, "exponential"))
    {
        t_sigma = p_fs_model->t_sigma;
        t_L = p_fs_model->t_L;
    }

    // Extracellular parameters
    sprintf_s(e_passive_mode, _MAX_PATH, p_fs_model->e_passive_mode);
    if (!strcmp(e_passive_mode, "exponential"))
    {
        e_sigma = p_fs_model->e_sigma;
        e_L = p_fs_model->e_L;
    }
    else
    {
        e_k_stiff = p_fs_model->e_k_stiff;
    }
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

    // This is used for force-balance
    original_x_vector = gsl_vector_alloc(hs_total_nodes);

    // Allocate space
    size_t nnz = size_t(5 * hs_total_nodes);
    sp_k_coo_bare = gsl_spmatrix_alloc_nzmax(hs_total_nodes, hs_total_nodes, nnz, GSL_SPMATRIX_COO);
    sp_k_csc_bare = gsl_spmatrix_alloc_nzmax(hs_total_nodes, hs_total_nodes, nnz, GSL_SPMATRIX_CSC);

    sp_F = gsl_vector_alloc(hs_total_nodes);
    sp_G = gsl_vector_alloc(hs_total_nodes);

    // Calculate the x positions

    // First calculate sp_ksc_bare
    calculate_sparse_k_bare();

    // Set the x_vector to an initial approximation
    gsl_linalg_solve_tridiag(tri_d_vector, tri_e_vector, tri_f_vector,
        f0_vector, x_vector);

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

    gsl_vector_free(original_x_vector);

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

    // Delete the transition vectors
    gsl_vector_free(transition_probs);

    // Recover space

    gsl_spmatrix_free(sp_k_coo_bare);
    gsl_spmatrix_free(sp_k_csc_bare);

    gsl_vector_free(sp_F);
    gsl_vector_free(sp_G);
}

// Functions

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

    for (int r = 0; r < n_rows; r++)
    {
        for (int c = 0; c < n_cols; c++)
        {
            p_mf[counter]->m_y = y;
            p_mf[counter]->m_z = z;

            y = y + (2.0 * y_distance);
            counter = counter + 1;
        }

        z = z + (3.0 * z_distance);
        
        if (GSL_IS_ODD(n_rows))
        {
            if (GSL_IS_ODD(r))
                y = 0.0;
            else
                y = y_distance;
        }
        else
        {
            if (GSL_IS_EVEN(r))
                y = -y_distance;
            else
                y = 0.0;
        }
    }

    // Now actins
    if (GSL_IS_ODD(n_rows))
    {
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
    else
    {
        counter = 0;
        y = -y_distance;
        z = (-1.0 * z_distance);
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

            switch (r % 4)
            {
                case 0:
                    z = z - z_distance;
                    break;

                case 1:
                    z = z + 4 * z_distance;
                    break;

                case 2:
                    z = z - z_distance;
                    break;

                case 3:
                    z = z + 4 * z_distance;
                    break;
            }

            if (((r % 4) == 0) || ((r % 4) == 1))
                y = 0;
            else
                y = -y_distance;
        }
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
    for (int i = 0; i < n_rows; i++)
        a_m_matrix[i] = new short int[n_cols];

    // Zero the array
    for (r = 0; r < n_rows; r++)
    {
        for (c = 0; c < n_cols; c++)
            a_m_matrix[r][c] = 0;
    }

    // Loop through m_rows filling in m_cols of myosins on alternate rows
    counter = 1;
    for (r = 0; r < m_rows; r++)
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

        for (c = 0; c < m_rows; c++)
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
    for (r = 1; r < (n_rows - 1); r++)
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
        a_m_matrix[n_rows - 1][c + 1] = a_m_matrix[0][c];
    }
    for (c = 1; c < (n_cols - 1); c++)
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
        for (r = (n_rows - 1); r >= 0; r--)
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
        for (r = 0; r < m_n; r++)
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

void half_sarcomere::handle_hs_variation(model_hs_variation* p_hs_variation)
{
    //! Code updates the model for this half-sarcomere

    // Variables
    string variable = p_hs_variation->model_variable;

    // Code
    if (variable == "m_filament_density")
    {
        p_fs_model->m_filament_density = p_fs_model->m_filament_density *
            gsl_vector_get(p_hs_variation->hs_multiplier, hs_id);
    }

    if (variable == "t_sigma")
    {
        p_fs_model->t_sigma = p_fs_model->t_sigma *
            gsl_vector_get(p_hs_variation->hs_multiplier, hs_id);
    }

    if (variable == "t_k_stiff")
    {
        p_fs_model->t_k_stiff = p_fs_model->t_k_stiff *
            gsl_vector_get(p_hs_variation->hs_multiplier, hs_id);
    }


    if (variable.rfind("m_kinetics", 0) == 0)
    {
        // Starts with m_kinetics
        int isotype;
        int state;
        int transition;
        int parameter_index;
        int no_of_digits = 4;
        int digits[4];

        double y;
        double new_y;

        for (int i = 0; i < no_of_digits; i++)
        {
            digits[i] = 0;
        }

        extract_digits(variable, digits, no_of_digits);

        isotype = digits[0] - 1;
        state = digits[1] - 1;
        transition = digits[2] - 1;
        parameter_index = digits[3] - 1;

        // Pull off value
        y = gsl_vector_get(
            p_fs_model->p_m_scheme[isotype]->p_m_states[state]->p_transitions[transition]->rate_parameters,
            parameter_index);

        new_y = y * gsl_vector_get(p_hs_variation->hs_multiplier, hs_id);

        // Set it
        gsl_vector_set(
            p_fs_model->p_m_scheme[isotype]->p_m_states[state]->p_transitions[transition]->rate_parameters,
            parameter_index,
            new_y);
    }
}

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

    // Some of the kinetics are faster than we need to update the positions,
    // so we can speed up the calculations by updating the kinetics with
    // small sub-steps
    int kinetic_repeats = 10;
    double local_time_step = time_step / (double)kinetic_repeats;

    for (int i = 0; i < kinetic_repeats; i++)
    {
        thin_filament_kinetics(local_time_step, pow(10, -pCa));
        thick_filament_kinetics(local_time_step);
    }

    // Branch depending on sim_mode
    if (sim_mode >= 0.0)
    {
        // Replace the delta_hsl with the delta_hsl found by
        // an iterative search
        delta_hsl = calculate_delta_hsl_for_force(sim_mode, time_step);
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
    hs_force = calculate_force(delta_hsl, time_step);

    // Calculate mean filament lengths
    calculate_mean_filament_lengths();

    // Hold state variables
    calculate_a_pops();
    calculate_m_pops();
    calculate_c_pops();

    // Return
    return x_solve_iterations;
}

void half_sarcomere::sarcomere_kinetics(double time_step, double set_pCa)
{
    //! Code updates the status of the thin, thick, and mybp-c molecules

    // Variables

    // Update the pCa
    pCa = set_pCa;

    // Update isoforms if required
    if (p_fs_model->p_m_iso_scheme != NULL)
        myosin_isotype_kinetics(time_step);
    if (p_fs_model->p_c_iso_scheme != NULL)
        mybpc_isotype_kinetics(time_step);

    // Map the filaments and run kinetics
    set_cb_nearest_a_n();
    set_pc_nearest_a_n();

    // Some of the kinetics are faster than we need to update the positions,
    // so we can speed up the calculations by updating the kinetics with
    // small sub-steps
    int kinetic_repeats = 10;
    double local_time_step = time_step / (double)kinetic_repeats;

    for (int i = 0; i < kinetic_repeats; i++)
    {
        thin_filament_kinetics(local_time_step, pow(10, -pCa));
        thick_filament_kinetics(local_time_step);
    }

    // Hold state variables
    calculate_a_pops();
    calculate_m_pops();
    calculate_c_pops();
}

size_t half_sarcomere::update_lattice(double time_step, double delta_hsl)
{
    //! Updates the positions of the nodes and the forces

    // Variables
    size_t x_solve_iterations;                  // The number of iterations to solve
                                                // the lattice positions

    // Code
    if (fabs(delta_hsl) > 0.0)
    {
        hs_length = hs_length + delta_hsl;
        update_f0_vector(delta_hsl);
    }

    // Calculate positions and deduce force
    x_solve_iterations = calculate_x_positions();

    unpack_x_vector();

    hs_force = calculate_force(delta_hsl, time_step);

    // Calculate mean filament lengths
    calculate_mean_filament_lengths();

    // Return
    return x_solve_iterations;
}

size_t half_sarcomere::update_lattice_for_force(double time_step, double target_force)
{
    //! Updates the hs_length for a target force
    
    // Variables
    size_t lattice_iterations;

    double delta_hsl;

    // Code
    delta_hsl = calculate_delta_hsl_for_force(target_force, time_step);

    lattice_iterations = update_lattice(time_step, delta_hsl);

    // Set the return value for use with threads
    thread_return_value = (double)lattice_iterations;

    // And the main return value
    return lattice_iterations;
}

double half_sarcomere::return_hs_length_for_force(double target_force, double time_step)
{
    //! Returns the half-sarcomere length required for a given force
    // Builds on half_sarcomere::calculate_delta_hsl_for_force

    // Variables
    double delta_hsl;

    // Code
    delta_hsl = calculate_delta_hsl_for_force(target_force, time_step);

    // Return
    return (hs_length + delta_hsl);
}

double half_sarcomere::calculate_delta_hsl_for_force(double target_force, double time_step)
{
    //! Returns the delta_hsl required for the half-sarcomere to generate the target force

    // Variables
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type* T;
    gsl_root_fsolver* s;
    double r = 0.0;
    double x_lo = GSL_MAX(-(hs_length - 10.0), -p_fs_options->hs_force_control_max_delta_hs_length);
    double x_hi = p_fs_options->hs_force_control_max_delta_hs_length;
    struct force_control_params params = { target_force, time_step, this };

    gsl_function F;
    F.function = &test_force_wrapper;
    F.params = &params;

    // Test
    
    double test_value;
    test_value = test_force_wrapper(x_lo, &params);
    //printf("x_lo: %g\t\ttest_value: %g\n", x_lo, test_value);

    test_value = test_force_wrapper(x_hi, &params);
    //printf("x_hi: %g\t\ttest_value: %g\n", x_hi, test_value);
    
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
        status = gsl_root_test_interval(x_lo, x_hi, 0.01, 0);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    return r;
}

double half_sarcomere::test_force_wrapper(double delta_hsl, void* params)
{
    struct force_control_params* p =
        (struct force_control_params*) params;
    half_sarcomere* p_hs = p->p_hs;

    double test_value;

    // Code
    test_value = p_hs->test_force_for_delta_hsl(delta_hsl, params);

    //printf("Wrapper: hs_l: %g\tdelta_hsl: %g\ttest_value: %g\n", p_hs->hs_length, delta_hsl, test_value);

    if (!gsl_finite(test_value))
    {
        if (delta_hsl > 0)
            test_value = DBL_MAX;
        else
            test_value = -DBL_MAX;
    }

    return test_value;
}

double half_sarcomere::test_force_for_delta_hsl(double delta_hsl, void *params)
{
    //! Returns the difference between the calculated force and the target force
    //  for a given delta_hsl

    // Variables
    struct force_control_params* p =
        (struct force_control_params*) params;
    double target_force = p->target_force;
    double time_step = p->time_step;
    double test_force;

    // Current state - need to return to this at the end
    double original_hs_length = hs_length;

    // Code

    // Save the x_vector which we need later
    gsl_vector_memcpy(original_x_vector, x_vector);

    // Apply the length length
    hs_length = hs_length + delta_hsl;
    update_f0_vector(delta_hsl);

    // Solve the x positions and calculate force
    calculate_x_positions();
    test_force = calculate_force(delta_hsl, time_step);

    // Restore the half_sarcomere
    hs_length = original_hs_length;
    gsl_vector_memcpy(x_vector, original_x_vector);
    update_f0_vector(-delta_hsl);

    // Return the difference
    return (test_force - target_force);
}

size_t half_sarcomere::calculate_x_positions()
{
    //! Calculates the x positions by solving kx = f

    // Variables
    size_t no_of_iterations;


    if (p_fs_options->calculate_x_mode == 1)
    {
        no_of_iterations = calculate_x_positions_sparse();
    }
    else
    {
        no_of_iterations = calculate_x_positions_ye_method();
    }

    return no_of_iterations;
}

size_t half_sarcomere::calculate_x_positions_sparse()
{
    //! Uses a sparse iterative approach to calculate x positions

    // Variables
    size_t iterations;

    gsl_vector* x_0;
    gsl_vector* x_worker;
    gsl_vector* x_total;
    gsl_vector* x_last;
    gsl_vector* x_diff;

    gsl_vector* f_rhs;
    gsl_vector* f_k0_x;

    double max_diff;
    double last_max_diff;
    double dx_scale;

    int keep_going;

    // Code

    // Allocate
    x_worker = gsl_vector_alloc(hs_total_nodes);
    x_0 = gsl_vector_alloc(hs_total_nodes);
    x_total = gsl_vector_alloc(hs_total_nodes);
    x_last = gsl_vector_alloc(hs_total_nodes);
    x_diff = gsl_vector_alloc(hs_total_nodes);
    f_rhs = gsl_vector_alloc(hs_total_nodes);
    f_k0_x = gsl_vector_alloc(hs_total_nodes);

    // Calculate x_0
    gsl_linalg_solve_tridiag(tri_d_vector, tri_e_vector, tri_f_vector, f0_vector, x_0);

    // Initialise
    gsl_vector_set_zero(x_worker);
    gsl_vector_memcpy(x_total, x_vector);
    gsl_vector_memcpy(x_last, x_vector);

    // Loop
    keep_going = 1;
    iterations = 0;

    // Set the dx scale
    dx_scale = 1;

    // And the last max_diff
    last_max_diff = GSL_POSINF;

    while (keep_going)
    {
        iterations = iterations + 1;

        calculate_sp_F_and_G(x_total);

        gsl_spblas_dgemv(CblasNoTrans, 1.0, sp_k_csc_bare, x_total, 0, f_k0_x);

        gsl_vector_memcpy(f_rhs, sp_F);
        gsl_vector_sub(f_rhs, f_k0_x);
        gsl_vector_sub(f_rhs, sp_G);

        gsl_linalg_solve_tridiag(tri_d_vector, tri_e_vector, tri_f_vector, f_rhs, x_worker);

        if (dx_scale < 1.0)
        {
            gsl_vector_scale(x_worker, dx_scale);
        }

        gsl_vector_memcpy(x_total, x_last);
        gsl_vector_add(x_total, x_worker);

        // Calculate the deviations
        gsl_vector_memcpy(x_diff, x_last);
        gsl_vector_sub(x_diff, x_total);

        // Calculate the max deviation
        max_diff = GSL_MAX(gsl_vector_max(x_diff), -gsl_vector_min(x_diff));

        // Decide whether to break out
        if (max_diff < p_fs_options->x_pos_rel_tol)
            keep_going = 0;

        if (iterations >= p_fs_options->x_vector_max_iterations)
            keep_going = 0;

        // Update the dx_scale
        if (max_diff > last_max_diff)
            dx_scale = dx_scale * 0.5;
        last_max_diff = max_diff;

        if (keep_going)
            gsl_vector_memcpy(x_last, x_total);
    }

    // End
    gsl_vector_memcpy(x_vector, x_total);

    // Free
    gsl_vector_free(x_0);
    gsl_vector_free(x_worker);
    gsl_vector_free(x_total);
    gsl_vector_free(x_last);
    gsl_vector_free(x_diff);
    gsl_vector_free(f_rhs);
    gsl_vector_free(f_k0_x);

    return iterations;
}

size_t half_sarcomere::calculate_x_positions_ye_method(void)
{
    //! Code calculates the x_positions using the method developed by
    //! Qiang Ye and used until ver 2.3
    //! This is fast but is not stable for simulations with very compliant
    //! filaments

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

    // Calculate df_vector
    calculate_df_vector(x_vector);

    // Calculate f_temp which is the part of f that does not depend on x
    // f_temp = f0_vector + df_vector
    gsl_vector_memcpy(f_temp, f0_vector);
    gsl_vector_add(f_temp, df_vector);

    // Prefill x_temp with our guess
    gsl_vector_memcpy(x_worker, x_vector);

    // Loop until convergence
    no_of_iterations = 0;
    keep_going = 1;
    while (keep_going)
    {
        // Calculate h = k0\f_temp
        gsl_linalg_solve_tridiag(tri_d_vector, tri_e_vector, tri_f_vector,
            f_temp, h_vector);

        // Now calculate x_worker which is -k0\g + h
        // First calculate g
        calculate_g_vector(x_worker);

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
        else if (no_of_iterations >= p_fs_options->x_vector_max_iterations)
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

void half_sarcomere::calculate_sparse_k_bare(void)
{
    //! Code adds elements to a sparse k matrix in triplet form
    //! that account for the thick and thin filaments only
    //! Then converts matrix to csc format

    // Variables
    int row_index;

    // Code

    // First zero matrix
    gsl_spmatrix_set_zero(sp_k_coo_bare);

    // Now cycle through the lattice

    row_index = 0;

    // Loop through the thin filaments
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        for (int node_counter = 0; node_counter < a_nodes_per_thin_filament; node_counter++)
        {
            if (node_counter == 0)
            {
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index, 2.0 * a_k_stiff);
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index + 1, -1.0 * a_k_stiff);
            }
            else if (node_counter == (a_nodes_per_thin_filament - 1))
            {
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index - 1, -1.0 * a_k_stiff);
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index, 1.0 * a_k_stiff);
            }
            else
            {
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index - 1, -1.0 * a_k_stiff);
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index, 2.0 * a_k_stiff);
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index + 1, -1.0 * a_k_stiff);
            }

            row_index = row_index + 1;
        }
    }

    // Loop through the thick filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int node_counter = 0; node_counter < m_nodes_per_thick_filament; node_counter++)
        {
            if (node_counter == 0)
            {
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index, 2.0 * m_k_stiff);
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index + 1, -1.0 * m_k_stiff);
            }
            else if (node_counter == (m_nodes_per_thick_filament - 1))
            {
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index - 1, -1.0 * m_k_stiff);
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index, 1.0 * m_k_stiff);
            }
            else
            {
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index - 1, -1.0 * m_k_stiff);
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index, 2.0 * m_k_stiff);
                gsl_spmatrix_set(sp_k_coo_bare, row_index, row_index + 1, -1.0 * m_k_stiff);
            }

            row_index = row_index + 1;
        }
    }

    // Now convert to csc format
    sp_k_csc_bare = gsl_spmatrix_ccs(sp_k_coo_bare);
}

void half_sarcomere::calculate_sp_F_and_G(gsl_vector* x)
{
    //! Calculates sp_F, a vector that holds the right-hand side of
    //! K_0 X + G(X) = F(X)
    //! and the corresponding G vector

    // Variables

    double temp;
    double f_temp;

    // Code

    // Start by copying across the f0 vector that hold sarcomere length
    gsl_vector_memcpy(sp_F, f0_vector);

    // Now zero the G vector
    gsl_vector_set_zero(sp_G);

    // Start with titin

    // Loop through the myosin filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        // Work out where titin attaches to the thick filament
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

            // Linear portion of G
            f_temp = t_k_stiff * (gsl_vector_get(x, thin_node_index) -
                gsl_vector_get(x, thick_node_index) + t_offset);

            temp = gsl_vector_get(sp_G, thin_node_index);
            gsl_vector_set(sp_G, thin_node_index, temp + f_temp);

            temp = gsl_vector_get(sp_G, thick_node_index);
            gsl_vector_set(sp_G, thick_node_index, temp - f_temp);

            if (!strcmp(t_passive_mode, "exponential"))
            {
                // Exponential portion of G
                double x_diff = gsl_vector_get(x, thick_node_index) -
                    gsl_vector_get(x, thin_node_index) - t_offset;

                f_temp = GSL_SIGN(x_diff) * t_sigma * (exp(fabs(x_diff) / t_L) - 1);

//                f_temp = t_sigma * exp(
//                    (gsl_vector_get(x, thick_node_index) - gsl_vector_get(x, thin_node_index)) / t_L);

                temp = gsl_vector_get(sp_G, thin_node_index);
                gsl_vector_set(sp_G, thin_node_index, temp - f_temp);

                temp = gsl_vector_get(sp_G, thick_node_index);
                gsl_vector_set(sp_G, thick_node_index, temp + f_temp);
            }
        }
    }

    // Now add in cross-bridges
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

                // Get the extension
                int cb_iso = gsl_vector_short_get(p_mf[m_counter]->cb_iso, cb_counter);
                int cb_state = gsl_vector_short_get(p_mf[m_counter]->cb_state, cb_counter);
                double ext = p_m_scheme[cb_iso - 1]->p_m_states[cb_state - 1]->extension;

                if (fabs(ext) > 0.0)
                {
                    // Need to adjust sp_F
                    f_temp = (m_k_cb * ext);

                    temp = gsl_vector_get(sp_F, thin_node_index);
                    gsl_vector_set(sp_F, thin_node_index, temp + f_temp);

                    temp = gsl_vector_get(sp_F, thick_node_index);
                    gsl_vector_set(sp_F, thick_node_index, temp - f_temp);
                }

                // Always need to adjust sp_G
                f_temp = m_k_cb * (gsl_vector_get(x, thin_node_index) -
                    gsl_vector_get(x, thick_node_index));

                temp = gsl_vector_get(sp_G, thin_node_index);
                gsl_vector_set(sp_G, thin_node_index, temp + f_temp);

                temp = gsl_vector_get(sp_G, thick_node_index);
                gsl_vector_set(sp_G, thick_node_index, temp - f_temp);
            }
        }
    }

    // Now add in myosin binding protein-C
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int pc_counter = 0; pc_counter < p_mf[m_counter]->c_no_of_pcs; pc_counter++)
        {
            // Check for a link
            if (gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_f, pc_counter) >= 0)
            {
                int thick_node_index = node_index('c', m_counter, pc_counter);
                int thin_node_index = node_index('a',
                    gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_f, pc_counter),
                    gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_n, pc_counter));

                // Get the extension
                int pc_iso = gsl_vector_short_get(p_mf[m_counter]->pc_iso, pc_counter);
                int pc_state = gsl_vector_short_get(p_mf[m_counter]->pc_state, pc_counter);
                double ext = p_c_scheme[pc_iso - 1]->p_m_states[pc_state - 1]->extension;

                if (fabs(ext) > 0.0)
                {
                    // Need to adjust sp_F

                    f_temp = (c_k_stiff * ext);

                    temp = gsl_vector_get(sp_F, thin_node_index);
                    gsl_vector_set(sp_F, thin_node_index, temp + f_temp);

                    temp = gsl_vector_get(sp_F, thick_node_index);
                    gsl_vector_set(sp_F, thick_node_index, temp - f_temp);
                }

                // Always need to adjust sp_G
                f_temp = c_k_stiff * (gsl_vector_get(x, thin_node_index) -
                    gsl_vector_get(x, thick_node_index));

                temp = gsl_vector_get(sp_G, thin_node_index);
                gsl_vector_set(sp_G, thin_node_index, temp + f_temp);

                temp = gsl_vector_get(sp_G, thick_node_index);
                gsl_vector_set(sp_G, thick_node_index, temp - f_temp);
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

    /*// Add in titin
    // Titins are added in g vector
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

            double df_adjustment = 0.0;

            gsl_vector_set(df_vector, thin_node_index,
                gsl_vector_get(df_vector, thin_node_index) - df_adjustment);

            gsl_vector_set(df_vector, thick_node_index,
                gsl_vector_get(df_vector, thick_node_index) + df_adjustment);
        }
    }*/

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
                double ext = p_m_scheme[cb_iso - 1]->p_m_states[cb_state - 1]->extension;

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

    // Add in mybpc
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int pc_counter = 0; pc_counter < c_no_of_pcs; pc_counter++)
        {
            // Check for a link
            if (gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_f, pc_counter) >= 0)
            {
                // Check whether there is an extension
                int pc_state = gsl_vector_short_get(p_mf[m_counter]->pc_state, pc_counter);
                int pc_iso = gsl_vector_short_get(p_mf[m_counter]->pc_iso, pc_counter);
                double ext = p_c_scheme[pc_iso - 1]->p_m_states[pc_state - 1]->extension;

                if (fabs(ext) > 0.0)
                {
                    int thick_node_index = node_index('c', m_counter, pc_counter);
                    int thin_node_index = node_index('a',
                        gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_f, pc_counter),
                        gsl_vector_short_get(p_mf[m_counter]->pc_bound_to_a_n, pc_counter));

                    double df_adjustment = c_k_stiff * ext;

                    gsl_vector_set(df_vector, thin_node_index,
                        gsl_vector_get(df_vector, thin_node_index) + df_adjustment);

                    gsl_vector_set(df_vector, thick_node_index,
                        gsl_vector_get(df_vector, thick_node_index) - df_adjustment);
                }
            }
        }
    }
};

void half_sarcomere::calculate_g_vector(gsl_vector* x_trial)
{
    //! Sets g_vector which is the part of f in kx = f that is
    //  cross-bridges, or mybpc and depends on x
    // and all the titin component because it is easier to put it in here

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

            double x_a = gsl_vector_get(x_trial, thin_node_index);
            double x_m = gsl_vector_get(x_trial, thick_node_index);

            // There is always a linear component
            g_adjustment = t_k_stiff * (x_a - x_m + t_offset);

            if (!strcmp(t_passive_mode, "exponential"))
            {
                g_adjustment = g_adjustment -
                    GSL_SIGN(x_m - x_a) * t_sigma * (exp(fabs(x_m - x_a - t_offset) / t_L) - 1);
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
                int thick_node_index = node_index('c', m_counter, pc_counter);

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

double half_sarcomere::calculate_force(double delta_hsl, double time_step)
{
    //! Calculate the force from the average strain in the last myosin node

    double holder = 0.0;

    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        holder = holder + (m_k_stiff *
            (hs_length - p_mf[m_counter]->m_lambda - m_inter_crown_rest_length -
                gsl_vector_get(x_vector, node_index('m', m_counter, 0))));
    }

    hs_titin_force = calculate_titin_force();
    hs_viscous_force = calculate_viscous_force(delta_hsl, time_step);
    hs_extracellular_force = calculate_extracellular_force();

    // Adjust for nm scale of filaments, proportion of non-fibrosis muscle,
    // proportion of myofilaments and density of thick filaments
    // Normalize to the number of thick filaments in the calculation
    // Return force in N m^-2
    holder = (holder * (1.0 - p_fs_model->prop_fibrosis) *
                    p_fs_model->prop_myofilaments * p_fs_model->m_filament_density *
                    1e-9 / (double)m_n) +
                hs_viscous_force +
                hs_extracellular_force;

    // Return
    return holder;
}

double half_sarcomere::calculate_viscous_force(double delta_hsl, double time_step)
{
    //! Calculates viscous force
   
    // Variables
    double vf = 0.0;

    // Code

    vf = viscosity * delta_hsl / time_step;

    // Return
    return vf;
}

double half_sarcomere::calculate_titin_force(void)
{
    //! Calculate the titin contribution to total force
    //! Updates the titin_force associated with each thin filament as it goes

    // Variables

    double t_strand_force = 0.0;            // force in a single strand of titin

    double holder = 0.0;

    // Code

    // First, zero all the titin forces
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        p_af[a_counter]->a_titin_force = 0.0;
    }

    // Now loop through thick filaments

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

            // There is always a linear force
            t_strand_force = t_k_stiff * (x_m - x_a - t_offset);

            if (!strcmp(t_passive_mode, "exponential"))
            {
                t_strand_force = t_strand_force +
                    GSL_SIGN(x_m - x_a - t_offset) * t_sigma * (exp(fabs(x_m - x_a - t_offset) / t_L) - 1);
            }

            // Update
            holder = holder + t_strand_force;

            // Adjust a_titin_force for the filament
            p_af[gsl_matrix_short_get(nearest_actin_matrix, m_counter, a_counter)]->a_titin_force =
                p_af[gsl_matrix_short_get(nearest_actin_matrix, m_counter, a_counter)]->a_titin_force + t_strand_force;
        }
    }

    // Adjust for nm scale of filaments
    // Normalize to the number of thick filaments in the calculation
    holder = holder * 1e-9 / (double)m_n;

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
    
    if (!strcmp(e_passive_mode, "linear"))
    {
        pas_force = p_fs_model->prop_fibrosis *
            e_k_stiff * (hs_length - e_slack_length);
    }

    return pas_force;
}

void half_sarcomere::update_f0_vector(double delta_hsl)
{
    //! Updates the f0_vector, the part of the right-hand side that holds information
    // on half-sarcomere length

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

    size_t nearest_thin_node;
    size_t bs_ind;

    double bs_angle;
    double temp_diff;

    // Variables
    gsl_vector* cb_to_thin_node_x = gsl_vector_alloc(a_nodes_per_thin_filament);
    gsl_vector* angle_differences = gsl_vector_alloc(a_bs_per_node);

    // Code

    // Loop through thick filaments
    for (size_t m_counter = 0; m_counter < m_n; m_counter++)
    {
        // Reset the matrix for the thick filament
        gsl_matrix_short_set_all(p_mf[m_counter]->cb_nearest_a_n, -1);

        for (int cb_counter = 0; cb_counter < p_mf[m_counter]->m_no_of_cbs; cb_counter++)
        {
            cb_nearest_a_f = gsl_vector_short_get(
                p_mf[m_counter]->cb_nearest_a_f, cb_counter);

            // Fill vector with distance from cb to nearest_node
            for (size_t node_counter = 0; node_counter < a_nodes_per_thin_filament; node_counter++)
            {
                double x1 = gsl_vector_get(p_mf[m_counter]->cb_x, cb_counter);
                double x2 = gsl_vector_get(x_vector,
                    ((size_t)cb_nearest_a_f * (size_t)a_nodes_per_thin_filament) +
                        node_counter);

                gsl_vector_set(cb_to_thin_node_x, node_counter, fabs(x1 - x2));
            }
            nearest_thin_node = gsl_vector_min_index(cb_to_thin_node_x);

            // Now loop throug the attachment span
            for (size_t span_i = 0; span_i < m_attachment_span; span_i++)
            {
                int thin_node_ind = (int)nearest_thin_node - adjacent_bs + (int)span_i;

                if ((thin_node_ind < 0) ||
                    (thin_node_ind >= a_nodes_per_thin_filament))
                {
                    // Node is not on thin filament
                    continue;
                }

                // There are a_bs_per_node binding sites at this node. Find the one pointing to the cb
                for (size_t bs_counter = 0; bs_counter < a_bs_per_node; bs_counter++)
                {
                    bs_ind = ((size_t)thin_node_ind * (size_t)a_bs_per_node) + bs_counter;
                    bs_angle = gsl_vector_get(p_af[cb_nearest_a_f]->bs_angle, bs_ind);

                    temp_diff = fabs(fmod(gsl_vector_get(p_mf[m_counter]->cb_angle, cb_counter) - bs_angle, 360.0));
                    temp_diff = GSL_MIN(temp_diff, 360.0 - temp_diff);

                    gsl_vector_set(angle_differences, bs_counter, temp_diff);
                }

                // Set the index as the bs that has the biggest angular difference from the cb
                gsl_matrix_short_set(p_mf[m_counter]->cb_nearest_a_n,
                    cb_counter, span_i,
                        (short)((thin_node_ind * a_bs_per_node) +
                            gsl_vector_max_index(angle_differences)));

                // Note the angular difference
                gsl_matrix_set(p_mf[m_counter]->cb_nearest_bs_angle_diff,
                    cb_counter, span_i,
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
            double c_jitter = gsl_ran_gaussian(rand_generator, p_fs_model->c_sd_angle_deviation);

            for (int i = 0; i < 6; i++)
            {
                // Calculate the angular difference, where difference between angles a and b is calculated as
                // norm_deg = mod(a-b,360)
                // diff_deg = min(norm_deg, 360-norm_deg)
                // Add a random jitter to the pc_angle
                double pc_angle = gsl_vector_get(p_mf[m_counter]->pc_angle, pc_counter) +
                    c_jitter;
                temp_diff = fabs(fmod(pc_angle - gsl_vector_get(thin_angles, i), 360.0));
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
    
    size_t nearest_thin_node;
    size_t bs_ind;

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
            for (size_t node_counter = 0; node_counter < a_nodes_per_thin_filament; node_counter++)
            {
                int pc_index = (a_n * a_nodes_per_thin_filament) +
                    (m_counter * m_nodes_per_thick_filament) +
                    gsl_vector_short_get(p_mf[m_counter]->pc_node_index, pc_counter);
                double x1 = gsl_vector_get(x_vector, pc_index);
                double x2 = gsl_vector_get(x_vector,
                    (((size_t)pc_nearest_a_f * (size_t)a_nodes_per_thin_filament) +
                        node_counter));

                gsl_vector_set(pc_to_thin_node_x, node_counter, fabs(x1 - x2));
            }
            nearest_thin_node = gsl_vector_min_index(pc_to_thin_node_x);

            // Now loop through the attachment span
            for (size_t span_i = 0; span_i < m_attachment_span; span_i++)
            {
                int thin_node_ind = (int)nearest_thin_node - adjacent_bs + (int)span_i;

                if ((thin_node_ind < 0) ||
                    (thin_node_ind >= a_nodes_per_thin_filament))
                {
                    // Node is not on filament
                    continue;
                }

                // There are a_bs_per_node_binding sites at this node. Find the one
                // pointing to the pc
                for (size_t bs_counter = 0; bs_counter < a_bs_per_node; bs_counter++)
                {
                    bs_ind = ((size_t)thin_node_ind * (size_t)a_bs_per_node) + bs_counter;
                    bs_angle = gsl_vector_get(p_af[pc_nearest_a_f]->bs_angle, bs_ind);

                    temp_diff = fabs(fmod(gsl_vector_get(p_mf[m_counter]->pc_angle, pc_counter) - bs_angle, 360.0));
                    temp_diff = GSL_MIN(temp_diff, 360.0 - temp_diff);

                    gsl_vector_set(angle_differences, bs_counter, temp_diff);
                }

                // Set the index as the bs that has the biggest angular difference the cb
                gsl_matrix_short_set(p_mf[m_counter]->pc_nearest_a_n,
                    pc_counter, span_i,
                    (short)((thin_node_ind * a_bs_per_node) +
                        gsl_vector_max_index(angle_differences)));
                
                // Set the angular difference
                gsl_matrix_set(p_mf[m_counter]->pc_nearest_bs_angle_diff,
                    pc_counter, span_i,
                        gsl_vector_max(angle_differences));
            }
        }
    }

    // Tidy up
    gsl_vector_free(pc_to_thin_node_x);
    gsl_vector_free(angle_differences);
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

    double coop_boost = 0.0;

    int thin_filament_repeat = 1;

    gsl_vector_short* bs_indices;

    // Thin filament kinetics are faster than myosin kinetics, so run this multiple times
    // with a short time-step

    double local_time_step = time_step / (double)thin_filament_repeat;

    // Code
    bs_indices = gsl_vector_short_alloc(a_bs_per_unit);

    // Update node forces
    for (int a_counter = 0; a_counter < a_n; a_counter++)
    {
        p_af[a_counter]->calculate_node_forces();
    }

    for (int rep_counter = 0; rep_counter < thin_filament_repeat; rep_counter++)
    {

        // Loop through thin filaments setting active neighbors
        for (int a_counter = 0; a_counter < a_n; a_counter++)
        {
            // Loop through strands, setting the active neighbors
            for (int str_counter = 0; str_counter < a_strands_per_filament; str_counter++)
            {
                // Loop through regulatory units
                for (int unit = 0; unit < a_regulatory_units_per_strand; unit++)
                {
                    int unit_counter = str_counter + (unit * a_strands_per_filament);

                    short int active_neigh = 0;

                    // Deduce the status of the neighbors
                    if (unit_counter >= a_strands_per_filament)
                    {
                        down_neighbor_status =
                            gsl_vector_short_get(p_af[a_counter]->unit_status,
                                ((size_t)unit_counter - (size_t)a_strands_per_filament));

                        if (down_neighbor_status == 2)
                            active_neigh = active_neigh + 1;
                    }

                    if (unit_counter <= (a_strands_per_filament * a_regulatory_units_per_strand - a_strands_per_filament - 1))
                    {
                        up_neighbor_status =
                            gsl_vector_short_get(p_af[a_counter]->unit_status,
                                ((size_t)unit_counter + (size_t)a_strands_per_filament));

                        if (up_neighbor_status == 2)
                            active_neigh = active_neigh + 1;
                    }

                    // Set active_neighbors value
                    gsl_vector_short_set(p_af[a_counter]->active_neighbors, unit_counter, active_neigh);
                }
            }

            // Loop through strands
            for (int str_counter = 0; str_counter < a_strands_per_filament; str_counter++)
            {
                // Loop through regulatory units
                for (int unit = 0; unit < a_regulatory_units_per_strand; unit++)
                {
                    int unit_counter = str_counter + (unit * a_strands_per_filament);

                    // Set the indices for the unit
                    p_af[a_counter]->set_regulatory_unit_indices(unit_counter, bs_indices);

                    if (gsl_vector_short_get(p_af[a_counter]->unit_status, unit_counter) == 1)
                    {
                        // Site is off and can turn on
                        coop_boost = a_gamma_coop *
                            (double)gsl_vector_short_get(p_af[a_counter]->active_neighbors, unit_counter);

                        rate = a_k_on * Ca_conc * (1.0 + coop_boost);

                        // Test event with a random number
                        rand_double = gsl_rng_uniform(rand_generator);

                        if (rand_double > exp(-rate * local_time_step))
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
                            coop_boost = a_gamma_coop *
                                (2.0 - (double)gsl_vector_short_get(p_af[a_counter]->active_neighbors, unit_counter));

                            rate = a_k_off * (1.0 + coop_boost);

                            // Test event with a random number
                            rand_double = gsl_rng_uniform(rand_generator);

                            if (rand_double > exp(-rate * local_time_step))
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

    }

    // Tidy up
    gsl_vector_short_free(bs_indices);
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

    mybpc_kinetics(time_step);
}

void half_sarcomere::myosin_kinetics(double time_step)
{
    //! Code implements myosin kinetics

    // Variables

    int cb_state;                       // cb_state (>=1)
    int cb_isotype;                     // cb_isotoype(>=0)

    int new_state;                      // New cb state

    char old_type;                      // Existing state type
    char new_type;                      // New state type

    bool s_allowed;                     // transition involving S-state allowed

    int transition_index;               // index to an m_transition

    m_state* p_m_state;                 // pointer to a myosin state

    int bs_state;                        // state of closest bs

    bool t_allowed;                     // transition involving odd heads and the SRX state

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

        // Reset the nearest BS state matrix for the thick filament
        gsl_matrix_short_set_all(p_mf[m_counter]->cb_nearest_a_n_states, -1);

        // Set the number of transitions to zero for each myosin
        //gsl_vector_int_set_zero(s_trans_count);

        for (int cb_counter = 0; cb_counter < m_cbs_per_thick_filament; cb_counter++)
        {            
            // For each CB, save the state of closest BS 
            // right before evaluating the transition events
            int a_f = gsl_vector_short_get(p_mf[m_counter]->cb_nearest_a_f, cb_counter);

            if (GSL_IS_EVEN(cb_counter))
                t_allowed = true; // prepare the flag for tracking SRX transitions

            for (size_t span_i = 0; span_i < m_attachment_span; span_i++)
            {
                int bs_ind = gsl_matrix_short_get(
                    p_mf[m_counter]->cb_nearest_a_n, cb_counter, span_i);

                if ((bs_ind < 0) || (bs_ind >= a_bs_per_thin_filament))
                {
                    bs_state = -1;          // binding site is not on filament
                }

                else if (gsl_vector_short_get(p_af[a_f]->bound_to_m_n, bs_ind) >= 0)
                {
                    bs_state = -1;          // binding site is occupied
                }

                else 
                {
                    bs_state = gsl_vector_short_get(p_af[a_f]->bs_state, bs_ind);
                                             // binding site state is OFF (1) or ON (2)
                }

                gsl_matrix_short_set(p_mf[m_counter]->cb_nearest_a_n_states,
                    cb_counter, span_i, bs_state); // update cb_nearest_a_n_states element
            }

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

                // Does it involve the S state
                if ((old_type == 'S') || (new_type == 'S'))
                {
                    s_allowed = true;

                    if (GSL_IS_ODD(cb_counter))
                        s_allowed = false;  
                    else 
                    {
                        if (cb_counter < (m_cbs_per_thick_filament - 1))
                        {
                            if (cb_state !=
                                gsl_vector_short_get(p_mf[m_counter]->cb_state,
                                    (size_t)cb_counter + 1))
                            {
                                // dimers have different states
                                s_allowed = false;
                            }
                        }
                    }

                    if (s_allowed == true)
                    {
                        // Implement transition for first head
                        handle_lattice_event(p_event[transition_index]);

                        // And the second
                        p_event[transition_index]->m_n = p_event[transition_index]->m_n + 1;
                        handle_lattice_event(p_event[transition_index]);

                        t_allowed = false; // this odd head cannot transition anymore within this time-step

                    }
                }
                else
                {
                    // Transition does not involve the S state

                    if (t_allowed == true) // odd head that already transitionned from/to SRX cannot transition a second time
                    {
                        handle_lattice_event(p_event[transition_index]);
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

    int bs_state;                       // state of closest bs

    // Code

    // Cycle through filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        // Reset the nearest BS state matrix for the thick filament
        gsl_matrix_short_set_all(p_mf[m_counter]->pc_nearest_a_n_states, -1);

        for (int pc_counter = 0; pc_counter < p_mf[m_counter]->c_no_of_pcs; pc_counter++)
        {
            
            // For each pc, save the state of closest BS 
            // right before evaluating the transition events

            int a_f = gsl_vector_short_get(p_mf[m_counter]->pc_nearest_a_f, pc_counter);

            for (size_t span_i = 0; span_i < m_attachment_span; span_i++)
            {
                int bs_ind = gsl_matrix_short_get(
                    p_mf[m_counter]->pc_nearest_a_n, pc_counter, span_i);

                if ((bs_ind < 0) || (bs_ind >= a_bs_per_thin_filament))
                {
                    bs_state = -1;           // binding site is not on filament
                }

                else if (gsl_vector_short_get(p_af[a_f]->bound_to_m_n, bs_ind) >= 0)
                {
                    bs_state = -1;
                }

                else
                {
                    bs_state = gsl_vector_short_get(p_af[a_f]->bs_state, bs_ind);
                    // binding site state is OFF (1) or ON (2)
                }

                gsl_matrix_short_set(p_mf[m_counter]->pc_nearest_a_n_states,
                    pc_counter, span_i, bs_state); // update pc_nearest_a_n_states element
            }
            
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
    // Code returns the transition for a partner head

    // Variables
    
    int cb_state;                       // cb_state (>=1)
    int cb_isotype;                     // cb_isotoype(>=1)

    int crown_index;                    // index of crown with cb

    int a_f;                            // relevant actin filament
    int a_n;                            // relevant binding site

    int bs_ind;                         // index of binding site

    int mybpc_state;                    // state number for MyBPC controlling cb
    int mybpc_iso;                      // isotype number for MyBPC controlling cb

    m_state* p_m_state;                 // pointer to a myosin state
    transition* p_trans;                // pointer to a transition

    int new_state;                      // new cb_state after transition

    double x;                           // distance between cn and relevant bs
    double x_ext;                       // cb state extension
    double node_f;                      // node_force

    short active_neigh;                 // number of active neighbors

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

    // Zero the transition vector
    // This vector is max_transitions * m_attachment_span
    // Binding events to different actin nodes are in different elements
    // Non-binding transitions are single elements

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
    for (int t_counter = 0; t_counter < max_m_transitions; t_counter++)
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

                    // Transition is possible
                    x = gsl_vector_get(p_mf[m_counter]->cb_x, cb_counter) -
                        gsl_vector_get(p_af[a_f]->bs_x, bs_ind);

                    angle = gsl_matrix_get(p_mf[m_counter]->cb_nearest_bs_angle_diff,
                        cb_counter, bs_counter);
                    alignment_factor = -cos(angle * M_PI / 180.0);

                    // Get cb extension
                    x_ext = p_m_state->extension;

                    prob = (1.0 - exp(-time_step * alignment_factor *
                        p_trans->calculate_rate(x, x_ext, node_f, mybpc_state, mybpc_iso, 0, this)));

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
                if (p_m_state->state_type == 'A')
                {
                    a_n = gsl_vector_short_get(p_mf[m_counter]->cb_bound_to_a_n, cb_counter);
                    x = gsl_vector_get(p_mf[m_counter]->cb_x, cb_counter) -
                        gsl_vector_get(p_af[a_f]->bs_x, a_n);

                    // Get the number of active neighbors
                    short bs_u = gsl_vector_short_get(p_af[a_f]->bs_unit, a_n);

                    active_neigh = gsl_vector_short_get(p_af[a_f]->active_neighbors, (size_t)(bs_u - 1));
                }
                else
                {
                    a_n = -1;
                    x = 0.0;
                    active_neigh = (short)0;
                }

                // Get cb extension
                x_ext = p_m_state->extension;

                prob = (1.0 - exp(-time_step *
                    p_trans->calculate_rate(x, x_ext, node_f, mybpc_state, mybpc_iso, active_neigh, this)));

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

    int prob_index;                 // index in the probability vector

    int event_index;                // -1 if no event occurs
                                    // otherwise, indices to a p_event

    double node_force;              // force at the node the mypbc is
                                    // associated with

    double pc_x;                    // position of mybpc

    double x;                       // distance between mybpc and actin bs

    double x_ext;                   // state extension

    double angle;                   // angle between mybpc and actin bs

    double alignment_factor;        // allows for angle between mybpc and bs

    double prob;                    // probability of transition

    m_state* p_c_state;             // pointer to a state
    transition* p_trans;            // pointer to a transition

    // Code

    // Set values
    c_state = gsl_vector_short_get(p_mf[m_counter]->pc_state, pc_counter);
    c_isotype = gsl_vector_short_get(p_mf[m_counter]->pc_iso, pc_counter);
    
    p_c_state = p_c_scheme[c_isotype - 1]->p_m_states[c_state - 1];

    // Zero the transition vector
    // This vector is max_transitions * m_attachment_span
    // Binding events to different actin nodes are in different elements
    // Non-binding transitions are single elements
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
    for (int t_counter = 0; t_counter < max_c_transitions; t_counter++)
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

                    // Get extension

                    x_ext = p_c_state->extension;

                    prob = (1.0 - exp(-time_step * alignment_factor *
                            p_trans->calculate_rate(x, x_ext, node_force, c_state, c_isotype, 0, this)));

                    if (gsl_isnan(prob))
                    {
                        printf("node_force: %g\n", node_force);
                        printf("isnan, stopping\n");
                        exit(1);
                    }

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

                // Get cb extension
                x_ext = p_c_state->extension;

                prob = (1.0 - exp(-time_step *
                    p_trans->calculate_rate(x, x_ext, node_force, c_state, c_isotype, 0, this)));

                if (gsl_isnan(prob))
                {
                    printf("isnan - single level, stopping\n");
                    exit(1);
                }

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

    // Return
    return event_index;
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

    gsl_vector* cum_prob;           // GSL vector holding the cumulative probability

    // Code

    cum_prob = gsl_vector_alloc(n);

    holder = 0.0;
    for (int i = 0; i < n; i++)
    {
        holder = holder + gsl_vector_get(prob, i);
        gsl_vector_set(cum_prob, i, holder);
    }

    // Scale if required (where the probabilities are very high)
    if ((holder > 1.0) && (rescaling_flag == 0))
    {
        // Flag
        printf("Probabilities are re-scaled\n");
        gsl_vector_scale(cum_prob, 1.0 / holder);
        rescaling_flag = 1;
    }

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

void half_sarcomere::handle_lattice_event(lattice_event* p_event)
{
    //! Handles lattice event

    // Variables
    int new_state;

    // Code

    if (p_event->mol_type == 'm')
    {
        // It's a myosin transition
        
        // Pull the states
        new_state = p_event->p_trans->new_state;

        // Set the new one
        gsl_vector_short_set(p_mf[p_event->m_f]->cb_state, p_event->m_n, new_state);

        // Handle lattice interactions
        if (p_event->p_trans->transition_type == 'a')
        {
            gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_type, p_event->a_n, 1);
            gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_f, p_event->a_n, p_event->m_f);
            gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_n, p_event->a_n, p_event->m_n);

            gsl_vector_short_set(p_mf[p_event->m_f]->cb_bound_to_a_f, p_event->m_n, p_event->a_f);
            gsl_vector_short_set(p_mf[p_event->m_f]->cb_bound_to_a_n, p_event->m_n, p_event->a_n);
        }

        if (p_event->p_trans->transition_type == 'd')
        {
            gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_type, p_event->a_n, 0);
            gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_f, p_event->a_n, -1);
            gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_n, p_event->a_n, -1);

            gsl_vector_short_set(p_mf[p_event->m_f]->cb_bound_to_a_f, p_event->m_n, -1);
            gsl_vector_short_set(p_mf[p_event->m_f]->cb_bound_to_a_n, p_event->m_n, -1);
        }
    }

    if (p_event->mol_type == 'c')
    {
        // It's a mybpc transition

        // Pull the states
        new_state = p_event->p_trans->new_state;

        // Set the new one
        gsl_vector_short_set(p_mf[p_event->m_f]->pc_state, p_event->m_n, new_state);

        // Handle lattice interactions
        if (p_event->p_trans->transition_type == 'a')
        {

            gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_type, p_event->a_n, 2);
            gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_f, p_event->a_n, p_event->m_f);
            gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_n, p_event->a_n, p_event->m_n);

            gsl_vector_short_set(p_mf[p_event->m_f]->pc_bound_to_a_f, p_event->m_n, p_event->a_f);
            gsl_vector_short_set(p_mf[p_event->m_f]->pc_bound_to_a_n, p_event->m_n, p_event->a_n);
        }

        if (p_event->p_trans->transition_type == 'd')
        {
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_type, p_event->a_n, 0);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_f, p_event->a_n, -1);
                gsl_vector_short_set(p_af[p_event->a_f]->bound_to_m_n, p_event->a_n, -1);

                gsl_vector_short_set(p_mf[p_event->m_f]->pc_bound_to_a_f, p_event->m_n, -1);
                gsl_vector_short_set(p_mf[p_event->m_f]->pc_bound_to_a_n, p_event->m_n, -1);
        }
    }
}

void half_sarcomere::myosin_isotype_kinetics(double time_step)
{
    //! Code deduces whether a myosin head will switch isotype

    // Variables
    iso_scheme* p_iso_scheme;           // pointer to the m_iso_scheme
    iso_transition* p_iso_trans;        // pointer to an isotype transition

    int m_iso;                          // isotype of the a myosin head
    int event_index;                    // index of isoform transition event

    double x;                           // the x position of each cross-bridge

    double prob;                        // probability of an event

    // Code
    p_iso_scheme = p_fs_model->p_m_iso_scheme;

    // Cycle through filaments
    for (int m_counter = 0; m_counter < m_n; m_counter++)
    {
        for (int cb_counter = 0; cb_counter < m_cbs_per_thick_filament; cb_counter++)
        {
            // Zero the transition vector
            gsl_vector_set_zero(p_iso_scheme->transition_probs);

            m_iso = gsl_vector_short_get(p_mf[m_counter]->cb_iso, cb_counter);

            for (int t_counter = 0; t_counter < p_iso_scheme->p_iso_types[m_iso - 1]->no_of_transitions;
                t_counter++)
            {
                // Calculate the probability of the event
                // Set the pointer
                p_iso_trans = p_iso_scheme->p_iso_types[m_iso - 1]->p_iso_transitions[t_counter];

                // And the cross-bridge position x
                x = gsl_vector_get(p_mf[m_counter]->cb_x, cb_counter);

                // Calculate it
                prob = (1.0 - exp(-time_step * p_iso_trans->calculate_rate(x, this)));

                // Assign
                gsl_vector_set(p_iso_scheme->transition_probs, t_counter, prob);
            }

            // Now use a random number to determine which event (if any) occurred
            event_index = return_event_index(p_iso_scheme->transition_probs);

            // Implement event if required
            if (event_index >= 0)
            {
                gsl_vector_short_set(p_mf[m_counter]->cb_iso, cb_counter,
                    (short)(p_iso_scheme->p_iso_types[m_iso - 1]->p_iso_transitions[event_index]->new_iso_type));
            }
        }
    }
}

void half_sarcomere::mybpc_isotype_kinetics(double time_step)
{
    //! Handles mybpc isotype switching

    // Variables

    // Code
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
            fprintf_s(output_file, "%.8f",
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

    gsl_vector_fprintf(output_file, p_vector, "%.8f");

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
        printf("Dump file: %s\ncould not be opened\n",
            output_file_string);
        exit(1);
    }

    fprintf_s(output_file, "{\n\"hs_data\": {\n");
    fprintf_s(output_file, "\t\"hs_id\": %i,\n", hs_id);
    fprintf_s(output_file, "\t\"time\": %g,\n", time_s);
    fprintf_s(output_file, "\t\"hs_length\": %.*g,\n", p_fs_options->dump_precision, hs_length);
    fprintf_s(output_file, "\t\"hs_force\": %.*g,\n", p_fs_options->dump_precision, hs_force);
    fprintf_s(output_file, "\t\"prop_fibrosis\": %.*g,\n", p_fs_options->dump_precision, p_fs_model->prop_fibrosis);
    fprintf_s(output_file, "\t\"prop_myofilaments\": %.*g,\n", p_fs_options->dump_precision, p_fs_model->prop_myofilaments);
    fprintf_s(output_file, "\t\"m_filament_density\": %.*g,\n", p_fs_options->dump_precision, p_fs_model->m_filament_density);
    fprintf_s(output_file, "\t\"temperature\": %.*g,\n", p_fs_options->dump_precision, p_fs_model->temperature);
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
                fprintf_s(output_file, "]\n");
            }
            else
            {
                fprintf_s(output_file, ", ");
            }
        }
    }
    fprintf_s(output_file, "},\n");

    // Titin parameters

    fprintf(output_file, "\"titin\": {\n");
    fprintf(output_file, "\t\"t_passive_mode\": \"%s\",\n", t_passive_mode);
    fprintf(output_file, "\t\"t_k_stiff\": %.*F,\n", p_fs_options->dump_precision, t_k_stiff);
    fprintf(output_file, "\t\"t_offset\": %.*F,\n", p_fs_options->dump_precision, t_offset);
    fprintf(output_file, "\t\"t_sigma\": %.*F,\n", p_fs_options->dump_precision, t_sigma);
    fprintf(output_file, "\t\"t_L\": %.*F,\n", p_fs_options->dump_precision, t_L);
    fprintf(output_file, "\t\"t_attach_a_node\": %i,\n", t_attach_a_node);
    fprintf(output_file, "\t\"t_attach_m_node\": %i\n", t_attach_m_node);
    fprintf(output_file, "},\n");
    
    fprintf(output_file, "\"thick\": [\n");

    for (int thick_counter = 0; thick_counter < m_n; thick_counter++)
    {

        fprintf_s(output_file, "{\n\t\"thick_id\": %i,\n", p_mf[thick_counter]->thick_id);
        fprintf_s(output_file, "\t\"m_y\": %.*F,\n", p_fs_options->dump_precision,
            p_mf[thick_counter]->m_y);
        fprintf_s(output_file, "\t\"m_z\": %.*F,\n", p_fs_options->dump_precision,
            p_mf[thick_counter]->m_z);
        fprintf_s(output_file, "\t\"m_no_of_cbs\": %i,\n", p_mf[thick_counter]->m_no_of_cbs);
        fprintf_s(output_file, "\t\"m_no_of_isotypes\": %i,\n", p_fs_model->m_no_of_isotypes);
        fprintf_s(output_file, "\t\"m_no_of_states\": %i,\n", p_fs_model->m_no_of_cb_states);
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
        fprintf_s(output_file, "\t\"c_no_of_isotypes\": %i,\n", p_fs_model->c_no_of_isotypes);
        fprintf_s(output_file, "\t\"c_no_of_states\": %i,\n", p_fs_model->c_no_of_pc_states);

        /** Whole "nearest actin filaments" matrix is now dumped in the log file 
        
        sprintf_s(temp_string, _MAX_PATH, "%s", "nearest_actin_filaments");
        JSON_functions::write_gsl_matrix_short_as_JSON_array(
            nearest_actin_matrix,
            output_file, temp_string,
            false, 1);

        */

        fprintf_s(output_file, "\t\"nearest_actin_filaments\": [");
        for (int i = 0; i < 6; i++)
        {
            fprintf_s(output_file, "%i", gsl_matrix_short_get(nearest_actin_matrix, thick_counter, i));
            if (i < 5)
                fprintf_s(output_file, ", ");
            else
                fprintf_s(output_file, "],\n");
        }
        
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

        sprintf_s(temp_string, _MAX_PATH, "cb_nearest_a_n_states");
        JSON_functions::write_gsl_matrix_short_as_JSON_array(
            p_mf[thick_counter]->cb_nearest_a_n_states,
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
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "pc_nearest_a_n_states");
        JSON_functions::write_gsl_matrix_short_as_JSON_array(
            p_mf[thick_counter]->pc_nearest_a_n_states,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "pc_nearest_bs_angle_diff");
        JSON_functions::write_gsl_matrix_as_JSON_array(
            p_mf[thick_counter]->pc_nearest_bs_angle_diff, output_file,
            temp_string, true, p_fs_options->dump_precision);



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
        fprintf_s(output_file, "\t\"m_no_of_isotypes\": %i,\n", p_fs_model->a_no_of_bs_isotypes);
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
            p_af[thin_counter]->bs_isotype,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "unit_status");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_af[thin_counter]->unit_status,
            output_file,
            temp_string, false);

        sprintf_s(temp_string, _MAX_PATH, "active_neighbors");
        JSON_functions::write_gsl_vector_short_as_JSON_array(
            p_af[thin_counter]->active_neighbors,
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
            temp_string, true);

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

    // Tidy up
    fclose(output_file);
}

void half_sarcomere::extract_digits(string test_string, int digits[], int no_of_digits)
{
    //! Function fills digits array with [0-9] extracted from string
    //! See https://en.cppreference.com/w/cpp/regex/regex_iterator

    // Variables
    int counter;

    // Code
    regex ex("[0-9]");

    auto digits_begin = sregex_iterator(test_string.begin(), test_string.end(), ex);
    auto digits_end = sregex_iterator();

    counter = 0;
    for (regex_iterator i = digits_begin; i != digits_end; i++)
    {
        smatch match = *i;
        digits[counter] = atoi((match.str().c_str()));
        counter = counter + 1;
    }
}
