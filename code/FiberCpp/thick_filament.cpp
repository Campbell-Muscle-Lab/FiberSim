/**
 * @file    thick_filament.cpp
 * @brief   Source file for the thick_filament class
 * @author  Ken Campbell
 */

#include <cstdio>

#include "thick_filament.h"
#include "half_sarcomere.h"
#include "FiberSim_model.h"
#include "FiberSim_options.h"

#include "gsl_vector.h"
#include "gsl_rng.h"


// Constructor
thick_filament::thick_filament(
    FiberSim_model* set_p_fs_model,
    FiberSim_options* set_p_fs_options,
    half_sarcomere* set_p_parent_hs,
    int set_thick_id)
{
    // Initialise

    // Set the pointers
    p_fs_model = set_p_fs_model;
    p_fs_options = set_p_fs_options;
    p_parent_hs = set_p_parent_hs;
    
    // Set the id
    thick_id = set_thick_id;

    // Log
    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file,
            "In constructor for thick_filament[%i] in half_sarcomere constructor[%i]\n",
                thick_id, p_parent_hs->hs_id);
    }

    // Set values from p_fs_model
    m_crowns_per_filament = p_fs_model->m_crowns_per_filament;
    m_hubs_per_crown = p_fs_model->m_hubs_per_crown;
    m_myosins_per_hub = p_fs_model->m_myosins_per_hub;

    m_inter_crown_rest_length = p_fs_model->m_inter_crown_rest_length;
    m_lambda = p_fs_model->m_lambda;

    m_starting_angle = p_fs_model->m_starting_angle;
    m_inter_crown_twist = p_fs_model->m_inter_crown_twist;
    m_within_hub_twist = p_fs_model->m_within_hub_twist;

    // Calculate number of cross-bridges
    m_no_of_cbs = m_crowns_per_filament * m_hubs_per_crown *
                    m_myosins_per_hub;

    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file,
            "%i cbs in thick filament[%i] in half_sarcomere[%i]\n",
            m_no_of_cbs, thick_id, p_parent_hs->hs_id);
    }

    // Allocate space for myosin arrays
    cb_x = gsl_vector_alloc(m_no_of_cbs);
    cb_angle = gsl_vector_alloc(m_no_of_cbs);

    cb_nearest_bs_angle_diff = gsl_vector_alloc(m_no_of_cbs);

    cb_state = new short int[m_no_of_cbs];
    cb_isoform = new short int[m_no_of_cbs];

    cb_bound_to_a_f = new short int[m_no_of_cbs];
    cb_bound_to_a_n = new short int[m_no_of_cbs];

    cb_nearest_a_f = new short int[m_no_of_cbs];
    cb_nearest_a_n = new short int[m_no_of_cbs];

    cb_controlling_pc_index = new short int[m_no_of_cbs];

    // Allocate space for node_forces
    node_forces = gsl_vector_alloc(m_crowns_per_filament);

    // Initialize arrays

    // Use special function for cb_x and cb_angle
    initalise_cb_x_and_cb_angle();

    // Initialise the random number generator
    const gsl_rng_type* T;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    rand_generator_iso = gsl_rng_alloc(T);
    gsl_rng_set(rand_generator_iso, unsigned long(100));


    // Other arrays are intialized with constants
    for (int i = 0; i < m_no_of_cbs; i++)
    {
        cb_state[i] = 1;

        // Now set cb_isoforms depending on the isoform proportion

        if (GSL_IS_EVEN(i)){ 

            if (p_fs_model->m_no_of_isoforms == 2) {

                double rand_iso = gsl_rng_uniform(rand_generator_iso);
                double first_pop = gsl_vector_get(p_fs_model->isoforms_props,0);

                if (rand_iso < first_pop) {
                    cb_isoform[i] = 0;
                    cb_isoform[i + 1] = 0; // Pairs of myosin heads are the same isoform
                }
                else {
                    cb_isoform[i] = 1;
                    cb_isoform[i+1] = 1; // Pairs of myosin heads are the same isoform
                }
            }

            else {
                double rand_iso = gsl_rng_uniform_int(rand_generator_iso, p_fs_model->m_no_of_isoforms); // random assignment between 3 or more isoforms
                cb_isoform[i] = rand_iso;
                cb_isoform[i + 1] = rand_iso; // Pairs of myosin heads are the same isoform
            }

        }

        cb_bound_to_a_f[i] = -1;
        cb_bound_to_a_n[i] = -1;
        cb_nearest_a_f[i] = -1;
        cb_nearest_a_n[i] = -1;
    }

    // Allocate space for MyBPC arrays
    c_thick_proximal_node = p_fs_model->c_thick_proximal_node;
    c_thick_stripes = p_fs_model->c_thick_stripes;
    c_thick_node_spacing = p_fs_model->c_thick_node_spacing;
    c_mols_per_node = p_fs_model->c_mols_per_node;
    c_starting_angle = p_fs_model->c_starting_angle;
    
    c_no_of_pcs = c_thick_stripes * c_mols_per_node;

    pc_angle = gsl_vector_alloc(c_no_of_pcs);

    pc_node_index = new short int[c_no_of_pcs];
    pc_state = new short int[c_no_of_pcs];
    pc_phos = new short int[c_no_of_pcs];

    pc_bound_to_a_f = new short int[c_no_of_pcs];
    pc_bound_to_a_n = new short int[c_no_of_pcs];

    pc_nearest_a_f = new short int[c_no_of_pcs];
    pc_nearest_a_n = new short int[c_no_of_pcs];

    // Use special function for cb_x and cb_angle
    initalise_pc_node_index_and_pc_angle();

    // Other arrays are initialised with constants
    for (int i = 0; i < c_no_of_pcs; i++)
    {
        pc_state[i] = 1;

        // Set the phoshorylation status

        double rand_iso = gsl_rng_uniform(rand_generator_iso);
        double phos_pop = p_fs_model->c_phos_frac; // get the fraction of phosphorylated MyBPC
        
        // Random distribution of (non-)phosphorylated MyBPC states

        if (rand_iso < phos_pop) { 
            pc_phos[i] = 1; 
        }
        else {
            pc_phos[i] = 0;
        }
        pc_bound_to_a_f[i] = -1;
        pc_bound_to_a_n[i] = -1;
        pc_nearest_a_f[i] = -1;
        pc_nearest_a_n[i] = -1;
    }

    // Initialise cb_controlling_pc_index
    initialise_cb_controlling_pc_index();

    // Log
    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file,
            "Finished constructor for thick_filament[%i] in half_sarcomere constructor[%i]\n",
            thick_id, p_parent_hs->hs_id);
    }
}

// Destructor
thick_filament::~thick_filament()
{
    // Tidy up
    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file,
            "In destructor for thick_filament[%i] in half_sarcomere constructor[%i]\n",
            thick_id, p_parent_hs->hs_id);
    }

    // Delete gsl_vectors
    gsl_vector_free(cb_x);
    gsl_vector_free(cb_angle);
    gsl_vector_free(cb_nearest_bs_angle_diff);

    // Delete the integer vectors
    delete [] cb_state;
    delete [] cb_isoform;
    delete [] cb_bound_to_a_f;
    delete [] cb_bound_to_a_n;
    delete [] cb_nearest_a_f;
    delete [] cb_nearest_a_n;
    delete [] cb_controlling_pc_index;

    // Delete MyBPC
    gsl_vector_free(pc_angle);

    delete [] pc_node_index;
    delete [] pc_state;
    delete [] pc_phos;
    delete [] pc_bound_to_a_f;
    delete [] pc_bound_to_a_n;
    delete [] pc_nearest_a_f;
    delete [] pc_nearest_a_n;

    // Delete node forces
    gsl_vector_free(node_forces);

    gsl_rng_free(rand_generator_iso);
}

void thick_filament::initalise_cb_x_and_cb_angle(void)
{
    /** Sets the cb_x and cb_angle values for each element in the arrays */

    // Variables
    int cb_counter = 0;

    double x = p_parent_hs->hs_length - m_lambda - m_inter_crown_rest_length;
    double angle = m_starting_angle;
    double base_angle = 0.0;
    double inter_hub_angle = 0.0;
    double inter_myosin_angle = m_within_hub_twist;

    // Some allocations
    if (m_hubs_per_crown > 0)
        inter_hub_angle = 360.0 / (double)m_hubs_per_crown;

    // Loop through crowns, hubs, and myosins in order
    for (int crown_counter = 0; crown_counter < m_crowns_per_filament; crown_counter++)
    {
        for (int hub_counter = 0; hub_counter < m_hubs_per_crown; hub_counter++)
        {
            for (int m_counter = 0; m_counter < m_myosins_per_hub; m_counter++)
            {
                // Set the cb_x value
                gsl_vector_set(cb_x, cb_counter, x);

                // Adjust the hug angle by the appropriate
                // inter-hub and inter-myosin twists
                // Then set the cb_angle value as the double remainder
                // of the angle / 360

                angle = base_angle +
                    (double)(hub_counter * inter_hub_angle) +
                    (double)(m_counter * inter_myosin_angle);

                gsl_vector_set(cb_angle, cb_counter, fmod(angle, 360.0));

                // Update counter
                cb_counter = cb_counter + 1;
            }
        }
        // Update x and base_angle for the next crown
        x = x - m_inter_crown_rest_length;
        base_angle = base_angle + m_inter_crown_twist;
    }
}

void thick_filament::initalise_pc_node_index_and_pc_angle(void)
{
    //! Sets the pc_node_index and pc_angle for each element in the array

    // Variables
    int ind;

    double inter_pc_angle = 360.0 / double(c_mols_per_node);

    // Code

    ind = 0;

    for (int stripe_counter = 0; stripe_counter < c_thick_stripes; stripe_counter++)
    {
        for (int pc_counter = 0; pc_counter < c_mols_per_node; pc_counter++)
        {
            // Set the index
            pc_node_index[ind] = c_thick_proximal_node +
                (stripe_counter * c_thick_node_spacing) - 1;

            // Now the angle
            double angle = c_starting_angle + (double)pc_counter * inter_pc_angle;
            gsl_vector_set(pc_angle, ind, angle);

            // Update index
            ind = ind + 1;
        }
    }
}

void thick_filament::initialise_cb_controlling_pc_index(void)
{
    //! Sets cb_controlling_pc_index
    // -1 if there is no MyBPC
    // index (>=0) of MyBPC if there is

    // Variables

    int cb_index;

    int proximal_crown_index;   // index of m_crown closest to M-line
                                // that is influenced by MyBPC

    int distal_crown_index;     // index of m_crown furthest from M-line
                                // that is influenced by MyBPC

    gsl_vector_int* pc_indices = gsl_vector_int_alloc(c_mols_per_node);
    gsl_vector* angle_differences = gsl_vector_alloc(c_mols_per_node);

    // Code

    // Set indices
    proximal_crown_index = c_thick_proximal_node - 1;
    distal_crown_index = proximal_crown_index + ((c_thick_stripes - 1) * c_thick_node_spacing);

    // Loop through cbs

    cb_index = 0;

    for (int crown_counter = 0; crown_counter < m_crowns_per_filament; crown_counter++)
    {
        for (int hub_counter = 0; hub_counter < m_hubs_per_crown; hub_counter++)
        {
            for (int m_counter = 0; m_counter < m_myosins_per_hub; m_counter++)
            {
                if ((crown_counter < proximal_crown_index) ||
                    (crown_counter > (distal_crown_index + c_thick_node_spacing - 1)))
                {
                    cb_controlling_pc_index[cb_index] = -1;
                }
                else
                {
                    // Find the controlling crown
                    int pc_crown_index = proximal_crown_index +
                        (((crown_counter - proximal_crown_index) / c_thick_node_spacing) *
                            c_mols_per_node);

                    // Now find the difference between 
                    int counter = 0;
                    int pc_index = 0;

                    gsl_vector_int_set_zero(pc_indices);
                    gsl_vector_set_zero(angle_differences);

                    int keep_going = 1;
                    while (keep_going)
                    {
                        if (pc_node_index[pc_index] == pc_crown_index)
                        {
                            gsl_vector_int_set(pc_indices, counter, pc_index);

                            double pc_ang = gsl_vector_get(pc_angle, pc_index);
                            double cb_ang = gsl_vector_get(cb_angle, cb_index);

                            double temp_diff = fabs(fmod(pc_ang - cb_ang, 360.0));
                            temp_diff = GSL_MIN(temp_diff, 360.0 - temp_diff);
                            
                            gsl_vector_set(angle_differences, counter, temp_diff);

                            counter = counter + 1;
                        }

                        pc_index = pc_index + 1;

                        if (counter == c_mols_per_node)
                            keep_going = 0;
                    }

                    cb_controlling_pc_index[cb_index] = 
                        gsl_vector_int_get(pc_indices, gsl_vector_min_index(angle_differences));
                }

                cb_index = cb_index + 1;
            }
        }
    }

    // Tidy up
    gsl_vector_int_free(pc_indices);
    gsl_vector_free(angle_differences);
}

void thick_filament::calculate_node_forces(void)
{
    //! Sets the node force

    // Variable
    
    int cb_index;
    int cb_index_spacing = m_hubs_per_crown * m_myosins_per_hub;

    double crown_force;

    // Code

    for (int crown_counter = 0; crown_counter < m_crowns_per_filament; crown_counter++)
    {
        cb_index = crown_counter * cb_index_spacing;

        if (crown_counter == 0)
        {
            crown_force = p_fs_model->m_k_stiff *
                (p_parent_hs->hs_length - m_lambda - m_inter_crown_rest_length -
                    gsl_vector_get(cb_x, cb_index));
        }
        else
        {
            crown_force = p_fs_model->m_k_stiff *
                (gsl_vector_get(cb_x, cb_index - cb_index_spacing) - 
                        gsl_vector_get(cb_x, cb_index) -
                        m_inter_crown_rest_length);
        }

        // Set
        gsl_vector_set(node_forces, crown_counter, crown_force);
    }
}