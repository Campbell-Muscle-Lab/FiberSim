/**
 * @file    thin_filament.cpp
 * @brief   Source file for the thin_filament class
 * @author  Ken Campbell
 */

#include <cstdio>

#include "thin_filament.h"
#include "half_sarcomere.h"
#include "FiberSim_model.h"
#include "FiberSim_options.h"

#include "gsl_vector.h"

// Constructor
thin_filament::thin_filament(
    FiberSim_model* set_p_fs_model,
    FiberSim_options* set_p_fs_options,
    half_sarcomere* set_p_parent_hs,
    int set_thin_id)
{
    // Initialise

    // Set the pointers
    p_fs_model = set_p_fs_model;
    p_fs_options = set_p_fs_options;
    p_parent_hs = set_p_parent_hs;
    
    // Set the id
    thin_id = set_thin_id;

    // Log
    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file,
            "In constructor for thin_filament[%i] in half_sarcomere constructor[%i]\n",
                thin_id, p_parent_hs->hs_id);
    }

    // Set values from the p_fs_model
    a_regulatory_units_per_strand = p_fs_model->a_regulatory_units_per_strand;
    a_bs_per_unit = p_fs_model->a_bs_per_unit;
    a_strands_per_filament = p_fs_model->a_strands_per_filament;
    a_bs_per_node = p_fs_model->a_strands_per_filament;

    a_inter_bs_rest_length = p_fs_model->a_inter_bs_rest_length;
    a_inter_bs_twist = p_fs_model->a_inter_bs_twist;

    // Set the total number of binding sites on the filament
    a_no_of_bs = a_strands_per_filament * a_regulatory_units_per_strand *
                    a_bs_per_unit;

    // Zero the titin force
    a_titin_force = 0.0;

    // Allocate space for arrays
    bs_x = gsl_vector_alloc(a_no_of_bs);
    bs_angle = gsl_vector_alloc(a_no_of_bs);

    bs_state = gsl_vector_short_alloc(a_no_of_bs);
    bs_isotype = gsl_vector_short_alloc(a_no_of_bs);
    bs_unit = gsl_vector_short_alloc(a_no_of_bs);

    unit_status = gsl_vector_short_alloc(a_regulatory_units_per_strand * a_strands_per_filament);

    bound_to_m_type = gsl_vector_short_alloc(a_no_of_bs);

    bound_to_m_f = gsl_vector_short_alloc(a_no_of_bs);
    bound_to_m_n = gsl_vector_short_alloc(a_no_of_bs);

    active_neighbors = gsl_vector_short_alloc(a_regulatory_units_per_strand * a_strands_per_filament);

    // Initialise arrays

    // Special for bs_x, bs_angle, and bs_unit
    initialise_bs_x_bs_angle_bs_unit();

    // Other arrays are initialized with constants
    gsl_vector_short_set_all(bs_state, 1);
    gsl_vector_short_set_all(bs_isotype, 1);
    gsl_vector_short_set_all(bound_to_m_type, 0);
    gsl_vector_short_set_all(bound_to_m_f, -1);
    gsl_vector_short_set_all(bound_to_m_n, -1);
    gsl_vector_short_set_all(unit_status, 1);
    gsl_vector_short_set_all(active_neighbors, 0);

    // Allocate space for node forces
    node_forces = gsl_vector_alloc(a_regulatory_units_per_strand * a_bs_per_unit);
    gsl_vector_set_zero(node_forces);

    // Allocate space for nearest_thick_filaments
    nearest_thick_filaments = gsl_vector_short_alloc(3);
    gsl_vector_short_set_zero(nearest_thick_filaments);

    // Log
    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file,
            "Finished constructor for thin_filament[%i] in half_sarcomere constructor[%i]\n",
            thin_id, p_parent_hs->hs_id);
    }
}

// Destructor
thin_filament::~thin_filament()
{
    // Tidy up
    if (p_fs_options->log_mode > 0)
    {
        fprintf(p_fs_options->log_file,
            "In destructor for thin_filament[%i] in half_sarcomere constructor[%i]\n",
            thin_id, p_parent_hs->hs_id);
    }

    // Delete gsl_vectors
    gsl_vector_free(bs_x);
    gsl_vector_free(bs_angle);

    // Delete the integer vectors
    gsl_vector_short_free(bs_state);
    gsl_vector_short_free(bs_isotype);
    gsl_vector_short_free(bs_unit);
    gsl_vector_short_free(unit_status);
    gsl_vector_short_free(active_neighbors);
    gsl_vector_short_free(bound_to_m_type);
    gsl_vector_short_free(bound_to_m_f);
    gsl_vector_short_free(bound_to_m_n);

    gsl_vector_free(node_forces);

    gsl_vector_short_free(nearest_thick_filaments);
}

// Functions

void thin_filament::initialise_bs_x_bs_angle_bs_unit(void)
{
    /** Sets the bs_x and bs_angle values for each element in the arrays */

    // Variables
    int bs_counter = 0;

    double x = a_inter_bs_rest_length;
    double base_angle = 0.0;
    double angle = 0.0;
    double inter_site_angle = 0.0;

    // Some allocation
    if (a_strands_per_filament > 0)
        inter_site_angle = 360.0 / (double)a_strands_per_filament;

    // Loop through
    for (int unit_counter = 0 ; unit_counter < a_regulatory_units_per_strand ;
            unit_counter++)
    {
        for (int site_counter = 0; site_counter < a_bs_per_unit; site_counter++)
        {
            for (int strand_counter = 0; strand_counter < a_strands_per_filament;
                strand_counter++)
            {
                // Set the bs_x
                gsl_vector_set(bs_x, bs_counter, x);

                // Adjust the angle by the appropriate inter_strand twist
                // Then set the bs_angle as the double remainder of the angle / 360
                angle = base_angle +
                    ((double)(strand_counter)*inter_site_angle);

                gsl_vector_set(bs_angle, bs_counter, fmod(angle, 360.0));

                // Set the bs_unit, which has values starting at 1
                gsl_vector_short_set(bs_unit, bs_counter,
                    (unit_counter * a_strands_per_filament) + strand_counter + 1);

                // Update counter
                bs_counter = bs_counter + 1;
            }

            // Update x and base_angle
            x = x + a_inter_bs_rest_length;
            base_angle = base_angle - a_inter_bs_twist; // thin fil has a right-handed helix structure
        }
    }
}

void thin_filament::set_unit_status(void)
{
    //! Code sets the status of each unit

    // Variables
    gsl_vector_short* bs_indices;

    // Code

    // Allocate
    bs_indices = gsl_vector_short_alloc(a_bs_per_unit);

    // Code
    for (int unit_counter = 0; unit_counter < (a_strands_per_filament * a_regulatory_units_per_strand);
        unit_counter++)
    {
        // Get the bs_indices for the unit
        set_regulatory_unit_indices(unit_counter, bs_indices);

        // Set the unit status from the first entry
        gsl_vector_short_set(unit_status, unit_counter,
            gsl_vector_short_get(bs_state, gsl_vector_short_get(bs_indices, 0)));
    }

    // Tidy up
    gsl_vector_short_free(bs_indices);
}


void thin_filament::set_regulatory_unit_indices(int unit_ind, gsl_vector_short* bs_indices)
{
    //! Fills the bs_indices array with the indices of the binding sites in the unit

    // Variables

    int units_along = unit_ind / a_strands_per_filament;
    int offset = (units_along * a_bs_per_unit * a_strands_per_filament) + 
        (unit_ind % a_strands_per_filament);

    // Code

    for (int i = 0; i < a_bs_per_unit; i++)
    {
        gsl_vector_short_set(bs_indices, i, offset + (i * a_strands_per_filament));
    }
}

void thin_filament::calculate_node_forces(void)
{
    //! Sets the node force

    // Variable

    int nodes_per_filament = a_regulatory_units_per_strand * a_bs_per_unit;

    double node_force;

    // Code

    for (int node_counter = 0; node_counter < nodes_per_filament; node_counter++)
    {
        int bs_index = node_counter * a_bs_per_node;

        if (node_counter == 0)
        {
            node_force = p_fs_model->a_k_stiff *
                (gsl_vector_get(bs_x, bs_index) - a_inter_bs_rest_length);
        }
        else
        {
            int a, b;
            a = bs_index;
            b = bs_index - a_bs_per_node;
            node_force = p_fs_model->a_k_stiff *
                (gsl_vector_get(bs_x, bs_index) -
                    gsl_vector_get(bs_x, bs_index - a_bs_per_node) -
                    a_inter_bs_rest_length);
        }

        // Set
        gsl_vector_set(node_forces, node_counter, node_force);
    }
}
