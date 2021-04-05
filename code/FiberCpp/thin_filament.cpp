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
    a_bs_per_node = p_parent_hs->a_bs_per_node;
    a_strands_per_filament = p_fs_model->a_strands_per_filament;

    a_inter_bs_rest_length = p_fs_model->a_inter_bs_rest_length;
    a_inter_bs_twist = p_fs_model->a_inter_bs_twist;

    // Set the total number of binding sites on the filament
    a_no_of_bs = a_strands_per_filament * a_regulatory_units_per_strand *
                    a_bs_per_unit;

    // Allocate space for arrays
    bs_x = gsl_vector_alloc(a_no_of_bs);
    bs_angle = gsl_vector_alloc(a_no_of_bs);

    bs_state = new short int[a_no_of_bs];
    bs_isoform = new short int[a_no_of_bs];
    bs_unit = new short int[a_no_of_bs];

    unit_status = new short int[a_regulatory_units_per_strand * a_strands_per_filament];

    bound_to_m_type = new short int[a_no_of_bs];

    bound_to_m_f = new short int[a_no_of_bs];
    bound_to_m_n = new short int[a_no_of_bs];

    nearest_m_f = new short int[a_no_of_bs];
    nearest_m_n = new short int[a_no_of_bs];

    // Initialise arrays

    // Special for bs_x, bs_angle, and bs_unit
    initialise_bs_x_bs_angle_bs_unit();

    // Other arrays are initialized with constants
    for (int i = 0; i < a_no_of_bs; i++)
    {
        bs_state[i] = 0;
        bs_isoform[i] = 0;

        bound_to_m_type[i] = 0;

        bound_to_m_f[i] = -1;
        bound_to_m_n[i] = -1;
        nearest_m_f[i] = -1;
        nearest_m_n[i] = -1;
    }

    for (int i = 0; i < (a_regulatory_units_per_strand * a_strands_per_filament); i++)
    {
        unit_status[i] = 0;
    }

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
    delete[] bs_state;
    delete[] bs_isoform;
    delete[] bs_unit;

    delete[] unit_status;

    delete[] bound_to_m_type;

    delete[] bound_to_m_f;
    delete[] bound_to_m_n;
    delete[] nearest_m_f;
    delete[] nearest_m_n;
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
                bs_unit[bs_counter] = (unit_counter * a_strands_per_filament) +
                                        strand_counter + 1;

                // Update counter
                bs_counter = bs_counter + 1;
            }

            // Update x and base_angle
            x = x + a_inter_bs_rest_length;
            base_angle = base_angle + a_inter_bs_twist;
        }
    }
}

void thin_filament::set_unit_status(void)
{
    //! Code sets the status of each unit

    // Variables
    int * bs_indices = new int[a_bs_per_unit];

    // Code
    for (int unit_counter = 0; unit_counter < (a_strands_per_filament * a_regulatory_units_per_strand);
        unit_counter++)
    {
        // Get the bs_indices for the unit
        set_regulatory_unit_indices(unit_counter, bs_indices);

        // Set the unit status from the first entry
        unit_status[unit_counter] = bs_state[bs_indices[0]];
    }

    // Tidy up
    delete bs_indices;
}


void thin_filament::set_regulatory_unit_indices(int unit_ind, int bs_indices[])
{
    //! Fills the bs_indices array with the indices of the binding sites in the unit

    // Variables

    int units_along = unit_ind / a_strands_per_filament;
    int offset = (units_along * a_bs_per_unit * a_strands_per_filament) + 
        (unit_ind % a_strands_per_filament);

    // Code

    for (int i = 0; i < a_bs_per_unit; i++)
    {
        bs_indices[i] = offset + (i * a_strands_per_filament);
    }
}
