/**
* @file		FiberSim_data.cpp
* @brief	Source file for the FiberSim_data class
* @author	ken Campbell
*/

#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "FiberSim_data.h"
#include "FiberSim_options.h"
#include "FiberSim_model.h"

#include "hs_data.h"
#include "kinetic_scheme.h"

#include "gsl_vector.h"
#include "gsl_math.h"

using namespace std::filesystem;

// Constructor
FiberSim_data::FiberSim_data(int set_no_of_time_points,
	FiberSim_options* set_p_fs_options = NULL,
	FiberSim_model* set_p_fs_model = NULL)
{
	// Initialise

	no_of_time_points = set_no_of_time_points;
	p_fs_options = set_p_fs_options;
	p_fs_model = set_p_fs_model;

	// Allocate space for data vectors
	fs_time = gsl_vector_alloc(no_of_time_points);
	fs_m_length = gsl_vector_alloc(no_of_time_points);
	fs_m_force = gsl_vector_alloc(no_of_time_points);

	fs_sc_extension = gsl_vector_alloc(no_of_time_points);
	fs_sc_force = gsl_vector_alloc(no_of_time_points);
	
	// Set to zero
	gsl_vector_set_zero(fs_time);
	gsl_vector_set_zero(fs_m_length);
	gsl_vector_set_zero(fs_m_force);

	gsl_vector_set_all(fs_sc_extension, GSL_NAN);
	gsl_vector_set_all(fs_sc_force, GSL_NAN);

	// Now create the data objects for each half-sarcomere
	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		p_hs_data[hs_counter] = new hs_data(no_of_time_points,
			this,
			p_fs_options,
			p_fs_model);
	}
}

// Destructor
FiberSim_data::~FiberSim_data(void)
{
	// Recover space

	// First the gsl_vectors
	gsl_vector_free(fs_time);
	gsl_vector_free(fs_m_length);
	gsl_vector_free(fs_m_force);
	gsl_vector_free(fs_sc_extension);
	gsl_vector_free(fs_sc_force);

	// Now the half-sarcomere data structures
	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		delete p_hs_data[hs_counter];
	}
}

// Functions

void FiberSim_data::write_data_to_delimited_file(char output_file_string[], char delimiter)
{
	//! Writes data to a delimited file

	// Variables
	FILE* output_file;

	// Code

	if (p_fs_options->log_mode > 0)
	{
		fprintf(p_fs_options->log_file, "Writing data to %s\n", output_file_string);
	}

	// Make sure results directory exists
	path output_file_path(output_file_string);

	if (!(is_directory(output_file_path.parent_path())))
	{
		if (create_directories(output_file_path.parent_path()))
		{
			std::cout << "\nCreating folder: " << output_file_path.string() << "\n";
		}
		else
		{
			std::cout << "\nError: Results folder could not be created: " <<
				output_file_path.parent_path().string() << "\n";
			exit(1);
		}
	}

	// Check file can be opened, abort if not
	errno_t err = fopen_s(&output_file, output_file_string, "w");
	if (err != 0)
	{
		printf("Results file: %s\ncould not be opened\n",
			output_file_string);
		exit(1);
	}

	// Write header
	fprintf_s(output_file, "time%c", delimiter);
	fprintf_s(output_file, "m_length%c", delimiter);
	fprintf_s(output_file, "m_force%c", delimiter);
	fprintf_s(output_file, "sc_extension%c", delimiter);
	fprintf_s(output_file, "sc_force%c", delimiter);

	// Now add in the headers for each half-sarcomere
	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		fprintf_s(output_file, "hs_%i_pCa%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_length%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_command_length%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_slack_length%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_a_length%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_m_length%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_force%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_titin_force%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_viscous_force%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_extracellular_force%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_a_k_on_t_force_factor%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_a_k_off_t_force_factor%c", hs_counter + 1, delimiter);
		fprintf_s(output_file, "hs_%i_a_k_coop_t_force_factor%c", hs_counter + 1, delimiter);

		// Build pops as loops
		for (int j = 0; j < p_fs_model->a_no_of_bs_states; j++)
		{
			fprintf_s(output_file, "hs_%i_a_pop_%i%c", hs_counter + 1, j + 1, delimiter);
		}

		for (int j = 0; j < p_fs_model->p_m_scheme[0]->no_of_states; j++)
		{
			fprintf_s(output_file, "hs_%i_m_pop_%i%c", hs_counter + 1, j + 1, delimiter);
		}

		for (int j = 0; j < p_fs_model->p_c_scheme[0]->no_of_states; j++)
		{
			fprintf_s(output_file, "hs_%i_c_pop_%i", hs_counter + 1, j + 1);

			if ((j == (p_fs_model->p_c_scheme[0]->no_of_states - 1)) &&
				(hs_counter == (p_fs_model->no_of_half_sarcomeres - 1)))
			{
				fprintf_s(output_file, "\n");
			}
			else
			{
				fprintf_s(output_file, "%c", delimiter);
			}
		}
	}

	// Loop through points
	for (int i = 0; i < no_of_time_points; i++)
	{
		fprintf_s(output_file, "%g%c%.3f%c%g%c%g%c%g%c",
			gsl_vector_get(fs_time, i), delimiter,
			gsl_vector_get(fs_m_length, i), delimiter,
			gsl_vector_get(fs_m_force, i), delimiter,
			gsl_vector_get(fs_sc_extension, i), delimiter,
			gsl_vector_get(fs_sc_force, i), delimiter);

		// Now add in the data for each  half-sarcomere
		for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
		{
			fprintf_s(output_file, "%.3g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_pCa, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_length, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_command_length, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_slack_length, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_a_length, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_m_length, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_force, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_titin_force, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_viscous_force, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_extracellular_force, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_a_k_on_t_force_factor, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_a_k_off_t_force_factor, i), delimiter);
			fprintf_s(output_file, "%g%c",
				gsl_vector_get(p_hs_data[hs_counter]->hs_a_k_coop_t_force_factor, i), delimiter);

			// Build pops as loops
			for (int j = 0; j < p_fs_model->a_no_of_bs_states; j++)
			{
				fprintf_s(output_file, "%g%c",
					gsl_matrix_get(p_hs_data[hs_counter]->hs_a_pops, i, j), delimiter);
			}

			for (int j = 0; j < p_fs_model->p_m_scheme[0]->no_of_states; j++)
			{
				fprintf_s(output_file, "%g%c",
					gsl_matrix_get(p_hs_data[hs_counter]->hs_m_pops, i, j), delimiter);
			}

			for (int j = 0; j < p_fs_model->p_c_scheme[0]->no_of_states; j++)
			{
				fprintf_s(output_file, "%g",
					gsl_matrix_get(p_hs_data[hs_counter]->hs_c_pops, i, j));

				if ((j == (p_fs_model->p_c_scheme[0]->no_of_states - 1)) &&
					(hs_counter == (p_fs_model->no_of_half_sarcomeres - 1)))
				{
					fprintf_s(output_file, "\n");
				}
				else
				{
					fprintf_s(output_file, "%c", delimiter);
				}
			}
		}
	}

	// Tidy up
	fclose(output_file);

	std::cout << "Finished writing data to: " << output_file_string << "\n";
}
