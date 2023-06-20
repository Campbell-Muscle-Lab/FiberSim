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
#include "kinetic_scheme.h"

#include "gsl_vector.h"

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
	fs_pCa = gsl_vector_alloc(no_of_time_points);
	fs_length = gsl_vector_alloc(no_of_time_points);
	fs_command_length = gsl_vector_alloc(no_of_time_points);
	fs_slack_length = gsl_vector_alloc(no_of_time_points);
	fs_force = gsl_vector_alloc(no_of_time_points);
	fs_titin_force = gsl_vector_alloc(no_of_time_points);
	fs_viscous_force = gsl_vector_alloc(no_of_time_points);
	fs_extracellular_force = gsl_vector_alloc(no_of_time_points);
	fs_a_length = gsl_vector_alloc(no_of_time_points);
	fs_m_length = gsl_vector_alloc(no_of_time_points);	
	
	// Allocate space for data matrices
	fs_a_pops = gsl_matrix_alloc(no_of_time_points, p_fs_model->a_no_of_bs_states);
	fs_m_pops = gsl_matrix_alloc(no_of_time_points, p_fs_model->p_m_scheme[0]->no_of_states);
	fs_c_pops = gsl_matrix_alloc(no_of_time_points, p_fs_model->p_c_scheme[0]->no_of_states);
	
	// Set to zero
	gsl_vector_set_zero(fs_time);
	gsl_vector_set_zero(fs_pCa);
	gsl_vector_set_zero(fs_length);
	gsl_vector_set_zero(fs_command_length);
	gsl_vector_set_zero(fs_slack_length);
	gsl_vector_set_zero(fs_force);
	gsl_vector_set_zero(fs_titin_force);
	gsl_vector_set_zero(fs_viscous_force);
	gsl_vector_set_zero(fs_extracellular_force);

	gsl_vector_set_zero(fs_a_length);
	gsl_vector_set_zero(fs_m_length);

	gsl_matrix_set_zero(fs_a_pops);
	gsl_matrix_set_zero(fs_m_pops);
	gsl_matrix_set_zero(fs_c_pops);
}

// Destructor
FiberSim_data::~FiberSim_data(void)
{
	// Recover space

	// First the gsl_vectors
	gsl_vector_free(fs_time);
	gsl_vector_free(fs_pCa);
	gsl_vector_free(fs_length);
	gsl_vector_free(fs_command_length);
	gsl_vector_free(fs_slack_length);
	gsl_vector_free(fs_force);
	gsl_vector_free(fs_titin_force);
	gsl_vector_free(fs_viscous_force);
	gsl_vector_free(fs_extracellular_force);

	gsl_vector_free(fs_a_length);
	gsl_vector_free(fs_m_length);

	// Then the matrices
	gsl_matrix_free(fs_a_pops);
	gsl_matrix_free(fs_m_pops);
	gsl_matrix_free(fs_c_pops);
}

// Functions

void FiberSim_data::write_progress_file(char progress_file_string[])
{
	// Writes a progress file for the FiberPyGUI

	// Variables
	FILE* progress_output_file;

	// Code
	path;

	errno_t err = fopen_s(&output_file, output_file_string, "w");
	if (err != 0)
	{
		printf("Results file: %s\ncould not be opened\n",
			output_file_string);
		exit(1);
	}
	
}

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
	fprintf_s(output_file, "pCa%c", delimiter);
	fprintf_s(output_file, "hs_length%c", delimiter);
	fprintf_s(output_file, "hs_command_length%c", delimiter);
	fprintf_s(output_file, "hs_slack_length%c", delimiter);
	fprintf_s(output_file, "force%c", delimiter);
	fprintf_s(output_file, "titin_force%c", delimiter);
	fprintf_s(output_file, "viscous_force%c", delimiter);
	fprintf_s(output_file, "extracellular_force%c", delimiter);
	fprintf_s(output_file, "a_fil_length%c", delimiter);
	fprintf_s(output_file, "m_fil_length%c", delimiter);

	// Build pops as loops
	for (int i = 0; i < p_fs_model->a_no_of_bs_states; i++)
	{
		fprintf_s(output_file, "a_pop_%i%c", i, delimiter);
	}
	for (int i = 0; i < p_fs_model->p_m_scheme[0]->no_of_states; i++)
	{
		fprintf_s(output_file, "m_pop_%i%c", i, delimiter);
	}
	for (int i = 0; i < p_fs_model->p_c_scheme[0]->no_of_states; i++)
	{
		fprintf_s(output_file, "c_pop_%i", i);
		if (i == (p_fs_model->p_c_scheme[0]->no_of_states - 1))
			
			fprintf_s(output_file, "\n");
		else
			fprintf_s(output_file, "%c", delimiter);
	}


	// Loop through points
	for (int i = 0; i < no_of_time_points; i++)
	{
		fprintf_s(output_file, "%g%c%.3f%c%g%c%g%c%g%c%g%c%g%c%g%c%g%c%g%c%g%c",
			gsl_vector_get(fs_time, i), delimiter,
			gsl_vector_get(fs_pCa, i), delimiter,
			gsl_vector_get(fs_length, i), delimiter,
			gsl_vector_get(fs_command_length, i), delimiter,
			gsl_vector_get(fs_slack_length, i), delimiter,
			gsl_vector_get(fs_force, i), delimiter,
			gsl_vector_get(fs_titin_force, i), delimiter,
			gsl_vector_get(fs_viscous_force, i), delimiter,
			gsl_vector_get(fs_extracellular_force, i), delimiter,
			gsl_vector_get(fs_a_length, i), delimiter,
			gsl_vector_get(fs_m_length, i), delimiter);

		// Build a pops and m_pops as loops
		for (int j = 0; j < p_fs_model->a_no_of_bs_states; j++)
		{
			fprintf_s(output_file, "%g%c",
				gsl_matrix_get(fs_a_pops, i, j), delimiter);
		}

		for (int j = 0; j < p_fs_model->p_m_scheme[0]->no_of_states; j++)
		{
			fprintf_s(output_file, "%g%c",
				gsl_matrix_get(fs_m_pops, i, j), delimiter);
		}

		for (int j = 0; j < p_fs_model->p_c_scheme[0]->no_of_states; j++)
		{
			fprintf_s(output_file, "%g", gsl_matrix_get(fs_c_pops, i, j));
			if (j == (p_fs_model->p_c_scheme[0]->no_of_states - 1))
				fprintf_s(output_file, "\n");
			else
				fprintf_s(output_file, "%c", delimiter);
		}
	}

	// Tidy up
	fclose(output_file);

	std::cout << "Finished writing data to: " << output_file_string << "\n";
}
