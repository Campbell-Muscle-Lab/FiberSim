/**
 * @file    muscle.cpp
 * @brief   Source file for the muscle class
 * @author  Ken Campbell
  */

#include <iostream>
#include <filesystem>

#include "muscle.h"

#include "FiberSim_model.h"
#include "FiberSim_options.h"
#include "FiberSim_protocol.h"
#include "FiberSim_data.h"

#include "half_sarcomere.h"
#include "kinetic_scheme.h"

#include "rapidjson\document.h"
#include "rapidjson\istreamwrapper.h"

using namespace std::filesystem;

// Constructor
muscle::muscle(char set_model_file_string[], char set_options_file_string[])
{
	// Set file_strings
	sprintf_s(model_file_string, _MAX_PATH, "%s", set_model_file_string);
	sprintf_s(options_file_string, _MAX_PATH, "%s", set_options_file_string);

	// Initialise the muscle id
	muscle_id = 0;

	// Load the options
	p_fs_options = new FiberSim_options(options_file_string);

	// Load the model
	p_fs_model = new FiberSim_model(model_file_string, p_fs_options);

	// Now create the half_sarcomeres
	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		p_hs[hs_counter] = new half_sarcomere(p_fs_model, p_fs_options, p_fs_protocol, this, hs_counter);
	}

	printf("Muscle created half-sarcomeres\n");
}

// Destructor
muscle::~muscle()
{
    // Tidy up

	// Delete the half-sarcomeres
	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		delete p_hs[hs_counter];
	}

	// Delete the FiberSim_model object
	delete p_fs_model;

	// Delete the FiberSim_options object
	delete p_fs_options;
}

// Functions

void muscle::implement_protocol(char set_protocol_file_string[], char set_results_file_string[])
{
	//! Code runs a muscle through a protocol

	// Variables

	// Code

		// Set the main results file
	sprintf_s(results_file_string, _MAX_PATH, "%s", set_results_file_string);

	// Update the protocol file_string
	sprintf_s(protocol_file_string, _MAX_PATH, "%s", set_protocol_file_string);


	// Load the protocol
	p_fs_protocol = new FiberSim_protocol(protocol_file_string);

	// Create a FiberSim_data object for the muscle
	p_fs_data = new FiberSim_data(p_fs_protocol->no_of_time_points, p_fs_options, p_fs_model);

	// And also for the half-sarcomeres
	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		p_hs[hs_counter]->p_fs_data = new FiberSim_data(p_fs_protocol->no_of_time_points,
										p_fs_options, p_fs_model);
	}

	// Implement the protocol
	for (int i = 0; i < p_fs_protocol->no_of_time_points; i++)
	{
		implement_time_step(i);
	}

	// Output main results file
	if (strlen(results_file_string) > 0)
	{
		printf("Muscle[%i]: Attempting to write results to: %s\n", muscle_id, results_file_string);
		p_fs_data->write_data_to_delimited_file(results_file_string);
	}

	// Delete the FiberSim_data object for the half-sarcomeres
	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		delete p_hs[hs_counter]->p_fs_data;
	}

	// And also for the muscle
	delete p_fs_data;

	// Delete the FiberSim_protocol object
	delete p_fs_protocol;
}

void muscle::implement_time_step(int protocol_index)
{
	//! Code implements a time-step

	// Variables
	size_t calculate_x_iterations;

	// Code

	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		calculate_x_iterations = 
			p_hs[hs_counter]->implement_time_step(
			gsl_vector_get(p_fs_protocol->dt, protocol_index),
			gsl_vector_get(p_fs_protocol->delta_hsl, protocol_index),
			gsl_vector_get(p_fs_protocol->sim_mode, protocol_index),
			gsl_vector_get(p_fs_protocol->pCa, protocol_index));

		if ((protocol_index % 100) == 0)
		{
			printf("muscle->hs[%i]->time_step: %i ->calculate_x_iterations: %i  force: %g\n",
				hs_counter, protocol_index, (int)calculate_x_iterations, p_hs[hs_counter]->hs_force);
		}
	}

	if (p_fs_model->no_of_half_sarcomeres > 1)
	{
		printf("Muscle::implement_time_step does not yet work with >1 half-sarcomere\n");
		exit(1);
	}

	// Update FiberSim_data
	gsl_vector_set(p_fs_data->fs_time, protocol_index, p_hs[0]->time_s);
	gsl_vector_set(p_fs_data->fs_length, protocol_index, p_hs[0]->hs_length);
	gsl_vector_set(p_fs_data->fs_force, protocol_index, p_hs[0]->hs_force);
	gsl_vector_set(p_fs_data->fs_titin_force, protocol_index, p_hs[0]->hs_titin_force);
	gsl_vector_set(p_fs_data->fs_extracellular_force, protocol_index,
		p_hs[0]->hs_extracellular_force);
	gsl_vector_set(p_fs_data->fs_pCa, protocol_index, p_hs[0]->pCa);

	gsl_vector_set(p_fs_data->fs_a_length, protocol_index, p_hs[0]->a_mean_fil_length);
	gsl_vector_set(p_fs_data->fs_m_length, protocol_index, p_hs[0]->m_mean_fil_length);

	// Update pops
	for (int i = 0; i < p_fs_model->a_no_of_bs_states; i++)
	{
		gsl_matrix_set(p_fs_data->fs_a_pops,
			protocol_index, i, gsl_vector_get(p_hs[0]->a_pops, i));
	}
	for (int i = 0; i < p_fs_model->p_m_scheme[0]->no_of_states; i++)
	{
		gsl_matrix_set(p_fs_data->fs_m_pops,
				protocol_index, i, gsl_vector_get(p_hs[0]->m_pops, i));
	}
	for (int i = 0; i < p_fs_model->p_c_scheme->no_of_states; i++)
	{
		gsl_matrix_set(p_fs_data->fs_c_pops,
			protocol_index, i, gsl_vector_get(p_hs[0]->c_pops, i));
	}

}