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

namespace fs = std::filesystem;

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

	// Get the model version

	sprintf_s(model_version, _MAX_PATH, p_fs_model->version);

	// Now create the half_sarcomeres
	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		p_hs[hs_counter] = new half_sarcomere(p_fs_model, p_fs_options, p_fs_protocol, this, hs_counter);
	}

	// Dump rate_functions to file
	if (strlen(p_fs_options->rate_file_string) > 0)
		p_hs[0]->p_m_scheme[0]->write_rate_functions_to_file(p_fs_options->rate_file_string);

	// Initialise_status_counter
	dump_status_counter = 1;

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
	size_t calculate_x_iterations;			// number of iterations required to solve
											// half-sarcomere force balance

	double sim_mode;						// value from protocol file

	double new_length;						// hs_length if muscle is slack
	double adjustment;						// length change to implose

	// Code

	for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres; hs_counter++)
	{
		// Update the hs_command_length
		p_hs[hs_counter]->hs_command_length = p_hs[hs_counter]->hs_command_length +
			gsl_vector_get(p_fs_protocol->delta_hsl, protocol_index);

		// Branch on control mode
		sim_mode = gsl_vector_get(p_fs_protocol->sim_mode, protocol_index);

		if (sim_mode == -1.0)
		{
			// Check slack length mode for ktr
			p_hs[hs_counter]->hs_slack_length =
				p_hs[hs_counter]->return_hs_length_for_force(0.0, gsl_vector_get(p_fs_protocol->dt, protocol_index));

			// The hs_length cannot be shorter than its slack length
			new_length = GSL_MAX(p_hs[hs_counter]->hs_slack_length,
				p_hs[hs_counter]->hs_command_length);

			// The increment might be smaller than delta_hsl if we are catching up
			// on slack
			adjustment = new_length - p_hs[hs_counter]->hs_length;

			// Make the adjustment
			calculate_x_iterations =
				p_hs[hs_counter]->implement_time_step(
					gsl_vector_get(p_fs_protocol->dt, protocol_index),
					adjustment,
					sim_mode,
					gsl_vector_get(p_fs_protocol->pCa, protocol_index));
		}
		else
		{
			// Over-write slack length
			p_hs[hs_counter]->hs_slack_length = GSL_NAN;

			// Normal operation
			calculate_x_iterations =
				p_hs[hs_counter]->implement_time_step(
					gsl_vector_get(p_fs_protocol->dt, protocol_index),
					gsl_vector_get(p_fs_protocol->delta_hsl, protocol_index),
					gsl_vector_get(p_fs_protocol->sim_mode, protocol_index),
					gsl_vector_get(p_fs_protocol->pCa, protocol_index));

			// Update command length with current hs_length
			// which will have changed in isotonic mode
			p_hs[hs_counter]->hs_command_length = p_hs[hs_counter]->hs_length;
		}

		if ((protocol_index % 100) == 0)
		{
			printf("muscle->hs[%i][%i] ->calculate_x_iterations: %i hsl: %.2f force: %g  a[0]: %g  m[0]: %g c[0]: %g\n",
				hs_counter, protocol_index, (int)calculate_x_iterations,
				p_hs[hs_counter]->hs_length,
				p_hs[hs_counter]->hs_force,
				gsl_vector_get(p_hs[hs_counter]->a_pops, 0),
				gsl_vector_get(p_hs[hs_counter]->m_pops, 0),
				gsl_vector_get(p_hs[hs_counter]->c_pops, 0));
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
	gsl_vector_set(p_fs_data->fs_command_length, protocol_index, p_hs[0]->hs_command_length);
	gsl_vector_set(p_fs_data->fs_slack_length, protocol_index, p_hs[0]->hs_slack_length);
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
	for (int i = 0; i < p_fs_model->p_c_scheme[0]->no_of_states; i++)
	{
		gsl_matrix_set(p_fs_data->fs_c_pops,
			protocol_index, i, gsl_vector_get(p_hs[0]->c_pops, i));
	}

		// Dump the hs_status files if required
		if (protocol_index >= (p_fs_options->start_status_time_step - 1))
	{
		if (protocol_index <= (p_fs_options->stop_status_time_step - 1))
		{
			if (dump_status_counter == 1)
			{
				// Dump status files for each half-sarcomere
				for (int hs_counter = 0; hs_counter < p_fs_model->no_of_half_sarcomeres;
					hs_counter++)
				{
					char hs_status_file_string[_MAX_PATH];
					sprintf_s(hs_status_file_string, _MAX_PATH, "%s/hs_%i_time_step_%i.json",
						p_fs_options->status_folder, hs_counter + 1, protocol_index + 1);
					p_hs[hs_counter]->write_hs_status_to_file(hs_status_file_string);
				}
			}
		}

		// Update dump_status_counter
		dump_status_counter++;

		if (dump_status_counter > p_fs_options->skip_status_time_step)
			dump_status_counter = 1;
	}
}