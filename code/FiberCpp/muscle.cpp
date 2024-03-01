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
#include "hs_data.h"

#include "half_sarcomere.h"
#include "kinetic_scheme.h"
#include "series_component.h"

#include "gsl_vector.h"
#include "gsl_multiroots.h"
#include "gsl_math.h"

#include <thread>

#include "rapidjson\document.h"
#include "rapidjson\istreamwrapper.h"

#include "BS_thread_pool.hpp"

#include "chrono"

namespace fs = std::filesystem;

// Structure used for root-finding for myofibril in length-or force control mode
struct m_control_params
{
	double time_step;
	muscle* p_m;
	double target_force;
};

// This is a function used by the root finding algorithm that handles the recasting of pointers
int wrapper_length_control_myofibril_with_series_compliance(const gsl_vector* x, void* params, gsl_vector* f);

struct force_control_params
{
	double target_force;
	double time_step;
	half_sarcomere* p_hs;
	double delta_hsl;
};

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

	// Check if we are doing an afterloaded contraction
	if (gsl_isnan(p_fs_options->afterload_load) == 0)
	{
		afterload_flag = 1;
		afterload_mode = 0;
		afterload_min_hs_length = GSL_POSINF;
	}
	else
	{
		afterload_flag = -1;
	}

	// Initialize the muscle length, and build it from the half-sarcomeres
	m_length = 0.0;

	// Load an initial copy of the model - we need this to deduce the number
	// of half-sarcomeres
	p_fs_model[0] = new FiberSim_model(model_file_string, p_fs_options);

	// Now create the half_sarcomeres
	for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
	{
		// Load another model if required
		if (hs_counter > 0)
		{
			p_fs_model[hs_counter] = new FiberSim_model(model_file_string, p_fs_options);
		}

		// Create the half-sarcomere
		p_hs[hs_counter] = new half_sarcomere(p_fs_model[hs_counter],
												p_fs_options, p_fs_protocol, this, hs_counter);

		m_length = m_length + p_hs[hs_counter]->hs_length;
	}

	// Get the model version
	sprintf_s(model_version, _MAX_PATH, p_fs_model[0]->version);

	// Set the muscle force to the force in the last half-sarcomere
	m_force = p_hs[0]->hs_force;

	// Make a series component if you need one
	if (!gsl_isnan(p_fs_model[0]->sc_k_stiff))
	{
		// Correct the series stiffness for different numbers of half-sarcomeres
		// This ensures that in myofibrils with different numbers of half-sarcomeres,
		// each half-sarcomere will shorten by the same amount
		p_fs_model[0]->sc_k_stiff = p_fs_model[0]->sc_k_stiff / (double)p_fs_model[0]->no_of_half_sarcomeres;

		p_sc = new series_component(p_fs_model[0], p_fs_options, this);
	}
	else
		p_sc = NULL;

	// Dump rate_functions to file
	if (strlen(p_fs_options->rate_file_string) > 0)
		write_rates_file();

	// Initialise_status_counter
	dump_status_counter = 1;

	printf("Muscle created %i half-sarcomeres\n", p_fs_model[0]->no_of_half_sarcomeres);

	// Initialized the clock
	last_time = std::chrono::high_resolution_clock::now();

	// Activate a worker pool if required
	if (p_fs_options->myofibril_multithreading > 0)
	{
		p_fs_options->no_of_worker_threads = p_fs_model[0]->no_of_half_sarcomeres;

		p_thread_pool = new BS::thread_pool(p_fs_options->no_of_worker_threads);
		printf("FiberSim has launched: %i threads\n", p_thread_pool->get_thread_count());
	}
	else
	{
		p_thread_pool = NULL;
		printf("FiberSim is running on a single thread\n");
	}
}

// Destructor
muscle::~muscle()
{
    // Tidy up

	// Delete the instances of the model and the half-sarcomeres
	for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
	{
		delete p_fs_model[hs_counter];
		delete p_hs[hs_counter];
	}

	// Delete the series component if it was created
	if (p_sc != NULL)
		delete p_sc;

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
	p_fs_data = new FiberSim_data(p_fs_protocol->no_of_time_points,
										p_fs_options, p_fs_model[0]);

	// And also for the half-sarcomeres
	for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
	{
		p_hs[hs_counter]->p_fs_data = new FiberSim_data(p_fs_protocol->no_of_time_points,
										p_fs_options, p_fs_model[hs_counter]);
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
	for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
	{
		delete p_hs[hs_counter]->p_fs_data;
	}

	// And also for the muscle
	delete p_fs_data;

	// Delete the FiberSim_protocol object
	delete p_fs_protocol;

	// Delete the thread_pool
	delete p_thread_pool;
}

void muscle::afterload_time_step(int protocol_index)
{
	//! Code implements a time-step with afterload control

	// Variables
	size_t lattice_iterations;					// number of iterations required to solve
	// x positions for lattice. If the code is
	// simulating a myofibril, this is the largest
	// number of iterations required

	double sim_mode;						// value from protocol file

	double time_step_s;
	double pCa;

	double new_length;						// hs_length if muscle is slack
	double adjustment;						// length change to implose

	int hs_counter;

	if ((p_fs_model[0]->no_of_half_sarcomeres == 1) && (p_sc == NULL))
	{

	}
	else
	{
		// We have a myofibril

		printf("afterload_flag: %i\n", afterload_flag);

		if (afterload_mode <= 0)
		{
			// Length control
			lattice_iterations = length_control_myofibril_with_series_compliance(protocol_index);

			if ((afterload_mode == 0) && (m_force >= p_fs_options->afterload_load))
				afterload_mode = 1;
		}

		if (afterload_mode == 1)
		{
			gsl_vector_set(p_fs_protocol->sim_mode, protocol_index, p_fs_options->afterload_load);

			// Force control
			lattice_iterations = force_control_myofibril_with_series_compliance(protocol_index);

			afterload_min_hs_length = GSL_MIN(afterload_min_hs_length,
				m_length / p_fs_model[0]->no_of_half_sarcomeres);

			if ((m_length - afterload_min_hs_length) > p_fs_options->afterload_break_delta_hs_length)
				afterload_mode = -1;
		}
	}

	printf("afterload_min_hs_length: %g\n", afterload_min_hs_length);

}

void muscle::implement_time_step(int protocol_index)
{
	//! Code implements a time-step

	// Variables
	size_t lattice_iterations;					// number of iterations required to solve
											// x positions for lattice. If the code is
											// simulating a myofibril, this is the largest
											// number of iterations required

	double sim_mode;						// value from protocol file

	double time_step_s;
	double pCa;

	double new_length;						// hs_length if muscle is slack
	double adjustment;						// length change to implose

	int hs_counter;

	// Code
	if (afterload_flag == 1)
	{
		afterload_time_step(protocol_index);
	}
	else
	{
		if ((p_fs_model[0]->no_of_half_sarcomeres == 1) && (p_sc == NULL))
		{
			// Simplest case of 1 half-sarcomere and no compliance
			hs_counter = 0;

			// Update the hs_command_length
			p_hs[hs_counter]->hs_command_length = p_hs[hs_counter]->hs_command_length +
				gsl_vector_get(p_fs_protocol->delta_hsl, protocol_index);

			// Update the kinetics
			time_step_s = gsl_vector_get(p_fs_protocol->dt, protocol_index);
			pCa = gsl_vector_get(p_fs_protocol->pCa, protocol_index);

			p_hs[hs_counter]->time_s = p_hs[hs_counter]->time_s + time_step_s;
			p_hs[hs_counter]->sarcomere_kinetics(time_step_s, pCa);

			if (p_hs[hs_counter]->hs_length > p_fs_options->min_hs_length)
			{
				// We are safe to continue

				// Branch on control mode
				sim_mode = gsl_vector_get(p_fs_protocol->sim_mode, protocol_index);

				// Semi-clever check for comparing sim_mode to -1.0
				if (gsl_fcmp(sim_mode, -1.0, 1e-3) == 0)
				{
					// Check slack length mode for ktr
					p_hs[hs_counter]->hs_slack_length =
						p_hs[hs_counter]->return_hs_length_for_force(0.0, gsl_vector_get(p_fs_protocol->dt, protocol_index));

					// The hs_length cannot be shorter than its slack length
					new_length = GSL_MAX(p_hs[hs_counter]->hs_slack_length,
						p_hs[hs_counter]->hs_command_length);

					adjustment = new_length - p_hs[hs_counter]->hs_length;
				}
				else
				{
					// Over-write slack length
					p_hs[hs_counter]->hs_slack_length = GSL_NAN;

					// Are we in force control
					if (sim_mode >= 0.0)
					{
						// Force control
						new_length = p_hs[hs_counter]->return_hs_length_for_force(sim_mode, time_step_s);

						adjustment = new_length - p_hs[hs_counter]->hs_length;

						// Update the command length which changes depending on force control
						p_hs[hs_counter]->hs_command_length = p_hs[hs_counter]->hs_length;
					}
					else
					{
						// Length control
						adjustment = gsl_vector_get(p_fs_protocol->delta_hsl, protocol_index);
					}
				}

				// Apply the adjustment
				lattice_iterations = p_hs[hs_counter]->update_lattice(time_step_s, adjustment);

				// Update the muscle length
				m_length = p_hs[hs_counter]->hs_command_length;

				// Update the muscle force
				m_force = p_hs[hs_counter]->hs_force;
			}
		}
		else
		{
			// We have a myofibril

			// Branch on control mode
			sim_mode = gsl_vector_get(p_fs_protocol->sim_mode, protocol_index);

			if (sim_mode >= 0.0)
			{
				// Force control
				lattice_iterations = force_control_myofibril_with_series_compliance(protocol_index);
			}
			else
			{
				// Length control
				lattice_iterations = length_control_myofibril_with_series_compliance(protocol_index);
			}
		}
	}

	if ((protocol_index % 100) == 0)
	{
		hs_counter = 0;

		printf("muscle->hs[%i][%i] ->lattice_iterations: %zi hsl: %.2f force: %g  a[0]: %g  m[0]: %g c[0]: %g\n",
			hs_counter, protocol_index, lattice_iterations,
			p_hs[hs_counter]->hs_length,
				p_hs[hs_counter]->hs_force,
				gsl_vector_get(p_hs[hs_counter]->a_pops, 0),
				gsl_vector_get(p_hs[hs_counter]->m_pops, 0),
				gsl_vector_get(p_hs[hs_counter]->c_pops, 0));

		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> fp_ms = t2 - last_time;

		std::cout << "Section took: " << fp_ms << "\n";

		last_time = std::chrono::high_resolution_clock::now();
	}

	// Update FiberSim_data

	// First the shared variables
	gsl_vector_set(p_fs_data->fs_time, protocol_index, p_hs[0]->time_s);
	gsl_vector_set(p_fs_data->fs_m_length, protocol_index, m_length);
	gsl_vector_set(p_fs_data->fs_m_force, protocol_index, m_force);

	if (p_sc != NULL)
	{
		gsl_vector_set(p_fs_data->fs_sc_extension, protocol_index, p_sc->sc_extension);
		gsl_vector_set(p_fs_data->fs_sc_force, protocol_index, p_sc->sc_force);
	}

	// Now the half-sarcomeres
	for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
	{
		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_pCa,
			protocol_index,
			p_hs[hs_counter]->pCa);

		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_length,
			protocol_index,
			p_hs[hs_counter]->hs_length);

		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_command_length,
			protocol_index,
			p_hs[hs_counter]->hs_command_length);

		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_slack_length,
			protocol_index,
			p_hs[hs_counter]->hs_slack_length);

		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_force,
			protocol_index,
			p_hs[hs_counter]->hs_force);

		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_titin_force,
			protocol_index,
			p_hs[hs_counter]->hs_titin_force);

		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_viscous_force,
			protocol_index,
			p_hs[hs_counter]->hs_viscous_force);

		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_extracellular_force,
			protocol_index,
			p_hs[hs_counter]->hs_extracellular_force);

		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_a_length,
			protocol_index,
			p_hs[hs_counter]->a_mean_fil_length);

		gsl_vector_set(p_fs_data->p_hs_data[hs_counter]->hs_m_length,
			protocol_index,
			p_hs[hs_counter]->m_mean_fil_length);

		// Update pops
		for (int i = 0; i < p_fs_model[hs_counter]->a_no_of_bs_states; i++)
		{
			gsl_matrix_set(p_fs_data->p_hs_data[hs_counter]->hs_a_pops,
				protocol_index, i, gsl_vector_get(p_hs[hs_counter]->a_pops, i));
		}
		for (int i = 0; i < p_fs_model[hs_counter]->p_m_scheme[0]->no_of_states; i++)
		{
			gsl_matrix_set(p_fs_data->p_hs_data[hs_counter]->hs_m_pops,
				protocol_index, i, gsl_vector_get(p_hs[hs_counter]->m_pops, i));
		}
		for (int i = 0; i < p_fs_model[hs_counter]->p_c_scheme[0]->no_of_states; i++)
		{
			gsl_matrix_set(p_fs_data->p_hs_data[hs_counter]->hs_c_pops,
				protocol_index, i, gsl_vector_get(p_hs[hs_counter]->c_pops, i));
		}
	}

	// Dump the hs_status files if required
	if (protocol_index >= (p_fs_options->start_status_time_step - 1))
	{
		if (protocol_index <= (p_fs_options->stop_status_time_step - 1))
		{
			if (dump_status_counter == 1)
			{
				// Dump status files for each half-sarcomere
				for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres;
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

void thread_sarcomere_kinetics(half_sarcomere* p_hs, double time_step, double pCa_value)
{
	// This is a non-class function that is called during multi-threading
	// This is an easy way to handle the pointers to half-sarcomeres

	p_hs->sarcomere_kinetics(time_step, pCa_value);
}

size_t muscle::length_control_myofibril_with_series_compliance(int protocol_index)
{
	//! Tries to impose length control on a system with a series compliance
	//! and/or 1 or more half-sarcomeres

	// Try to find a vector x such that the forces in each half-sarcomere
	// and the force in the series elastic element (which is the length
	// of the muscle - the combined length of the half-sarcomeres)
	// are all equal. x is filled with an initial guess.

	// Variables

	// Stuff to do with the root finding
	int status;
	int myofibril_iterations;
	size_t x_length;
	gsl_vector* x;

	// Other stuff
	double time_step_s;
	double holder_hs_length;

	size_t max_lattice_iterations;

	// Code

	// First adjust the muscle length
	m_length = m_length + ((double)p_fs_model[0]->no_of_half_sarcomeres *
		gsl_vector_get(p_fs_protocol->delta_hsl, protocol_index));


	// Now run the kinetics on each half-sarcomere
	time_step_s = gsl_vector_get(p_fs_protocol->dt, protocol_index);

	const auto t1 = std::chrono::high_resolution_clock::now();

	if (p_thread_pool != NULL)
	{
		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			p_hs[hs_counter]->time_s = p_hs[hs_counter]->time_s + time_step_s;

			double pCa_temp = gsl_vector_get(p_fs_protocol->pCa, protocol_index);

			p_thread_pool->push_task(&half_sarcomere::sarcomere_kinetics, p_hs[hs_counter],
				time_step_s, gsl_vector_get(p_fs_protocol->pCa, protocol_index));

		}
		p_thread_pool->wait_for_tasks();
	}
	else
	{
		// Serial processing
		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			p_hs[hs_counter]->time_s = p_hs[hs_counter]->time_s + time_step_s;

			p_hs[hs_counter]->sarcomere_kinetics(time_step_s, gsl_vector_get(p_fs_protocol->pCa, protocol_index));
		}
	}

	auto t_kinetics = std::chrono::high_resolution_clock::now();
	
	const std::chrono::duration<double, std::milli> dur_kinetics = t_kinetics - t1;

	//std::cout << "Kinetics took " << dur_kinetics << "\n";

	// Now deduce the length of x
	x_length = p_fs_model[0]->no_of_half_sarcomeres + 1;

	// Allocate the vector
	x = gsl_vector_alloc(x_length);

	// The x-vector has all the half-sarcomere lengths followed by the force in the first_half_sarcomere
	for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
	{
		gsl_vector_set(x, hs_counter, p_hs[hs_counter]->hs_length);
	}
	gsl_vector_set(x, x_length - 1, p_sc->sc_force);

	// Do the root finding
	const gsl_multiroot_fsolver_type* T;
	gsl_multiroot_fsolver* s;
	const size_t calculation_size = x_length;

	m_control_params* par = new m_control_params;
	par->p_m = this;
	par->time_step = gsl_vector_get(p_fs_protocol->dt, protocol_index);

	gsl_multiroot_function f = { &wrapper_length_control_myofibril_with_series_compliance, calculation_size, par };

	//T = gsl_multiroot_fsolver_hybrid;
	T = gsl_multiroot_fsolver_hybrid;
	s = gsl_multiroot_fsolver_alloc(T, calculation_size);
	gsl_multiroot_fsolver_set(s, &f, x);

	myofibril_iterations = 0;

	do
	{
		gsl_vector* y = gsl_vector_alloc(x_length);
		
		status = gsl_multiroot_fsolver_iterate(s);

		myofibril_iterations++;

		if (status)
		{
			printf("Myofibril multiroot solver break - Status: %i\t", status);

			if (status == GSL_EBADFUNC)
			{
				printf("Bad function value\n");
			}

			if (status == GSL_ENOPROG)
			{
				printf("Not making progress\n");
			}

			if (status == GSL_ENOPROGJ)
			{
				printf("Jacobian evaluations are not helping\n");
			}
		}

		status = gsl_multiroot_test_delta(s->dx, s->x, p_fs_options->myofibril_force_tolerance, 0);

		gsl_vector_free(y);
	}
	while ((status == GSL_CONTINUE) && (myofibril_iterations < p_fs_options->myofibril_max_iterations));

	// At this point, the s->x vector contains the lengths of the n half-sarcomeres
	// followed by the force in the series element

	// Implement the change
	max_lattice_iterations = 0;

	holder_hs_length = 0.0;
	
	for (int hs_counter = 0 ; hs_counter < p_fs_model[0]->no_of_half_sarcomeres ; hs_counter++)
	{
		double new_hs_length = gsl_vector_get(s->x, hs_counter);
		double delta_hsl = new_hs_length - p_hs[hs_counter]->hs_length;

		size_t lattice_iterations;

		lattice_iterations = p_hs[hs_counter]->update_lattice(time_step_s, delta_hsl);
		max_lattice_iterations = GSL_MAX(lattice_iterations, max_lattice_iterations);

		holder_hs_length = holder_hs_length + new_hs_length;
	}

	p_sc->sc_extension = (m_length - holder_hs_length);
	p_sc->sc_force = p_sc->return_series_force(p_sc->sc_extension);

	// Update muscle force
	m_force = p_sc->sc_force;

	// Tidy up
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);

	delete par;

	// Return the max number of lattice iterations
	return max_lattice_iterations;
}

int wrapper_length_control_myofibril_with_series_compliance(const gsl_vector* x, void* p, gsl_vector* f)
{
	//! This is a wrapper around muscle::check_residuals_for_myofibril_length_control()
	//! that handles the re-casting of pointers

	// Variables
	int f_return_value;

	struct m_control_params* params =
		(struct m_control_params*)p;

	// Code

	muscle* p_m = params->p_m;

	f_return_value = (int)p_m->worker_length_control_myofibril_with_series_compliance(x, params, f);

	return f_return_value;
}

void thread_test_force_wrapper(half_sarcomere* p_hs, force_control_params* fp, gsl_vector* f, int i)
{
	//! This is a helper function that can be launched as a seperate thread to calculate
	//! the residual for force balance

	double y;
	y= p_hs->test_force_wrapper(fp->delta_hsl, fp);

	gsl_vector_set(f, i, y);
}

size_t muscle::worker_length_control_myofibril_with_series_compliance(
	const gsl_vector* x, void* p, gsl_vector* f)
{
	//! This code calculates the f_vector that is minimized towards 0 to
	//! ensure that the force in the myofibril and its length are constrained

	// Variables
	struct m_control_params* params =
		(struct m_control_params*)p;

	double delta_hsl;
	double cum_hs_length;
	double test_se_length;
	double force_diff;

	// Code

	// The f-vector is the difference between the force in each half-sarcomere and
	// the force in the series component
	// The x vector has a series of lengths followed by a muscle force
	// Calculate the force in each half-sarcomere and compare it to the force
	// Store up the half-sarcomere lengths as you go, and use that to calculate the length
	// of the series component

	// We need the force-control params for the calculation

	// Decide whether we are doing this in parallel with threads
	if (p_thread_pool != NULL)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		// We need an array of force-control points, one for each thread
		force_control_params* f_p[MAX_NO_OF_HALF_SARCOMERES];

		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			f_p[hs_counter] = new force_control_params();
			f_p[hs_counter]->target_force = gsl_vector_get(x, x->size - 1);
			f_p[hs_counter]->time_step = params->time_step;
			f_p[hs_counter]->p_hs = p_hs[hs_counter];
		}

		//auto t2 = std::chrono::high_resolution_clock::now();
		//const std::chrono::duration<double, std::milli> make_fp = t2 - t1;
		//printf("make_f_p: %f\n", make_fp);

		// Now cycle through the half-sarcomeres, checking the residual for each
		// half sarcomere in a separate thread

		cum_hs_length = 0;

		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			cum_hs_length = cum_hs_length + gsl_vector_get(x, hs_counter);

			delta_hsl = gsl_vector_get(x, hs_counter) - p_hs[hs_counter]->hs_length;

			// Constrain delta_hsl to a plausible range
			if (delta_hsl > p_fs_options->myofibril_max_delta_hs_length)
				delta_hsl = p_fs_options->myofibril_max_delta_hs_length;
			if (delta_hsl < -p_fs_options->myofibril_max_delta_hs_length)
				delta_hsl = -p_fs_options->myofibril_max_delta_hs_length;

			f_p[hs_counter]->delta_hsl = delta_hsl;

			p_thread_pool->push_task(thread_test_force_wrapper, p_hs[hs_counter], f_p[hs_counter], f, hs_counter);
		}
		p_thread_pool->wait_for_tasks();

		//t2 = std::chrono::high_resolution_clock::now();
		//const std::chrono::duration<double, std::milli> adjust_fp = t2 - t1;
		//printf("adjust_f_p: %f\n", adjust_fp);

		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
			delete f_p[hs_counter];

		const std::chrono::high_resolution_clock::time_point t2 =
			std::chrono::high_resolution_clock::now();

		const std::chrono::duration<double, std::milli> total = t2 - t_start;
//		printf("total_parallel: %f\n", total);
	}
	else
	{
		// Serial operation

		// Set the force control parameters
		force_control_params* fp = new force_control_params;
		fp->target_force = gsl_vector_get(x, x->size - 1);
		fp->time_step = params->time_step;

		cum_hs_length = 0;

		auto t1 = std::chrono::high_resolution_clock::now();

		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			// Update the pointer
			fp->p_hs = p_hs[hs_counter];

			// Update the length
			cum_hs_length = cum_hs_length + gsl_vector_get(x, hs_counter);

			delta_hsl = gsl_vector_get(x, hs_counter) - p_hs[hs_counter]->hs_length;

			// Constrain delta_hsl to a plausible range
			if (delta_hsl > p_fs_options->myofibril_max_delta_hs_length)
				delta_hsl = p_fs_options->myofibril_max_delta_hs_length;
			if (delta_hsl < -p_fs_options->myofibril_max_delta_hs_length)
				delta_hsl = -p_fs_options->myofibril_max_delta_hs_length;

			force_diff = p_hs[hs_counter]->test_force_wrapper(delta_hsl, fp);

			gsl_vector_set(f, hs_counter, force_diff);
		}

		delete fp;

		auto t2 = std::chrono::high_resolution_clock::now();
		const std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
	}

	// Now deduce the series elastic force
	test_se_length = m_length - cum_hs_length;
	force_diff = p_sc->return_series_force(test_se_length) - gsl_vector_get(x, x->size - 1);

	gsl_vector_set(f, f->size - 1, force_diff);

	return GSL_SUCCESS;
}

size_t muscle::force_control_myofibril_with_series_compliance(int protocol_index)
{
	//! Impose force control on a system with a series compliance
	//! and/or 1 or more half-sarcomeres

	// Variables

	// Other stuff
	double time_step_s;
	double target_force;

	size_t lattice_iterations;
	size_t max_lattice_iterations = 0;

	// Run the kinetics on each half-sarcomere
	time_step_s = gsl_vector_get(p_fs_protocol->dt, protocol_index);

	const auto t1 = std::chrono::high_resolution_clock::now();

	if (p_thread_pool != NULL)
	{
		// We are using multi-threading
		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			p_hs[hs_counter]->time_s = p_hs[hs_counter]->time_s + time_step_s;

			double pCa_temp = gsl_vector_get(p_fs_protocol->pCa, protocol_index);

			p_thread_pool->push_task(&half_sarcomere::sarcomere_kinetics, p_hs[hs_counter],
				time_step_s, gsl_vector_get(p_fs_protocol->pCa, protocol_index));
		}
		p_thread_pool->wait_for_tasks();
	}
	else
	{
		// Serial processing
		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			p_hs[hs_counter]->time_s = p_hs[hs_counter]->time_s + time_step_s;

			p_hs[hs_counter]->sarcomere_kinetics(time_step_s, gsl_vector_get(p_fs_protocol->pCa, protocol_index));
		}
	}

	auto t_kinetics = std::chrono::high_resolution_clock::now();
	
	const std::chrono::duration<double, std::milli> dur_kinetics = t_kinetics - t1;

	//std::cout << "Kinetics took: " << dur_kinetics << "\n";

	// Now implement the force control on each half-sarcomere

	// Get the target_force
	target_force = gsl_vector_get(p_fs_protocol->sim_mode, protocol_index);

	if (p_thread_pool != NULL)
	{
		// We are using multi-threading
		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			p_thread_pool->push_task(&half_sarcomere::update_lattice_for_force, p_hs[hs_counter],
				time_step_s, target_force);
		}
		p_thread_pool->wait_for_tasks();

		max_lattice_iterations = 0;
		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			if ((size_t)(p_hs[hs_counter]->thread_return_value) > max_lattice_iterations)
				max_lattice_iterations = (size_t)(p_hs[hs_counter]->thread_return_value);
		}
	}
	else
	{
		// Serial processing
		max_lattice_iterations = 0;
		
		for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
		{
			lattice_iterations = p_hs[hs_counter]->update_lattice_for_force(time_step_s, target_force);

			if (lattice_iterations > max_lattice_iterations)
				max_lattice_iterations = lattice_iterations;
		}
	}

	// Update the series component
	p_sc->sc_extension = p_sc->return_series_extension(target_force);
	p_sc->sc_force = target_force;

	// Update the muscle length
	m_length = p_sc->sc_extension;
	for (int hs_counter = 0; hs_counter < p_fs_model[0]->no_of_half_sarcomeres; hs_counter++)
	{
		m_length = m_length + p_hs[hs_counter]->hs_length;
	}

	// Update muscle force
	m_force = p_sc->sc_force;
	
	// Return the max number of lattice iterations
	return max_lattice_iterations;
}

void muscle::write_rates_file()
{
	//! Function writes the m and c rate functions to file in JSON format

	// Variables
	int isotype_counter;					// isotype counter

	char file_write_mode[_MAX_PATH];		// mode for opening file

	char JSON_append_string[_MAX_PATH];		// written after scheme to keep JSON
											// structure, should be , if other entries follow
											// otherwise ""

	FILE* output_file;						// pointer for output file

	// Make sure directory exists
	path output_file_path(p_fs_options->rate_file_string);

	if (!(is_directory(output_file_path.parent_path())))
	{
		if (create_directories(output_file_path.parent_path()))
		{
			std::cout << "\nCreating folder: " << output_file_path.string() << "\n";
		}
		else
		{
			std::cout << "\nError: folder for rates file could not be created: " <<
				output_file_path.parent_path().string() << "\n";
			exit(1);
		}
	}

	// Check file can be opened in write mode, abort if not
	errno_t err = fopen_s(&output_file, p_fs_options->rate_file_string, "w");

	if (err != 0)
	{
		printf("muscle::write_rates_file(): %s\ncould not be opened\n",
			p_fs_options->rate_file_string);
		exit(1);
	}

	// Start JSON structure
	fprintf_s(output_file, "{\n\t\"FiberSim_rates\":\n\t{\n");
	fprintf_s(output_file, "\t\t\"myosin\":\n");
	fprintf_s(output_file, "\t\t[\n");
	fclose(output_file);

	// Set the file write mode
	sprintf_s(file_write_mode, _MAX_PATH, "a");

	// Now cycle through the m isotypes
	for (isotype_counter = 0; isotype_counter < p_fs_model[0]->m_no_of_isotypes; isotype_counter++)
	{
		// Set the append string
		if (isotype_counter < (p_fs_model[0]->m_no_of_isotypes - 1))
		{
			sprintf_s(JSON_append_string, _MAX_PATH, ",");
		}
		else
		{
			sprintf_s(JSON_append_string, _MAX_PATH, "");
		}

		p_hs[0]->p_m_scheme[isotype_counter]->write_rate_functions_to_file(
			p_fs_options->rate_file_string, file_write_mode,
			JSON_append_string, p_hs[0]);
	}

	// Re-open the file, close the myosin array, and prep for the c array
	fopen_s(&output_file, p_fs_options->rate_file_string, "a");
	fprintf_s(output_file, "\t\t],\n");
	fprintf_s(output_file, "\t\t\"mybpc\":\n\t\t[\n");
	fclose(output_file);

	// Now through the c_isotypes
	for (isotype_counter = 0; isotype_counter < p_fs_model[0]->c_no_of_isotypes; isotype_counter++)
	{
		// Set the append string
		if (isotype_counter < (p_fs_model[0]->c_no_of_isotypes - 1))
		{
			sprintf_s(JSON_append_string, _MAX_PATH, ",");
		}
		else
		{
			sprintf_s(JSON_append_string, _MAX_PATH, "");
		}

		p_hs[0]->p_c_scheme[isotype_counter]->write_rate_functions_to_file(
			p_fs_options->rate_file_string, file_write_mode,
			JSON_append_string);
	}

	// Now tidy up the rates file
	// Re-open the file, close the cc array
	fopen_s(&output_file, p_fs_options->rate_file_string, "a");
	fprintf_s(output_file, "\t\t]\n\t}\n}\n");
	fclose(output_file);
}



