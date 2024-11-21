/**
* @file		FiberSim_data.cpp
* @brief	Source file for the FiberSim_data class
* @author	ken Campbell
*/

#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "FiberSim_data.h"
#include "hs_data.h"
#include "FiberSim_options.h"
#include "FiberSim_model.h"
#include "kinetic_scheme.h"

#include "gsl_vector.h"

using namespace std::filesystem;

// Constructor
hs_data::hs_data(int set_no_of_time_points,
	FiberSim_data* set_p_fs_data = NULL,
	FiberSim_options* set_p_fs_options = NULL,
	FiberSim_model* set_p_fs_model = NULL)
{
	// Initialise

	p_fs_data = set_p_fs_data;

	no_of_time_points = set_no_of_time_points;
	p_fs_options = set_p_fs_options;
	p_fs_model = set_p_fs_model;

	// Allocate space for data vectors
	hs_time = gsl_vector_alloc(no_of_time_points);
	hs_pCa = gsl_vector_alloc(no_of_time_points);
	hs_length = gsl_vector_alloc(no_of_time_points);
	hs_command_length = gsl_vector_alloc(no_of_time_points);
	hs_slack_length = gsl_vector_alloc(no_of_time_points);
	hs_force = gsl_vector_alloc(no_of_time_points);
	hs_titin_force = gsl_vector_alloc(no_of_time_points);
	hs_viscous_force = gsl_vector_alloc(no_of_time_points);
	hs_extracellular_force = gsl_vector_alloc(no_of_time_points);
	hs_a_length = gsl_vector_alloc(no_of_time_points);
	hs_m_length = gsl_vector_alloc(no_of_time_points);

	hs_inter_hs_titin_force_effect = gsl_vector_alloc(no_of_time_points);
	
	// Allocate space for data matrices
	hs_a_pops = gsl_matrix_alloc(no_of_time_points, p_fs_model->a_no_of_bs_states);
	hs_m_pops = gsl_matrix_alloc(no_of_time_points, p_fs_model->p_m_scheme[0]->no_of_states);
	hs_c_pops = gsl_matrix_alloc(no_of_time_points, p_fs_model->p_c_scheme[0]->no_of_states);
	
	// Set to zero
	gsl_vector_set_zero(hs_time);
	gsl_vector_set_zero(hs_pCa);
	gsl_vector_set_zero(hs_length);
	gsl_vector_set_zero(hs_command_length);
	gsl_vector_set_zero(hs_slack_length);
	gsl_vector_set_zero(hs_force);
	gsl_vector_set_zero(hs_titin_force);
	gsl_vector_set_zero(hs_viscous_force);
	gsl_vector_set_zero(hs_extracellular_force);

	gsl_vector_set_zero(hs_a_length);
	gsl_vector_set_zero(hs_m_length);

	gsl_vector_set_zero(hs_inter_hs_titin_force_effect);

	gsl_matrix_set_zero(hs_a_pops);
	gsl_matrix_set_zero(hs_m_pops);
	gsl_matrix_set_zero(hs_c_pops);
}

// Destructor
hs_data::~hs_data(void)
{
	// Recover space

	// First the gsl_vectors
	gsl_vector_free(hs_time);
	gsl_vector_free(hs_pCa);
	gsl_vector_free(hs_length);
	gsl_vector_free(hs_command_length);
	gsl_vector_free(hs_slack_length);
	gsl_vector_free(hs_force);
	gsl_vector_free(hs_titin_force);
	gsl_vector_free(hs_viscous_force);
	gsl_vector_free(hs_extracellular_force);

	gsl_vector_free(hs_a_length);
	gsl_vector_free(hs_m_length);

	gsl_vector_free(hs_inter_hs_titin_force_effect);

	// Then the matrices
	gsl_matrix_free(hs_a_pops);
	gsl_matrix_free(hs_m_pops);
	gsl_matrix_free(hs_c_pops);
}
