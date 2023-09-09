/**
 * @file    half_sarcomere.cpp
 * @brief   Source file for the half_sarcomere class
 * @author  Ken Campbell
 */

#include <cstdio>
#include <chrono>

#include "series_component.h"
#include "muscle.h"
#include "FiberSim_model.h"
#include "FiberSim_options.h"

#include "gsl_math.h"


// Constructor
series_component::series_component(
	FiberSim_model* set_p_fs_model,
	FiberSim_options* set_p_fs_options,
	muscle* set_p_parent_m)
{
	// Initialise

	// Set the pointers
	p_fs_model = set_p_fs_model;
	p_fs_options = set_p_fs_options;
	p_parent_m = set_p_parent_m;

	// Set the stiffness
	sc_k_stiff = p_fs_model->sc_k_stiff;

	// Set the initial length to 0
	sc_extension = 0.0;

	// And the force to 0
	sc_force = 0.0;
}

// Destructor
series_component::~series_component()
{
	printf("\n\nsc_k_stiff: %g\n\n", sc_k_stiff);

	// Tidy up
}

// Other functions

double series_component::return_series_extension(double muscle_force)
{
	//! Returns the extension of the series component for a given force
	
	// Variables
	double ext;

	// Code
	if (gsl_isnan(sc_k_stiff))
	{
		printf("Error in series_component::return_series_extension, sc_k_stiff is NAN");
		exit(1);
	}

	ext = muscle_force / sc_k_stiff;

	return ext;
}

double series_component::return_series_force(double series_extension)
{
	//! Returns the force in the series component for a given extension

	// Variables
	double series_force;

	// Code
	series_force = series_extension * sc_k_stiff;

	return series_force;
}

