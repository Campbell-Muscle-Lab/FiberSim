/**
* @file		transition.cpp
* @brief	Source file for the transition class
* @author	Ken Campbell
*/

#include <cstdio>
#include <math.h>

#include "transition.h"
#include "m_state.h"
#include "kinetic_scheme.h"
#include "FiberSim_model.h"
#include "global_definitions.h"
#include "JSON_functions.h"

#include "rapidjson\document.h"

#include "gsl_vector.h"
#include "gsl_math.h"


// Constructor
transition::transition(const rapidjson::Value& tr, m_state* set_p_parent_m_state)
{
	// Set p_parent_m_state
	p_parent_m_state = set_p_parent_m_state;

	// Set transition_type to unknown - will be set later on
	transition_type = 'x';

	JSON_functions::check_JSON_member_int(tr, "new_state");
	new_state = tr["new_state"].GetInt();

	JSON_functions::check_JSON_member_string(tr, "rate_type");
	sprintf_s(rate_type, _MAX_PATH, tr["rate_type"].GetString());

	// Read in parameters
	JSON_functions::check_JSON_member_array(tr, "rate_parameters");
	const rapidjson::Value& rp = tr["rate_parameters"];

	rate_parameters = gsl_vector_alloc(MAX_NO_OF_RATE_PARAMETERS);
	gsl_vector_set_all(rate_parameters, GSL_NAN);

	for (int i = 0; i < (int)rp.Size(); i++)
	{
		gsl_vector_set(rate_parameters, i, rp[i].GetDouble());
	}
}

transition::transition()
{
	// Default constructor - used if there is no defined transition
	new_state = 0;
	transition_type = 'x';
	sprintf_s(rate_type, _MAX_PATH, "");
	rate_parameters = gsl_vector_alloc(MAX_NO_OF_RATE_PARAMETERS);
	gsl_vector_set_all(rate_parameters, GSL_NAN);
}

// Destructor
transition::~transition(void)
{
	// Tidy up
	gsl_vector_free(rate_parameters);
}

// Functions

double transition::calculate_rate(double x, double x_ext, double node_force, int mybpc_state, int mybpc_iso)
{
	//! Returns the rate for a transition with a given x
	//! 
	// printf("%s \n ", rate_type);

	// Variables
	double rate = 0.0;

	// Code

	// Constant
	if (!strcmp(rate_type, "constant"))
	{
		rate = gsl_vector_get(rate_parameters, 0);
	}

	if (!strcmp(rate_type, "MyBPC_dependent"))
	{
		rate = gsl_vector_get(rate_parameters, 0);
		
		if (mybpc_state == 1)
		{
			rate = rate * gsl_vector_get(rate_parameters, 1);
		}
	}

	// Force-dependent
	if (!strcmp(rate_type, "force_dependent"))
	{
		rate = gsl_vector_get(rate_parameters, 0) *
			(1.0 + (gsl_max(node_force,0.0) * gsl_vector_get(rate_parameters, 1)));
	}

	// Force and MyBPC-dependent
	if (!strcmp(rate_type, "force_and_MyBPC_dependent"))
	{
		// rate = k_base * modifier_base * (1 + node_force * k_force * modifier_force)
		// where the modifier depends on the isotype
		// rate_parameters is a vector
		// [k_base, k_force,
		//		modifier_base[isotype=1,isostate=1], modifier_force[isotype=1, isostate=1]
		//		modifier_base[isotype=1,isostate=2], modifier_force[isotype=1, isostate=2]
		//		... and so on for different isostates
		//		modifier_base[isotype=2,isostate=1], modifier_force[isotype=2, isostate=1]
		//		modifier_base[isotype=2,isostate=2], modifier_force[isotype=2, isostate=2]
		//		... and so on for different isotypes
		// So rate_parameters=[1, 10, 0.5, 1, 1, 1, 1, 1, 2, 1]
		// describes a situation with 2 isostypes, each with 2 states
		// where isotype 1, state 1 stabilizes SRX (modifer_base * 0.5)
		// and isotype 2, state 2, destabilizes SRX (modifier_base * 2)
		// 
		// CURRENTLY, ONLY WORKS FOR ISOTYPES THAT ALL HAVE 2 POSSIBLE STATES
		// GOING FURTHER MIGHT NEED PASSING IN INFORMATION ABOUT THE ENTIRE SCHEME

		double k_base = gsl_vector_get(rate_parameters, 0);
		double k_force = gsl_vector_get(rate_parameters, 1);

		double modifier_base = 1.0;
		double modifier_force = 1.0;

		if (mybpc_state > 0)
		{
			// Work out how far to step into the rate_parameters vector
			int base_index = 2 + ((mybpc_iso - 1) * 4) + ((mybpc_state - 1) * 2);
			int force_index = 2 + ((mybpc_iso - 1) * 4) + ((mybpc_state - 1) * 2) + 1;

			modifier_base = gsl_vector_get(rate_parameters, base_index);
			modifier_force = gsl_vector_get(rate_parameters, force_index);
		}

		rate = k_base * modifier_base * (1.0 + (gsl_max(node_force,0.0) * k_force * modifier_force));
	}

	// Gaussian
	if (!strcmp(rate_type, "gaussian"))
	{
		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k_cb = p_model->m_k_cb;

		rate = gsl_vector_get(rate_parameters, 0) *
			exp(-(0.5 * k_cb * gsl_pow_int(x, 2)) /
			(1e18 * 1.38e-23 * 310.0));
	}

	// Force-dependent gaussian
	if (!strcmp(rate_type, "force_dependent_gaussian"))
	{
		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k_cb = p_model->m_k_cb;

		rate = gsl_vector_get(rate_parameters, 0) *
			(1.0 + (gsl_max(node_force, 0.0) * gsl_vector_get(rate_parameters, 1))) *
				exp(-(0.5 * k_cb * gsl_pow_int(x, 2)) /
					(1e18 * 1.38e-23 * 310.0));
	}


	// Gaussian MyBP-C 

	if (!strcmp(rate_type, "gaussian_pc"))
	{
		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k_pc = p_model->c_k_stiff; // use MyBPC stiffness

		rate = gsl_vector_get(rate_parameters, 0) *
			exp(-(0.5 * k_pc * gsl_pow_int(x, 2)) /
				(1e18 * 1.38e-23 * 310.0));
	}

	// Force-dependent Gaussian

	if (!strcmp(rate_type, "force_dependent_gaussian"))
	{

		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k_cb = p_model->m_k_cb;

		rate = gsl_vector_get(rate_parameters, 0) *
			exp(-(0.5 * k_cb * gsl_pow_int(x, 2)) /
				(1e18 * 1.38e-23 * 310.0)) *
			(1.0 + gsl_max(node_force, 0.0) * gsl_vector_get(rate_parameters, 1));
	}

	// Force-dependent Gaussian for MyBPC

	if (!strcmp(rate_type, "force_dependent_gaussian_pc"))
	{

		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k_pc = p_model->c_k_stiff;

		rate = gsl_vector_get(rate_parameters, 0) *
			exp(-(0.5 * k_pc * gsl_pow_int(x, 2)) /
				(1e18 * 1.38e-23 * 310.0)) *
			(1.0 + gsl_max(node_force, 0.0) * gsl_vector_get(rate_parameters, 1));
	}

	// Poly
	if (!strcmp(rate_type, "poly"))
	{
		double x_center = gsl_vector_get(rate_parameters, 3); // optional parameter defining the zero of the polynomial

		if (gsl_isnan(x_center)) { // optional parameter is not specified, use the state extension instead
			x_center = x_ext;
		}	

		rate = gsl_vector_get(rate_parameters, 0) +
				(gsl_vector_get(rate_parameters, 1) *
					gsl_pow_int(x + x_center, (int)gsl_vector_get(rate_parameters, 2)));

	}

	// Poly_asymmetric
	if (!strcmp(rate_type, "poly_asym"))
	{
		double x_center = gsl_vector_get(rate_parameters, 5); // optional parameter defining the zero of the polynomial

		if (gsl_isnan(x_center)) { // optional parameter is not specified, use the state extension instead
			x_center = x_ext;
		}

		if (x > x_center)
			rate = gsl_vector_get(rate_parameters, 0) +
				(gsl_vector_get(rate_parameters, 1) *
					gsl_pow_int(x + x_center, (int)gsl_vector_get(rate_parameters, 3)));
		else
			rate = gsl_vector_get(rate_parameters, 0) +
			(gsl_vector_get(rate_parameters, 2) *
				gsl_pow_int(x + x_center, (int)gsl_vector_get(rate_parameters, 4)));
	}

	// Curtail at max rate
	FiberSim_options* p_options = p_parent_m_state->p_parent_scheme->p_fs_options;
	if (rate > (p_options->max_rate))
		rate = p_options->max_rate;

	if (rate < 0.0)
		rate = 0.0;

	// Return
	return rate;
}

