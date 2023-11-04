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
#include "half_sarcomere.h"
#include "muscle.h"
#include "global_definitions.h"
#include "JSON_functions.h"

#include "rapidjson\document.h"

#include "gsl_vector.h"
#include "gsl_math.h"
#include "gsl_const_mksa.h"


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

double transition::calculate_rate(double x, double x_ext, double node_force,
	int mybpc_state, int mybpc_iso, short int active_neigh, half_sarcomere* p_hs)
{
	//! Returns the rate for a transition with a given x
	//! 
	// printf("%s \n ", rate_type);

	// Variables
	double rate = 0.0;

	FiberSim_options* p_options;

	// Code

	// Set options
	p_options = p_parent_m_state->p_parent_scheme->p_fs_options;

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
			(1.0 + (gsl_max(node_force, 0.0) * gsl_vector_get(rate_parameters, 1)));
	}

	// Force and adjacent hs dependent
	if (!strcmp(rate_type, "force_and_adjacent_hs_dependent"))
	{
		// Assume that first half-sarcomere starts with it's Z-line on the left hand side

		int hs_ind;
		int hs_across_Z;
		int hs_across_M;

		half_sarcomere* p_hs_across_Z = NULL;
		half_sarcomere* p_hs_across_M = NULL;

		muscle* p_parent_m;

		FiberSim_model* p_fs_model;

		double factor_across_Z = 0;
		double factor_across_M = 0;

		// Calculate node_force factor
		double node_force_factor = gsl_max(node_force, 0.0) * gsl_vector_get(rate_parameters, 1);

		// Set the pointers
		hs_ind = p_hs->hs_id;
		p_parent_m = p_hs->p_parent_m;
		p_fs_model = p_hs->p_fs_model;

		// Need to decide if this is a half-sarcomere with Z-disk or M-line on left
		if (GSL_IS_EVEN(hs_ind))
		{
			// Half-sarcomere has Z-disk on low (left) side
			if (hs_ind > 1)
			{
				hs_across_Z = hs_ind - 1;
				p_hs_across_Z = p_parent_m->p_hs[hs_across_Z];
				factor_across_Z = gsl_vector_get(rate_parameters, 2) *
					(p_hs->hs_titin_force - p_hs_across_Z->hs_titin_force);
			}

			if (hs_ind < (p_fs_model->no_of_half_sarcomeres - 1))
			{
				hs_across_M = hs_ind + 1;
				p_hs_across_M = p_parent_m->p_hs[hs_across_M];
				factor_across_M = gsl_vector_get(rate_parameters, 3) *
					(p_hs->hs_titin_force - p_hs_across_M->hs_titin_force);
			}
		}
		else
		{
			// Half-sarcomere has M-line on low (left) side
			if (hs_ind > 0)
			{
				hs_across_M = hs_ind - 1;
				p_hs_across_M = p_parent_m->p_hs[hs_across_M];
				factor_across_M = gsl_vector_get(rate_parameters, 3) *
					(p_hs->hs_titin_force - p_hs_across_M->hs_titin_force);
			}

			if (hs_ind < (p_fs_model->no_of_half_sarcomeres - 1))
			{
				hs_across_Z = hs_ind + 1;
				p_hs_across_Z = p_parent_m->p_hs[hs_across_Z];
				factor_across_Z = gsl_vector_get(rate_parameters, 2) *
					(p_hs->hs_titin_force - p_hs_across_Z->hs_titin_force);
			}
		}

		rate = gsl_vector_get(rate_parameters, 0) *
			(1.0 + (node_force_factor + gsl_max(factor_across_Z, 0.0) + gsl_max(factor_across_M, 0.0)));

		//printf("hs_index: %i  titin_f: %g\t\tnode_f: %g\t\tf_across_Z: %g\t\tf_across_M: %g\n",
			//hs_ind, p_hs->hs_titin_force, node_force_factor, factor_across_Z, factor_across_M);
	}

	// Force and MyBPC-dependent
	if (!strcmp(rate_type, "force_and_mybpc_dependent"))
	{
		int no_of_c_states;

		double k_base = gsl_vector_get(rate_parameters, 0);
		double k_force = gsl_vector_get(rate_parameters, 1);

		double modifier_base = 1.0;
		double modifier_force = 1.0;

		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;

		// Deduce the no of c_states
		no_of_c_states = p_model->p_c_scheme[0]->no_of_states;

		// rate = k_base * modifier_base * (1 + node_force * k_force * modifier_force)
		// where the modifier depends on the isotype of MyBP-C
		// rate_parameters is a vector
		// [k_base, k_force,
		//		modifier_base[isotype=1,isostate=1], modifier_force[isotype=1, isostate=1]
		//		modifier_base[isotype=1,isostate=2], modifier_force[isotype=1, isostate=2]
		//		modifier_base[isotyp=1, isostate=3], modifier_force[isotype=1, isostate=3]
		//		... and so on for different isostates
		//		modifier_base[isotype=2,isostate=1], modifier_force[isotype=2, isostate=1]
		//		modifier_base[isotype=2,isostate=2], modifier_force[isotype=2, isostate=2]
		//		... and so on for different isotypes
		// So rate_parameters=[1, 10, 0.5, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1]
		// describes a situation with 2 isostypes, each with 3 states
		// where isotype 1, state 1 stabilizes SRX (modifer_base * 0.5)
		// and isotype 2, state 1, destabilizes SRX (modifier_base * 2)

		if (mybpc_state > 0)
		{
			// Work out how far to step into the rate_parameters vector
			int base_index = 2 + ((mybpc_iso - 1) * 2 * no_of_c_states) + ((mybpc_state - 1) * 2);
			int force_index = base_index + 1;

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
				(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature));
	}

	// Gaussian_hsl
	if (!strcmp(rate_type, "gaussian_hsl"))
	{
		// Distance between surface of thick and thin filaments is
		// (2/3)*d_1,0 - r_thin - r_thick
		// Assume d_1_0 at hsl = 1100 nm is 37 nm, r_thin = 5.5 nm, t_thick = 7.5 nm
		// d at hsl = x is (2/3) * (37 / sqrt(x/1100)) - 5/5 - 7.5
		// first passage time to position y is t = y^2 / (2*D)
		// rate is proportional to 1/t
		// rate at hsl == x is ref_rate * (y_ref / y_x)^2
		// See PMID 35450825 and first passage in
		// Mechanics of motor proteins and the cytoskeleton, Joe Howard book

		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k_cb = p_model->m_k_cb;

		double hs_length;
		double y_ref;		// distance between filaments at 1100 nm
		double y_actual;	// distance between filaments at current hsl
		double r_thick = 7.5;
		double r_thin = 5.5;

		y_ref = ((2.0 / 3.0) * 37.0) - r_thick - r_thin;

		if (p_hs == NULL)
			hs_length = 1100.0;
		else
			hs_length = p_hs->hs_length;

		y_actual = (2.0 / 3.0) * (37.0 / sqrt(hs_length / 1100.0)) - r_thick - r_thin;

		rate = gsl_vector_get(rate_parameters, 0) *
			exp(-(0.5 * k_cb * gsl_pow_int(x, 2)) /
				(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature));

		rate = rate * gsl_pow_2(y_ref / y_actual);
	}

	// Gaussian MyBP-C 

	if (!strcmp(rate_type, "gaussian_pc"))
	{
		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k_pc = p_model->c_k_stiff; // use MyBPC stiffness as default

		double temp = gsl_vector_get(rate_parameters, 1); // optional parameter which sets mybpc stiffness

		if (!gsl_isnan(temp)) { // optional parameter is not specified, use the state extension instead
			k_pc = temp;
		}

		rate = gsl_vector_get(rate_parameters, 0) *
			exp(-(0.5 * k_pc * gsl_pow_int(x, 2)) /
					(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature));
	}

	// Force-dependent Gaussian

	if (!strcmp(rate_type, "force_dependent_gaussian"))
	{

		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k_cb = p_model->m_k_cb;

		rate = gsl_vector_get(rate_parameters, 0) *
			exp(-(0.5 * k_cb * gsl_pow_int(x, 2)) /
					(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature)) *
			(1.0 + gsl_max(node_force, 0.0) * gsl_vector_get(rate_parameters, 1));
	}

	// Force-dependent Gaussian for MyBPC

	if (!strcmp(rate_type, "force_dependent_gaussian_pc"))
	{

		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k_pc = p_model->c_k_stiff;

		rate = gsl_vector_get(rate_parameters, 0) *
			exp(-(0.5 * k_pc * gsl_pow_int(x, 2)) /
					(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature)) *
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

	// Decreasing load-dependent exponential
	if (!strcmp(rate_type, "exp"))
	{
		double A = gsl_vector_get(rate_parameters, 0);
		double B = gsl_vector_get(rate_parameters, 1);
		double C = gsl_vector_get(rate_parameters, 2);
		double x_center = gsl_vector_get(rate_parameters, 3);
		double x_wall = gsl_vector_get(rate_parameters, 4);

		if (x < x_wall)
			rate = A + B * exp(-C * (x + x_center));
		else
			rate = p_options->max_rate;
	}

	if (!strcmp(rate_type, "exp_wall"))
	{
		// Variables

		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k0 = gsl_vector_get(rate_parameters, 0);
		double F = p_model->m_k_cb * (x + x_ext);
		double d = gsl_vector_get(rate_parameters, 1);
		double x_wall = gsl_vector_get(rate_parameters, 2);
		double x_smooth = gsl_vector_get(rate_parameters, 3);

		// Code
		rate = k0 * exp(-(F * d) /
				(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature));

		rate = rate + p_options->max_rate * (1 /
			(1 + exp(-x_smooth * (x - x_wall))));
	}

	if (!strcmp(rate_type, "exp_wall_sweep"))
	{
		// Variables

		FiberSim_model* p_model = p_parent_m_state->p_parent_scheme->p_fs_model;
		double k0 = gsl_vector_get(rate_parameters, 0);
		double F = p_model->m_k_cb * (x + x_ext);
		double d = gsl_vector_get(rate_parameters, 1);
		double x_wall = gsl_vector_get(rate_parameters, 2);
		double x_smooth = gsl_vector_get(rate_parameters, 3);
		double sweep = gsl_vector_get(rate_parameters, 4);

		// Code
		rate = k0 * exp(-(F * d) /
			(1e18 * GSL_CONST_MKSA_BOLTZMANN * p_model->temperature));

		rate = rate + p_options->max_rate * (1 /
			(1 + exp(-x_smooth * (x - x_wall))));

		// Add in the sweep, whereby neighboring units in the off state
		// increase the detachment rate
		rate = rate + sweep * (2.0 - (double)active_neigh);

		//printf("active_neigh: %g  rate: %f\n", (double)active_neigh, rate);
	}

	// Curtail at max rate

	if (rate > (p_options->max_rate))
		rate = p_options->max_rate;

	if (rate < 0.0)
		rate = 0.0;
	
	// Return
	return rate;
}

