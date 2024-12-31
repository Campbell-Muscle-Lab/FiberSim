/**
* @file		iso_transition.cpp
* @brief	Source file for the iso_transition class
* @author	Ken Campbell
*/

#include <cstdio>
#include <math.h>

#include "FiberSim_options.h"

#include "FiberSim_model.h"
#include "half_sarcomere.h"
#include "muscle.h"

#include "iso_transition.h"
#include "iso_type.h"
#include "iso_scheme.h"

#include "global_definitions.h"
#include "JSON_functions.h"

#include "rapidjson\document.h"

#include "gsl_vector.h"
#include "gsl_math.h"
#include "gsl_const_mksa.h"


// Constructor
iso_transition::iso_transition(const rapidjson::Value& iso_tr, iso_type* set_parent_isotype)
{
	// Set p_parent_m_state
	p_parent_iso_type = set_parent_isotype;

	JSON_functions::check_JSON_member_int(iso_tr, "new_number");
	new_iso_type = iso_tr["new_number"].GetInt();

	JSON_functions::check_JSON_member_string(iso_tr, "rate_type");
	sprintf_s(rate_type, _MAX_PATH, iso_tr["rate_type"].GetString());

	// Read in parameters
	JSON_functions::check_JSON_member_array(iso_tr, "rate_parameters");
	const rapidjson::Value& rp = iso_tr["rate_parameters"];

	rate_parameters = gsl_vector_alloc(MAX_NO_OF_RATE_PARAMETERS);
	gsl_vector_set_all(rate_parameters, GSL_NAN);

	for (int i = 0; i < (int)rp.Size(); i++)
	{
		gsl_vector_set(rate_parameters, i, rp[i].GetDouble());
	}
}

iso_transition::iso_transition()
{
	// Default constructor - used if there is no defined transition
	p_parent_iso_type = 0;
	new_iso_type = 0;
	sprintf_s(rate_type, _MAX_PATH, "");
	rate_parameters = gsl_vector_alloc(MAX_NO_OF_RATE_PARAMETERS);
	gsl_vector_set_all(rate_parameters, GSL_NAN);
}

// Destructor
iso_transition::~iso_transition(void)
{
	// Tidy up
	gsl_vector_free(rate_parameters);
}

// Functions
double iso_transition::calculate_rate(double x, half_sarcomere* p_hs)
{
	//! Returns the rate for a transition where the cross-bridge is in a given position x

	// Variables
	double rate = 0.0;						// transition rate

	FiberSim_options* p_options;			// pointer to options

	// Code

	// Set options
	p_options = p_parent_iso_type->p_parent_iso_scheme->p_fs_options;

	// Constant
	if (!strcmp(rate_type, "constant"))
	{
		double amp = gsl_vector_get(rate_parameters, 0);

		rate = amp;
	}

	if (!strcmp(rate_type, "exp_from_Z_line"))
	{
		double amp = gsl_vector_get(rate_parameters, 0);
		double k = gsl_vector_get(rate_parameters, 1);

		rate = amp * exp(-k * x);
	}

	if (!strcmp(rate_type, "wall_from_Z_line"))
	{
		double amp = gsl_vector_get(rate_parameters, 0);
		double wall = gsl_vector_get(rate_parameters, 1);

		if (fabs(x) <= wall)
			rate = amp;
		else
			rate = 0;
	}

	// Curtail at max rate

	if (rate > (p_options->max_rate))
		rate = p_options->max_rate;

	if (rate < 0.0)
		rate = 0.0;
	
	// Return
	return rate;
}
