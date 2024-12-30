#pragma once

/**
* @file		transition.h
* @brief	Header file for the transition class
* @author	Ken Campbell
*/

#include "half_sarcomere.h"
#include "rapidjson/document.h"
#include "JSON_functions.h"
#include "FiberSim_options.h"

#include "gsl_vector.h"

// Forward declaration
class half_sarcomere;

class iso_transition
{
public:

	// Variables

	int parent_isotype;				/**< integer defining the parent isotype */

	int new_isotype;				/**< integer defining the new isotype */

	char rate_type[_MAX_PATH];		/**< char array defining the transition type */

	gsl_vector* rate_parameters;	/**< gsl_vector holding parameter variables */

	// Functions

	// Constructor

	/**
	* Normal constructor called with an entry into a JSON document
	*/
	iso_transition(const rapidjson::Value& iso_tr, int parent_isotype);

	/**
	* Constructor that sets default values
	* Called when the transition is empty
	*/
	iso_transition();

	/**
	* Destructor
	*/
	~iso_transition(void);

	// Functions

	/**
	* double calculate_rate(double x)
	* @param x double defining the cb_x position - the bs_x position
	* @return the rate in units of s^-1
	*/
/*	double calculate_rate(double x, double x_ext, double node_force,
							int mybpc_state, int mybpc_iso,
							short int active_neigh,
							half_sarcomere* p_hs);
*/
};
