#pragma once

/**
* @file		kinetic_scheme.h
* @brief	header file for the kinetic_scheme class
* @author	Ken Campbell
*/

#include "rapidjson/document.h"
#include "JSON_functions.h"

#include "gsl_vector.h"

#include "FiberSim_model.h"
#include "FiberSim_options.h"

#include "global_definitions.h"

// Forward declarations
class m_state;
class half_sarcomere;

class kinetic_scheme
{
public:

	// Variables

	FiberSim_model* p_fs_model;			/**< pointer to a FiberSim_model object */

	FiberSim_options* p_fs_options;		/**< pointer to a FiberSim_options objects */

	int no_of_states;					/**< int defining the number of different states */

	int max_no_of_transitions;			/**< int defining the maximum number of transitions from a state */

	m_state* p_m_states[MAX_NO_OF_KINETIC_STATES];
										/**< pointer to an array of m_state objects */

	// Functions

	/**
	* Constructor
	* takes a FiberSim_model and parses it to give the kinetic scheme
	*/
	kinetic_scheme(const rapidjson::Value& m_ks, FiberSim_model* set_p_fs_model,
		FiberSim_options* set_p_fs_options);

	/**
	* Destructor
	*/
	~kinetic_scheme(void);

	/**
	* void set_transition_types(void)
	* loops through transitions from each state and sets the type, 'a' for attach,
	* 'd' for detach, and 'n' for neither
	* needs to be run after all the states are known so can't be included in the
	* transition constructor as these are built in sequence with each state
	* @return void
	*/
	void set_transition_types(void);

	/**
	* void write_kinetic_scheme_to_file(char output_file_string)
	* writes kinetic_scheme to specified file in JSON format
	* @return void
	*/
	void write_kinetic_scheme_to_file(char output_file_string[]);

	/**
	* void write_myosin_rate_functions_to_file(char output_file_string)
	* as a tab-delimited file
	* @return void
	*/
	void write_rate_functions_to_file(char output_file_string[], char file_write_mode[],
										char JSON_append_string[],
										half_sarcomere* p_hs = NULL);

};