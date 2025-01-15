#pragma once

/**
* @file		iso_type.h
* @brief	Header file for the iso_type class
* @author	Ken Campbell
*/

#include "gsl_vector.h"
#include "global_definitions.h"

#include "JSON_functions.h"
#include "rapidjson\document.h"


class iso_scheme;
class iso_transition;

class iso_type
{
public:

	// Variables

	iso_scheme* p_parent_iso_scheme;
									/**< pointer to the parent iso scheme */

	int iso_number;					/**< integer defining the iso number */

	int no_of_transitions;			/**< integer defining the number of transitions
										 from the iso_typee */

	iso_transition* p_iso_transitions[MAX_NO_OF_TRANSITIONS];

	/**
	* Constructor
	*/
	iso_type(const rapidjson::Value& iso_ty, iso_scheme* set_p_parent_scheme);

	/**
	* Destuctor
	*/
	~iso_type(void);

};

