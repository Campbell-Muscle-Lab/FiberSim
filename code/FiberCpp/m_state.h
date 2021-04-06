#pragma once

/**
* @file		m_state.h
* @brief	Header file for the m_state class
* @author	Ken Campbell
*/

#include "gsl_vector.h"
#include "global_definitions.h"

#include "JSON_functions.h"
#include "rapidjson\document.h"


class kinetic_scheme;
class transition;

class m_state
{
public:

	// Variables

	kinetic_scheme* p_parent_scheme;
									/**< pointer to the parent kinetic scheme */

	int state_number;				/**< integer defining the state number */

	char state_type;				/**< char defining the state type
										 'S', super-relaxed
										 'D', disorded-relaxed
										 'A', attached */

	double extension;				/**< double defing the link extension in nm */

	transition* p_transitions[MAX_NO_OF_TRANSITIONS];

	/**
	* Constructor
	*/
	m_state(const rapidjson::Value& m_st, kinetic_scheme* set_p_parent_scheme);

	/**
	* Destuctor
	*/
	~m_state(void);

};

