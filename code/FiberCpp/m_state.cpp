/**
* @file		m_state.cpp
* @brief	Source file for the m_state class
* @author	Ken Campbell
*/

#include <cstdio>

#include "m_state.h"
#include "kinetic_scheme.h"
#include "transition.h"

#include "JSON_functions.h"
#include "rapidjson\document.h"

// Constructor
m_state::m_state(const rapidjson::Value& m_st, kinetic_scheme* set_p_parent_scheme)
{
	
	char temp_string[_MAX_PATH];

	p_parent_scheme = set_p_parent_scheme;

	JSON_functions::check_JSON_member_int(m_st, "number");
	state_number = m_st["number"].GetInt();

	JSON_functions::check_JSON_member_string(m_st, "type");
	sprintf_s(temp_string, _MAX_PATH, m_st["type"].GetString());
	state_type = temp_string[0];

	JSON_functions::check_JSON_member_number(m_st, "extension");
	extension = m_st["extension"].GetDouble();

	// Pull array of transitions
	JSON_functions::check_JSON_member_array(m_st, "transition");
	const rapidjson::Value& trans = m_st["transition"];

	for (int i = 0; i < (int)trans.Size(); i++)
	{
		p_transitions[i] = new transition(trans[i],this);
	}
	// Fill in gaps if not all transitions are set
	if (trans.Size() < (size_t)(p_parent_scheme->max_no_of_transitions))
	{
		for (int i = trans.Size(); i < p_parent_scheme->max_no_of_transitions; i++)
		{
			p_transitions[i] = new transition();
		}
	}
}

// Destructor
m_state::~m_state(void)
{
	printf("in m_state destructor\n");

	// Tidy up
	for (int i = 0; i < p_parent_scheme->max_no_of_transitions; i++)
		delete p_transitions[i];

}