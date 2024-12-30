/**
* @file		kinetic_scheme.cpp
* @brief	Source file for the kinetic_scheme class
* @author	Ken Campbell
*/

#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "FiberSim_model.h"
#include "FiberSim_options.h"
#include "iso_scheme.h"
#include "iso_type.h"
#include "iso_transition.h"

#include "JSON_functions.h"
#include "rapidjson\document.h"

using namespace std::filesystem;

// Constructor
iso_scheme::iso_scheme(const rapidjson::Value& iso_sch,
	FiberSim_model* set_p_fs_model, FiberSim_options * set_p_fs_options)
{
	// Initialise

	// Set the pointer to the model
	p_fs_model = set_p_fs_model;

	// Set the pointer to the options
	p_fs_options = set_p_fs_options;

	// Check the isotypes
	JSON_functions::check_JSON_member_array(iso_sch, "type");
	const rapidjson::Value& types = iso_sch["type"];

	// Deduce the max number of transitions - we need this to make the types
	max_no_of_transitions = 0;
	for (rapidjson::SizeType i = 0; i < types.Size(); i++)
	{
		const rapidjson::Value& trans = types[i]["transition"];
		max_no_of_transitions = GSL_MAX(max_no_of_transitions, (int)trans.Size());
	}

	// Now cycle through making the isotypes
	for (rapidjson::SizeType i = 0; i < types.Size(); i++)
	{
		p_iso_types[i] = new iso_type(types[i], this);
	}

	no_of_isotypes = types.Size();

	printf("\nNo_of_isotypes: %i\n\n\n", no_of_isotypes);
	exit(1);

	// Now that we know the state properties, set the transition types
//	set_transition_types();
}

// Destructor
iso_scheme::~iso_scheme(void)
{
	//! Destructor

	// Tidy up

	// Delete the iso_types
	for (int type_counter = 0; type_counter < no_of_isotypes; type_counter++)
	{
		delete p_iso_types[type_counter];
	}
}

// Functions
/*
void kinetic_scheme::set_transition_types(void)
{
	//! Cycle through the states, identifying the transition type for each one

	// Code
	for (int state_counter = 0; state_counter < no_of_states; state_counter++)
	{
		char current_state_type = p_m_states[state_counter]->state_type;

		for (int t_counter = 0; t_counter < max_no_of_transitions; t_counter++)
		{
			int new_state = p_m_states[state_counter]->p_transitions[t_counter]->new_state;

			if (new_state == 0)
			{
				// Transition is not allowed - skip out
				continue;
			}

			char new_state_type = p_m_states[new_state - 1]->state_type;

			if (current_state_type == 'S')
			{
				if (new_state_type == 'S')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
				if (new_state_type == 'D')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
				if (new_state_type == 'A')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'a';
				}
			}

			if (current_state_type == 'D')
			{
				if (new_state_type == 'S')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
				if (new_state_type == 'D')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
				if (new_state_type == 'A')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'a';
				}
			}

			if (current_state_type == 'A')
			{
				if (new_state_type == 'S')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'd';
				}
				if (new_state_type == 'D')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'd';
				}
				if (new_state_type == 'A')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
			}
		}
	}
}

*/
