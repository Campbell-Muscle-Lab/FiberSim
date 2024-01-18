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
#include "kinetic_scheme.h"
#include "m_state.h"
#include "transition.h"

#include "JSON_functions.h"
#include "rapidjson\document.h"

using namespace std::filesystem;

// Constructor
kinetic_scheme::kinetic_scheme(const rapidjson::Value& m_ks,
	FiberSim_model* set_p_fs_model, FiberSim_options * set_p_fs_options)
{
	// Initialise

	// Set the pointer to the model
	p_fs_model = set_p_fs_model;

	// Set the pointer to the options
	p_fs_options = set_p_fs_options;

	// Pull array
	JSON_functions::check_JSON_member_array(m_ks, "state");
	const rapidjson::Value& state = m_ks["state"];

	// Deduce the max number of transitions - we need this to make the states
	max_no_of_transitions = 0;
	for (rapidjson::SizeType i = 0; i < state.Size(); i++)
	{
		const rapidjson::Value& trans = state[i]["transition"];
		max_no_of_transitions = GSL_MAX(max_no_of_transitions, (int)trans.Size());
	}

	// Now cycle through making states
	for (rapidjson::SizeType i = 0; i < state.Size(); i++)
	{
		p_m_states[i] = new m_state(state[i], this);
	}

	no_of_states = state.Size();

	// Now that we know the state properties, set the transition types
	set_transition_types();
}

// Destructor
kinetic_scheme::~kinetic_scheme(void)
{
	//! Destructor

	// Tidy up

	// Delete the kinetic states
	for (int state_counter = 0; state_counter < no_of_states;
		state_counter++)
	{
		delete p_m_states[state_counter];
	}
}

// Functions
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

void kinetic_scheme::write_rate_functions_to_file(char output_file_string[], char file_write_mode[],
	char JSON_append_string[], half_sarcomere* p_hs)
{
	//! Writes rate functions to output file

	// Variables
	int counter = 0;

	double x_limit = 10;
	double x_bin = 0.5;

	m_state* p_m_state;			// pointer to an m_state
	transition* p_trans;		// pointer to a transition

	int new_state;				// integer for new state

	FILE* output_file;

	// Code

	// Make sure directory exists
	path output_file_path(output_file_string);

	if (!(is_directory(output_file_path.parent_path())))
	{
		if (create_directories(output_file_path.parent_path()))
		{
			std::cout << "\nCreating folder: " << output_file_path.string() << "\n";
		}
		else
		{
			std::cout << "\nError: Folder for rates file could not be created: " <<
				output_file_path.parent_path().string() << "\n";
			exit(1);
		}
	}

	// Check file can be opened, abort if not
	// file write status should be "w" for new, or "a" for append
	errno_t err = fopen_s(&output_file, output_file_string, file_write_mode);

	if (err != 0)
	{
		printf("write_rate_functions_to_file(): %s\ncould not be opened\n",
			output_file_string);
		exit(1);
	}

	// Write the JSON bracket
	fprintf_s(output_file, "\t\t\t{\n\"scheme\":\n\"\n");

	// Cycle through transitions and rates writing the column headers
	for (int state_counter = 0; state_counter < no_of_states; state_counter++)
	{
		p_m_state = p_m_states[state_counter];

		for (int t_counter = 0; t_counter < max_no_of_transitions; t_counter++)
		{
			p_trans = p_m_state->p_transitions[t_counter];
			new_state = p_trans->new_state;

			if (new_state > 0)
			{
				// It's a transition
				counter = counter + 1;
				if (counter == 1)
					fprintf_s(output_file, "x\tr_%i", counter);
				else
					fprintf_s(output_file, "\tr_%i", counter);
			}
		}
	}
	fprintf_s(output_file, "\n");

	// Cycle through x values and bins
	for (double x = -x_limit; x <= x_limit; x = x + x_bin)
	{
		fprintf_s(output_file, "%8g", x);

		for (int state_counter = 0; state_counter < no_of_states; state_counter++)
		{
			p_m_state = p_m_states[state_counter];

			for (int t_counter = 0; t_counter < max_no_of_transitions; t_counter++)
			{
				p_trans = p_m_state->p_transitions[t_counter];
				new_state = p_trans->new_state;

				if (new_state > 0)
				{
					// It's a transition
					double x_ext = p_m_state->extension;
					double rate = p_trans->calculate_rate(x, x_ext, 0, 0, 0, 0, p_hs);

					fprintf_s(output_file, "\t%8g", rate);
				}
			}
		}
		fprintf_s(output_file, "\n");
	}

	// Close the JSON bracket
	fprintf_s(output_file, "\"\n\t\t\t}%s\n", JSON_append_string);

	fclose(output_file);
}

void kinetic_scheme::write_kinetic_scheme_to_file(char output_file_string[])
{
	//! Writes kinetics scheme to output file

	// Variables
	FILE* output_file;

	// Code

	// Make sure irectory exists
	path output_file_path(output_file_string);

	if (!(is_directory(output_file_path.parent_path())))
	{
		if (create_directories(output_file_path.parent_path()))
		{
			std::cout << "\nCreating folder: " << output_file_path.string() << "\n";
		}
		else
		{
			std::cout << "\nError: Results folder could not be created: " <<
				output_file_path.parent_path().string() << "\n";
			exit(1);
		}
	}

	// Check file can be opened, abort if not
	errno_t err = fopen_s(&output_file, output_file_string, "w");
	if (err != 0)
	{
		printf("write_kinetic_scheme_to_file(): %s\ncould not be opened\n",
			output_file_string);
		exit(1);
	}

	// Kinetic scheme information
	fprintf_s(output_file, "{\n");
	fprintf_s(output_file, "\t\"m_kinetics\": {\n");
	fprintf_s(output_file, "\t\t\"no_of_states\": %i,\n", no_of_states);
	fprintf_s(output_file, "\t\t\"max_no_of_transitions\": %i,\n", max_no_of_transitions);
	fprintf_s(output_file, "\t\t\"scheme\": [\n");

	for (int state_counter = 0; state_counter < no_of_states; state_counter++)
	{
		fprintf_s(output_file, "\t\t{\n");
		fprintf_s(output_file, "\t\t\t\"number\": %i,\n", p_m_states[state_counter]->state_number);
		fprintf_s(output_file, "\t\t\t\"type\": \"%c\",\n", p_m_states[state_counter]->state_type);
		fprintf_s(output_file, "\t\t\t\"extension\": \"%g\",\n", p_m_states[state_counter]->extension);
		fprintf_s(output_file, "\t\t\t\"transition\":\n");
		fprintf_s(output_file, "\t\t\t[\n");

		for (int t_counter = 0; t_counter < max_no_of_transitions; t_counter++)
		{
			fprintf_s(output_file, "\t\t\t\t{\n");
			fprintf_s(output_file, "\t\t\t\t\t\"new_state\": %i,\n",
				p_m_states[state_counter]->p_transitions[t_counter]->new_state);
			fprintf_s(output_file, "\t\t\t\t\t\"transition_type\": \"%c\",\n",
				p_m_states[state_counter]->p_transitions[t_counter]->transition_type);
			fprintf_s(output_file, "\t\t\t\t\t\"rate_type\": \"%s\",\n",
				p_m_states[state_counter]->p_transitions[t_counter]->rate_type);
			fprintf_s(output_file, "\t\t\t\t\t\"rate_parameters\": [");

			for (int p_counter = 0; p_counter < MAX_NO_OF_RATE_PARAMETERS; p_counter++)
			{
				fprintf_s(output_file, "%g",
					gsl_vector_get(p_m_states[state_counter]->p_transitions[t_counter]->rate_parameters,
						p_counter));
				if (p_counter == (MAX_NO_OF_RATE_PARAMETERS - 1))
					fprintf_s(output_file, "]\n");
				else
					fprintf_s(output_file, ", ");
			}
			fprintf_s(output_file, "\t\t\t\t}");

			if (t_counter == (max_no_of_transitions - 1))
				fprintf_s(output_file, "\n");
			else
				fprintf_s(output_file, ",\n");
		}
		fprintf_s(output_file, "\t\t\t]\n");

		fprintf_s(output_file, "\t\t}");
		if (state_counter == (no_of_states - 1))
			fprintf_s(output_file, "\n");
		else
			fprintf_s(output_file, ",\n");
	}
	fprintf_s(output_file, "\t}\n");
	fprintf_s(output_file, "}\n");

	// Tidy up
	fclose(output_file);
}
