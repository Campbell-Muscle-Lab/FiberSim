/**
 * @file    FiberSim.cpp
 * @brief   Core file for the FiberSim model
 * @author  Ken Campbell
 */

#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "muscle.h"
#include "FiberSim_version.h"

#include "gsl_math.h"
#include <math.h>

// Pointers
muscle* p_muscle;       /**< pointer to a muscle */

using namespace std;

void version_comp(string a, string b);

int main(int argc, char* argv[])
{
/**
 Main function
 + the entry point for FiberSim
 */

    // Variables
    char model_file_string[_MAX_PATH];
    char options_file_string[_MAX_PATH];
    char protocol_file_string[_MAX_PATH];
    char results_file_string[_MAX_PATH];

    // Unpack inputs
    sprintf_s(model_file_string, _MAX_PATH, "%s", argv[1]);
    sprintf_s(options_file_string, _MAX_PATH, "%s", argv[2]);
    sprintf_s(protocol_file_string, _MAX_PATH, "%s", argv[3]);
    sprintf_s(results_file_string, _MAX_PATH, "%s", argv[4]);

    // Display
    printf("FiberSim: starting\n");

    // Make a muscle
    p_muscle = new muscle(model_file_string, options_file_string);

	// Make sure the model version and the code version are compatible
	version_comp(p_muscle->model_version, code_version);

    // Run the protocol
    p_muscle->implement_protocol(protocol_file_string, results_file_string);

    // Tidy up
    delete p_muscle;

    // Display
    printf("FiberSim: closing correctly\n");
}

void version_comp(string a, string b)

{
	// Vector of string to save tokens
	vector <string> tokens;
	vector <string> tokens2;

	// stringstream class check1
	stringstream check1(a);

	// stringstream class check1
	stringstream check2(b);

	string intermediate;

	// Tokenizing w.r.t. space ' '
	while (getline(check1, intermediate, '.'))
	{
		tokens.push_back(intermediate);
	}

	// Tokenizing w.r.t. space ' '
	while (getline(check2, intermediate, '.'))
	{
		tokens2.push_back(intermediate);
	}

	if (tokens[0] < tokens2[0])
	{
		printf("NOT COMPATIBLE VERSIONS\n");
		exit(1);
	}

}
