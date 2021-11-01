/**
 * @file    FiberSim.cpp
 * @brief   Core file for the FiberSim model
 * @author  Ken Campbell
 */

#include <stdio.h>
#include <iostream>
#include <filesystem>
#include <sstream>

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
	// Function variables

	bool incompatible_version = false;

	// Vector of string to save tokens
	vector <string> token_m;
	vector <string> token_c;

	// stringstream class
	stringstream model_version(a);
	stringstream code_version(b);
	string intermediate;

	// Tokenizing 
	while (getline(model_version, intermediate, '.'))
	{
		token_m.push_back(intermediate);
	}

	while (getline(code_version, intermediate, '.'))
	{
		token_c.push_back(intermediate);
	}

	// Check compatibility 

	if (token_c[0] > token_m[0])
	{
		incompatible_version = true;
	}
	else if(token_c[1] < token_m[1])
	{
		incompatible_version = true;
	}

	if (incompatible_version)
	{
		cout << "FiberSim version problem" << "\n";
		cout << "Code version: " << b << "\n";
		cout << "Model version: " << a << "\n";
		exit(1);
	}

}
