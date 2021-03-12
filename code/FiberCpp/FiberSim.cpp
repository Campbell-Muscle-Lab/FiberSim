/**
 * @file    FiberSim.cpp
 * @brief   Core file for the FiberSim model
 * @author  Ken Campbell
 */

#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "muscle.h"

// Pointers
muscle* p_muscle;       /**< pointer to a muscle */

using namespace std::filesystem;

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

    // Make a muscle
    p_muscle = new muscle(model_file_string, options_file_string);

    /*
    // Create a path to the results folder
    path results_path(results_folder);
    */

    // Run the protocol
    p_muscle->implement_protocol(protocol_file_string, results_file_string);

    // Tidy up
    delete p_muscle;
}
