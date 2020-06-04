/**
 * @file    FiberSim.cpp
 * @brief   Core file for the FiberSim model
 * @author  Ken Campbell
 */

#include <iostream>

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

#include "muscle.h"

// Pointers
muscle* p_muscle;       /**< pointer to a muscle */

int main(int argc, char* argv[])
{
/**
 Main function
 + the entry point for FiberSim
 */

    std::cout << "Hello World!\n";
    
    double x = 5.0;

/**
 a test of GSL
 */
    double y = gsl_sf_bessel_J0(x);
    printf("J0(%g) = %.18e\n", x, y);

    // Variables
    char model_file_string[_MAX_PATH];
    char options_file_string[_MAX_PATH];
    char protocol_file_string[_MAX_PATH];
    char results_folder[_MAX_PATH];

    // Unpack inputs
    sprintf_s(model_file_string, _MAX_PATH, "%s", argv[1]);
    sprintf_s(options_file_string, _MAX_PATH, "%s", argv[2]);
    sprintf_s(protocol_file_string, _MAX_PATH, "%s", argv[3]);
    sprintf_s(results_folder, _MAX_PATH, "%s", argv[4]);

    // Make a muscle
    p_muscle = new muscle(model_file_string, options_file_string);

    // Run the protocol
    p_muscle->implement_protocol(protocol_file_string, results_folder);

    // Tidy up
    delete p_muscle;
}
