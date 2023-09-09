#pragma once

/**
* @file		FiberSim_data.h
* @brief	header file for the FiberSim_data class
* @uathor	Ken Campbell
*/

#include "gsl_vector.h"
#include "gsl_matrix.h"

#include "global_definitions.h"

class FiberSim_options;
class FiberSim_model;
class hs_data;

class FiberSim_data
{
public:
	/**
	* Constructor
	* param integer number of time-points
	*/
	FiberSim_data(int no_of_time_points,
		FiberSim_options* set_p_fs_options,
		FiberSim_model* set_p_fs_model);

	/**
	* Destructor
	*/
	~FiberSim_data(void);

	/**
	* void write_data_to_delimited_file(char output_file_string[], char delimiter)

	* @param output_file_string[] a character array holding the output file name
	* @param delimiter a character holding the delimiter
	* @return void
	*/
	void write_data_to_delimited_file(char output_file_string[], char delimiter = '\t');

	// Variables
	int no_of_time_points;		/**< integer number of time-points in the simulation */

	FiberSim_options* p_fs_options;
								/**< pointer to a FiberSim options object */

	FiberSim_model* p_fs_model;
								/**< pointer to a FiberSim_model object
									 used to set the size of the matrices holding
									 the bs and cb state variables */

	gsl_vector* fs_time;		/**< gsl_vector holding time in second at
									 each time-point */

	gsl_vector* fs_m_length;	/**< gsl_vector holding muscle length (nm) for each time-point */

	gsl_vector* fs_m_force;		/**< gsl_vector holding muscle force (N m^-2) for each time-point */

	gsl_vector* fs_sc_extension;
								/**< gsl_vector holding series component extension (nm)
									 for each time-step */

	gsl_vector* fs_sc_force;
								/**< gsl_vector holding series component force
									 in N m^-2 for each time-step */

	hs_data* p_hs_data [MAX_NO_OF_HALF_SARCOMERES];
								/**< pointer to an array of data objects for
									 individual half-sarcomeres */
};
