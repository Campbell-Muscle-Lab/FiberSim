#pragma once

/**
* @file		FiberSim_data.h
* @brief	header file for the FiberSim_data class
* @uathor	Ken Campbell
*/

#include "gsl_vector.h"
#include "gsl_matrix.h"

class FiberSim_options;
class FiberSim_model;

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
	
	gsl_vector* fs_pCa;			/**< gsl_vector holding pCa for each time-point */

	gsl_vector* fs_length;		/**< gsl_vector holding hs_length (nm) for each time-point */

	gsl_vector* fs_command_length;
								/**< gsl_vector holding command_length for each time point */

	gsl_vector* fs_slack_length;
								/**< gsl_vector holding slack length for each time-point */

	gsl_vector* fs_a_length;	/**< gsl_vector holder thin_filament length (nm) for
										 each time-point */

	gsl_vector * fs_m_length;	/**< gsl_vector holder thick_filament length (nm) for
									 each time-point */

	gsl_vector* fs_force;		/**< gsl_vector holding hs_force for each time-point */

	gsl_vector* fs_titin_force;	/**< gsl_vector holding hs_titin_force for each time-point */

	gsl_vector* fs_extracellular_force;
								/**< gsl_vector holding hs_extracellular_force for each time-point */

	gsl_matrix* fs_a_pops;		/**< gsl_matrix holding the proportion of binding sites
									 in each state at each time-point */

	gsl_matrix* fs_m_pops;		/**< gsl_matrix holding the proportion of cross-bridges
									 in each state at each time-point */

	gsl_matrix* fs_c_pops;		/**< gsl_matrix holding the proportion of MyBPC
									 in each state at each time-point */
};
