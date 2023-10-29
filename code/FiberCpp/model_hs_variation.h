#pragma once

/**
* @file		model_hs_variation.h
* @brief	header file for the model_hs_variation class
* @author	Ken Campbell
*/

#include <stdio.h>

#include "rapidjson/document.h"

#include "gsl_vector.h"
#include "gsl_matrix.h"

class FiberSim_model;

class model_hs_variation
{
public:
	/**
	* Constructor
	* param integer number of time-points
	*/
	model_hs_variation(FiberSim_model* set_p_fs_model, const rapidjson::Value& hsv);

	/**
	* Destructor
	*/
	~model_hs_variation(void);

	// Variables

	FiberSim_model* p_fs_model;			/**< pointer to parent model */

	char model_variable[_MAX_PATH];		/**< character array for the model parameter that
												changes across half-sarcomeres */
	
	gsl_vector* hs_multiplier;			/**< gsl_vector holding multipliers for each
												half-sarcomere */
};