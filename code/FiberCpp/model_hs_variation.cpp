/**
* @file		model_hs_variation.cpp
* @brief	Source file for the model_hs_variation class
* @author	ken Campbell
*/

#include <iostream>
#include <filesystem>

#include "model_hs_variation.h"
#include "FiberSim_model.h"

#include "rapidjson/document.h"
#include "rapidjson\filereadstream.h"
#include "JSON_functions.h"

#include "gsl_vector.h"

using namespace std::filesystem;

// Constructor
model_hs_variation::model_hs_variation(
	FiberSim_model* set_p_fs_model,
	const rapidjson::Value& hsv)
{
	// Initialise
	p_fs_model = set_p_fs_model;

	// Now parse the structure
	JSON_functions::check_JSON_member_string(hsv, "variable");
	sprintf_s(model_variable, _MAX_PATH, hsv["variable"].GetString());

	JSON_functions::check_JSON_member_array(hsv, "multiplier");
	const rapidjson::Value& a = hsv["multiplier"];

	hs_multiplier = gsl_vector_alloc(a.Size());
	for (int i = 0; i < (int)a.Size(); i++)
	{
		gsl_vector_set(hs_multiplier, i, a[i].GetDouble());
	}
}

// Destructor
model_hs_variation::~model_hs_variation(void)
{
	// Recover space

	// First the gsl_vectors
	gsl_vector_free(hs_multiplier);
}
