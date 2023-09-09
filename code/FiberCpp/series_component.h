#pragma once

/** 
 * @file    series_component.h
 * @brief   header file for the Series_component class
 * @author  Ken Campbell
 */

#include <filesystem>

#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

#include "global_definitions.h"

class FiberSim_model;
class FiberSim_options;
class FiberSim_data;

class muscle;

using namespace std::filesystem;

class series_component
{
public:

    // Variables

    FiberSim_model* p_fs_model;                 /**< pointer to a FiberSim_model */

    FiberSim_options* p_fs_options;             /**< pointer to a FiberSim options */

    muscle* p_parent_m;                         /**< pointer to the parent muscle */

    double sc_extension;                        /**< double holding sc length in nm */

    double sc_k_stiff;                          /**< double holdering sc stiffness in N m^-1 */

    // Functions

    /**
     * Constructor
     * + initialized with
         + pointer to a model
         + pointer to options
         + pointer to the parent muscle
     */
    series_component(
        FiberSim_model* set_p_fs_model,
        FiberSim_options* set_pfs_options,
        muscle* p_parent_m
    );

    /**
     * Destructor
     */
    ~series_component(void);

    /**
    * double return_series_extension(double muscle_force)
    * returns the series extension for a given force
    * @param force the muscle force
    * @return double
    */
    double return_series_extension(double muscle_force);

    /**
    * double return_series_force(double series_extension)
    * returns the series force for a given extension
    * @param series_extension the extension of the series component
    * @return double
    */
    double return_series_force(double series_extension);

};