#pragma once

/** 
 * @file    muscle.h
 * @brief   header file for the Muscle class
 * @author  Ken Campbell
 */

#include <filesystem>

#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

#include "global_definitions.h"

class FiberSim_model;
class FiberSim_protocol;
class FiberSim_options;
class FiberSim_data;
class half_sarcomere;

using namespace std::filesystem;

class muscle
{
public:
    /**
     * Constructor
     * @param set_model_file_string a char array with the file_string for the model description
     *          in JSON format
     * @param set_options_file_string a char array with the file_string for the simulation options
                in JSON format
     */
    muscle(char set_model_file_string[], char set_options_file_string[]);
    /*
    muscle::muscle(
        FiberSim_model* set_FiberSim_model,
        FiberSim_protocol* set_FiberSim_protocol,
        FiberSim_options* set_FiberSim_options,
        int set_thread_number);
        */


    /**
     * Destructor
     */
    ~muscle(void);

    /**
    * void implement_protocol(char set_protocol_file_string[])
    * runs a muscle through an experimental protocol
    * @param set_protocol_file_string[] a char array with the file string for the protocol
    * @return void
    */
    void implement_protocol(char set_protocol_file_string[], char set_results_file_string[]);

    /**
    * implement_time_step(int protocol_index)
    * @ param protocol_index, an integer defining the index in the protocol arrays
    */
    void implement_time_step(int protocol_index);

    // Variables

    //static inline std::string const version_number = "2.1.1";

    int muscle_id;                      /**< integer labeling the muscle */

    int dump_status_counter;            /**< Integer used to track which status files
                                             to dump */

    char model_file_string[_MAX_PATH];  /**< character array for the model file string */

    char options_file_string[_MAX_PATH];
                                        /**< character array for the simulation options file string */

    char protocol_file_string[_MAX_PATH];
                                        /**< character array for the simulation protocol file string */

    char results_file_string[_MAX_PATH];
                                        /**< character array for the simulation results file string */

    FiberSim_model* p_fs_model;         /**< pointer to FiberSim_model */

    FiberSim_options* p_fs_options;     /**< pointer to FiberSim_options */

    FiberSim_protocol* p_fs_protocol;   /**< Pointer to FiberSim_protocol */

    FiberSim_data* p_fs_data;           /**< Pointer to FiberSim_data */

    half_sarcomere * p_hs [MAX_NO_OF_HALF_SARCOMERES];
                                        /**< pointer to an array of half-sarcomere objects */

    // Functions
};