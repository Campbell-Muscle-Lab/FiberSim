#pragma once

/** 
 * @file    muscle.h
 * @brief   header file for the Muscle class
 * @author  Ken Campbell
 */

#include <filesystem>

#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

#include "gsl_vector.h"

#include "global_definitions.h"

#include "BS_thread_pool.hpp"

class FiberSim_model;
class FiberSim_protocol;
class FiberSim_options;
class FiberSim_data;
class half_sarcomere;
class series_component;

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

    /**
    * void write_rates_file()
    * writes the m and c rate functions to file
    * @ return void
    */
    void write_rates_file();

    /**
    */
    void force_control_muscle_system();

    size_t length_control_myofibril_with_series_compliance(int protocol_index);

    //int wrapper_length_control_myofibril_with_series_compliance(const gsl_vector* x, void* p, gsl_vector* f);

    size_t worker_length_control_myofibril_with_series_compliance(const gsl_vector* x, void* p, gsl_vector* f);

    //int length_control_wrapper(const gsl_vector* x, void* params, gsl_vector* f);

    // Variables

    char model_version[_MAX_PATH];      /**< FiberSim version from the model file */

    int muscle_id;                      /**< integer labeling the muscle */

    double m_length;                    /**< double defining the length of the muscle in nm */

    double m_force;                     /**< double defining the stress in the muscle in N m^-2 */

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

    series_component* p_sc;             /**< Pointer to a series elastic component */

    BS::thread_pool pool;                /**< Pointer to a worker pool */

    // Functions
};