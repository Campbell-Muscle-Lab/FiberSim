#pragma once

/** 
 * @file    FiberSim_options.h
 * @brief   header file for FiberSim_options
 * @author  Ken Campbell
 */

// Definitions for JSON parsing
#ifndef _RAPIDJSON_DOCUMENT
#define _RAPIDJSON_DOCUMENT
#include "rapidjson/document.h"
#endif

#include "stdio.h"

class FiberSim_options
{
public:
    /**
     * Constructor
     */
    FiberSim_options(char JSON_options_file_string[]);

    /**
     * Destructor
     */
    ~FiberSim_options(void);

    // Variables
   
    double max_rate;                    /**< double defining the maximum rate allowed in calculations */

    int adjacent_bs;                    /**< default value of adjacent binding sites
                                                 available for myosin or mybpc attachment
                                                 0 restricts to nearest site */

    double x_pos_rel_tol;               /**< double defining the relative tolerance for calculating
                                             x positions */

    double lambda_jitter;               /**< double defining lambda jitter. The first myosin crown on
                                             each thick filament will be at
                                             x = lambda + rand()*lambda_jitter
                                             Set to 0 for no jitter (aligned filaments). See
                                             https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1635681/
                                             for justification */

    bool rand_jitter;                   /**< bool defining whether the randon number generator is
                                             initialized based on the system clock or not
                                             defaults to false if not specified in the  options file */
   
    int dump_precision;                 /**< integer defining the precision of doubles dumped in the half-
                                             sarcomere status files */

    char log_relative_to[_MAX_PATH];    /**< char array used to direct paths */

    char log_folder[_MAX_PATH];         /**< Folder to hold files about the program status:
                                             intended for debugging and testing purposes
                                            "none" means do not write files */

    char log_file_string[_MAX_PATH];    /**< Log file holding information about the program execution
                                             intended mainly for debugging and testing purposes */

    FILE* log_file;                     /**< Pointer to a log file */

    int log_mode;                       /**< Integer defining whether or not to dump program status
                                             0 means do not log
                                             1 means log
                                             Value is set by the presence of log_folder */

    char status_relative_to[_MAX_PATH]; /**< char array used to direct paths */

    char status_folder[_MAX_PATH];      /**< Folder to hold status files
                                             "none" means do not write files */

    char rate_relative_to[_MAX_PATH];  /**< char array used to direct paths */

    char rate_file_string[_MAX_PATH];   /**< file_string for rate functions.*/

    char time_steps_string[_MAX_PATH];  /**< String definning which time_steps to dump */

    int start_status_time_step;         /**< Integer of first time-step to dump status file */

    int stop_status_time_step;          /**< Integer of last time-step to dump status file */

    int skip_status_time_step;          /**< Integer of skips between dump of status file */

    // Functions

    /**
     * a function that sets FiberSim_options parameters from a JSON file
     * @param json_file_string the filename for the JSON file
     */
    void set_FiberSim_options_from_JSON_file_string(char JSON_file_string[]);


    /**
     * a function that writes model options to file
     */
    void write_FiberSim_options_to_file(void);
};