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

    int dump_hs_status;                 /**< Integer defining whether or not to dump hs_status for
                                             each time-step in log/hs_status folder */

    char hs_status_folder[_MAX_PATH];
                                        /**< Folder to hold hs_status files. These are used for
                                             debugging and more advanced analyses */

    int no_of_repeats;                  /** integer defining the number of repeats to perform */

    int multithreading;                 /**< integer
                                            + 0 means no multithreading
                                            + 1 means multithreading */

    double x_pos_rel_tol;               /**< double defining the relative tolerance for calculating
                                             x positions */

    double max_rate;                    /**< double defining the maximum rate allowed in calculations */
    
    int dump_precision;                 /**< integer defining the precision of doubles dumped in the half-
                                         sarcomere status files */

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