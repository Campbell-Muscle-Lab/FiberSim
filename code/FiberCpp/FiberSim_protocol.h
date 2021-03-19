#pragma once

/**
 * @file    FiberSim_protocol.h
 * @brief   header file for the FiberSim_protocol
 * @author  Ken Campbell
 */

//#include "rapidjson/document.h"

#include "global_definitions.h"
#include "gsl_vector.h"

class FiberSim_protocol
{
public:
    /**
    * Constructor - creates a protocol from a tab-delimited text file
    * param protocol_file_string a char array defining the protocol
    */
    FiberSim_protocol(char protocol_file_string[]);

    /**
     * Destructor
     */
    ~FiberSim_protocol(void);


    // Variables
    int no_of_time_points;          /**< integer number of time-points */

    gsl_vector* dt;                 /**< gsl_vector holding dt for
                                         each time-point */
    gsl_vector* delta_hsl;          /**< gsl_vector holding delta-hsl for
                                         each time-point */
    gsl_vector* sim_mode;           /**< gsl_vector holding sim_mode for
                                         each time-point */
    gsl_vector* pCa;                /**< gsl_vector holding pCa for
                                         each time-point */
};