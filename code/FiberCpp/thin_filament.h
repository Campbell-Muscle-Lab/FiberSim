#pragma once

/**
 * @file    thin_filament.h
 * @brief   header file for the thin_filament class
 * @author  Ken Campbell
 */

#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

#include "gsl_vector.h"

class FiberSim_model;
class FiberSim_options;
class half_sarcomere;

class thin_filament
{
public:
    // Variables
    FiberSim_model* p_fs_model;     /**< Pointer to the FiberSim_model */

    FiberSim_options* p_fs_options; /**< Pointer to the FiberSim_model */

    half_sarcomere* p_parent_hs;    /**< Pointer to the parent half-sarcomere */

    int thin_id;                    /**< integer identifier for the thick_filament */

    int a_regulatory_units_per_strand;
                                    /**< integer defining the number of regulatory
                                         units along a single trand */
    int a_bs_per_unit;              /**< integer defining the number of binding sites
                                         per strand */
    int a_strands_per_filament;     /**< integer defining the number of strands
                                         per filament */
    int a_no_of_bs;                 /**< integer defining the total number of binding
                                         sites on a thin filament */

    int a_bs_per_node;              /**< integer defining the number of binding sites
                                         per node */

    double a_inter_bs_rest_length;  /**< double defining the rest-length of the
                                         spring between binding-sites in nm */

    double a_inter_bs_twist;        /**< double defining the inter-bs twist
                                         along the thin filament in degrees */

    double a_y;                     /**< the y coordinate of the filament in nm */
    double a_z;                     /**< the z coordinate of the filament in nm */

    // Arrays
    gsl_vector* bs_x;               /**< pointer to a gsl_vector holding
                                         bs positions */
    gsl_vector* bs_angle;           /**< pointer to a gsl_vector holding
                                         bs angles */

    gsl_vector_short* bs_state;     /**< pointer to a gsl array of signed short
                                         integers indicating which state the
                                         bs is in
                                         1 means off 
                                         2 means on */

    gsl_vector_short* bs_isoform;   /**< pointer to a gsl array of signed short
                                         integers indicating the isoform of
                                         the bs */

    gsl_vector_short* bs_unit;      /**< pointer to a gsl array of signed short
                                         integers indicating the regulatory
                                         unit the binding site is part of */

    gsl_vector_short* unit_status;  /**< pointer to a gsl array of signed short
                                         integers indicating the status of
                                         the unit
                                         1 means off
                                         2 means on */

    gsl_vector_short* active_neighbors;
                                    /**< pointer to a gsl_array of signed short
                                         integers, indicating the number of
                                         adjacent regulatory units in the active
                                         state */

    gsl_vector_short* bound_to_m_type;
                                    /**< pointer to a gsl array of signed short
                                         integers indicating what the site is
                                         bound to:
                                         0 = nothing
                                         1 = myosin
                                         2 = mybpc
                                         */

    gsl_vector_short* bound_to_m_f; /**< pointer to a gsl array of signed short
                                         integers indicating which m_f the bs
                                         is bound to */
    gsl_vector_short* bound_to_m_n; /**< pointer to a gsl array of signed short
                                         integers indicating which m_n the bs
                                         is bound to */

    gsl_vector* node_forces;        /**< pointer to a gsl_vector holding force
                                         at each node */

    /**
    * Constructor
    * + initialized with:
    + pointer to a model
    + pointer to options
    + pointer to the parent half_sarcomere
    + an id number
    */
    thin_filament(
        FiberSim_model* set_p_fs_model,
        FiberSim_options* set_p_fs_options,
        half_sarcomere* p_parent_hs,
        int set_thin_id);

    /**
    * Destructor
    */
    ~thin_filament(void);

    /**
     * initialise_bs_x_bs_angle_bs_unit()
     * + for each bs on the filament, sets
     *   + the bs_x value
     *   + the bs_angle value
     *   + which unit the binding site is a member of
     */
    void initialise_bs_x_bs_angle_bs_unit(void);

    /**
    * set unit status
    * return void
    */
    void set_unit_status(void);

    /**
    * void set_regulatory_unit_bs_indices(int unit, int bs_indices[])
    * fills an integer array with the indices of the specified regulatory unit
    * @return void
    */
    void set_regulatory_unit_indices(int unit_ind, gsl_vector_short*  bs_indices);

    /**
    * void calculate node forces(void)
    * calculates the forces at each node
    * @return(void)
    */
    void calculate_node_forces(void);

};
