#pragma once

/**
 * @file    thick_filament.h
 * @brief   header file for the thick_filament class
 * @author  Ken Campbell
 */

#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

#include "gsl_vector.h"
#include "gsl_rng.h"

class FiberSim_model;
class FiberSim_options;
class half_sarcomere;

class thick_filament
{
public:
    /**
     * Constructor
     * + initialized with:
         + pointer to a model
         + pointer to options
         + pointer to the parent half_sarcomere
         + an id number
     */
    thick_filament(
        FiberSim_model* set_p_fs_model,
        FiberSim_options* set_p_fs_options,
        half_sarcomere* p_parent_hs,
        int set_thick_id);

    /**
    * Destructor
    */
    ~thick_filament(void);

    /**
     * initialise_cb_x_and_cb_angle()
     * + for each cb on the filament, sets
     *   + the cb_x value
     *   + the cb_angle value
     */
    void initalise_cb_x_and_cb_angle(void);

    /**
    * initialise_pc_node_index_and_pc_angle()
    * + for each pc on the filament, sets
    *   + the node_index, the node where the pc is
    *   + the pc_angle value
    */
    void initalise_pc_node_index_and_pc_angle(void);

    /**
    * initialise_cb_controlling_pc_index()
    * sets cb_controlling_pc_index
    * + -1 if there is no MyBPC
    * + index (>=0) of controlling MyBPC
    * @return void
    */
    void initialise_cb_controlling_pc_index(void);

    /**
    * void calculate_node_forces(void)
    * calculates the force at each node
    * @return void
    */
    void calculate_node_forces(void);

    // Variables
    FiberSim_model* p_fs_model;         /**< Pointer to the FiberSim_model */

    FiberSim_options* p_fs_options;     /**< Pointer to the FiberSim_model */

    half_sarcomere* p_parent_hs;        /**< Pointer to the parent half-sarcomere */

    int thick_id;                       /**< integer identifier for the thick_filament */

    int m_crowns_per_filament;          /**< integer defining the number of crowns per thick filament */
    int m_hubs_per_crown;               /**< integer defining number of hubs per crown */
    int m_myosins_per_hub;              /**< integer defining the number of myosins per hub */
    int m_no_of_cbs;                    /**< integer defining the total number of myosins per filament */

    double m_inter_crown_rest_length;   /**< double defining the rest-length of the spring
                                             between crowns in nm */
    double m_lambda;                    /**< double defining the barezone length in nm */

    double m_starting_angle;            /**< double defining the starting angle for the first crown */

    double m_inter_crown_twist;         /**< double defining the inter-crown twist between
                                             crowns  in degrees*/
    double m_within_hub_twist;          /**< double definiting the angle between myosins
                                             within a hub */

    double m_y;                         /**< the y coordinate of the filament in nm */
    double m_z;                         /**< the z coordinate of the filament in nm */

    // Arrays
    gsl_vector* cb_x;                   /**< pointer to gsl_vector holding cb_x positions */
    gsl_vector* cb_angle;               /**< pointer to gsl_vector holding cb_angles */

    gsl_vector* cb_nearest_bs_angle_diff;
                                        /**< pointer to gsl_vector holding angular difference
                                             between cb and nearest bs */

    signed short int* cb_state;         /**< pointer to an array of signed short integers
                                             indicating which state the cb is in */
    
    signed short int* cb_iso;           /**< pointer to an array of signed short integers
                                             indicating which isoform type the cb is */

    signed short int* cb_bound_to_a_f;  /**< pointer to an array of signed short integers
                                             indicating which a_f, the cb is bound to */
    signed short int* cb_bound_to_a_n;  /**< pointer to an array of signed short integers
                                             indicating which a_n on the corresponding a_f
                                             the cb is bound to */

    signed short int* cb_nearest_a_f;   /**< pointer to an array of signed short integers
                                             indicating which a_f, the cb is nearest to */
    signed short int* cb_nearest_a_n;   /**< pointer to an array of signed short integers
                                             indicating which a_n on the corresponding a_f
                                             the cb is nearest to */

    signed short int* cb_controlling_pc_index;
                                        /**< pointer to an array of signed short integers
                                             indicating the index of the MyBP-C that
                                             influences the myosin head */

    // MyBPC

    // Variables

    int c_thick_proximal_node;          /**< integer defining the first node for MyBPC */
    int c_thick_stripes;                /**< integer defining the number of MyBPC stripes */
    int c_thick_node_spacing;           /**< integer defining the node spacing for MyBPC */
    int c_mols_per_node;                /**< integer defining the number of MyBPC per node */

    int c_no_of_pcs;                    /**< integer defining the number of MyBPC
                                             molecues */

    double c_starting_angle;            /**< double defining the angle of the first MyBPC
                                             on each crown */

    // Arrays

    gsl_vector* pc_angle;               /**< pointer to gsl_vector holding pc_x positions */

    signed short int* pc_node_index;    /**< pointer to an array of signed short integers
                                             indicating which node the pc is associated
                                             with */

    signed short int* pc_state;         /**< pointer to an array of signed short integers
                                             indicating which state the pc is in */

    signed short int* pc_iso;          /**< pointer to an array of signed short integers
                                             indicating where pc is phosphorylated */

    signed short int* pc_bound_to_a_f;  /**< pointer to an array of signed short integers
                                             indicating which a_f the pc is bound to */
    signed short int* pc_bound_to_a_n;  /**< pointer to an array of signed short integers
                                             indicating which a_n the pc is bound to */

    signed short int* pc_nearest_a_f;   /**< pointer to an array of signed short integers
                                             indicating which a_f the pc is nearest to */
    signed short int* pc_nearest_a_n;   /**< pointer to an array of signed short integers
                                             indicating which a_n the pc is nearest to */

    // Nodes
    gsl_vector* node_forces;            /**< pointer to gsl_vector holding force at each node */

    // Random numbers
    gsl_rng* rand_generator_iso;            /**< pointer to a random number generator */
};

