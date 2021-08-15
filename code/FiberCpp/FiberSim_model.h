#pragma once

/** 
 * @file    FiberSim_model.h
 * @brief   header file for the FiberSim_model
 * @author  Ken Campbell
 */

// Definitions for JSON parsing
#ifndef _RAPIDJSON_DOCUMENT
#define _RAPIDJSON_DOCUMENT
#include "rapidjson/document.h"
#endif

#include "global_definitions.h"
#include "gsl_vector.h"

class FiberSim_options;
class kinetic_scheme;

class FiberSim_model
{
public:
    /**
     * Constructor
     */
    FiberSim_model(char JSON_model_file_string[], FiberSim_options* set_p_fs_options);

    /**
     * Destructor
     */
    ~FiberSim_model(void);

    // Variables

    FiberSim_options* p_fs_options;     /**< Pointer to FiberSim_options object
                                             used for logging */

    // Muscle
    int no_of_half_sarcomeres;          /**< Number of half-sarcomeres in model */
    int no_of_myofibrils;               /**< Number of myofibrils in model */

    double initial_hs_length;           /**< double defining the initial hs_length  in nm */

    double prop_fibrosis;               /**< double defining the proportion of the
                                             cross-section occupied by fibrosis. This
                                             contributes extracellular passive tension */

    double prop_myofilaments;           /**< double defining the proportion of the
                                             non-fibrotic cross-section occupied by
                                             myofilaments. This contributes titin and
                                             cross-bridged mediated force */

    double m_filament_density;          /**< double defining the number of thick filaments
                                             per square meter of cross-section */
        
    // Filaments

    // Thick structure
    int m_n;                            /**< Number of thick filaments per half-sarcomere */

    int m_crowns_per_filament;          /**< Number of crowns per thick filament */
    int m_hubs_per_crown;               /**< Number of start-points per crown */
    int m_myosins_per_hub;              /**< Number of myosins per hub */

    double m_inter_crown_rest_length;   /**< double defining the rest-length of the spring
                                             between crowns in nm */
    double m_lambda;                    /**< double defining the barezone length in nm */

    double m_starting_angle;            /**< double defining the starting angle for the first crown */

    double m_inter_crown_twist;         /**< double defining the inter-crown twist between
                                             crowns  in degrees*/
    double m_within_hub_twist;          /**< double defining the twist between myosins
                                             in the same hub */

    // Thick parameters

    kinetic_scheme* p_m_scheme[MAX_NO_OF_ISOTYPES];
                                        /**< pointer to a kinetic scheme array for myosin */
    
    int m_no_of_cb_states;              /**< integer defining the number of states a
                                             myosin head can transition between */

    double m_k_stiff;                   /**< double defining the stiffness of a
                                             myosin filament spring in N m^-1 */

    int m_no_of_isotypes;               /**< Number of myosin isotypes */

    gsl_vector* m_isotype_props;	    /**< gsl_vector holding the myosin isotype proportions */

    // Myosin parameters

    double m_k_cb;                      /**< double defining the stiffness of a
                                             cross-bridge link in N m^-1 */

    // Thin structure
    int a_regulatory_units_per_strand;
                                        /**< integer defining the number of regulatory
                                             units along a single trand */
    int a_bs_per_unit;                  /**< integer defining the number of binding sites
                                             per strand */
    int a_strands_per_filament;         /**< integer defining the number of strands
                                             per filament */
    int a_no_of_bs;                     /**< integer defining the total number of binding
                                             sites on a thin filament */

    double a_inter_bs_rest_length;      /**< double defining the rest-length of the
                                             spring between binding-sites in nm */

    double a_inter_bs_twist;            /**< double defining the inter-bs twist
                                             along the thin filament in degrees */

    // Thin parameters

    int a_no_of_bs_states;              /**< integer defining the number of states that
                                             a binding site can tranisiton between */

    double a_k_stiff;                   /**< double defining the stiffness of an
                                             actin spring in N m^-1 */

    double a_k_on;                      /**< double defining second order rate constant
                                             for binding site activation in units of
                                             s^-1 M^-1 of Ca2+ */

    double a_k_off;                      /**< double defining rate constant
                                              for binding site de-activation in units
                                              of s^-1 */

    double a_k_coop;                     /**< double defining cooperativity
                                              dimensionless */

    // Titin structure
    int t_attach_a_node;                /**< int defining the thin node at which titin
                                             attaches */

    int t_attach_m_node;                /**< int defining the thick node at which titin
                                             attaches */


    // Titin parameters
    char t_passive_mode[_MAX_PATH];     /**< char array defining the passive_mode */

    double t_k_stiff;                   /**< double definining the stiffness of titin
                                             in N m^-1 */

    double t_slack_length;              /**< double defining the slack length of titin
                                             in nm */

    // Extracellular parameters
    char e_passive_mode[_MAX_PATH];     /**< char array defining the extracellular
                                             passive mode */

    double e_sigma;                     /**< double defining passive scaling factor */

    double e_L;                         /**< double defining passive curvature */

    double e_slack_length;              /**< double defining slack length of extracellular
                                             component */

    // MyBPC structure

    int c_thick_proximal_node;          /**< integer defining the first node for MyBPC */

    int c_thick_stripes;                /**< integer defining the number of MyBPC stripes */

    int c_thick_node_spacing;           /**< integer defining the node spacing for MyBPC */

    int c_mols_per_node;                /**< integer defining the number of MyBPC per node */

    double c_starting_angle;            /**< double defining the angle of the first MyBPC
                                             on each crown */

    double c_inter_stripe_twist;        /**< double defining the inter-stripe twist between
                                             MyBPC stripes in degrees */

    double c_k_stiff;                   /**< double defining the stiffness of a MyBPC link */

    int c_no_of_isotypes;               /**< Number of C-protein isotypes */

    gsl_vector* c_isotype_props;	    /**< gsl_vector holding the C-protein isotypes proportions */

    kinetic_scheme* p_c_scheme[MAX_NO_OF_ISOTYPES];        /**< pointer to a kinetic scheme array for MyBPC */

    // Functions

    /**
     * a function that writes FiberSim properties to file
     */
    void write_FiberSim_model_to_file(void);

    /**
     * a function that sets FiberSim_model parameters from a JSON file
     * @param json_file_string a char[] holding the filename for the JSON file
     */
    void set_FiberSim_model_parameters_from_JSON_file_string(char JSON_file_string[]);

    /**
    * a function that creates a kinetic scheme
    * @param
    */
    kinetic_scheme* create_kinetic_scheme(const rapidjson::Value& ks);
};