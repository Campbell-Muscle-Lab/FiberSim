/**
 * @file    FiberSim_model.cpp
 * @brief   Source file for the FiberSim_model class
 * @author  Ken Campbell
 */

#include <cstdio>
#include <string>

#include "FiberSim_model.h"
#include "FiberSim_options.h"

#include "kinetic_scheme.h"
#include "m_state.h"
#include "transition.h"
#include "model_hs_variation.h"

#include "iso_scheme.h"
#include "iso_type.h"
#include "iso_transition.h"

#include "gsl_math.h"

#include "JSON_functions.h"

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"

#include "global_definitions.h"

// Structure for thin_kinetics
struct a_kinetics
{
    double a_k_on;
    double a_k_off;
    double a_k_coop;
    double a_k_on_titin_effect;
    double a_k_coop_titin_effect;
};

// Constructor
FiberSim_model::FiberSim_model(char JSON_model_file_string[],
    FiberSim_options * set_p_fs_options)
{
    // Initialise

    // Set the pointer
    p_fs_options = set_p_fs_options;

    // And null some vectors that might not be needed
    m_isotype_ints = NULL;
    c_isotype_ints = NULL;

    // Zero the number of a_isotypes
    a_no_of_bs_iso_types = 0;

    // Allocate vectors for inter-hs force weights
    inter_hs_t_force_effects = gsl_vector_alloc(MAX_NO_OF_RATE_PARAMETERS);

    // Set to NaN
    gsl_vector_set_all(inter_hs_t_force_effects, GSL_NAN);

    // Null some pointers that might not be used
    p_m_iso_scheme = NULL;
    p_c_iso_scheme = NULL;
    p_thin_iso_scheme = NULL;

    // Log
    if (p_fs_options->log_mode > 0)
    {
        fprintf_s(p_fs_options->log_file, "In FiberSim_model constructor\n");
        fprintf_s(p_fs_options->log_file, "JSON_model_file_string: %s\n", JSON_model_file_string);
    }

    set_FiberSim_model_parameters_from_JSON_file_string(JSON_model_file_string);

    // Check whether we need to adjust m_filament_density
    if (!gsl_isnan(p_fs_options->m_filament_density_ref_hs_length))
    {
        m_filament_density = m_filament_density * (initial_hs_length /
            p_fs_options->m_filament_density_ref_hs_length);
    }

    if (strlen(p_fs_options->log_folder) > 0)
    {
        // Dumps model file
        //write_FiberSim_model_to_file();
        //sprintf_s(p_fs_options->log_folder, _MAX_PATH, "c:\\temp");

        // Dumps the kinetic scheme separarately
        char model_JSON_file_string[_MAX_PATH];
        sprintf_s(model_JSON_file_string, _MAX_PATH, "%s\\kinetic_scheme.json",
                    p_fs_options->log_folder);
        for (int i = 0; i < m_no_of_isotypes; i ++)
        {            p_m_scheme[i]->write_kinetic_scheme_to_file(model_JSON_file_string);
        }
        
    }
}

// Destructor
FiberSim_model::~FiberSim_model()
{
    // Tidy up
    if (p_fs_options->log_mode > 0)
    {
        fprintf_s(p_fs_options->log_file, "In FiberSim_model destructor\n");
    }

    // Delete thick filaments
    for (int i = 0; i < m_no_of_isotypes; i++)
    {
        delete p_m_scheme[i];
    }

    for (int i = 0; i < c_no_of_isotypes; i++)
    {
        delete p_c_scheme[i];
    }

    for (int i = 0; i < no_of_model_hs_variations; i++)
    {
        delete p_model_hs_variation[i];
    }

    // Delete gsl_vector
    gsl_vector_free(m_isotype_props);
    gsl_vector_free(c_isotype_props);

    gsl_vector_free(inter_hs_t_force_effects);

    // Delete a_kinetics
    for (int i = 0; i < a_no_of_bs_iso_types; i++)
    {
        delete p_a_kinetics[i];
    }

    // Delete arrays if necessary
    if (m_isotype_ints != NULL)
        gsl_vector_short_free(m_isotype_ints);

    if (c_isotype_ints != NULL)
        gsl_vector_short_free(c_isotype_ints);

    // Delete iso_schemes if necessary
    if (p_m_iso_scheme != NULL)
        delete p_m_iso_scheme;
    if (p_c_iso_scheme != NULL)
        delete p_c_iso_scheme;
    if (p_thin_iso_scheme != NULL)
        delete p_thin_iso_scheme;
}

// Functions

void FiberSim_model::set_FiberSim_model_parameters_from_JSON_file_string(char JSON_file_string[])
{
    // Load the model from the JSON_file_string
    
    if (p_fs_options->log_mode > 0)
    {
        fprintf_s(p_fs_options->log_file, "In set_FiberSim_model_parameters_from_JSON_file_string\n");
    }

    errno_t file_error;

    FILE *fp;
    file_error = fopen_s(&fp, JSON_file_string, "rb");
    if (file_error != 0)
    {
        printf("Error opening JSON model file: %s\n", JSON_file_string);
        exit(1);
    }

    char readBuffer[65536];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

    rapidjson::Document doc;
    doc.ParseStream(is);

    fclose(fp);

    if (p_fs_options->log_mode > 0)
    {
        fprintf_s(p_fs_options->log_file, "Finished reading JSON model stream\n");
        fprintf_s(p_fs_options->log_file, "Now setting model data\n");
    }

    // Load the FiberSim version number

    JSON_functions::check_JSON_member_object(doc, "FiberSim");
    const rapidjson::Value& fs = doc["FiberSim"];

    JSON_functions::check_JSON_member_string(fs, "version");
    sprintf_s(version, _MAX_PATH, fs["version"].GetString());

    // Load the muscle parameters
    JSON_functions::check_JSON_member_object(doc, "muscle");
    const rapidjson::Value& mus = doc["muscle"];

    JSON_functions::check_JSON_member_int(mus, "no_of_half_sarcomeres");
    no_of_half_sarcomeres = mus["no_of_half_sarcomeres"].GetInt();

    JSON_functions::check_JSON_member_int(mus, "no_of_myofibrils");
    no_of_myofibrils = mus["no_of_myofibrils"].GetInt();

    JSON_functions::check_JSON_member_number(mus, "initial_hs_length");
    initial_hs_length = mus["initial_hs_length"].GetDouble();

    JSON_functions::check_JSON_member_number(mus, "prop_fibrosis");
    prop_fibrosis = mus["prop_fibrosis"].GetDouble();

    JSON_functions::check_JSON_member_number(mus, "prop_myofilaments");
    prop_myofilaments = mus["prop_myofilaments"].GetDouble();

    JSON_functions::check_JSON_member_number(mus, "m_filament_density");
    m_filament_density = mus["m_filament_density"].GetDouble();

    // Check if the series elastic component is specified
    if (JSON_functions::check_JSON_member_exists(mus, "sc_k_stiff"))
    {
        // Is it a number
        if (JSON_functions::valid_JSON_member_number(mus, "sc_k_stiff"))
        {
            sc_k_stiff = mus["sc_k_stiff"].GetDouble();
        }
        else
        {
            sc_k_stiff = GSL_NAN;
        }
    }
    else
        sc_k_stiff = GSL_NAN;

    // Check if temperature is specified - This part of the code ensures compatibility with FiberSim V2.0.2

    if (JSON_functions::check_JSON_member_exists(mus, "temperature"))
        temperature = mus["temperature"].GetDouble();
    else
        temperature = 310;

    // Load the lattice_parameters

    // Check if lattice_parameters is specified - This part of the code ensures compatibility with FiberSim V2.0.2

    if (JSON_functions::check_JSON_member_exists(doc, "lattice_parameters"))
    { 
        const rapidjson::Value& lp = doc["lattice_parameters"];
        JSON_functions::check_JSON_member_number(lp, "viscosity");
        viscosity = lp["viscosity"].GetDouble();
    }
    else
    {
        viscosity = 0.0;
    }

    // Load the thin_structure variables
    JSON_functions::check_JSON_member_object(doc, "thin_structure");
    const rapidjson::Value& thin_structure = doc["thin_structure"];

    JSON_functions::check_JSON_member_int(thin_structure, "a_strands_per_filament");
    a_strands_per_filament = thin_structure["a_strands_per_filament"].GetInt();

    JSON_functions::check_JSON_member_int(thin_structure, "a_regulatory_units_per_strand");
    a_regulatory_units_per_strand = thin_structure["a_regulatory_units_per_strand"].GetInt();

    JSON_functions::check_JSON_member_int(thin_structure, "a_bs_per_unit");
    a_bs_per_unit = thin_structure["a_bs_per_unit"].GetInt();

    JSON_functions::check_JSON_member_number(thin_structure, "a_inter_bs_rest_length");
    a_inter_bs_rest_length = thin_structure["a_inter_bs_rest_length"].GetDouble();

    JSON_functions::check_JSON_member_number(thin_structure, "a_inter_bs_twist");
    a_inter_bs_twist = thin_structure["a_inter_bs_twist"].GetDouble();

    // Load the thick_structure variables
    JSON_functions::check_JSON_member_object(doc, "thick_structure");
    const rapidjson::Value& thick_structure = doc["thick_structure"];

    JSON_functions::check_JSON_member_int(thick_structure, "m_n");
    m_n = thick_structure["m_n"].GetInt();

    JSON_functions::check_JSON_member_int(thick_structure, "m_crowns_per_filament");
    m_crowns_per_filament = thick_structure["m_crowns_per_filament"].GetInt();

    JSON_functions::check_JSON_member_int(thick_structure, "m_hubs_per_crown");
    m_hubs_per_crown = thick_structure["m_hubs_per_crown"].GetInt();

    JSON_functions::check_JSON_member_int(thick_structure, "m_myosins_per_hub");
    m_myosins_per_hub = thick_structure["m_myosins_per_hub"].GetInt();

    JSON_functions::check_JSON_member_number(thick_structure, "m_inter_crown_rest_length");
    m_inter_crown_rest_length = thick_structure["m_inter_crown_rest_length"].GetDouble();

    JSON_functions::check_JSON_member_number(thick_structure, "m_lambda");
    m_lambda = thick_structure["m_lambda"].GetDouble();

    JSON_functions::check_JSON_member_number(thick_structure, "m_starting_angle");
    m_starting_angle = thick_structure["m_starting_angle"].GetDouble();

    JSON_functions::check_JSON_member_number(thick_structure, "m_inter_crown_twist");
    m_inter_crown_twist = thick_structure["m_inter_crown_twist"].GetDouble();

    JSON_functions::check_JSON_member_number(thick_structure, "m_within_hub_twist");
    m_within_hub_twist = thick_structure["m_within_hub_twist"].GetDouble();

    // Load the titin_structure variables
    JSON_functions::check_JSON_member_object(doc, "titin_structure");
    const rapidjson::Value& titin_structure = doc["titin_structure"];

    JSON_functions::check_JSON_member_int(titin_structure, "t_attach_a_node");
    t_attach_a_node = titin_structure["t_attach_a_node"].GetInt();
    
    JSON_functions::check_JSON_member_int(titin_structure, "t_attach_m_node");
    t_attach_m_node = titin_structure["t_attach_m_node"].GetInt();

    // Load the MyBPC structure variables
    JSON_functions::check_JSON_member_object(doc, "mybpc_structure");
    const rapidjson::Value& mybpc_structure = doc["mybpc_structure"];

    JSON_functions::check_JSON_member_int(mybpc_structure, "c_thick_proximal_node");
    c_thick_proximal_node = mybpc_structure["c_thick_proximal_node"].GetInt();

    JSON_functions::check_JSON_member_int(mybpc_structure, "c_thick_stripes");
    c_thick_stripes = mybpc_structure["c_thick_stripes"].GetInt();

    JSON_functions::check_JSON_member_int(mybpc_structure, "c_thick_node_spacing");
    c_thick_node_spacing = mybpc_structure["c_thick_node_spacing"].GetInt();

    JSON_functions::check_JSON_member_int(mybpc_structure, "c_mols_per_node");
    c_mols_per_node = mybpc_structure["c_mols_per_node"].GetInt();

    JSON_functions::check_JSON_member_number(mybpc_structure, "c_starting_angle");
    c_starting_angle = mybpc_structure["c_starting_angle"].GetDouble();

    // Check if c_inter_stripe_twist is specified
    // This part of the code ensures compatibility with FiberSim V1.1.2
    if (JSON_functions::check_JSON_member_exists(mybpc_structure, "c_inter_stripe_twist"))
        c_inter_stripe_twist = mybpc_structure["c_inter_stripe_twist"].GetDouble();
    else
        c_inter_stripe_twist = 10.0;

    // Check if c_sd_angle_deviation is specified, set to 0 if not
    if (JSON_functions::check_JSON_member_exists(mybpc_structure, "c_sd_angle_deviation"))
        c_sd_angle_deviation = mybpc_structure["c_sd_angle_deviation"].GetDouble();
    else
        c_sd_angle_deviation = 0.0;

    // Load the thin_parameters
    JSON_functions::check_JSON_member_object(doc, "thin_parameters");
    const rapidjson::Value& thin_parameters = doc["thin_parameters"];

    JSON_functions::check_JSON_member_number(thin_parameters, "a_k_stiff");
    a_k_stiff = thin_parameters["a_k_stiff"].GetDouble();

    JSON_functions::check_JSON_member_number(thin_parameters, "a_no_of_bs_states");
    a_no_of_bs_states = thin_parameters["a_no_of_bs_states"].GetInt();

    // Check whether thin_kinetics is specified
    if (JSON_functions::check_JSON_member_exists(doc, "thin_kinetics"))
    {
        // It does, create arrays of structures
        JSON_functions::check_JSON_member_array(doc, "thin_kinetics");
        const rapidjson::Value& thin_kin = doc["thin_kinetics"].GetArray();

        a_no_of_bs_iso_types = thin_kin.Size();
        for (rapidjson::SizeType i = 0; i < thin_kin.Size(); i++)
        {
            p_a_kinetics[i] = new a_kinetics();
            p_a_kinetics[i]->a_k_on = thin_kin[i]["a_k_on"].GetDouble();
            p_a_kinetics[i]->a_k_off = thin_kin[i]["a_k_off"].GetDouble();
            p_a_kinetics[i]->a_k_coop = thin_kin[i]["a_k_coop"].GetDouble();

            if (JSON_functions::check_JSON_member_exists(thin_kin[i], "a_k_on_titin_effect"))
            {
                p_a_kinetics[i]->a_k_on_titin_effect = thin_kin[i]["a_k_on_titin_effect"].GetDouble();
            }
            else
            {
                p_a_kinetics[i]->a_k_on_titin_effect = 0.0;
            }

            if (JSON_functions::check_JSON_member_exists(thin_kin[i], "a_k_coop_titin_effect"))
            {
                p_a_kinetics[i]->a_k_coop_titin_effect = thin_kin[i]["a_k_coop_titin_effect"].GetDouble();
            }
            else
            {
                p_a_kinetics[i]->a_k_coop_titin_effect = 0.0;
            }

        }
    }
    else
    {
        // Thin kinetics is not specified. Create a single struct and fill it
        // from thin_parameters
        a_no_of_bs_iso_types = 1;
        p_a_kinetics[0] = new a_kinetics();
        p_a_kinetics[0]->a_k_on = thin_parameters["a_k_on"].GetDouble();
        p_a_kinetics[0]->a_k_off = thin_parameters["a_k_off"].GetDouble();
        p_a_kinetics[0]->a_k_coop = thin_parameters["a_k_coop"].GetDouble();
        p_a_kinetics[0]->a_k_on_titin_effect = 0.0;
        p_a_kinetics[0]->a_k_coop_titin_effect = 0.0;
    }

    // Load the thick_parameters
    JSON_functions::check_JSON_member_object(doc, "thick_parameters");
    const rapidjson::Value& thick_parameters = doc["thick_parameters"];

    JSON_functions::check_JSON_member_number(thick_parameters, "m_k_stiff");
    m_k_stiff = thick_parameters["m_k_stiff"].GetDouble();

    // Load the titin_parameters
    JSON_functions::check_JSON_member_object(doc, "titin_parameters");
    const rapidjson::Value& titin_parameters = doc["titin_parameters"];

    //JSON_functions::check_JSON_member_string(titin_parameters, "t_passive_mode");
    //sprintf_s(t_passive_mode, _MAX_PATH, titin_parameters["t_passive_mode"].GetString());

    check_and_assign_double(titin_parameters, "t_k_stiff", &t_k_stiff, 0.0);
    check_and_assign_double(titin_parameters, "t_sigma", &t_sigma, 0.0);
    check_and_assign_double(titin_parameters, "t_L", &t_L, 1e6);
    check_and_assign_double(titin_parameters, "t_offset", &t_offset, 0.0);

    // Load the extracellular_parameters
    JSON_functions::check_JSON_member_object(doc, "extracellular_parameters");
    const rapidjson::Value& extracellular_parameters = doc["extracellular_parameters"];

    check_and_assign_double(extracellular_parameters, "e_slack_length", &e_slack_length, 1000.0);
    check_and_assign_double(extracellular_parameters, "e_sigma", &e_sigma, 0.0);
    check_and_assign_double(extracellular_parameters, "e_L", &e_L, 1e6);
    check_and_assign_double(extracellular_parameters, "e_k_stiff", &e_k_stiff, 0.0);

    // Inter_half-sarcomere effects
    if (JSON_functions::is_JSON_member(doc, "inter_half_sarcomere_parameters"))
    {
        const rapidjson::Value& ihs = doc["inter_half_sarcomere_parameters"];

        if (JSON_functions::is_JSON_member(ihs, "titin_force_effects"))
        {
            const rapidjson::Value& tfe = ihs["titin_force_effects"];

            for (int i = 0; i < (int)tfe.Size(); i++)
            {
                gsl_vector_set(inter_hs_t_force_effects, i, tfe[i].GetDouble());
            }
        }
    }

    // Load the m_parameters
    JSON_functions::check_JSON_member_object(doc, "m_parameters");
    const rapidjson::Value& m_parameters = doc["m_parameters"];

    JSON_functions::check_JSON_member_number(m_parameters, "m_k_cb");
    m_k_cb = m_parameters["m_k_cb"].GetDouble();

    if (JSON_functions::is_JSON_member(m_parameters, "m_isotype_proportions"))
    {
        const rapidjson::Value& mip = m_parameters["m_isotype_proportions"];

        m_no_of_isotypes = mip.Size();
        m_isotype_props = gsl_vector_alloc(MAX_NO_OF_ISOTYPES);
        gsl_vector_set_zero(m_isotype_props);

        for (int i = 0; i < (int)mip.Size(); i++)
        {
            gsl_vector_set(m_isotype_props, i, mip[i].GetDouble());
        }
    }
    
    // If there is an array of isotypes, save it to the appropriate short vector
    if (JSON_functions::is_JSON_member(m_parameters, "m_isotype_ints"))
    {
        const rapidjson::Value& mi_ints = m_parameters["m_isotype_ints"];

        int model_cb_n = (int)mi_ints.Size();

        m_isotype_ints = gsl_vector_short_alloc(model_cb_n);
        gsl_vector_short_set_zero(m_isotype_ints);

        for (int i = 0; i < model_cb_n; i++)
        {
            gsl_vector_short_set(m_isotype_ints, i, (short)mi_ints[i].GetInt());
        }
    }

    // Kinetic scheme for myosin - this is complicated so it's done in a different file
    JSON_functions::check_JSON_member_array(doc, "m_kinetics");
    const rapidjson::Value& m_ks = doc["m_kinetics"].GetArray();

    // Set the kinetic scheme for each isotype
    m_no_of_cb_states = 0;
    for (rapidjson::SizeType i = 0; i < m_ks.Size(); i++)
    {
        p_m_scheme[i] = create_kinetic_scheme(m_ks[i]);
        m_no_of_cb_states = GSL_MAX(m_no_of_cb_states, p_m_scheme[i]->no_of_states);
    }

    m_no_of_isotypes = m_ks.Size();

    // Load the MyBPC parameters variables
    JSON_functions::check_JSON_member_object(doc, "mybpc_parameters");
    const rapidjson::Value& mybpc_parameters = doc["mybpc_parameters"];

    JSON_functions::check_JSON_member_number(mybpc_parameters, "c_k_stiff");
    c_k_stiff = mybpc_parameters["c_k_stiff"].GetDouble();

    if (JSON_functions::is_JSON_member(mybpc_parameters, "c_isotype_proportions"))
    {
        const rapidjson::Value& cip = mybpc_parameters["c_isotype_proportions"];

        c_no_of_isotypes = cip.Size();

        c_isotype_props = gsl_vector_alloc(MAX_NO_OF_ISOTYPES);
        gsl_vector_set_zero(c_isotype_props);

        for (int i = 0; i < (int)cip.Size(); i++)
        {
            gsl_vector_set(c_isotype_props, i, cip[i].GetDouble());
        }
    }

    // If there is an array of isotypes, save it to the appropriate short vector
    if (JSON_functions::is_JSON_member(mybpc_parameters, "c_isotype_ints"))
    {
        const rapidjson::Value& ci_ints = mybpc_parameters["c_isotype_ints"];

        int model_mybpc_n = (int)ci_ints.Size();

        c_isotype_ints = gsl_vector_short_alloc(model_mybpc_n);
        gsl_vector_short_set_zero(c_isotype_ints);

        for (int i = 0; i < model_mybpc_n; i++)
        {
            gsl_vector_short_set(c_isotype_ints, i, (short)ci_ints[i].GetInt());
        }
    }
    

    // Kinetic scheme for MyBPC
    JSON_functions::check_JSON_member_array(doc, "c_kinetics");
    const rapidjson::Value& c_ks = doc["c_kinetics"].GetArray();

    c_no_of_pc_states = 0;
    for (rapidjson::SizeType i = 0; i < c_ks.Size(); i++)
    {
        p_c_scheme[i] = create_kinetic_scheme(c_ks[i]);
        c_no_of_pc_states = GSL_MAX(c_no_of_pc_states, p_c_scheme[i]->no_of_states);
    }

    c_no_of_isotypes = c_ks.Size();

    // Try to load the half-sarcomere variation
    if (JSON_functions::check_JSON_member_exists(doc, "half_sarcomere_variation"))
    {
        JSON_functions::check_JSON_member_array(doc, "half_sarcomere_variation");
        const rapidjson::Value& hsv = doc["half_sarcomere_variation"].GetArray();

        no_of_model_hs_variations = hsv.Size();

        for (rapidjson::SizeType i = 0; i < hsv.Size(); i++)
        {
            p_model_hs_variation[i] = new model_hs_variation(this, hsv[i]);
        }
    }

    // Check for isotype switches
    if (JSON_functions::check_JSON_member_exists(doc, "m_iso_switching"))
    {
        const rapidjson::Value& m_iso = doc["m_iso_switching"];
        
        p_m_iso_scheme = new iso_scheme(m_iso, this, p_fs_options);
    }

    // Check for isotype switches
    if (JSON_functions::check_JSON_member_exists(doc, "c_iso_switching"))
    {
        const rapidjson::Value& c_iso = doc["c_iso_switching"];

        p_c_iso_scheme = new iso_scheme(c_iso, this, p_fs_options);
    }

    // Check for isotype switches
    if (JSON_functions::check_JSON_member_exists(doc, "thin_iso_switching"))
    {
        const rapidjson::Value& thin_iso = doc["thin_iso_switching"];

        p_thin_iso_scheme = new iso_scheme(thin_iso, this, p_fs_options);
    }

    if (p_fs_options->log_mode > 0)
    {
        fprintf_s(p_fs_options->log_file, "Finished setting model data\n");
    }
}

void FiberSim_model::check_and_assign_double(const rapidjson::Value& doc, string tag, double* p_double, double default_value)
{
    // Variables
    char tag_string[_MAX_PATH];

    // Code
    sprintf_s(tag_string, _MAX_PATH, "%s", tag.c_str());

    if (JSON_functions::check_JSON_member_exists(doc, tag_string))
    {
        *p_double = doc[tag_string].GetDouble();
    }
    else
    {
        *p_double = default_value;
    }

    printf("%s: %g\n", tag_string, *p_double);
}

kinetic_scheme* FiberSim_model::create_kinetic_scheme(const rapidjson::Value& ks)
{
    //! Loads kinetic scheme

    // Variables
    kinetic_scheme* p_scheme;

    // Create the kinetic scheme
    p_scheme = new kinetic_scheme(ks, this, p_fs_options);

    // Return the pointer
    return p_scheme;
}
