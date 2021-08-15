/**
 * @file    FiberSim_model.cpp
 * @brief   Source file for the FiberSim_model class
 * @author  Ken Campbell
 */

#include <cstdio>

#include "FiberSim_model.h"
#include "FiberSim_options.h"

#include "kinetic_scheme.h"
#include "m_state.h"
#include "transition.h"

#include "kinetic_scheme.h"
#include "JSON_functions.h"

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"

#include "global_definitions.h"

// Constructor
FiberSim_model::FiberSim_model(char JSON_model_file_string[],
    FiberSim_options * set_p_fs_options)
{
    // Initialise

    // Set the pointer
    p_fs_options = set_p_fs_options;

    // Set pointer to kinetic scheme to NULL
    // p_m_scheme[MAX_NO_OF_ISOTYPES] = {NULL};

    // Log
    if (p_fs_options->log_mode > 0)
    {
        fprintf_s(p_fs_options->log_file, "In FiberSim_model constructor\n");
        fprintf_s(p_fs_options->log_file, "JSON_model_file_string: %s\n", JSON_model_file_string);
    }

    set_FiberSim_model_parameters_from_JSON_file_string(JSON_model_file_string);

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

    // Delete gsl_vector
    gsl_vector_free(m_isotype_props);
    gsl_vector_free(c_isotype_props);
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

    // Check if c_inter_stripe_twist is specified - This part of the code ensures compatibility with FiberSim V1.1.2
    int c_inter_stripe_flag = JSON_functions::is_JSON_member(mybpc_structure, "c_inter_stripe_twist");

    if (c_inter_stripe_flag == 1) // User specified inter_stripe_twist
    {
        JSON_functions::check_JSON_member_number(mybpc_structure, "c_inter_stripe_twist");
        c_inter_stripe_twist = mybpc_structure["c_inter_stripe_twist"].GetDouble();
    }

    else if (c_inter_stripe_flag == 0) // User did not specified inter_stripe_twist, use default value
    {    
        c_inter_stripe_twist = 0.0;
    }

    // Load the MyBPC parameters variables
    JSON_functions::check_JSON_member_object(doc, "mybpc_parameters");
    const rapidjson::Value& mybpc_parameters = doc["mybpc_parameters"];

    JSON_functions::check_JSON_member_number(mybpc_parameters, "c_k_stiff");
    c_k_stiff = mybpc_parameters["c_k_stiff"].GetDouble();

    JSON_functions::check_JSON_member_array(mybpc_parameters, "c_isotype_proportions");
    const rapidjson::Value& cip = mybpc_parameters["c_isotype_proportions"];

    c_no_of_isotypes = cip.Size();

    c_isotype_props = gsl_vector_alloc(MAX_NO_OF_ISOTYPES);
    gsl_vector_set_zero(c_isotype_props);

    for (int i = 0; i < (int)cip.Size(); i++)
    {
        gsl_vector_set(c_isotype_props, i, cip[i].GetDouble());
    }

    // Load the thin_parameters
    JSON_functions::check_JSON_member_object(doc, "thin_parameters");
    const rapidjson::Value& thin_parameters = doc["thin_parameters"];

    JSON_functions::check_JSON_member_number(thin_parameters, "a_k_stiff");
    a_k_stiff = thin_parameters["a_k_stiff"].GetDouble();

    JSON_functions::check_JSON_member_number(thin_parameters, "a_no_of_bs_states");
    a_no_of_bs_states = thin_parameters["a_no_of_bs_states"].GetInt();

    JSON_functions::check_JSON_member_number(thin_parameters, "a_k_on");
    a_k_on = thin_parameters["a_k_on"].GetDouble();

    JSON_functions::check_JSON_member_number(thin_parameters, "a_k_off");
    a_k_off = thin_parameters["a_k_off"].GetDouble();

    JSON_functions::check_JSON_member_number(thin_parameters, "a_k_coop");
    a_k_coop = thin_parameters["a_k_coop"].GetDouble();

    // Load the thick_parameters
    JSON_functions::check_JSON_member_object(doc, "thick_parameters");
    const rapidjson::Value& thick_parameters = doc["thick_parameters"];

    JSON_functions::check_JSON_member_number(thick_parameters, "m_k_stiff");
    m_k_stiff = thick_parameters["m_k_stiff"].GetDouble();

    // Load the m_parameters
    JSON_functions::check_JSON_member_object(doc, "m_parameters");
    const rapidjson::Value& m_parameters = doc["m_parameters"];

    JSON_functions::check_JSON_member_number(m_parameters, "m_k_cb");
    m_k_cb = m_parameters["m_k_cb"].GetDouble();

    JSON_functions::check_JSON_member_array(m_parameters, "m_isotype_proportions");
    const rapidjson::Value& mip = m_parameters["m_isotype_proportions"];

    m_no_of_isotypes = mip.Size();
    m_isotype_props = gsl_vector_alloc(MAX_NO_OF_ISOTYPES);
    gsl_vector_set_zero(m_isotype_props);

    for (int i = 0; i < (int)mip.Size(); i++)
    {
        gsl_vector_set(m_isotype_props, i, mip[i].GetDouble());
    }

    // Load the titin_parameters
    JSON_functions::check_JSON_member_object(doc, "titin_parameters");
    const rapidjson::Value& titin_parameters = doc["titin_parameters"];

    JSON_functions::check_JSON_member_string(titin_parameters, "t_passive_mode");
    sprintf_s(t_passive_mode, _MAX_PATH, titin_parameters["t_passive_mode"].GetString());

    JSON_functions::check_JSON_member_number(titin_parameters, "t_k_stiff");
    t_k_stiff = titin_parameters["t_k_stiff"].GetDouble();

    JSON_functions::check_JSON_member_number(titin_parameters, "t_slack_length");
    t_slack_length = titin_parameters["t_slack_length"].GetDouble();

    // Load the extracellular_parameters
    JSON_functions::check_JSON_member_object(doc, "extracellular_parameters");
    const rapidjson::Value& extracellular_parameters = doc["extracellular_parameters"];

    JSON_functions::check_JSON_member_string(extracellular_parameters, "e_passive_mode");
    sprintf_s(e_passive_mode, _MAX_PATH, extracellular_parameters["e_passive_mode"].GetString());

    JSON_functions::check_JSON_member_number(extracellular_parameters, "e_sigma");
    e_sigma = extracellular_parameters["e_sigma"].GetDouble();

    JSON_functions::check_JSON_member_number(extracellular_parameters, "e_L");
    e_L = extracellular_parameters["e_L"].GetDouble();

    JSON_functions::check_JSON_member_number(extracellular_parameters, "e_slack_length");
    e_slack_length = extracellular_parameters["e_slack_length"].GetDouble();

    // Kinetic scheme for myosin - this is complicated so it's done in a different file
    JSON_functions::check_JSON_member_array(doc, "m_kinetics");
    const rapidjson::Value& m_ks = doc["m_kinetics"].GetArray();  

    for (rapidjson::SizeType i = 0; i < m_ks.Size(); i++)
    {
        p_m_scheme[i] = create_kinetic_scheme(m_ks[i]);
    }

    // Kinetic scheme for MyBPC
    JSON_functions::check_JSON_member_array(doc, "c_kinetics");
    const rapidjson::Value& c_ks = doc["c_kinetics"].GetArray();

    for (rapidjson::SizeType i = 0; i < c_ks.Size(); i++)
    {
        p_c_scheme[i] = create_kinetic_scheme(c_ks[i]);
    }

    if (p_fs_options->log_mode > 0)
    {
        fprintf_s(p_fs_options->log_file, "Finished setting model data\n");
    }
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

void FiberSim_model::write_FiberSim_model_to_file(void)
{
    // Code writes FiberSim_model parameters to file

    // Variables
    char output_file_string[_MAX_PATH];
    FILE* output_file;

    // Code
    sprintf_s(output_file_string, _MAX_PATH, "%s\\%s",
        p_fs_options->log_folder, "FiberSim_model.json");

    errno_t err = fopen_s(&output_file, output_file_string, "w");
    if (err != 0)
    {
        printf("Options log file file: %s\ncould not be opened\n",
            output_file_string);
        exit(1);
    }

    fprintf_s(output_file, "\"muscle\":{\n");
    fprintf_s(output_file, "\t\"no_of_half_sarcomeres\": %i,\n", no_of_half_sarcomeres);
    fprintf_s(output_file, "\t\"no_of_myofibrils\": %i},\n", no_of_myofibrils);
    fprintf_s(output_file, "\t\"initial_hs_length\": %g},\n", initial_hs_length);
    fprintf_s(output_file, "\t\"m_filament_density\": %g},\n", m_filament_density);
    
    fprintf_s(output_file, "\"thick_structure\":{\n");
    fprintf_s(output_file, "\t\"m_n\": %i,\n", m_n);
    fprintf_s(output_file, "\t\"m_crowns_per_filament\": %i,\n", m_crowns_per_filament);
    fprintf_s(output_file, "\t\"m_hubs_per_crown\": %i,\n", m_hubs_per_crown);
    fprintf_s(output_file, "\t\"m_myosins_per_hub\": %i,\n", m_myosins_per_hub);
    fprintf_s(output_file, "\t\"m_inter_crown_rest_length\": %g,\n", m_inter_crown_rest_length);
    fprintf_s(output_file, "\t\"m_lambda\": %g,\n", m_lambda);
    fprintf_s(output_file, "\t\"m_inter_crown_twist\": %g,\n", m_inter_crown_twist);
    fprintf_s(output_file, "\t\"m_within_hub_twist\": %g},\n", m_within_hub_twist);
    
    fprintf_s(output_file, "\"thin_structure\":{\n");
    fprintf_s(output_file, "\t\"a_regulatory_units_per_strand\": %i,\n", a_regulatory_units_per_strand);
    fprintf_s(output_file, "\t\"a_bs_per_unit\": %i,\n", a_bs_per_unit);
    fprintf_s(output_file, "\t\"a_strands_per_filament\": %i,\n", a_strands_per_filament);
    fprintf_s(output_file, "\t\"a_inter_bs_rest_length\": %g,\n", a_inter_bs_rest_length);
    fprintf_s(output_file, "\t\"a_inter_bs_twist\": %g},\n", a_inter_bs_twist);

    fprintf_s(output_file, "\"titin_structure\":{\n");
    fprintf_s(output_file, "\t\"t_attach_a_node\": %i,\n", t_attach_a_node);
    fprintf_s(output_file, "\t\"t_attach_m_node\": %i},\n", t_attach_m_node);

    fprintf_s(output_file, "\"mybpc_structure\":{\n");
    fprintf_s(output_file, "\t\"c_thick_proximal_node\": %i,\n", c_thick_proximal_node);
    fprintf_s(output_file, "\t\"c_thick_stripes\": %i,\n", c_thick_stripes);
    fprintf_s(output_file, "\t\"c_thick_node_spacing\": %i,\n", c_thick_node_spacing);
    fprintf_s(output_file, "\t\"c_mols_per_node\": %i},\n", c_mols_per_node);

    fprintf_s(output_file, "\"mybpc_parameters\":{\n");
    fprintf_s(output_file, "\t\"c_k_stiff\": %g},\n", c_k_stiff);
    
    fprintf_s(output_file, "\"thick_parameters\":{\n");
    fprintf_s(output_file, "\t\"m_k_stiff\": %g},\n", m_k_stiff);

    fprintf_s(output_file, "\"m_parameters\":{\n");
    fprintf_s(output_file, "\t\"m_k_cb\": %g},\n", m_k_cb);

    fprintf_s(output_file, "\"thin_parameters\":{\n");
    fprintf_s(output_file, "\t\"a_no_of_bs_states\": %i,", a_no_of_bs_states);
    fprintf_s(output_file, "\t\"a_k_stiff\": %g,\n", a_k_stiff);
    fprintf_s(output_file, "\t\"a_k_on\": %g,\n", a_k_on);
    fprintf_s(output_file, "\t\"a_k_off\": %g,\n", a_k_off);
    fprintf_s(output_file, "\t\"a_k_coop\": %g},\n", a_k_coop);

    fprintf_s(output_file, "\"titin_parameters\":{\n");
    fprintf_s(output_file, "\t\"t_k_stiff\": %g,\n", t_k_stiff);
    fprintf_s(output_file, "\t\"t_slack_length\": %g,\n", t_slack_length);
    
    // Tidy up
    fclose(output_file);
}

