/**
 * @file    FiberSim_options.cpp
 * @brief   Source file for the FiberSim_options class
 * @author  Ken Campbell
 */

#include <cstdio>
#include <filesystem>

#include "FiberSim_options.h"
#include "JSON_functions.h"

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"

using namespace std::filesystem;

// Constructor
FiberSim_options::FiberSim_options(char JSON_options_file_string[])
{
    // Initialise

    printf("In FiberSim_options constructor\n");
    printf("JSON_options_file_string: %s\n", JSON_options_file_string);

    // Set default values
    x_pos_rel_tol = 0.001;                  /**< default value of relative tolerance for
                                                 calculating x positions */

    max_rate = 5000.0;                      /**< default value for maximum rate allowed
                                                 in calculations */
    
    set_FiberSim_options_from_JSON_file_string(JSON_options_file_string);

    // Do some processing on the options
    if (strlen(log_folder) > 0)
    {
        // Set the log mode
        log_mode = 1;
        
        // Clear log direcory
        printf("Clearing log folder\n");

        if (is_directory(log_folder))
        {
            // This removes the directory as well as the files
            std::uintmax_t n = remove_all(log_folder);
            printf("Deleted %i file(s) from %s\n", n, log_folder);
        }

        // Now create the log directory
        if (create_directory(log_folder))
        {
            printf("Log folder created at: %s\n", log_folder);
        }
        else
        {
            printf("Log folder could not be created: %s\n", log_folder);
            exit(1);
        }

        // Create the log file
        sprintf_s(log_file_string, _MAX_PATH, "%s/log_file.log", log_folder);
        errno_t err = fopen_s(&log_file, log_file_string, "w");
        if (err != 0)
        {
            printf("log file: %s\ncould not be opened\n", log_file_string);
            exit(1);
        }
        else
        {
            printf("log file opened: %s\n", log_file_string);
            fprintf_s(log_file, "Log file opened\n");
        }

        // Check whether we need another folder for the hs_status files
        if (dump_hs_status)
        {
            sprintf_s(hs_status_folder, _MAX_PATH, "%s/%s", log_folder, "hs_status");
            if (create_directory(hs_status_folder))
            {
                printf("HS status folder created at :%s\n", hs_status_folder);
            }
            else
            {
                printf("HS status folder could not be created: %s\n", hs_status_folder);
                exit(1);
            }
        }
        else
        {
            // Set the folder as an emptry string
            sprintf_s(hs_status_folder, _MAX_PATH, "");
        }

        //write_FiberSim_options_to_file();
    }
    else
    {
        log_mode = 0;
    }
}

// Destructor
FiberSim_options::~FiberSim_options()
{
    // Tidy up
    printf("In FiberSim_options destructor\n");

    if (log_mode > 0)
    {
        // Close log file
        fclose(log_file);
    }
}

void FiberSim_options::set_FiberSim_options_from_JSON_file_string(char JSON_file_string[])
{
    printf("In set_FiberSim_options_from_JSON_file_string\n");

    errno_t file_error;

    FILE *fp;
    file_error = fopen_s(&fp, JSON_file_string, "rb");
    if (file_error != 0)
    {
        printf("Error opening JSON options file: %s\n", JSON_file_string);
        exit(1);
    }

    char readBuffer[65536];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

    rapidjson::Document doc;
    doc.ParseStream(is);

    fclose(fp);
    printf("Finished reading stream\n");

    // Check we have options
    JSON_functions::check_JSON_member_object(doc, "options");
    const rapidjson::Value& options = doc["options"];

    // Now check for the log folder
    if (JSON_functions::is_JSON_member(options, "log_folder"))
    {
        JSON_functions::check_JSON_member_string(options, "log_folder");
        sprintf_s(log_folder, _MAX_PATH, "%s", options["log_folder"].GetString());

        // If we have a log folder, we need to check the dump_hs_status
        if (JSON_functions::is_JSON_member(options, "dump_hs_status"))
        {
            JSON_functions::check_JSON_member_int(options, "dump_hs_status");
            dump_hs_status = options["dump_hs_status"].GetInt();
        }
    }
    else
    {
        log_mode = 0;
        sprintf_s(log_folder, _MAX_PATH, "");
        dump_hs_status = 0;
    }


    // Check we have entries and set them

    if (JSON_functions::is_JSON_member(options, "x_pos_rel_tol"))
    {
        JSON_functions::check_JSON_member_number(options, "x_pos_rel_tol");
        x_pos_rel_tol = options["x_pos_rel_tol"].GetDouble();
    }


    if (JSON_functions::is_JSON_member(options, "max_rate"))
    {
        JSON_functions::check_JSON_member_number(options, "max_rate");
        max_rate = options["max_rate"].GetDouble();
    }


    printf("log folder: %s\n", log_folder);

    // Set status mode
    if (strcmp(log_folder, "none") == 0)
    {
        log_mode = 0;
    }
    else
    {
        log_mode = 1;
    }

}

void FiberSim_options::write_FiberSim_options_to_file(void)
{
    // Code writes FiberSim_options to file

    // Variables
    char output_file_string[_MAX_PATH];
    FILE* output_file;

    // Code
    sprintf_s(output_file_string, _MAX_PATH, "%s\\%s",
        log_folder, "FiberSim_options.log");

    errno_t err = fopen_s(&output_file, output_file_string, "w");
    if (err != 0)
    {
        printf("Options log file file: %s\ncould not be opened\n",
            output_file_string);
        exit(1);
    }

    fprintf_s(output_file, "log_folder: %s\n", log_folder);
    fprintf_s(output_file, "log_file_string: %s\n", log_file_string);
    fprintf_s(output_file, "log_mode: %i\n", log_mode);
    fprintf_s(output_file, "dump_hs_status: %i\n", dump_hs_status);
    fprintf_s(output_file, "hs_status_folder: %s\n", hs_status_folder);
    fprintf_s(output_file, "no_of_repeats: %i\n", no_of_repeats);
    fprintf_s(output_file, "multithreading: %i\n", multithreading);
    fprintf_s(output_file, "max_rate: %g\n", max_rate);

    // Tidy up
    fclose(output_file);
}
