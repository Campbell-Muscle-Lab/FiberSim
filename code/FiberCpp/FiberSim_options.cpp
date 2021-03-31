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

namespace fs = std::filesystem;

// Constructor
FiberSim_options::FiberSim_options(char JSON_options_file_string[])
{
    // Initialise

    // Set default values
    x_pos_rel_tol = 0.001;                  /**< default value of relative tolerance for
                                                 calculating x positions */

    max_rate = 5000.0;                      /**< default value for maximum rate allowed
                                                 in calculations */

    dump_precision = 8;                     /**< default value for dump precision */


    // Update values from log file
    set_FiberSim_options_from_JSON_file_string(JSON_options_file_string);

    // Do some processing on the options
    if (strlen(status_folder) > 0)
    {
        if (!strcmp(status_relative_to, "this_file"))
        {
            fs::path options_file = JSON_options_file_string;
            fs::path options_path = options_file.parent_path();
            fs::path status_path = options_path / status_folder;

            // Make sure the status folder exists
            if (fs::is_directory(status_path))
            {
                // Clean the directory
                int n = (int) fs::remove_all(status_path);
                printf("Deleting %i files from status_folder: %s",
                    n, status_path.string().c_str());
            }

            // Now create the directory
            if (fs::create_directory(status_path))
            {
                printf("Status folder created at: %s\n", status_path.string().c_str());
            }
            else
            {
                printf("Status folder could not be created: %s\n", status_path.string().c_str());
                exit(1);
            }

            // Set the status folder
            sprintf_s(status_folder, _MAX_PATH, "%s", status_path.string().c_str());
        }

        // Parse the time_steps string
        std::string ts_string = time_steps_string;

        size_t first_sep = ts_string.find_first_of(":");
        size_t last_sep = ts_string.find_last_of(":");

        start_status_time_step = (int)std::stoi(ts_string.substr(0, first_sep));
        skip_status_time_step = (int)std::stoi(ts_string.substr((first_sep+1), last_sep));
        stop_status_time_step = (int)std::stoi(ts_string.substr(last_sep+1));
    }

/*

        }
        // Set the log mode
        log_mode = 1;
        
        // Clear log direcory
        printf("Clearing log folder\n");

        if (is_directory(log_folder))
        {
            // This removes the directory as well as the files
            std::uintmax_t n = remove_all(log_folder);
            printf("Deleted %i file(s) from %s\n", (int)n, log_folder);
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
*/
}

// Destructor
FiberSim_options::~FiberSim_options()
{
    // Tidy up

    if (log_mode > 0)
    {
        // Close log file
        fclose(log_file);
    }
}

void FiberSim_options::set_FiberSim_options_from_JSON_file_string(char JSON_file_string[])
{
    // Variables
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

    // Check we have options
    JSON_functions::check_JSON_member_object(doc, "options");
    const rapidjson::Value& options = doc["options"];

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

    // Check if the dump precision was specified.
    if (JSON_functions::is_JSON_member(options, "dump_precision"))
    {
        JSON_functions::check_JSON_member_int(options, "dump_precision");
        dump_precision = options["dump_precision"].GetInt();
    }

    // Now check for logging
    if (JSON_functions::is_JSON_member(options, "logging"))
    {
        const rapidjson::Value& logging = options["logging"];

        JSON_functions::check_JSON_member_string(logging, "relative_to");
        sprintf_s(log_relative_to, _MAX_PATH, "%s", logging["relative_to"].GetString());

        JSON_functions::check_JSON_member_string(logging, "log_folder");
        sprintf_s(log_folder, _MAX_PATH, "%s", logging["log_folder"].GetString());
    }

    // Now check for status files
    if (JSON_functions::is_JSON_member(options, "status_files"))
    {
        const rapidjson::Value& status_files = options["status_files"];

        JSON_functions::check_JSON_member_string(status_files, "relative_to");
        sprintf_s(status_relative_to, _MAX_PATH, "%s", status_files["relative_to"].GetString());

        JSON_functions::check_JSON_member_string(status_files, "status_folder");
        sprintf_s(status_folder, _MAX_PATH, "%s", status_files["status_folder"].GetString());

        JSON_functions::check_JSON_member_string(status_files, "time_steps");
        sprintf_s(time_steps_string, _MAX_PATH, "%s", status_files["time_steps"].GetString());
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
    fprintf_s(output_file, "max_rate: %g\n", max_rate);
    fprintf_s(output_file, "dump_precision: %i\n", dump_precision);

    // Tidy up
    fclose(output_file);
}
