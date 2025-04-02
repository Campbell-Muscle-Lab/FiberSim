/**
 * @file    FiberSim_options.cpp
 * @brief   Source file for the FiberSim_options class
 * @author  Ken Campbell
 */

#include <cstdio>
#include <filesystem>
#include <string>

#include "FiberSim_options.h"
#include "JSON_functions.h"

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"

#include "gsl_math.h"

namespace fs = std::filesystem;

// Constructor
FiberSim_options::FiberSim_options(char JSON_options_file_string[])
{
    // Initialise

    // Set default values
    x_pos_rel_tol = 0.001;                  /**< default value of absolute tolerance for
                                                 calculating x positions */

    x_vector_max_iterations = 100;          /**< default value of max iterations for
                                                 x vector calculation */

    hs_force_control_max_delta_hs_length = 100;
                                            /**< default value of max_delta_hs_length for
                                                 force-control for a lattice */

    adjacent_bs = 0;                        /**< default value of adjacent binding sites
                                                 available for myosin or mybpc attachment
                                                 0 restricts to nearest site */

    max_rate = 1000.0;                      /**< default value for maximum rate allowed
                                                 in calculations */

    thin_filament_sub_steps = 10;            /**< default value for number of sub_steps
                                                 used for thin filament kinetics */

    lambda_jitter = 0.0;                    /**< default value for lambda jitter */

    m_filament_density_ref_hs_length = GSL_NAN;
                                            /**< default value for scaling length
                                                 GSL_NAN implies no adjustment */

    dump_precision = 6;                     /**< default value for dump precision */

    // Set default values for strings
    sprintf_s(rate_file_string, _MAX_PATH, "");
    sprintf_s(status_folder, _MAX_PATH, "");
    sprintf_s(log_folder, _MAX_PATH, "");
    sprintf_s(time_steps_string, _MAX_PATH, "");

    myofibril_force_tolerance = 0.001;      /**< default value for the force tolerance
                                                 for myofibril multiroot calculations */

    myofibril_max_iterations = 100;         /**< default value for the maximum number of
                                                 iterations in myofibril multiroot
                                                 calculations */

    min_hs_length = 500.0;                  /**< default value for the half-sarcomere length
                                                 at which the simulation collapses and gives
                                                 up trying to keep the calculations going */

    no_of_worker_threads = 0;               /**< default value for number of worker threads */

    calculate_x_mode = 1;                   /**< default value for calculate_x_mode */

    afterload_load = GSL_NAN;               /**< default value for afterload load */
    afterload_break_delta_hs_length = GSL_NAN;
                                            /**< default value for breakout length */
    afterload_post_break_wait_s = GSL_POSINF;
                                            /**< default value for post break wait s */
    afterload_restretch_vel = 0;            /**< default value for restretch velocity */
    afterload_factor_s = GSL_POSINF;        /**< default value for afterload factor s */
    afterload_factor_multiplier = 1.0;      /**< default value for afterload multiplier */

    start_status_time_step = 0;             /**< default value */
    stop_status_time_step = 0;              /**< default value */
    skip_status_time_step = 0;              /**< default value */

    // Update values from log file
    set_FiberSim_options_from_JSON_file_string(JSON_options_file_string);

    // Do some processing on the options

    if (strlen(rate_file_string) > 0)
    {
        if (!strcmp(rate_relative_to, "this_file"))
        {
            fs::path options_file = JSON_options_file_string;
            fs::path options_path = options_file.parent_path();
            fs::path rate_path = options_path / rate_file_string;

            sprintf_s(rate_file_string, _MAX_PATH, "%s", rate_path.string().c_str());
        }
    }

    if (strlen(status_folder) > 0)
    {
        if (!strcmp(status_relative_to, "this_file"))
        {
            fs::path options_file = JSON_options_file_string;
            fs::path options_path = options_file.parent_path();
            fs::path status_path = fs::absolute(options_path / status_folder);

            // Make sure the status folder exists
            if (fs::exists(status_path))
            {
                // Clean the directory
                printf("Cleaning status_folder: %s\n", status_path.string().c_str());
                for (auto const& dir_entry : fs::recursive_directory_iterator(status_path))
                {
                    fs::remove(dir_entry);
                }
            }
            else
            {
                // Create the directory
                if (fs::create_directories(status_path))
                {
                    printf("Status folder created at: %s\n", status_path.string().c_str());
                }
                else
                {
                    printf("Status folder could not be created: %s\n", status_path.string().c_str());
                    exit(1);
                }
                
            }
            // Set the status folder
            sprintf_s(status_folder, _MAX_PATH, "%s", status_path.string().c_str());
        }

        /*
        // Parse the time_steps string
        std::string ts_string = time_steps_string;

        size_t first_sep = ts_string.find_first_of(":");
        size_t last_sep = ts_string.find_last_of(":");

        start_status_time_step = (int)std::stoi(ts_string.substr(0, first_sep));
        skip_status_time_step = (int)std::stoi(ts_string.substr((first_sep+1), last_sep));
        stop_status_time_step = (int)std::stoi(ts_string.substr(last_sep+1));
        */
    }

    if (strlen(log_folder) > 0)
    {
        log_mode = 1;
        if (!strcmp(log_relative_to, "this_file"))
        {
            fs::path options_file = JSON_options_file_string;
            fs::path options_path = options_file.parent_path();
            fs::path log_path = options_path / log_folder;

            // Make sure the status folder exists
            if (fs::is_directory(log_path))
            {
                // Clean the directory
                int n = (int)fs::remove_all(log_path);
                printf("Deleting %i files from status_folder: %s\n",
                    n, log_path.string().c_str());
            }

            // Now create the directory
            if (fs::create_directories(log_path))
            {
                printf("Log folder created at: %s\n", log_path.string().c_str());
            }
            else
            {
                printf("Log folder could not be created: %s\n", log_path.string().c_str());
                exit(1);
            }

            // Set the log folder
            sprintf_s(log_folder, _MAX_PATH, "%s", log_path.string().c_str());

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

            write_FiberSim_options_to_file();
        }
    }
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

    if (JSON_functions::is_JSON_member(options, "max_rate"))
    {
        JSON_functions::check_JSON_member_number(options, "max_rate");
        max_rate = options["max_rate"].GetDouble();
    }

    if (JSON_functions::is_JSON_member(options, "thin_filament_sub_steps"))
    {
        JSON_functions::check_JSON_member_int(options, "thin_filament_sub_steps");
        thin_filament_sub_steps = options["thin_filament_sub_steps"].GetInt();
    }

    if (JSON_functions::is_JSON_member(options, "adjacent_bs"))
    {
        JSON_functions::check_JSON_member_int(options, "adjacent_bs");
        adjacent_bs = options["adjacent_bs"].GetInt();
    }

    if (JSON_functions::is_JSON_member(options, "x_pos_rel_tol"))
    {
        JSON_functions::check_JSON_member_number(options, "x_pos_rel_tol");
        x_pos_rel_tol = options["x_pos_rel_tol"].GetDouble();
    }

    // In FiberSim versions > 1.2.0, x_pos_tol notation is preferred over x_pos_rel_tol

    if (JSON_functions::is_JSON_member(options, "x_pos_tol"))
    {
        JSON_functions::check_JSON_member_number(options, "x_pos_tol");
        x_pos_rel_tol = options["x_pos_tol"].GetDouble();
    }

    if (JSON_functions::is_JSON_member(options, "x_vector_max_iterations"))
    {
        JSON_functions::check_JSON_member_int(options, "x_vector_max_iterations");
        x_vector_max_iterations = options["x_vector_max_iterations"].GetInt();
    }

    if (JSON_functions::is_JSON_member(options, "hs_force_control_max_delta_hs_length"))
    {
        JSON_functions::check_JSON_member_number(options, "hs_force_control_max_delta_hs_length");
        hs_force_control_max_delta_hs_length = options["hs_force_control_max_delta_hs_length"].GetDouble();
    }

    if (JSON_functions::is_JSON_member(options, "min_hs_length"))
    {
        JSON_functions::check_JSON_member_number(options, "min_hs_length");
        min_hs_length = options["min_hs_length"].GetDouble();
    }

    // Check if lambda_jitter was specified.
    if (JSON_functions::is_JSON_member(options, "lambda_jitter"))
    {
        JSON_functions::check_JSON_member_number(options, "lambda_jitter");
        lambda_jitter = options["lambda_jitter"].GetDouble();
    }

    // Check for m_filament_density scaling
    if (JSON_functions::is_JSON_member(options, "m_filament_density_ref_hs_length"))
    {
        JSON_functions::check_JSON_member_number(options, "m_filament_density_ref_hs_length");
        m_filament_density_ref_hs_length = options["m_filament_density_ref_hs_length"].GetDouble();
    }

    // Check for calculate_x_mode
    if (JSON_functions::is_JSON_member(options, "calculate_x_mode"))
    {
        JSON_functions::check_JSON_member_int(options, "calculate_x_mode");
        calculate_x_mode = options["calculate_x_mode"].GetInt();
    }

    // Check for rand_seed - set to empty string if missing
    if (JSON_functions::is_JSON_member(options, "rand_seed"))
    {
        JSON_functions::check_JSON_member_string(options, "rand_seed");
        sprintf_s(rand_seed, _MAX_PATH, "%s", options["rand_seed"].GetString());
    }
    else
    {
        sprintf_s(rand_seed, _MAX_PATH, "");
    }

    // Check for myofibrils
    if (JSON_functions::check_JSON_member_exists(options, "myofibrils"))
    {
        const rapidjson::Value& myofibrils = options["myofibrils"];

        JSON_functions::check_JSON_member_number(myofibrils, "force_tolerance");
        myofibril_force_tolerance = myofibrils["force_tolerance"].GetDouble();

        JSON_functions::check_JSON_member_int(myofibrils, "max_iterations");
        myofibril_max_iterations = myofibrils["max_iterations"].GetInt();

        JSON_functions::check_JSON_member_number(myofibrils, "max_delta_hs_length");
        myofibril_max_delta_hs_length = myofibrils["max_delta_hs_length"].GetDouble();

        JSON_functions::check_JSON_member_int(myofibrils, "multithreading");
        myofibril_multithreading = myofibrils["multithreading"].GetInt();
    }

    // Check for afterload
    if (JSON_functions::check_JSON_member_exists(options, "afterload"))
    {
        const rapidjson::Value& afterload = options["afterload"];

        JSON_functions::check_JSON_member_number(afterload, "load");
        afterload_load = afterload["load"].GetDouble();

        JSON_functions::check_JSON_member_number(afterload, "break_delta_hs_length");
        afterload_break_delta_hs_length = afterload["break_delta_hs_length"].GetDouble();

        if (JSON_functions::check_JSON_member_exists(afterload, "post_break_wait_s"))
        {
            afterload_post_break_wait_s = afterload["post_break_wait_s"].GetDouble();
        }

        if (JSON_functions::check_JSON_member_exists(afterload, "restretch_vel"))
        {
            afterload_restretch_vel = afterload["restretch_vel"].GetDouble();
        }

        if (JSON_functions::check_JSON_member_exists(afterload, "factor_s"))
        {
            afterload_factor_s = afterload["factor_s"].GetDouble();
        }

        if (JSON_functions::check_JSON_member_exists(afterload, "factor_multiplier"))
        {
            afterload_factor_multiplier = afterload["factor_multiplier"].GetDouble();
        }
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

    // Check for rate logging
    if (JSON_functions::is_JSON_member(options, "rate_files"))
    {
        const rapidjson::Value& rate_files = options["rate_files"];

        JSON_functions::check_JSON_member_string(rate_files, "relative_to");
        sprintf_s(rate_relative_to, _MAX_PATH, "%s", rate_files["relative_to"].GetString());

        JSON_functions::check_JSON_member_string(rate_files, "file");
        sprintf_s(rate_file_string, _MAX_PATH, "%s", rate_files["file"].GetString());
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

        // Check if the dump precision was specified.
        if (JSON_functions::is_JSON_member(status_files, "dump_precision"))
        {
            JSON_functions::check_JSON_member_int(status_files, "dump_precision");
            dump_precision = status_files["dump_precision"].GetInt();
        }

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

void FiberSim_options::parse_status_time_steps_string(int no_of_time_points)
{
    // Sets the time-points values

    // Code
    if (strcmp(time_steps_string, "last") == 0)
    {
        start_status_time_step = no_of_time_points - 1;
        skip_status_time_step = 1;
        stop_status_time_step = start_status_time_step;
    }
    else
    {
        // Parse the time_steps string
        std::string ts_string = time_steps_string;

        size_t first_sep = ts_string.find_first_of(":");
        size_t last_sep = ts_string.find_last_of(":");

        start_status_time_step = (int)std::stoi(ts_string.substr(0, first_sep));
        skip_status_time_step = (int)std::stoi(ts_string.substr((first_sep + 1), last_sep));
        stop_status_time_step = (int)std::stoi(ts_string.substr(last_sep + 1));
    }
}