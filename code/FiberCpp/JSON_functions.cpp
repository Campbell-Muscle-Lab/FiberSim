/**
 * @file    JSON_functions.cpp
 * @brief   Source file for functions used to handle JSON data
 * @author  Ken Campbell
 */

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"

#include "JSON_functions.h"

namespace JSON_functions {

    //! Verifies that JSON Value has member and it is an object.
    void check_JSON_member_object(const rapidjson::Value& doc,
        const char mem_name[])
    {
        char temp_string[_MAX_PATH];
        try {
            if (!doc.HasMember(mem_name)) {
                sprintf_s(temp_string, "\"%s\" not specified.", mem_name);
                throw std::runtime_error(temp_string);
            }
            else if (!doc[mem_name].IsObject()) {
                sprintf_s(temp_string, "\"%s\" must be an object.", mem_name);
                throw std::runtime_error(temp_string);
            }
        }
        catch (std::exception & e) {
            printf("Exception: %s", e.what());
            exit(1);
        }
    }
    
    //! Verifies that JSON Value has a member and it is an integer
    void check_JSON_member_int(const rapidjson::Value& doc,
        const char mem_name[])
    {
        char temp_string[_MAX_PATH];
        try {
            if (!doc.HasMember(mem_name)) {
                sprintf_s(temp_string, "\"%s\" not specified.", mem_name);
                throw std::runtime_error(temp_string);
            }
            else if (!doc[mem_name].IsInt()) {
                sprintf_s(temp_string, "\"%s\" must be an integer.", mem_name);
                throw std::runtime_error(temp_string);
            }
        }
        catch (std::exception & e) {
            printf("Exception: %s", e.what());
            exit(1);
        }
    }

    //! Verifies that JSON Value has member and it is a number.
    void check_JSON_member_number(const rapidjson::Value& doc,
        const char mem_name[])
    {
        char temp_string[_MAX_PATH];
        try {
            if (!doc.HasMember(mem_name)) {
                sprintf_s(temp_string, "\"%s\" not specified.", mem_name);
                throw std::runtime_error(temp_string);
            }
            else if (!doc[mem_name].IsNumber()) {
                sprintf_s(temp_string, "\"%s\" must be a number.", mem_name);
                throw std::runtime_error(temp_string);
            }
        }
        catch (std::exception & e) {
            printf("Exception: %s", e.what());
            exit(1);
        }
    }

    //! Verifies that JSON Value has member and it is a string.
    void check_JSON_member_string(const rapidjson::Value& doc,
        const char mem_name[])
    {
        char temp_string[_MAX_PATH];
        try {
            if (!doc.HasMember(mem_name)) {
                sprintf_s(temp_string, "\"%s\" not specified.", mem_name);
                throw std::runtime_error(temp_string);
            }
            else if (!doc[mem_name].IsString()) {
                sprintf_s(temp_string, "\"%s\" must be a string.", mem_name);
                throw std::runtime_error(temp_string);
            }
        }
        catch (std::exception & e) {
            printf("Exception: %s", e.what());
            exit(1);
        }
    }

    //! Verifies that JSON Value has member and it is an array.
    void check_JSON_member_array(const rapidjson::Value& doc, const char mem_name[])
    {
        char temp_string[_MAX_PATH];
        try {
            if (!doc.HasMember(mem_name)) {
                sprintf_s(temp_string, "\"%s\" not specified.", mem_name);
                throw std::runtime_error(temp_string);
            }
            else if (!doc[mem_name].IsArray()) {
                sprintf_s(temp_string, "\"%s\" must be an array.", mem_name);
                throw std::runtime_error(temp_string);
            }
        }
        catch (std::exception & e) {
            printf("Exception: %s", e.what());
            exit(1);
        }
    }
        
    //! Checks for existance of JSON member number
    int is_JSON_member(const rapidjson::Value& doc, const char mem_name[])
    {
        if (doc.HasMember(mem_name))
            return 1;
        else
            return 0;
    }

    //! Writes gsl_vector to JSON file
    void write_gsl_vector_as_JSON_array(gsl_vector* p_v, FILE* output_file,
        char label_string[], bool is_last_entry, int precision)
    {
        // Code writes a gsl_vector to the file which must be open

        // Write the label
        fprintf_s(output_file, "\t\"%s\": [", label_string);

        // Write the values
        for (int i = 0; i < (p_v->size); i++)
        {
            fprintf(output_file, "%.*F", precision, gsl_vector_get(p_v, i));
            if (i == ((p_v->size) - 1))
            {
                // It's the last entry in the array
                fprintf_s(output_file, "]");
                if (!is_last_entry)
                {
                    fprintf_s(output_file, ",");
                }
                fprintf_s(output_file, "\n");
            }
            else
                fprintf_s(output_file, ", ");
        }
    }

    //! Writes gsl_vector_short to JSON file
    void write_gsl_vector_short_as_JSON_array(gsl_vector_short* p_v, FILE* output_file,
        char label_string[], bool is_last_entry)
    {
        // Code writes a gsl_vector to the file which must be open

        // Write the label
        fprintf_s(output_file, "\t\"%s\": [", label_string);

        // Write the values
        for (int i = 0; i < (p_v->size); i++)
        {
            fprintf(output_file, "%i", gsl_vector_short_get(p_v, i));
            if (i == ((p_v->size) - 1))
            {
                // It's the last entry in the array
                fprintf_s(output_file, "]");
                if (!is_last_entry)
                {
                    fprintf_s(output_file, ",");
                }
                fprintf_s(output_file, "\n");
            }
            else
                fprintf_s(output_file, ", ");
        }
    }

    //! Writes gsl_matrix to JSON file
    void write_gsl_matrix_as_JSON_array(gsl_matrix* p_v, FILE* output_file,
        char label_string[], bool is_last_entry, int precision)
    {
        // Code writes a gsl_matrix to the file which must be open

        // Write the label for the first col
        for (int col = 0; col < p_v->size2; col++)
        {
            fprintf_s(output_file, "\t\"%s[x_%i]\": [", label_string, col);
            for (int row = 0; row < p_v->size1 ; row++)
            {
                fprintf_s(output_file, "%.*F", precision, gsl_matrix_get(p_v, row, col));
                if (row == (p_v->size1 - 1))
                {
                    // It's the last entry in the array
                    fprintf_s(output_file, "]");
                    if ((is_last_entry) && (col == (p_v->size2-1)))
                    {
                        fprintf_s(output_file, "\n");
                    }
                    else
                    {
                        fprintf_s(output_file, ",\n");
                    }
                }
                else
                    fprintf_s(output_file, ", ");
            }
        }
    }

    //! Writes gsl_matrix to JSON file
    void write_gsl_matrix_short_as_JSON_array(gsl_matrix_short* p_v, FILE* output_file,
        char label_string[], bool is_last_entry, int transpose)
    {
        // Code writes a gsl_matrix to the file which must be open

        // Write the label for the first col
        for (int col = 0; col < p_v->size2; col++)
        {
            if (transpose==0)
                fprintf_s(output_file, "\t\"%s[x_%i]\": [", label_string, col);
            else
                fprintf_s(output_file, "\t\"%s[%i_x]\": [", label_string, col);

            for (int row = 0; row < p_v->size1; row++)
            {
                fprintf_s(output_file, "%i", gsl_matrix_short_get(p_v, row, col));
                if (row == (p_v->size1 - 1))
                {
                    // It's the last entry in the array
                    fprintf_s(output_file, "]");
                    if (!(is_last_entry) || (col<(p_v->size2-1)))
                    {
                        fprintf_s(output_file, ",");
                    }
                    fprintf_s(output_file, "\n");
                }
                else
                    fprintf_s(output_file, ", ");
            }
        }
    }

    void write_short_int_array_as_JSON_array(
        short int p_v[], int n_entries,
        FILE* output_file, char label_string[], bool is_last_entry)
    {
        // Code writes a gsl_vector to the file which must be open

        // Write the label
        fprintf_s(output_file, "\t\"%s\": [", label_string);

        // Write the values
        for (int i = 0 ; i < n_entries ; i++)
        {
            fprintf(output_file, "%i", p_v[i]);
            if (i == (n_entries - 1))
            {
                // It's the last entry in the array
                fprintf_s(output_file, "]");
                if (!is_last_entry)
                {
                    fprintf_s(output_file, ",");
                }
                fprintf_s(output_file, "\n");
            }
            else
                fprintf_s(output_file, ", ");
        }
    }
}

