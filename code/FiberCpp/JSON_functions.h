#pragma once

/** 
 * @file    muscle.h
 * @brief   header file for the Muscle class
 * @author  Ken Campbell
 */

#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"

#include "gsl_vector.h"
#include "gsl_matrix.h"

namespace JSON_functions {

    /**
    * a function that checks whether an object exists in a JSON doc
    * @param doc pointer to a rapidjson::Document
    * @param mem_name pointer to a member_name char array
    * @return void
    */
    void check_JSON_member_object(const rapidjson::Value& doc, const char mem_name[]);

    /**
     * a function that checks whether mem_name is an integer member of doc
     * @param doc a pointer to a rapidjson::Document
     * @param mem_name a char array
     * @return void
     */
    void check_JSON_member_int(const rapidjson::Value& doc, const char mem_name[]);

    /**
     * a function that checks whether mem_name is a number member of doc
     * @param doc a pointer to a rapidjson::Document
     * @param mem_name a char array
     * @return void
     */
    void check_JSON_member_number(const rapidjson::Value& doc, const char mem_name[]);

    /**
     * a function that checks whether mem_name is a string in doc
     * @param doc a pointer to a rapidjson::Document
     * @param mem_name a char array
     * @return void
     */
    void check_JSON_member_string(const rapidjson::Value& doc, const char mem_name[]);

    /**
    *a function that checks whether mem_name is an array
    * @param doc a pointer to a rapidjson::Document
    * @param mem_name an array
    * @return void
    */
    void check_JSON_member_array(const rapidjson::Value& doc, const char mem_name[]);

    /**
    * a function that checks whether mem_name exists
    * @param doc a pointer to a rapidjson::Document
    * @param memname[] a character array
    * @return int, 0 for does not exist, 1 for does exist
    */
    int is_JSON_member(const rapidjson::Value& doc, const char mem_name[]);

    /**
    * a function that adds a gsl_vector to a file as a line in JSON format
    * @param p_v a pointer to a gsl_vector
    * @param output_file a pointer to the file stream (which must be open)
    * @param label_string a char array for the label  (for example, "cb_x")
    * @param is_last_entry a boolean saying whether this is the last entry
    *        + if it is not, add a , to the end of the line
    * @param precision the precision that doubles are output with
    */
    void write_gsl_vector_as_JSON_array(gsl_vector* p_v, FILE* output_file,
        char label_string[], bool is_last_entry, int precision);

    /**
    * a function that adds a gsl_vector_short to a file as a line in JSON format
    * @param p_v a pointer to a gsl_vector_short
    * @param output_file a pointer to the file stream (which must be open)
    * @param label_string a char array for the label  (for example, "cb_x")
    * @param is_last_entry a boolean saying whether this is the last entry
    *        + if it is not, add a , to the end of the line
    * @param precision the precision that doubles are output with
    */
    void write_gsl_vector_short_as_JSON_array(gsl_vector_short* p_v, FILE* output_file,
        char label_string[], bool is_last_entry);

    /**
    * a function that writes a short inte array to a file as a line in JSON format
    * @param p_v a pointer to a short int array
    * @param n_entries an integer with the number of entries to write
    * @param output_file a pointer to the file stream (which must be open)
    * @param label_string a char array for the label  (for example, "cb_x")
    * @param is_last_entry a boolean saying whether this is the last entry
    *        + if it is not, add a , to the end of the line
    */
    void write_short_int_array_as_JSON_array(short int p_v[], int n_entries,
        FILE* output_file, char label_string[], bool is_last_entry);

    /**
    * a function that adds a gsl_natrix to a file as a line in JSON format
    * @param p_v a pointer to a gsl_matrix
    * @param output_file a pointer to the file stream(which must be open)
    * @param label_string a char array for the label(for example, "cb_x")
    * @param is_last_entry a boolean saying whether this is the last entry
    * +if it is not, add a, to the end of the line
    * @param precision the precision that doubles are output with
    */
    void write_gsl_matrix_as_JSON_array(gsl_matrix* p_v, FILE* output_file,
        char label_string[], bool is_last_entry, int precision);

    /**
   * a function that adds a gsl_natrix short to a file as a line in JSON format
   * @param p_v a pointer to a gsl_matrix
   * @param output_file a pointer to the file stream(which must be open)
   * @param label_string a char array for the label(for example, "cb_x")
   * @param is_last_entry a boolean saying whether this is the last entry
   * +if it is not, add a, to the end of the line
   */
    void write_gsl_matrix_short_as_JSON_array(gsl_matrix_short* p_v, FILE* output_file,
        char label_string[], bool is_last_entry);

};