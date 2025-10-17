/**
 * @file preprocessor.h
 * @brief Command-line parsing and preprocessing utilities
 */

#pragma once

#include "common/dtypes.h"
#include "elphC.h"

#define PH_X_READ_BUF_SIZE 512
#define ELPH_MAX_ENV_SIZE 64
#define PH_SAVE_DIR_NAME_DEFAULT "ph_save"

/**
 * @brief Parse command-line arguments
 * @param argc Number of command-line arguments
 * @param argv Array of argument strings
 * @param calc_info Calculation details structure (out)
 * @return void
 * @details Parses options like -i, -w, --help, --version and sets calc_info accordingly
 */
void ELPH_cli_parser(int argc, char* argv[], struct calc_details* calc_info);

/**
 * @brief Create ph_save directory for Quantum ESPRESSO %phonon output
 * @param inp_file Quantum ESPRESSO input file path
 * @return void
 * @details Extracts outdir/prefix and creates ph_save subdirectory structure
 */
void create_ph_save_dir_pp_qe(const char* inp_file);

/**
 * @brief Print program version and copyright information
 * @return void
 */
void ELPH_print_version(void);

/**
 * @brief Print help message with usage information
 * @return void
 */
void ELPH_print_help(void);
