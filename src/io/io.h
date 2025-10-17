/**
 * @file io.h
 * @brief Input/output operations for DFT data and results
 */

#pragma once

#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <stddef.h>

#include "common/dtypes.h"
#include "elphC.h"

#define NC4_DEFAULT_CHUCK_KB 2048
// default chunking for large nc varaibles (in Kilobytes)

/**
 * @brief Read and allocate DFT save data from Quantum ESPRESSO output
 * @param SAVEdir Path to DFT save directory
 * @param Comm MPI communicator structure
 * @param start_band Starting band index
 * @param end_band Ending band index
 * @param wfcs Wave function structure (out)
 * @param ph_save_dir \%Phonon save directory path
 * @param lattice Lattice structure (out)
 * @param pseudo Pseudopotential data (out)
 * @param phonon \%Phonon data (out)
 * @param dft_code DFT code type
 * @return void
 */
void read_and_alloc_save_data(char* SAVEdir, const struct ELPH_MPI_Comms* Comm,
                              ND_int start_band, ND_int end_band,
                              struct WFC** wfcs, char* ph_save_dir,
                              struct Lattice* lattice, struct Pseudo* pseudo,
                              struct Phonon* phonon,
                              enum ELPH_dft_code dft_code);

/**
 * @brief Free \%phonon data structures
 * @param phonon \%Phonon data to free
 * @return void
 */
void free_phonon_data(struct Phonon* phonon);

/**
 * @brief Free all DFT save data structures
 * @param wfcs Wave function data
 * @param lattice Lattice data
 * @param pseudo Pseudopotential data
 * @param phonon \%Phonon data
 * @return void
 */
void free_save_data(struct WFC* wfcs, struct Lattice* lattice,
                    struct Pseudo* pseudo, struct Phonon* phonon);

/**
 * @brief Parse UPF pseudopotential file for local potential
 * @param filename UPF file path
 * @param loc_pseudo Local pseudopotential structure (out)
 * @return void
 */
void parse_upf(const char* filename, struct local_pseudo* loc_pseudo);

/**
 * @brief Extract element symbol and valence charge from UPF file
 * @param filename UPF file path
 * @param atomic_sym Element symbol (out)
 * @param Zval Valence charge (out)
 * @return void
 */
void get_upf_element(const char* filename, char* atomic_sym, ELPH_float* Zval);

/**
 * @brief Initialize user input configuration structure
 * @param input Input structure pointer (out)
 * @return void
 */
void init_usr_input(struct usr_input** input);

/**
 * @brief Free user input configuration structure
 * @param input Input structure to free
 * @return void
 */
void free_usr_input(struct usr_input* input);

/**
 * @brief Read and parse configuration file
 * @param input_file Configuration file path
 * @param input_data Input data structure (out)
 * @param MPI_world_comm MPI communicator
 * @return void
 * @details Broadcasts parsed data to all MPI processes
 */
void read_input_file(const char* input_file, struct usr_input** input_data,
                     MPI_Comm MPI_world_comm);

/**
 * @brief Define NetCDF variable with optional chunking
 * @param ncid NetCDF file ID
 * @param varid Variable ID (out)
 * @param rank Variable rank (number of dimensions)
 * @param xtype NetCDF data type
 * @param dims Dimension sizes array
 * @param var_name Variable name
 * @param dim_names Dimension names array
 * @param chunksize Chunk sizes for each dimension
 * @return void
 */
void def_ncVar(const int ncid, int* varid, ND_int rank, nc_type xtype,
               ND_int* dims, const char* var_name, char** dim_names,
               size_t* chunksize);

/**
 * @brief Write basic metadata to NetCDF file
 * @param ncid NetCDF file ID
 * @param lattice Lattice structure
 * @param phonon \%Phonon data
 * @param kernel_str Kernel type string
 * @param convention_str Convention string
 * @return void
 */
void write_basic_data(const int ncid, struct Lattice* lattice,
                      struct Phonon* phonon, const char* kernel_str,
                      const char* convention_str);
