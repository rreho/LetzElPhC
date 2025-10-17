#pragma once
#include <mpi.h>
#include <stdbool.h>

#include "common/dtypes.h"
#include "elphC.h"

/**
 * @brief Main driver for electron-%phonon coupling calculations
 * @param ELPH_input_file Path to configuration file
 * @param dft_code DFT code type (QE, ABINIT, etc.)
 * @param comm_world MPI communicator for all processes
 * @return void
 * @details Coordinates the full electron-%phonon calculation workflow including
 * DFT data loading, phonon matrix computations, and output generation
 */
void elph_driver(const char* ELPH_input_file, enum ELPH_dft_code dft_code,
                 MPI_Comm comm_world);

/**
 * @brief Compute and write electron-%phonon coupling matrices for a q-point
 * @param wfcs Wave function data structure
 * @param lattice Lattice/cell information
 * @param pseudo Pseudopotential data
 * @param phonon %Phonon mode information
 * @param iqpt q-point index
 * @param eigVec Phonon eigenvectors
 * @param dVscfq Derivative of local potential in q-space
 * @param ncid_elph NetCDF file ID for coupling matrices
 * @param varid_elph NetCDF variable ID for coupling
 * @param ncid_dmat NetCDF file ID for deformation potentials
 * @param varid_dmat NetCDF variable ID for deformation
 * @param non_loc Include non-local contributions flag
 * @param kminusq Use k-q symmetry flag
 * @param Comm MPI communicator structure
 * @return void
 */
void compute_and_write_elphq(struct WFC* wfcs, struct Lattice* lattice,
                             struct Pseudo* pseudo, struct Phonon* phonon,
                             const ND_int iqpt, ELPH_cmplx* eigVec,
                             ELPH_cmplx* dVscfq, const int ncid_elph,
                             const int varid_elph, const int ncid_dmat,
                             const int varid_dmat, const bool non_loc,
                             const bool kminusq,
                             const struct ELPH_MPI_Comms* Comm);

/**
 * @brief Compute and write deformation potential matrices
 * @param file_name Output NetCDF file path
 * @param wfcs Wave function data
 * @param lattice Lattice information
 * @param nph_sym Number of %phonon symmetry operations
 * @param sym_data Symmetry operation matrices and translations
 * @param Comm MPI communicator structure
 * @return void
 */
void compute_and_write_dmats(const char* file_name, const struct WFC* wfcs,
                             const struct Lattice* lattice,
                             const ND_int nph_sym,
                             const struct symmetry* sym_data,
                             const struct ELPH_MPI_Comms* Comm);

/**
 * @brief Initialize kernel with default settings
 * @param kernel Kernel info structure to initialize
 * @return void
 */
void init_kernel(struct kernel_info* kernel);

/**
 * @brief Set kernel type and options from string
 * @param kernel_str Kernel name string (bare_loc, dfpt, etc.)
 * @param kernel Kernel info structure to configure
 * @return void
 */
void set_kernel(const char* kernel_str, struct kernel_info* kernel);
