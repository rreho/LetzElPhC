#pragma once

#include <mpi.h>

#include "common/dtypes.h"
#include "elphC.h"

/**
 * @brief Main driver for Wannier interpolation of electron-%phonon coupling
 * @param ELPH_input_file Path to configuration file
 * @param dft_code DFT code type (QE, ABINIT, etc.)
 * @param comm_world MPI communicator for all processes
 * @return void
 * @details Interpolates electron-%phonon coupling from DFT q-points to arbitrary q-points using Wannier functions
 */
void interpolation_driver(const char* ELPH_input_file,
                          enum ELPH_dft_code dft_code, MPI_Comm comm_world);
