/**
 * @file wfc.h
 * @brief Wave function manipulation and symmetry operations
 */

#pragma once
#include <mpi.h>
#include <stdbool.h>

#include "elphC.h"

/**
 * @brief Rotate G-vectors under symmetry operation
 * @param Gvecs Input G-vector coordinates
 * @param sym Symmetry rotation matrix
 * @param ngvecs Number of G-vectors
 * @param lat_vec Lattice vectors
 * @param inverse Apply inverse symmetry flag
 * @param crystal Use crystal coordinates flag
 * @param G0 Output rotated G-vectors
 * @param Gvec_out Output G-vector list
 * @return void
 */
void rotateGvecs(const ELPH_float* Gvecs, const ELPH_float* sym,
                 const ND_int ngvecs, const ELPH_float* lat_vec,
                 const bool inverse, const bool crystal, ELPH_float* G0,
                 ELPH_float* Gvec_out);

/**
 * @brief Apply SU(2) transformation to wave functions with spin-orbit coupling
 * @param nspinor Number of spinor components (1 or 2)
 * @param npw Number of plane waves
 * @param nsets Number of bands
 * @param su2mat SU(2) transformation matrix
 * @param wfc Wave function coefficients (in/out)
 * @return void
 */
void su2rotate(const int nspinor, const ND_int npw, const ND_int nsets,
               const ELPH_cmplx* su2mat, ELPH_cmplx* wfc);

/**
 * @brief Apply translational symmetry phase factor to wave functions
 * @param trans_vec Translation vector
 * @param kvec k-point coordinate
 * @param nsets Number of bands
 * @param npw Number of plane waves
 * @param gvecs G-vector list
 * @param wfc_G Wave function in G-space (in/out)
 * @param conjugate Apply conjugate operation flag
 * @return void
 */
void apply_trans_wfc(const ELPH_float* trans_vec, const ELPH_float* kvec,
                     const ND_int nsets, const ND_int npw,
                     const ELPH_float* gvecs, ELPH_cmplx* wfc_G,
                     const bool conjugate);

/**
 * @brief Generate FFT box from energy cutoff
 * @param EcutRy Energy cutoff in Rydberg
 * @param blat Reciprocal lattice vectors
 * @param fft_box Output FFT box dimensions [Nx, Ny, Nz]
 * @param commK MPI communicator for k-points
 * @return void
 */
void get_fft_box(const ELPH_float EcutRy, const ELPH_float* blat,
                 ND_int* fft_box, MPI_Comm commK);

/**
 * @brief Construct SU(2) matrix for spin-orbit rotation
 * @param sym_in Symmetry rotation matrix (3x3)
 * @param nspinor Number of spinor components
 * @param invert_sym Apply inverse operation flag
 * @param time_rev Apply time-reversal symmetry flag
 * @param su2mat Output SU(2) matrix
 * @return void
 */
void SU2mat(const ELPH_float* sym_in, const ND_int nspinor,
            const bool invert_sym, const bool time_rev, ELPH_cmplx* su2mat);

/**
 * @brief Sort and redistribute plane waves across MPI processes
 * @param npw_tot Total number of plane waves
 * @param npw_loc Local number of plane waves per process
 * @param fft_dims FFT grid dimensions
 * @param gvec_in Input G-vector list
 * @param wfc_in Input wave function coefficients
 * @param nsets Number of bands
 * @param npw_loc_out Output local plane wave count
 * @param nG_xy Output G_xy dimension
 * @param gvec_out Output redistributed G-vector list
 * @param wfc_out Output redistributed wave functions
 * @param mpi_comm MPI communicator
 * @return void
 */
void Sort_pw(const ND_int npw_tot, const ND_int npw_loc, const ND_int* fft_dims,
             const ELPH_float* gvec_in, const ELPH_cmplx* wfc_in,
             const ND_int nsets, ND_int* npw_loc_out, ND_int* nG_xy,
             int** gvec_out, ELPH_cmplx** wfc_out, MPI_Comm mpi_comm);
