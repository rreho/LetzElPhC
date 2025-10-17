/**
 * @file Vnonloc.h
 * @brief Non-local pseudopotential electron-\%phonon coupling
 */

#pragma once
#include "common/dtypes.h"
#include "elphC.h"

/**
 * @brief Sum non-local potential contribution across k and k+q points
 * @param K_ptr Projector coefficients at k
 * @param Kp_ptr Projector coefficients at k+q
 * @param fcoeff Non-local \%phonon coupling coefficient
 * @param nspin Number of spin states
 * @param nbnd Number of bands
 * @param nspinor Spinor components (1 or 2)
 * @param out Output non-local coupling matrix element
 * @return void
 */
void sum_VNL_KKp(ELPH_cmplx* K_ptr, ELPH_cmplx* Kp_ptr, ELPH_cmplx* fcoeff,
                 ND_int nspin, ND_int nbnd, ND_int nspinor, ELPH_cmplx* out);

/**
 * @brief Add non-local pseudopotential electron-\%phonon coupling
 * @param wfcs Wave function data
 * @param lattice Lattice structure
 * @param pseudo Pseudopotential data including projectors
 * @param ikq k+q point index
 * @param ik k-point index
 * @param kqsym Symmetry index for k+q
 * @param ksym Symmetry index for k
 * @param eigVec \%Phonon eigenvector
 * @param elph_kq_mn Non-local coupling matrix (in/out)
 * @param Comm MPI communicator structure
 * @return void
 */
void add_elphNonLocal(struct WFC* wfcs, struct Lattice* lattice,
                      struct Pseudo* pseudo, int ikq, int ik, int kqsym,
                      int ksym, ELPH_cmplx* eigVec, ELPH_cmplx* elph_kq_mn,
                      const struct ELPH_MPI_Comms* Comm);
