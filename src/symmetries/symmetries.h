#pragma once
#include <stdbool.h>

#include "common/dtypes.h"
#include "elphC.h"

/**
 * @brief Expand irreducible Brillouin zone k-points to full Brillouin zone
 * @param Nibz Number of irreducible k-points
 * @param Nsym Number of symmetry operations
 * @param ibz_kpts Irreducible k-point coordinates
 * @param symms Symmetry operations (rotation + translation)
 * @param lat_vec Lattice vectors
 * @param kpoints Output full BZ k-points
 * @param kstar Symmetry star indices
 * @param kmap k-point mapping array
 * @return Total number of k-points in full BZ
 */
ND_int bz_expand(const ND_int Nibz, const ND_int Nsym,
                 const ELPH_float* ibz_kpts, const struct symmetry* symms,
                 const ELPH_float* lat_vec, ELPH_float* kpoints, ND_int* kstar,
                 int* kmap);

/**
 * @brief Generate irreducible Brillouin zone k-points from grid
 * @param kgrid K-point mesh dimensions [nk1, nk2, nk3]
 * @param Nsym Number of symmetry operations
 * @param symms Symmetry operations (rotation + translation)
 * @param lat_vec Lattice vectors
 * @param blat Reciprocal lattice vectors
 * @param ibz_kpts Output irreducible k-point coordinates
 * @param crystal Use crystal coordinates flag
 * @return Number of irreducible k-points
 */
ND_int generate_iBZ_kpts(const ND_int* kgrid, const ND_int Nsym,
                         const struct symmetry* symms,
                         const ELPH_float* lat_vec, const ELPH_float* blat,
                         ELPH_float* ibz_kpts, const bool crystal);

/**
 * @brief Generate electronic band representations under symmetry
 * @param wfcs Wave function data
 * @param lattice Lattice structure
 * @param Rsym_mat Rotation matrix of symmetry operation
 * @param tauR Translation of symmetry operation
 * @param tim_revR Time-reversal symmetry flag
 * @param ikBZ k-point index in full BZ
 * @param Dkmn_rep Output representation matrix
 * @param Comm MPI communicator structure
 * @return void
 */
void electronic_reps(const struct WFC* wfcs, const struct Lattice* lattice,
                     const ELPH_float* Rsym_mat, const ELPH_float* tauR,
                     const bool tim_revR, const ND_int ikBZ,
                     ELPH_cmplx* Dkmn_rep, const struct ELPH_MPI_Comms* Comm);

/**
 * @brief Rotate electron-%phonon coupling matrix under symmetry
 * @param Dmats_l Left deformation potential matrix
 * @param elph_mat_q Electron-%phonon coupling at q
 * @param Dmats_r Right deformation potential matrix
 * @param lattice Lattice structure
 * @param tim_rev Time-reversal symmetry flag
 * @param elph_mat_Sq Output rotated coupling matrix at Sq
 * @return void
 */
void elph_q_rotate(const ELPH_cmplx* Dmats_l, const ELPH_cmplx* elph_mat_q,
                   const ELPH_cmplx* Dmats_r, const struct Lattice* lattice,
                   const bool tim_rev, ELPH_cmplx* elph_mat_Sq);

/**
 * @brief Rotate %phonon eigenvectors under symmetry
 * @param sym Symmetry operation
 * @param lattice Lattice structure
 * @param qpt q-point coordinate
 * @param eig_q Input eigenvectors at q
 * @param eig_Sq Output eigenvectors at rotated point Sq
 * @return void
 */
void rotate_eig_vecs(struct symmetry* sym, const struct Lattice* lattice,
                     const ELPH_float* qpt, const ELPH_cmplx* eig_q,
                     ELPH_cmplx* eig_Sq);

/**
 * @brief Rotate derivative of local potential under symmetry
 * @param dvscf_in Input derivative potential
 * @param sym Symmetry operation
 * @param lattice Lattice structure
 * @param composite_form Composite wave vector representation flag
 * @param dvscf_out Output rotated potential
 * @param commK MPI communicator for k-points
 * @return void
 */
void rotate_dvscf(const ELPH_cmplx* dvscf_in, struct symmetry* sym,
                  const struct Lattice* lattice, const bool composite_form,
                  ELPH_cmplx* restrict dvscf_out, MPI_Comm commK);
