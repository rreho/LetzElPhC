/**
 * @file dvloc.h
 * @brief Derivative of local potential computations
 */

#pragma once
#include <stdbool.h>

#include "common/dtypes.h"
#include "elphC.h"

/**
 * @brief Compute local potential in G-space from radial grid
 * @param work_arr Working array for interpolation
 * @param cutoff Cutoff type flag
 * @param Gnorm G-vector norm
 * @param ngrid Number of radial grid points
 * @param Vloc_atomic Atomic local potential on radial grid
 * @param r_grid Radial grid points
 * @param rab_grid Radial grid spacing weights
 * @param Zval Valence charge
 * @param eta Gaussian spreading parameter
 * @param volume Unit cell volume
 * @return Local potential value in G-space
 */
ELPH_float Vloc_Gspace(ELPH_float* work_arr, const char cutoff,
                       const ELPH_float Gnorm, const ND_int ngrid,
                       const ELPH_float* Vloc_atomic, const ELPH_float* r_grid,
                       const ELPH_float* rab_grid, const ELPH_float Zval,
                       const ELPH_float eta, const ELPH_float volume);

/**
 * @brief Add local potential derivative from QE format
 * @param dVscf Derivative scattering potential (in/out)
 * @param dVloc Derivative of local potential
 * @param lattice Lattice structure
 * @return void
 */
void add_dvscf_qe(ELPH_cmplx* dVscf, const ELPH_cmplx* dVloc,
                  const struct Lattice* lattice);

/**
 * @brief Compute local electron-%phonon coupling matrix element
 * @param qpt q-point coordinate
 * @param wfcs Wave function data
 * @param lattice Lattice structure
 * @param ikq k+q point index
 * @param ik k-point index
 * @param kqsym Symmetry index for k+q
 * @param ksym Symmetry index for k
 * @param dVlocr Derivative potential real-space
 * @param Comm MPI communicator structure
 * @param elph_kq Output coupling matrix element
 * @return void
 */
void elphLocal(const ELPH_float* qpt, struct WFC* wfcs, struct Lattice* lattice,
               int ikq, int ik, int kqsym, int ksym, ELPH_cmplx* dVlocr,
               const struct ELPH_MPI_Comms* Comm, ELPH_cmplx* elph_kq);

/**
 * @brief Compute local potential derivative for q-point
 * @param qpt q-point coordinate
 * @param lattice Lattice structure
 * @param pseudo Pseudopotential data
 * @param eigVec %Phonon eigenvector
 * @param Vlocr Output local potential derivative real-space
 * @param commK MPI communicator for k-points
 * @return void
 */
void dVlocq(const ELPH_float* qpt, struct Lattice* lattice,
            struct Pseudo* pseudo, const ELPH_cmplx* eigVec, ELPH_cmplx* Vlocr,
            MPI_Comm commK);

/**
 * @brief Compute long-range interaction electron-%phonon coupling
 * @param qpt q-point coordinate
 * @param gvecs G-vector list
 * @param npw_loc Number of local plane waves
 * @param Zvals Atomic charge values
 * @param epslion Dielectric constant tensor
 * @param Zeu Born effective charge tensor
 * @param Qpole Quadrupole tensor
 * @param natom Number of atoms
 * @param atom_pos Atomic positions
 * @param diminsion Dimensionality (3D/2D)
 * @param volume Unit cell volume
 * @param zlat Z-lattice parameter (for 2D)
 * @param EcutRy Energy cutoff
 * @param elph_lr_out Output long-range coupling
 * @return void
 */
void dVlong_range_kernel(const ELPH_float* qpt, const ELPH_float* gvecs,
                         const ND_int npw_loc, const ELPH_float* Zvals,
                         const ELPH_float* epslion, const ELPH_float* Zeu,
                         const ELPH_float* Qpole, const ND_int natom,
                         const ELPH_float* atom_pos, const char diminsion,
                         const ELPH_float volume, const ELPH_float zlat,
                         const ELPH_float EcutRy, ELPH_cmplx* elph_lr_out);

/**
 * @brief Add long-range correction to potential derivative
 * @param qpt q-point coordinate
 * @param lattice Lattice structure
 * @param phonon %Phonon data
 * @param Zvals Atomic charge values
 * @param eigVec %Phonon eigenvector
 * @param dVscf Derivative potential (in/out)
 * @param sign Sign for correction (1 or -1)
 * @param only_induced_part Compute only induced part flag
 * @param EcutRy Energy cutoff
 * @param nmags_add Additional magnetic contribution flags
 * @param commK MPI communicator for k-points
 * @return void
 */
void dV_add_longrange(const ELPH_float* qpt, struct Lattice* lattice,
                      struct Phonon* phonon, const ELPH_float* Zvals,
                      const ELPH_cmplx* eigVec, ELPH_cmplx* dVscf,
                      const ND_int sign, const bool only_induced_part,
                      const ELPH_float EcutRy, const bool* nmags_add,
                      MPI_Comm commK);

/**
 * @brief Change basis representation of potential derivative
 * @param dvscf Derivative potential (in/out)
 * @param rot_vecs Rotation vectors for basis change
 * @param nsets Number of modes
 * @param nmodes Number of phonon modes
 * @param nmag Number of magnetic sites
 * @param Nx FFT X dimension
 * @param Ny FFT Y dimension
 * @param Nz FFT Z dimension
 * @param blas_char BLAS character flag
 * @return void
 */
void dVscf_change_basis(ELPH_cmplx* dvscf, const ELPH_cmplx* rot_vecs,
                        const ND_int nsets, const ND_int nmodes,
                        const ND_int nmag, const ND_int Nx, const ND_int Ny,
                        const ND_int Nz, const char blas_char);

/**
 * @brief Multiply potential by phase factor exp(ikr)
 * @param pot_grid Potential on real-space grid (in/out)
 * @param qpt_crys q-point in crystal coordinates
 * @param lattice Lattice structure
 * @param nsets Number of modes
 * @param sign Phase factor sign (+1 or -1)
 * @return void
 */
void multiply_eikr(ELPH_cmplx* pot_grid, const ELPH_float* qpt_crys,
                   const struct Lattice* lattice, const ND_int nsets,
                   const ND_int sign);

/**
 * @brief Create lookup table for local potential values
 * @param lattice Lattice structure
 * @param pseudo Pseudopotential data
 * @param Comm MPI communicator structure
 * @return void
 */
void create_vlocg_table(const struct Lattice* lattice, struct Pseudo* pseudo,
                        const struct ELPH_MPI_Comms* Comm);

/**
 * @brief Free local potential lookup table
 * @param vloc_table Lookup table to free
 * @return void
 */
void free_vlocg_table(struct Vloc_table* vloc_table);
