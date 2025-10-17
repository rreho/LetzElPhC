/**
 * @file phonon.h
 * @brief Phonon mode handling and dynamical matrix operations
 */

#pragma once
#include "common/dtypes.h"
#include "elphC.h"

/**
 * @brief Mass-normalize %phonon polarization vectors
 * @param atomic_masses Atomic masses
 * @param nsets Number of modes
 * @param natoms Number of atoms
 * @param power Normalization power (typically -0.5)
 * @param pol_vecs Polarization vectors (in/out)
 * @return void
 */
void mass_normalize_pol_vecs(const ELPH_float* atomic_masses,
                             const ND_int nsets, const ND_int natoms,
                             const ELPH_float power, ELPH_cmplx* pol_vecs);

/**
 * @brief Convert %phonon polarization vectors to dynamical matrix
 * @param omega %Phonon frequencies
 * @param natom Number of atoms
 * @param atomic_masses Atomic masses
 * @param pol_vecs Polarization vectors (in/out as dynamical matrix)
 * @return void
 */
void pol_vecs_to_dyn(const ELPH_float* omega, const ND_int natom,
                     const ELPH_float* atomic_masses, ELPH_cmplx* pol_vecs);

/**
 * @brief Add long-range interaction to %phonon dynamical matrix
 * @param qpt q-point coordinate
 * @param lattice Lattice structure
 * @param phonon Phonon data with Born effective charges
 * @param Ggrid FFT grid for long-range corrections
 * @param sign Sign for long-range kernel (1 or -1)
 * @param atomic_masses Atomic masses
 * @param dyn_mat_asr ASR-corrected dynamical matrix
 * @param dyn_mat Output dynamical matrix with long-range (in/out)
 * @return void
 */
void add_ph_dyn_long_range(const ELPH_float* qpt, struct Lattice* lattice,
                           struct Phonon* phonon, const ND_int* Ggrid,
                           const ND_int sign, const ELPH_float* atomic_masses,
                           ELPH_cmplx* dyn_mat_asr, ELPH_cmplx* dyn_mat);

/**
 * @brief Compute long-range and acoustic sum rule correction
 * @param lattice Lattice structure
 * @param phonon %Phonon data with Born effective charges
 * @param Ggrid FFT grid for long-range corrections
 * @param atomic_masses Atomic masses
 * @param dyn_mat_asr Output ASR-corrected dynamical matrix
 * @return void
 */
void compute_dyn_lr_asr_correction(struct Lattice* lattice,
                                   struct Phonon* phonon, const ND_int* Ggrid,
                                   const ELPH_float* atomic_masses,
                                   ELPH_cmplx* dyn_mat_asr);
