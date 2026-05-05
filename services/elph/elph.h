#pragma once
#include <mpi.h>
#include <stdbool.h>

#include "common/dtypes.h"
#include "elphC.h"

/*
 * Per-(iq_BZ, ik_BZ) fill callback for yambo COLL integration.
 * Called once per BZ (q,k) pair on commK_rank==0 processes.
 *   iq_BZ, ik_BZ : 0-based BZ indices
 *   data          : ELPH_cmplx buffer, C row-major (nmodes, nspin, nbnds, nbnds)
 *   nq..nb_start  : full-BZ dimensions (constant across calls); nb_start is 1-based
 */
typedef void (*elph_fill_fn)(int iq_BZ, int ik_BZ,
                              const void* data,
                              int nq, int nk, int nmodes, int nspin,
                              int nbnds, int nb_start);

void elph_driver(const char* ELPH_input_file, enum ELPH_dft_code dft_code,
                 MPI_Comm comm_world);

/* Callback-enabled variant: skips ndb.elph; calls fill_fn per (q,k) instead. */
void elph_driver_cb(const char* ELPH_input_file, enum ELPH_dft_code dft_code,
                    MPI_Comm comm_world, elph_fill_fn fill_fn);

void compute_and_write_elphq(struct WFC* wfcs, struct Lattice* lattice,
                             struct Pseudo* pseudo, struct Phonon* phonon,
                             const ND_int iqpt, ELPH_cmplx* eigVec,
                             ELPH_cmplx* dVscfq, const int ncid_elph,
                             const int varid_elph, const int ncid_dmat,
                             const int varid_dmat, const bool non_loc,
                             const bool kminusq,
                             const struct ELPH_MPI_Comms* Comm,
                             elph_fill_fn fill_fn);

void compute_and_write_dmats(const char* file_name, const struct WFC* wfcs,
                             const struct Lattice* lattice,
                             const ND_int nph_sym,
                             const struct symmetry* sym_data,
                             const struct ELPH_MPI_Comms* Comm);

void init_kernel(struct kernel_info* kernel);

void set_kernel(const char* kernel_str, struct kernel_info* kernel);
