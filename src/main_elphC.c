/**
 * @file main_elphC.c
 * @brief Main program entry point for LetzElPhC
 * @author Muralidhar Nalabothula
 * 
 * This file contains the entry point of the LetzElPhC program. It handles:
 * - MPI initialization and communication setup
 * - Command-line argument parsing
 * - Routing to appropriate calculation modules (preprocessing, electron-phonon coupling, interpolation)
 * - Final cleanup and MPI finalization
 * 
 * @section usage Usage
 * 
 * The program can be invoked with various commands:
 * - `lelphc --help`: Show help message
 * - `lelphc --version`: Show version information
 * - `lelphc -F <input_file>`: Run main electron-phonon calculation
 * - `lelphc -code <input_file>`: DFT code. Currently only qe is supported
 * - `lelphc -i `: Interpolation mode.

 * 
 * @section mpi MPI Execution
 * 
 * The program uses MPI for parallelization. Run with mpirun/mpiexec:
 * ```
 * mpirun -n <num_procs> ./lelphc -F input.in
 * ```
 * 
 * OpenMP is also supported for hybrid parallelization when compiled with
 * the ELPH_OMP_PARALLEL_BUILD flag.
 */

#include "elphC.h"
#include "preprocessor/preprocessor.h"
#include "interpolation/interpolation.h"
#include "elph/elph.h"
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>

/**
 * @brief Main program entry point
 * 
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * 
 * @return Exit status (0 on success, non-zero on error)
 * 
 * @details
 * 
 * The main function performs the following steps:
 * 1. Initialize MPI communication (with thread support if compiled with OpenMP)
 * 2. Parse command-line arguments on the master process (rank 0)
 * 3. Broadcast configuration to all processes
 * 4. Route execution to appropriate module:
 *    - elph_driver(): Compute electron-phonon matrix elements
 *    - interpolation_driver(): Perform Wannier interpolation
 *    - Preprocessor operations: Prepare input files
 * 5. Finalize MPI and return
 */
int main(int argc, char* argv[])
{
#if defined(ELPH_OMP_PARALLEL_BUILD)
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
    MPI_Init(&argc, &argv);
#endif

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    struct calc_details* calc_info = malloc(sizeof(struct calc_details));

    ND_int dft_code_tmp;
    // this is used to broad cast the enum
    bool run_elph_calc = false;
    bool run_interpolation = false;
    if (my_rank == 0)
    {
        ELPH_cli_parser(argc, argv, calc_info);

        dft_code_tmp = calc_info->code;

        if (calc_info->calc == CALC_ELPH)
        {
            run_elph_calc = true;
        }
        else if (calc_info->calc == CALC_INTERPOLATION)
        {
            run_interpolation = true;
        }
    }

    MPI_Bcast(&dft_code_tmp, 1, ELPH_MPI_ND_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&run_elph_calc, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&run_interpolation, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(calc_info->input_file, sizeof(calc_info->input_file), MPI_CHAR, 0, MPI_COMM_WORLD);

    calc_info->code = dft_code_tmp;

    if (run_elph_calc)
    {
        elph_driver(calc_info->input_file, calc_info->code, MPI_COMM_WORLD);
    }
    else if (run_interpolation)
    {
        interpolation_driver(calc_info->input_file, calc_info->code, MPI_COMM_WORLD);
    }
    else
    {
        // run preprocessor
        if (my_rank == 0)
        {
            if (calc_info->calc == CALC_HELP)
            {
                ELPH_print_help();
            }
            else if (calc_info->calc == CALC_VERSION)
            {
                ELPH_print_version();
            }
            else if (calc_info->calc == CALC_PH_SAVE_CREATE)
            {
                if (calc_info->code == DFT_CODE_QE)
                {
                    create_ph_save_dir_pp_qe(calc_info->input_file);
                }
            }
        }
    }

    free(calc_info);

    MPI_Finalize();

    return 0;
}
