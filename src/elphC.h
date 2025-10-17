/**
 * @file elphC.h
 * @brief Core type definitions and configuration for LetzElPhC
 * @author Muralidhar Nalabothula
 * 
 * This header file contains fundamental type definitions and compile-time configuration
 * macros for the LetzElPhC code. It ensures C99 compliance and sets up the floating-point
 * precision (single or double) used throughout the code.
 */

#pragma once

#if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#ifdef __STDC_NO_COMPLEX__
#error Your compiler does not support C99 complex numbers, Please use a supported compiler.
#endif
#else
#error Your compiler does not support C99 standard.
#endif

/**
 * @typedef ND_int
 * @brief Large integer type for array indexing and counters
 * 
 * Always set to at least 64-bit integer to handle large arrays
 * and ensure compatibility with MPI communications.
 */
typedef long long int ND_int;

/**
 * @def ELPH_MPI_ND_INT
 * @brief MPI datatype corresponding to ND_int
 */
#define ELPH_MPI_ND_INT MPI_LONG_LONG_INT

/**
 * @defgroup Precision Floating-Point Precision Configuration
 * @brief Configure single or double precision for all calculations
 * 
 * The floating-point precision can be selected at compile time by defining
 * `COMPILE_ELPH_DOUBLE` flag. By default, double precision is used.
 * @{
 */

#if defined(COMPILE_ELPH_DOUBLE)

/**
 * @typedef ELPH_float
 * @brief Main floating-point type for the code (double precision)
 */
typedef double ELPH_float;

/**
 * @typedef ELPH_cmplx
 * @brief Complex floating-point type (double precision)
 */
typedef double _Complex ELPH_cmplx;

/**
 * @def ELPH_MPI_float
 * @brief MPI datatype for ELPH_float (double precision)
 */
#define ELPH_MPI_float MPI_DOUBLE

/**
 * @def ELPH_MPI_cmplx
 * @brief MPI datatype for ELPH_cmplx (double complex)
 */
#define ELPH_MPI_cmplx MPI_C_DOUBLE_COMPLEX

/**
 * @def ELPH_NC4_IO_FLOAT
 * @brief NetCDF4 output datatype (double precision)
 */
#define ELPH_NC4_IO_FLOAT NC_DOUBLE

#else

/**
 * @typedef ELPH_float
 * @brief Main floating-point type for the code (single precision)
 */
typedef float ELPH_float;

/**
 * @typedef ELPH_cmplx
 * @brief Complex floating-point type (single precision)
 */
typedef float _Complex ELPH_cmplx;

/**
 * @def ELPH_MPI_float
 * @brief MPI datatype for ELPH_float (single precision)
 */
#define ELPH_MPI_float MPI_FLOAT

/**
 * @def ELPH_MPI_cmplx
 * @brief MPI datatype for ELPH_cmplx (single complex)
 */
#define ELPH_MPI_cmplx MPI_C_FLOAT_COMPLEX

/**
 * @def ELPH_NC4_IO_FLOAT
 * @brief NetCDF4 output datatype (single precision)
 */
#define ELPH_NC4_IO_FLOAT NC_FLOAT

#endif
/** @} */
