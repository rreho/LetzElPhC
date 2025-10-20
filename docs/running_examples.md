# Running LetzElPhC Examples {#running_examples}

## Overview

LetzElPhC includes example calculations that demonstrate how to use the code for electron-%phonon coupling calculations.
These examples use Quantum ESPRESSO as the DFT backend and some python fortran scripts. Scripts are provided in the folder `scripts` and the necessary conventional libraries must be pre-installed. For the Fortran files, the same dependencies as the ones employed for QE should also be installed. The following sections provide detailed instructions for running these examples.

## Available Examples

### Silicon Bulk Calculation

Located in: `examples/qe/silicon/`

This is a simple example computing electron-%phonon coupling for bulk silicon at the Gamma point.

#### Files in the Example

- `Si.scf.in`: Self-consistent field (SCF) calculation input for pw.x
- `Si.ph.in`: %Phonon calculation input for ph.x
- `Si.elph.in`: LetzElPhC configuration file for electron-phonon coupling calculation
- `%README.md`: Detailed instructions for this example

#### Running the Silicon Example

**Step 1: Run the DFT Calculation**

First, execute the Quantum ESPRESSO SCF calculation:

```bash
cd examples/qe/silicon/
pw.x < Si.scf.in > Si.scf.out
```

This generates the ground state wave functions and density.

**Step 2: Run the %Phonon Calculation**

Next, compute the %phonon frequencies and polarization vectors:

```bash
ph.x < Si.ph.in > Si.ph.out
```

This creates the dynamical matrix files needed by LetzElPhC.

**Step 3: Run LetzElPhC**

Finally, compute the electron-phonon coupling matrices:

```bash
../../src/elphC Si.elph.in > Si.elph.out
```

Or with MPI parallelization:

```bash
mpirun -np 4 ../../src/elphC Si.elph.in > Si.elph.out
```

#### Output Files

The calculation produces NetCDF files containing:
- Electron-%phonon coupling matrix elements
- Band structure information
- %Phonon properties
- Symmetry operations applied

These files can be processed by tools like EPW or other electron-%phonon post-processing codes.

## Input File Format

The LetzElPhC configuration file (`Si.elph.in`) contains settings for:

### DFT Input
```ini
[DFT]
# Path to Quantum ESPRESSO output directory
dft_dir = ./

# DFT code type (currently supports QE)
dft_type = QE
```

### System Parameters
```ini
[System]
# Number of bands to include
nbnd = 20

# Number of q-points for %phonon calculations
nq = 4
```

### Interpolation Settings
```ini
[Interpolation]
# FFT grid for interpolation
fft_ngrid = 20 20 20

# Use Wannier functions for interpolation
use_wannier = false
```

### Output Options
```ini
[Output]
# Output file name
outdir = results/

# Output format (netcdf, hdf5, etc.)
output_format = netcdf
```

For complete parameter documentation, see the sample configuration file in `sample_config/`.

## Advanced Usage

### MPI Parallelization

Run with multiple processes:

```bash
mpirun -np 16 ./src/elphC config.in
```

### OpenMP Threading

Combine MPI with OpenMP for hybrid parallelization:

```bash
OMP_NUM_THREADS=4 mpirun -np 4 ./src/elphC config.in
```

### Custom Grid and k-point Sampling

Modify the input configuration to use different grids:

```ini
[Interpolation]
fft_ngrid = 40 40 40  # Finer interpolation grid
nkpts = 50            # More k-points
nqpts = 25            # More q-points
```

## Post-processing Results

LetzElPhC outputs NetCDF files containing electron-%phonon data. These can be processed with:

### Python Analysis

```python
import netCDF4
import numpy as np

# Load the output file
nc = netCDF4.Dataset('results/elph_coupling.nc', 'r')

# Access electron-%phonon matrix elements
g_matrices = nc.variables['g_matrix'][:]
frequencies = nc.variables['frequencies'][:]

nc.close()
```

### Further Calculations

The output matrices can be used with:
- **EPW**: Superconductivity, transport properties
- **BoltzWann**: Boltzmann transport
- **Perturbo**: Ultrafast dynamics
- Custom analysis scripts

## Troubleshooting

### Common Issues

**Problem**: "DFT output directory not found"
- **Solution**: Ensure `dft_dir` path in config file is correct and contains all necessary DFT output files

**Problem**: "Mismatch in number of k-points"
- **Solution**: Check that the k-point grid in DFT calculation matches the configuration file

**Problem**: "Symmetry operations inconsistent"
- **Solution**: Use the same pseudo-potentials and DFT parameters as the DFT calculation

### Getting Help

If you encounter issues:
1. Check the output file for error messages
2. Review the example configuration files
3. Check the [main documentation](../../docs/main.pdf)
4. Open an issue on GitLab: https://gitlab.com/lumen-code/LetzElPhC

## Example Workflow Summary

```bash
# Prepare working directory
mkdir -p my_calculation
cd my_calculation

# Copy example files
cp -r ../examples/qe/silicon/* .

# Build LetzElPhC (if not done)
cd ../../src
make clean && make -j4
cd ../examples/qe/silicon

# Run complete workflow
pw.x < Si.scf.in > Si.scf.out           # SCF
ph.x < Si.ph.in > Si.ph.out             # Phonons
mpirun -np 4 ../../src/elphC Si.elph.in # Electron-phonon
```

After successful completion, analyze the output NetCDF files with your post-processing tools.

## Next Steps

- Modify the example for your system of interest
- Explore different interpolation grids and k-point meshes
- Combine with EPW for transport properties
- Develop custom analysis workflows