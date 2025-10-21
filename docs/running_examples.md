# Running LetzElPhC Examples {#running_examples}

## Overview

LetzElPhC includes example calculations that demonstrate how to use the code for electron-%phonon coupling calculations.
These examples use Quantum ESPRESSO as the DFT backend and some python fortran scripts. Scripts are provided in the folder `scripts` and the necessary conventional libraries must be pre-installed. For the Fortran files, the same dependencies as the ones employed for QE should also be installed. The following sections provide detailed instructions for running these examples.

## Available Examples

- Si: `examples/qe/silicon`
- hBN:`examples/qe/hBN`
- MoS2:`examples/qe/mos2`

The DFT ground states calculation can be run all at once with the script `run.sh` available in the relative folders.
You can run in parallel changing the `NCORES` available in `run.sh`.
Post-processing plots can be generated with python scripts, available in the `scripts` folder.

### hBN Bulk Calculation


This is a simple example computing electron-%phonon coupling for bulk hBN.

#### Files in the Example

- `scf/scf.in`: Self-consistent field (SCF) calculation input for pw.x
- `nscf/nscf.in`: %Phonon calculation input for ph.x
- `ph/ph.in`: LetzElPhC configuration file for electron-phonon coupling calculation

#### Running the Silicon Example

**Step 1: Run the DFT Calculation**

First, execute the Quantum ESPRESSO SCF and nscf calculation:

```bash
cd examples/qe/silicon/scf/
pw.x < scf.in > scf.out
cd ../nscf/
pw.x < nscf.in > nscf.out
```

This generates the ground state wave functions and density.

**Step 2: Run the %Phonon Calculation**

Next, compute the %phonon frequencies and polarization vectors:

```bash
cd examples/qe/hBN/ph/
ph.x < ph.in > ph.out
```

This creates the dynamical matrix files needed by LetzElPhC.

**STEP 3: Initialize Yambo databases (optional)**
Initialize the `yambo` database and create the `SAVE` dir:
This step is needed only if you want to post-process the data to compute excitons-related observables.

```bash
cd examples/qe/hBN/nscf/hBN.save/
p2y
yambo
## now compute the elph matrix elements
cd elph
$LELPH -F elph.in
cd ..
cd bse
mpirun -np $NCORES $YAMBO -F bse.in -J BSE -C BSE -I ../nscf/hBN.save/SAVE
cd ..

**Step 4: Run LetzElPhC**

Finally, compute the electron-phonon coupling matrices:
First run the pre-processing tool of `LetzElPhC` to generate `ph_save` directory in the `ph folder`  
```bash
cd ph
$LELPH -pp --code=qe -F ph.in
cd ..
## now compute the elph matrix elements
cd elph
$LELPH -F elph.in
cd examples/qe/hBN/elph/
../../src/lelphc elph.in > elph.out
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

The LetzElPhC configuration file (`elph.in`) contains settings for:

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

LetzElPhC outputs NetCDF files containing electron-%phonon data. Those data can be interfaced with the output data of other codes such as `Yambo`.
As examples, we use **Python** scripts to:

- plot the exicton-phonon matrix elements on the Brillouin Zone;
- compute phonon-assisted luminescence;
- compute one- and two- phonon Raman spectra.

### Python Analysis


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

## Next Steps

- Modify the example for your system of interest
- Explore different interpolation grids and k-point meshes
- Combine with EPW for transport properties
- Develop custom analysis workflows
