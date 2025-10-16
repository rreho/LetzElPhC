# LetzElPhC Documentation {#mainpage}

## Introduction

**LetzElPhC** (Lëtzebuerg Electron-Phonon Code) is a high-performance computational physics software
designed for calculating **electron-phonon coupling matrix elements** from first-principles calculations.

The code processes output from Quantum ESPRESSO (QE) DFT calculations and computes the electron-phonon
interaction matrices needed for studying various electron-phonon phenomena in condensed matter systems.

## Key Features

- Direct calculation of electron-phonon matrix elements
- Parallelized with MPI and OpenMP for large-scale calculations
- Efficient Wannier interpolation of electron-phonon matrices
- Currently supports Quantum ESPRESSO
- Exploits crystal symmetries for computational efficiency
- Supports spin-orbit coupling and magnetic systems

## Building and Installation

See the main [README.md](../../README.md) for detailed installation instructions.


## API Documentation

- Browse the [File List](files.html) to explore the source code
- Visit the [Data Structures](annotated.html) for structure definitions
- Check the [Functions](functions.html) for available function calls
- See [Running Examples](pages.html) for usage tutorials

## Authors and Acknowledgments

**Main Developers:**
- Muralidhar Nalabothula (Code development)
- Ludger Wirtz (Supervision)
- Riccardo Reho (Documentation and extensions)

**Acknowledgments:**
- Fulvio Paleari (Testing and Yambopy interface)
- University of Luxembourg (Funding)
- HPC @ Uni.lu (Computing resources)

## License

LetzElPhC is distributed under the MIT License.
See LICENSE file for details.

## Contributing

Contributions to the code are always welcomed. Please follow the code style guidelines
and submit pull requests with clear descriptions of your changes.

For bug reports and feature requests, please open an issue at:
https://gitlab.com/lumen-code/LetzElPhC

## Important Notes

This code computes **only electron-phonon matrix elements**. For calculating derived 
physical properties (mobilities, superconducting properties, etc.), use complementary 
tools such as the [EPW](https://epw-code.org/) code.

The code is currently under active development. Please use at your own risk and 
report any issues encountered.

## References

For methodological details, see:
- Giustino, F., Cohen, M. L., & Louie, S. G. (2007). PhysRevB.76.165108
- Poncé, S., Margine, E. R., Verdi, C., & Giustino, F. (2018). COMP.MATER.SCI.140:1-22
- Nalabothula, M., et al. (2024). [Publication details]