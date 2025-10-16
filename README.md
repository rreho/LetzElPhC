<img src="./logo.png" alt="LetzElPhC Logo" width="1200">

LetzElPhC abbreviates to _"Lëtzebuerg Electron Phonon Code"_.
_"Lëtzebuerg"_ is the Luxembourgish name for the Luxembourg Country.

LetzElPhC is distributed under the MIT license.

## Table of Contents

- [Building & Installation](#building--installation)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [Examples](#examples)
- [Development](#development)
- [Acknowledgments](#acknowledgments)

## Building & Installation

Please refer to [**Installation & Usage Guide**](../../main.pdf) (PDF) for detailed **installation and usage instructions**. 

A simple example can be found in ```examples/qe/silicon``` folder.

## Quick Start

1. **Configure the build**: Copy one of the sample configuration files and adapt it to your system:
   ```bash
   cp sample_config/make_*.inc src/make.inc
   # Edit src/make.inc with your compiler and library paths
   ```

2. **Build the code**:
   ```bash
   cd src
   make clean
   make -j4
   ```

3. **Run a calculation**: See the examples directory for sample input files and configurations.

## Documentation

### Generating API Documentation

LetzElPhC includes comprehensive API documentation built with Doxygen. 

#### Prerequisites

Install Doxygen and Graphviz:
```bash
# macOS
brew install doxygen graphviz

# Ubuntu/Debian
sudo apt-get install doxygen graphviz

# CentOS/RHEL
sudo yum install doxygen graphviz
```

#### Generate Documentation

From the `docs/` directory:

```bash
# Full documentation with graphs
cd docs
./generate_docs.sh

# Quick generation without graphs (faster)
./generate_docs.sh -q

# Clean previous build and regenerate with verbose output
./generate_docs.sh -c -v

# Generate and automatically open in browser
./generate_docs.sh -o
```

#### View Documentation

After generation, open the documentation in your browser. Run the script with `-o` flag to auto-open:
```bash
./generate_docs.sh -o
```

Or manually open `docs/doxygen_output/html/index.html` with your web browser.

#### Documentation Options

The generate script supports several flags:
- `-h, --help`:   Show help message
- `-q, --quick`:  Quick generation (disables graphs for faster build)
- `-o, --open`:   Open documentation in browser after generation
- `-c, --clean`:  Clean previous documentation before generating
- `-v, --verbose`: Show all Doxygen messages

## Examples

Simple example calculations are provided in the `examples/` directory:
- `examples/qe/silicon/`: Silicon bulk calculation with Quantum ESPRESSO

## Development

### Status

**This code is currently under development. Please use at your own risk.**

### Contributing

Contributions to the code are always welcomed. Thank you very much if you are willing to contribute.

### Reporting Issues

In case of any bugs or other issues, please open a GitLab issue at:
https://gitlab.com/lumen-code/LetzElPhC

## Important Notes

**Scope**: This code computes **only electron-phonon matrix elements**. For calculating physical properties 
from electron-phonon coupling (such as mobilities, carrier scattering rates, superconducting properties, etc.),
use complementary tools such as:
- **EPW**: Electron-Phonon Wannier code
- **BoltzWann**: Boltzmann transport with Wannier functions
- **Perturbo**: Ultrafast dynamics from first principles

# Acknowledgments

- **Ludger Wirtz** (Supervisor)  
- **Fulvio Paleari** (Testing and Yambopy interface)  
- **University of Luxembourg** (Funding)  
- **HPC @ Uni.lu** (Computing resources)  
- **Henry Fried** (Logo)

# Roadmap

- [x] Support XML format for dynamical matrices
- [x] Image parallelization support for ph.x (preprocessor)
- [x] Automatic test cases
- [ ] Improve OpenMP efficiency
- [ ] Implement acoustic sum rule
- [x] Multiple kernel options for users
- [ ] Fröhlich polaron interpolation
- [ ] DFT + U support



