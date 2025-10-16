# Doxygen Documentation Setup for LetzElPhC

## Summary of Changes

This document outlines all the improvements made to set up professional Doxygen documentation for the LetzElPhC project.

### 1. Doxyfile Improvements

**File**: `docs/Doxyfile`

#### Issues Fixed:
- ✅ **C Language Optimization**: Changed `OPTIMIZE_FOR_FORTRAN = YES` to `NO`
- ✅ **Extension Mapping**: Removed incorrect Fortran mapping for C files
- ✅ **Duplicate Configuration**: Removed duplicate markdown/config sections
- ✅ **Output Directory**: Changed from absolute path `/doxygen_output` to relative `./doxygen_output`
- ✅ **Path Stripping**: Updated `STRIP_FROM_PATH` to use relative paths
- ✅ **Removed Non-Existent Paths**: Cleaned up references to non-existent example/image directories
- ✅ **CSS Stylesheet**: Updated path to `./custom_dark.css`
- ✅ **Dark Theme**: Set `HTML_COLORSTYLE = DARK` for better aesthetics
- ✅ **Warning Log**: Updated path to `./doxygen_warnings.log`

#### New Configuration Features:
- ✅ **Project Logo**: Added project logo reference
- ✅ **Math Support**: Enabled MathJax for mathematical formulas
- ✅ **Custom Aliases**: Added documentation aliases for authors (@mn, @lw, @rr) and custom sections (@formula, @algorithm, @complexity, @physics, @reference)
- ✅ **Call Graphs**: Enabled call and caller graphs for better code understanding
- ✅ **SVG Output**: SVG format for diagrams (better quality and interactive)
- ✅ **Dynamic Features**: Code folding, clipboard copy, and interactive SVG

### 2. Documentation Generation Script

**File**: `docs/generate_docs.sh`

#### Improvements:
- ✅ **Path Consistency**: Updated all paths to match Doxyfile configuration
- ✅ **Cleaned Paths**: Removed `doc/` prefix from all output directories
- ✅ **Better Error Messages**: Clear feedback on dependency checks
- ✅ **Statistics Report**: Shows number of HTML files, SVG graphs, and total size
- ✅ **Cross-Platform Support**: Works on macOS, Linux with automatic browser opening
- ✅ **Better User Feedback**: Color-coded output messages

#### Usage Examples:
```bash
# Full documentation with graphs
./generate_docs.sh

# Quick mode (no graphs, faster)
./generate_docs.sh -q

# Clean and regenerate with verbose output
./generate_docs.sh -c -v

# Generate and open in browser
./generate_docs.sh -o

# Show help
./generate_docs.sh -h
```

### 3. Custom Dark Theme

**File**: `docs/doxygen/custom_dark.css`

Features (already present, no changes needed):
- ✅ Elegant dark background (#0e0e11)
- ✅ Light text for readability
- ✅ Blue accent colors (#66b3ff)
- ✅ Smooth scrolling
- ✅ Professional typography
- ✅ Dark navigation panel

### 4. Main Documentation Page

**File**: `docs/mainpage.md` (NEW)

Created comprehensive main documentation page with:
- ✅ Project introduction
- ✅ Key features overview
- ✅ Code modules description
- ✅ Quick start guide
- ✅ Physics background with LaTeX formulas
- ✅ Authors and acknowledgments
- ✅ Important notes about code scope
- ✅ Contributing guidelines
- ✅ References section

### 5. README.md Enhancements

**File**: `README.md`

Major improvements:
- ✅ **Better Organization**: Added table of contents
- ✅ **Build Instructions**: Clear step-by-step guide
- ✅ **Documentation Section**: Comprehensive guide to generating and viewing docs
- ✅ **Code Structure**: Documented all source modules
- ✅ **Examples**: Listed available example calculations
- ✅ **Contributing Guidelines**: How to contribute to the project
- ✅ **Roadmap**: Converted TODO list to checkbox roadmap
- ✅ **Tool Recommendations**: Listed complementary codes for physical properties calculations

### 6. Source Code Documentation

#### Files Enhanced:

##### `src/elphC.h`
- ✅ **File-Level Documentation**: Description and authors
- ✅ **Type Definitions**: Documented all typedefs with Doxygen commands
- ✅ **Precision Group**: Created documentation group for floating-point precision settings
- ✅ **MPI Types**: Documented MPI datatype mappings
- ✅ **NetCDF Types**: Documented I/O datatype mappings

##### `src/main_elphC.c`
- ✅ **File Documentation**: Overview of program entry point
- ✅ **Main Function**: Comprehensive documentation with:
  - Purpose description
  - Parameter documentation
  - Return value documentation
  - Detailed algorithm steps
  - Usage examples
  - MPI execution guide

### 7. Documentation Features to Note

#### For Developers:
- **Custom Aliases** in Doxyfile allow consistent documentation:
  ```c
  /// \mn       // Mark as Muralidhar's code
  /// \lw       // Mark as Ludger's code
  /// \rr       // Mark as Riccardo's code
  /// \formula  // Document mathematical formulas
  /// \algorithm// Document algorithms
  /// \physics  // Add physics background
  ```

#### For Users:
- **Beautiful Dark Theme**: Professional appearance with dark background
- **Interactive Diagrams**: Call graphs and include hierarchies with SVG
- **Search Functionality**: Full-text search built-in
- **Mathematical Formulas**: LaTeX support via MathJax
- **Cross-References**: Automatic linking between related code sections
- **Source Browsing**: View actual source code with syntax highlighting

#### Performance Options:
- **Quick Mode**: Generates without graphs for faster builds
- **Verbose Mode**: Shows detailed Doxygen output
- **Clean Build**: Remove old documentation before regenerating

## Next Steps: Adding More Documentation

To add documentation to more C source files, use this template:

```c
/**
 * @file filename.c
 * @brief One-line description of the file
 * @author Author Name
 * 
 * Detailed description of what this file contains and does.
 * Explain the purpose, key functions, and algorithms.
 */

/**
 * @brief Brief description of function
 * 
 * @param param1 Description of param1
 * @param param2 Description of param2
 * 
 * @return Description of return value
 * 
 * @details Extended description with algorithm details, notes, or warnings
 * 
 * @note Any special notes
 * @warning Any warnings
 * @see Related functions or structures
 */
void function_name(int param1, float param2);
```

## Generation Instructions for Users

### Prerequisites
```bash
# macOS
brew install doxygen graphviz

# Ubuntu/Debian
sudo apt-get install doxygen graphviz

# CentOS/RHEL
sudo yum install doxygen graphviz
```

### Generate Documentation
```bash
cd docs
./generate_docs.sh -o  # Generate and open in browser
```

### View Documentation
```bash
# Open in your browser
open doxygen_output/html/index.html  # macOS
xdg-open doxygen_output/html/index.html  # Linux
```

## Troubleshooting

### Issue: Graphviz not installed
**Solution**: Script automatically falls back to quick mode (no graphs)

### Issue: Documentation not generating
**Solution**: Run with verbose mode:
```bash
./generate_docs.sh -v
```

### Issue: Permission denied on generate_docs.sh
**Solution**: Make it executable:
```bash
chmod +x docs/generate_docs.sh
```

## Files Modified

1. `docs/Doxyfile` - Configuration file (improved)
2. `docs/generate_docs.sh` - Generation script (path fixes)
3. `docs/doxygen/custom_dark.css` - Custom theme (no changes needed)
4. `docs/mainpage.md` - Main documentation page (NEW)
5. `README.md` - Project README (enhanced)
6. `src/elphC.h` - Core header (added Doxygen comments)
7. `src/main_elphC.c` - Main program (added Doxygen comments)

## Future Enhancements

Recommended additions for complete documentation:

1. **Add Comments to All Modules**:
   - `src/elph/elph.h` - Electron-phonon coupling API
   - `src/io/io.h` - Input/output interfaces
   - `src/interpolation/interpolation.h` - Interpolation functions
   - `src/symmetries/symmetries.h` - Symmetry operations
   - And other key header files

2. **Create Architecture Documentation**:
   - Create `docs/architecture.md`
   - Explain data flow and module dependencies
   - Add workflow diagrams

3. **Add Tutorial Pages**:
   - Create `docs/tutorials.md`
   - Step-by-step examples with code snippets
   - Common use cases and best practices

4. **Performance Guidelines**:
   - Create `docs/performance.md`
   - MPI/OpenMP tuning recommendations
   - Profiling and optimization tips

## Maintenance

Keep documentation updated when:
- Adding new functions
- Changing function signatures
- Modifying important algorithms
- Adding new features

Run periodic checks:
```bash
cd docs
./generate_docs.sh -c -v  # Clean rebuild with warnings
```

Review the warnings log to catch undocumented code:
```bash
cat docs/doxygen_warnings.log
```