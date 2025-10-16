#!/bin/bash
################################################################################
# LetzElPhC Documentation Generator
# 
# This script generates comprehensive Doxygen documentation for the LetzElPhC code
#
# Usage: ./generate_docs.sh [options]
#
# Options:
#   -h, --help       Show this help message
#   -q, --quick      Quick generation (no graphs)
#   -o, --open       Open documentation in browser after generation
#   -c, --clean      Clean previous documentation before generating
#   -v, --verbose    Verbose output
################################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Default options
QUICK_MODE=false
OPEN_BROWSER=false
CLEAN_FIRST=false
VERBOSE=false

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

################################################################################
# Functions
################################################################################

print_header() {
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC}  🌊 LetzElPhC Documentation Generator                          ${BLUE}║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

show_help() {
    cat << EOF
LetzElPhC Documentation Generator

Usage: ./generate_docs.sh [options]

Options:
  -h, --help       Show this help message
  -q, --quick      Quick generation (disable graphs for faster build)
  -o, --open       Open documentation in browser after generation
  -c, --clean      Clean previous documentation before generating
  -v, --verbose    Verbose output (show all Doxygen messages)

Examples:
  ./generate_docs.sh                    # Generate full documentation
  ./generate_docs.sh -q -o              # Quick generation and open in browser
  ./generate_docs.sh -c -v              # Clean and generate with verbose output

EOF
}

check_dependencies() {
    print_info "Checking dependencies..."
    
    # Check for doxygen
    if ! command -v doxygen &> /dev/null; then
        print_error "Doxygen is not installed!"
        echo "Please install Doxygen:"
        echo "  macOS:   brew install doxygen"
        echo "  Ubuntu:  sudo apt-get install doxygen"
        exit 1
    fi
    
    DOXYGEN_VERSION=$(doxygen --version)
    print_success "Doxygen found (version $DOXYGEN_VERSION)"
    
    # Check for dot (graphviz)
    if ! command -v dot &> /dev/null; then
        print_warning "Graphviz (dot) is not installed - graphs will not be generated"
        print_info "To install Graphviz:"
        echo "  macOS:   brew install graphviz"
        echo "  Ubuntu:  sudo apt-get install graphviz"
        QUICK_MODE=true
    else
        DOT_VERSION=$(dot -V 2>&1 | head -n1)
        print_success "Graphviz found ($DOT_VERSION)"
    fi
    
    echo ""
}

clean_documentation() {
    print_info "Cleaning previous documentation..."
    
    if [ -d "doxygen_output" ]; then
        rm -rf doxygen_output
        print_success "Removed doxygen_output/"
    fi
    
    if [ -f "doxygen_warnings.log" ]; then
        rm -f doxygen_warnings.log
        print_success "Removed doxygen_warnings.log"
    fi
    
    echo ""
}

prepare_doxyfile() {
    print_info "Preparing Doxyfile..."
    
    # Check if Doxyfile exists
    if [ ! -f "Doxyfile" ]; then
        print_error "Doxyfile not found in current directory!"
        exit 1
    fi
    
    # Create temporary Doxyfile if quick mode
    if [ "$QUICK_MODE" = true ]; then
        print_info "Quick mode: Disabling graphs for faster generation"
        cp Doxyfile Doxyfile.tmp
        sed -i.bak 's/HAVE_DOT.*=.*YES/HAVE_DOT = NO/' Doxyfile.tmp
        sed -i.bak 's/CALL_GRAPH.*=.*YES/CALL_GRAPH = NO/' Doxyfile.tmp
        sed -i.bak 's/CALLER_GRAPH.*=.*YES/CALLER_GRAPH = NO/' Doxyfile.tmp
        rm -f Doxyfile.tmp.bak
        DOXYFILE_TO_USE="Doxyfile.tmp"
    else
        DOXYFILE_TO_USE="Doxyfile"
    fi
    
    echo ""
}

generate_documentation() {
    print_info "Generating documentation..."
    echo ""
    
    START_TIME=$(date +%s)
    
    if [ "$VERBOSE" = true ]; then
        doxygen "$DOXYFILE_TO_USE"
    else
        doxygen "$DOXYFILE_TO_USE" 2>&1 | grep -E "(Generating|Parsing|Building|Warning|Error)" || true
    fi
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    
    echo ""
    print_success "Documentation generated in ${DURATION} seconds"
    
    # Clean up temporary file
    if [ "$QUICK_MODE" = true ]; then
        rm -f Doxyfile.tmp
    fi
}

show_statistics() {
    print_info "Documentation statistics:"
    
    OUTPUT_DIR="doxygen_output/html"
    
    if [ -d "$OUTPUT_DIR" ]; then
        NUM_HTML_FILES=$(find "$OUTPUT_DIR" -name "*.html" | wc -l | tr -d ' ')
        NUM_SVG_FILES=$(find "$OUTPUT_DIR" -name "*.svg" | wc -l | tr -d ' ')
        TOTAL_SIZE=$(du -sh "$OUTPUT_DIR" | cut -f1)
        
        echo "  📄 HTML files: $NUM_HTML_FILES"
        echo "  📊 SVG graphs: $NUM_SVG_FILES"
        echo "  💾 Total size: $TOTAL_SIZE"
    fi
    
    # Check for warnings
    if [ -f "doxygen_warnings.log" ]; then
        NUM_WARNINGS=$(wc -l < doxygen_warnings.log | tr -d ' ')
        if [ "$NUM_WARNINGS" -gt 0 ]; then
            print_warning "$NUM_WARNINGS warnings found (see doxygen_warnings.log)"
        else
            print_success "No warnings!"
        fi
    fi
    
    echo ""
}

open_documentation() {
    print_info "Opening documentation in browser..."
    
    INDEX_FILE="doxygen_output/html/index.html"
    
    if [ ! -f "$INDEX_FILE" ]; then
        print_error "Documentation index file not found!"
        exit 1
    fi
    
    # Detect OS and open browser
    case "$(uname -s)" in
        Darwin*)
            open "$INDEX_FILE"
            ;;
        Linux*)
            if command -v xdg-open &> /dev/null; then
                xdg-open "$INDEX_FILE"
            elif command -v firefox &> /dev/null; then
                firefox "$INDEX_FILE" &
            else
                print_warning "Could not detect browser. Please open manually:"
                echo "  file://$INDEX_FILE"
            fi
            ;;
        *)
            print_warning "Unknown OS. Please open manually:"
            echo "  file://$INDEX_FILE"
            ;;
    esac
    
    print_success "Documentation opened in browser"
}

################################################################################
# Parse command line arguments
################################################################################

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -q|--quick)
            QUICK_MODE=true
            shift
            ;;
        -o|--open)
            OPEN_BROWSER=true
            shift
            ;;
        -c|--clean)
            CLEAN_FIRST=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        *)
            print_error "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

################################################################################
# Main execution
################################################################################

print_header

# Check dependencies
check_dependencies

# Clean if requested
if [ "$CLEAN_FIRST" = true ]; then
    clean_documentation
fi

# Prepare Doxyfile
prepare_doxyfile

# Generate documentation
generate_documentation

# Show statistics
show_statistics

# Open in browser if requested
if [ "$OPEN_BROWSER" = true ]; then
    open_documentation
fi

# Final message
echo -e "${GREEN}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║${NC}  ✅ Documentation generation complete!                         ${GREEN}║${NC}"
echo -e "${GREEN}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "📖 View documentation at: doxygen_output/html/index.html"
echo ""

exit 0