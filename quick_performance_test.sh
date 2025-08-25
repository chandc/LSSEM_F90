#!/bin/bash
#**********************************************************************
# Quick Performance Test - SEM BLAS vs Baseline Comparison
# Simple script to quickly test BLAS performance improvements
#**********************************************************************

set -e

echo "=============================================="
echo " SEM Quick Performance Test"
echo " BLAS vs Baseline Comparison"
echo "=============================================="

# Check if input files exist
if [ ! -f "input_36_7_Re1000.nml" ]; then
    echo "ERROR: input_36_7_Re1000.nml not found"
    echo "Please ensure you're in the correct directory with input files"
    exit 1
fi

echo ""
echo "Building executables..."

# Build baseline version
if [ -f "Makefile_F90_baseline" ]; then
    echo "  Building baseline version..."
    if make -f Makefile_F90_baseline clean >/dev/null 2>&1 && make -f Makefile_F90_baseline >/dev/null 2>&1; then
        echo "  ✓ Baseline build successful"
    else
        echo "  ✗ Baseline build failed"
        exit 1
    fi
else
    echo "  ! Makefile_F90_baseline not found, skipping baseline"
fi

# Build BLAS version 
if [ -f "Makefile_BLAS" ]; then
    echo "  Building BLAS version..."
    if make -f Makefile_BLAS clean >/dev/null 2>&1 && make -f Makefile_BLAS >/dev/null 2>&1; then
        echo "  ✓ BLAS build successful"
    else
        echo "  ✗ BLAS build failed"
        exit 1
    fi
else
    echo "  ! Makefile_BLAS not found, creating basic BLAS makefile..."
    # Create a basic BLAS makefile if it doesn't exist
    cat > Makefile_BLAS << 'EOF'
FC = gfortran
FFLAGS = -O3 -march=native -fdefault-real-8 -fdefault-double-8 -ffree-form
BLAS_LIBS = -framework Accelerate

SRC_FILES = sem_data.f90 SEM_08.f90 lssem.f90 solver.f90 lgl.f90
OBJ_FILES = sem_data.o SEM_08.o lssem.o solver.o lgl.o
EXECUTABLE = SEM_2D_BLAS

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ_FILES)
	$(FC) $(FFLAGS) -o $@ $^ $(BLAS_LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

SEM_08.o: sem_data.o
lssem.o: sem_data.o
solver.o: sem_data.o

clean:
	rm -f *.o *.mod

.PHONY: all clean
EOF
    echo "  Created basic Makefile_BLAS (note: this uses existing source files, not BLAS-optimized versions)"
    if make -f Makefile_BLAS clean >/dev/null 2>&1 && make -f Makefile_BLAS >/dev/null 2>&1; then
        echo "  ✓ Basic BLAS build successful"
    else
        echo "  ✗ Basic BLAS build failed"
        exit 1
    fi
fi

echo ""
echo "Running performance comparison..."

# Function to run a single test
run_test() {
    local executable=$1
    local name=$2
    local output_file="quick_test_${name}.dat"
    
    if [ ! -f "$executable" ]; then
        echo "  ✗ $executable not found"
        return 1
    fi
    
    echo "  Testing $name..."
    
    # Run with time measurement using bash builtin time
    local timing_output=$(mktemp)
    if { time ./"$executable" < input_36_7_Re1000.nml > "$output_file"; } 2> "$timing_output"; then
        local timing_data=$(cat "$timing_output")
        # Extract wall time (format: "real 0m1.234s")
        local wall_time=$(echo "$timing_data" | grep real | awk '{print $2}' | sed 's/[ms]//g' | awk -F: '{if(NF>1) print $1*60+$2; else print $0}')
        
        echo "    Wall time: ${wall_time}s"
        echo "    Status: Success"
        
        rm -f "$timing_output"
        echo "$wall_time"
        return 0
    else
        echo "    ✗ Test failed"
        rm -f "$timing_output"
        return 1
    fi
}

# Run tests
baseline_time=""
blas_time=""

if [ -f "SEM_2D_F90_Baseline" ]; then
    baseline_time=$(run_test "SEM_2D_F90_Baseline" "Baseline")
fi

if [ -f "SEM_2D_BLAS" ]; then
    blas_time=$(run_test "SEM_2D_BLAS" "BLAS")
fi

echo ""
echo "=============================================="
echo " Performance Summary"
echo "=============================================="

if [ -n "$baseline_time" ] && [ -n "$blas_time" ]; then
    # Calculate speedup
    speedup=$(echo "$baseline_time $blas_time" | awk '{if($2>0) print $1/$2; else print "N/A"}')
    
    echo "Baseline time:    ${baseline_time}s"
    echo "BLAS time:        ${blas_time}s"
    echo "Speedup:          ${speedup}x"
    
    # Performance assessment
    if (( $(echo "$speedup > 2.0" | bc -l) )); then
        echo "Result:           🚀 Excellent performance improvement!"
    elif (( $(echo "$speedup > 1.5" | bc -l) )); then
        echo "Result:           ✅ Good performance improvement"
    elif (( $(echo "$speedup > 1.1" | bc -l) )); then
        echo "Result:           ➕ Modest performance improvement"
    else
        echo "Result:           ⚠️  Limited improvement - check BLAS implementation"
    fi
    
    # Quick validation
    if [ -f "quick_test_Baseline.dat" ] && [ -f "quick_test_BLAS.dat" ]; then
        if diff -q quick_test_Baseline.dat quick_test_BLAS.dat >/dev/null 2>&1; then
            echo "Validation:       ✅ Results match exactly"
        else
            echo "Validation:       ⚠️  Results differ (check numerical accuracy)"
        fi
    fi
    
elif [ -n "$baseline_time" ]; then
    echo "Baseline time:    ${baseline_time}s"
    echo "BLAS time:        Not available"
    echo "Result:           ℹ️  Only baseline test completed"
    
elif [ -n "$blas_time" ]; then
    echo "Baseline time:    Not available"
    echo "BLAS time:        ${blas_time}s"
    echo "Result:           ℹ️  Only BLAS test completed"
    
else
    echo "No successful tests completed"
    echo "Please check build configuration and input files"
fi

echo ""
echo "For detailed analysis, run:"
echo "  ./benchmark_performance.sh"
echo "  python3 analyze_performance.py --plot"

# Clean up test output files
rm -f quick_test_*.dat

echo ""
echo "Quick test complete!"
