# SEM Performance Comparison Scripts

This directory contains comprehensive performance benchmarking tools for comparing different optimizations of the SEM (Spectral Element Method) solver.

## Scripts Overview

### 1. `benchmark_performance.sh` (Main Benchmark Script)
A comprehensive bash script that automatically builds, runs, and compares different versions of the SEM solver.

**Features:**
- Automated building of baseline, BLAS-optimized, and OpenMP versions
- System information collection (CPU, memory, OS, compiler versions)
- Multiple test runs for statistical accuracy
- Automatic timing measurement using `/usr/bin/time`
- Results validation and numerical comparison
- Optional profiling with gprof and valgrind
- CSV and text report generation

**Usage:**
```bash
# Basic benchmark (builds all versions and runs comparison)
./benchmark_performance.sh

# Include profiling analysis
./benchmark_performance.sh --profile

# Show help
./benchmark_performance.sh --help
```

### 2. `analyze_performance.py` (Advanced Analysis Tool)
A Python script for detailed performance analysis with statistical processing and visualization.

**Features:**
- Statistical analysis (mean, median, standard deviation)
- Performance scaling analysis
- Speedup calculations
- Automated plotting (if matplotlib available)
- JSON/CSV export for further analysis
- Load and reanalyze existing results

**Usage:**
```bash
# Run complete benchmark with 10 iterations per test
python3 analyze_performance.py --runs 10 --plot

# Generate report from existing results
python3 analyze_performance.py --report-only --load-results benchmark_results/benchmark_results_20240817_143022.json

# Create plots only
python3 analyze_performance.py --plot --load-results results.json
```

## Expected Performance Improvements

Based on the BLAS and OpenMP optimizations identified:

### BLAS Optimizations:
- **Derivative computations**: 2-5x speedup from O(n³) → O(n²) DGEMV operations
- **BiCGSTAB vector operations**: 1.5-3x speedup from vectorized DDOT, DAXPY, DNRM2
- **Memory access patterns**: 10-20% improvement from cache-optimized BLAS routines

### OpenMP Optimizations:
- **Element-level parallelism**: 2-4x speedup (scales with cores)
- **Vector operations**: Additional 1.5-2x speedup in BiCGSTAB loops

### Combined Expected Speedup:
- **Single-threaded BLAS**: 2-4x overall speedup
- **Multi-threaded BLAS + OpenMP**: 4-16x total speedup on modern multi-core systems

## Output Files

The benchmark scripts generate several output files in the `benchmark_results/` directory:

### Timing Data:
- `timing_data_TIMESTAMP.csv` - Raw timing data for all runs
- `benchmark_results_TIMESTAMP.json` - Complete results in JSON format
- `benchmark_TIMESTAMP.log` - Detailed execution log

### Analysis Reports:
- `performance_report_TIMESTAMP.txt` - Human-readable performance comparison
- `performance_comparison_TIMESTAMP.png` - Performance visualization plots

### Profiling Data (if enabled):
- `profile_TIMESTAMP.txt` - gprof function-level profiling
- `memory_profile_TIMESTAMP.txt` - Memory usage analysis
- `massif_TIMESTAMP.out` - Valgrind massif memory profiling data

## Validation

The scripts automatically validate numerical results by comparing output files between different versions:

1. **Exact comparison**: Checks if results are identical
2. **Numerical tolerance**: Uses relative difference analysis (if Python/NumPy available)
3. **Acceptable thresholds**: 
   - Machine precision: < 1e-12 relative difference
   - Acceptable tolerance: < 1e-6 relative difference

## System Requirements

### Required:
- Bash shell (for main benchmark script)
- `gfortran` or compatible Fortran compiler
- `/usr/bin/time` command (for accurate timing)
- BLAS/LAPACK library (Accelerate Framework on macOS, OpenBLAS/MKL on Linux)

### Optional (for enhanced analysis):
- Python 3.6+ with NumPy, matplotlib, pandas
- `gprof` (for function profiling)
- `valgrind` (for memory analysis - Linux only)

## Platform-Specific Notes

### macOS:
- Uses Apple Accelerate Framework for BLAS (highly optimized)
- System information collected via `sysctl` commands
- Memory profiling limited (no valgrind)

### Linux:
- Supports multiple BLAS implementations (OpenBLAS, MKL, ATLAS)
- Full profiling support including valgrind
- System information via `/proc` filesystem

## Troubleshooting

### Common Issues:

1. **BLAS library not found:**
   ```bash
   # macOS: Ensure Xcode command line tools installed
   xcode-select --install
   
   # Linux: Install BLAS/LAPACK
   sudo apt-get install libblas-dev liblapack-dev  # Ubuntu/Debian
   sudo yum install blas-devel lapack-devel        # RHEL/CentOS
   ```

2. **Compilation errors:**
   - Check Fortran compiler version and flags in Makefile
   - Verify all source files are present
   - Check for proper module dependencies

3. **Timing inconsistencies:**
   - Ensure system is not under heavy load during benchmarking
   - Disable CPU frequency scaling if possible
   - Run multiple iterations to account for variance

### Performance Validation:

Expected timing relationships:
```
Baseline > BLAS > OpenMP > BLAS+OpenMP (fastest)
```

If this relationship doesn't hold, check:
- Compiler optimization flags (-O3)
- BLAS library linking
- OpenMP thread count (should match core count)
- System thermal throttling

## Example Output

```
================================================
 SEM Performance Benchmark Suite
 Mon Aug 17 14:30:22 PDT 2025
================================================

System Information
----------------------------------------
  OS: macOS 14.0
  CPU: Apple M1 Pro
  Cores: 8
  Memory: 16 GB

Benchmarking Baseline
----------------------------------------
Test Case: Re1000 (input_36_7_Re1000.nml)
  Run 1/5... Wall: 12.456s, User: 12.234s, Memory: 45678KB
  ...
  Results:
    Average Wall Time: 12.123s
    Speedup: 1.00x (baseline)

Benchmarking BLAS-Optimized
----------------------------------------
Test Case: Re1000 (input_36_7_Re1000.nml)
  Run 1/5... Wall: 4.123s, User: 4.001s, Memory: 48234KB
  ...
  Results:
    Average Wall Time: 4.234s
    Speedup: 2.86x

Performance Analysis
----------------------------------------
BLAS-Optimized:
  Wall Time: 4.234s (vs 12.123s baseline)
  Speedup: 2.86x
  Memory: 48234KB (ratio: 1.06x)
```

This comprehensive benchmarking framework will help you quantify the performance improvements from BLAS optimization and guide further optimization efforts.
