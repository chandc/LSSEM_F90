#!/bin/bash

# Script to run OpenMP scaling benchmark, record execution times, and verify output.

EXECUTABLE="./SEM_2D_F90_OMP"
INPUT_FILE="input_36_7_Re1000.nml"
OUTPUT_CSV="scaling_results.csv"
MAX_THREADS=16
BASELINE_OUTPUT_FILE="cavity_run_Re1000_1_threads.dat"

# Check if the executable exists before starting
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable '$EXECUTABLE' not found."
    echo "Please compile the OpenMP version first using 'make -f Makefile_OMP'."
    exit 1
fi

# Print the header to the CSV file
echo "threads,real_time_s" > "$OUTPUT_CSV"

echo "Starting OpenMP scaling benchmark for '$EXECUTABLE'..."
echo "Results will be saved to '$OUTPUT_CSV'."

# Loop from 1 to the maximum number of threads
for (( threads=1; threads<=MAX_THREADS; threads++ ))
do
    echo "  Running with $threads thread(s)..."
    export OMP_NUM_THREADS=$threads
    
    OUTPUT_FILE="cavity_run_Re1000_${threads}_threads.dat"
    
    # Use /usr/bin/time -p for stable, parsable output.
    # Redirect stderr (where time output goes) to a temporary file.
    TIME_LOG="time.log"
    /usr/bin/time -p "$EXECUTABLE" < "$INPUT_FILE" > "$OUTPUT_FILE" 2> "$TIME_LOG"
    
    # Extract the 'real' time value from the log file
    REAL_TIME_S=$(grep real "$TIME_LOG" | awk '{print $2}')
    
    # Append the thread count and the total seconds to our CSV file
    echo "$threads,$REAL_TIME_S" >> "$OUTPUT_CSV"
    
    rm "$TIME_LOG"
done

echo "Benchmark finished successfully."
echo ""
echo "Verifying output files against the 1-thread baseline..."

# Verification step
for (( threads=2; threads<=MAX_THREADS; threads++ ))
do
    CURRENT_OUTPUT_FILE="cavity_run_Re1000_${threads}_threads.dat"
    echo "  Comparing: $BASELINE_OUTPUT_FILE vs $CURRENT_OUTPUT_FILE"
    
    # Use diff to compare. The -q flag makes it quiet unless there's a difference.
    diff -q "$BASELINE_OUTPUT_FILE" "$CURRENT_OUTPUT_FILE"
    
    if [ $? -eq 0 ]; then
        # If diff returns 0, files are the same. We can remove the redundant file.
        rm "$CURRENT_OUTPUT_FILE"
    else
        # If files are different, report it.
        echo "    !! WARNING: Output for $threads threads differs from baseline."
    fi
done

# If all diffs were successful, only the baseline output file will remain.
if [ -f "cavity_run_Re1000_2_threads.dat" ]; then
     echo "Verification failed. Some output files were different."
else
     echo "Verification successful! All multi-threaded runs produced identical output to the single-thread run."
     # Rename the baseline for clarity
     mv "$BASELINE_OUTPUT_FILE" "cavity_run_36_7_Re1000_OMP_verified.dat"
     echo "Final verified output saved as: cavity_run_36_7_Re1000_OMP_verified.dat"
fi

