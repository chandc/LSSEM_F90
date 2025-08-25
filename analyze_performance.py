#!/usr/bin/env python3
"""
SEM Performance Analysis Tool
Advanced performance comparison and visualization for SEM solver optimizations
"""

import os
import sys
import subprocess
import time
import json
import csv
from pathlib import Path
import argparse
from typing import Dict, List, Tuple, Optional
import statistics

try:
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    print("Warning: matplotlib/pandas not available. Plotting disabled.")

class SEMBenchmark:
    def __init__(self, results_dir: str = "benchmark_results"):
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(exist_ok=True)
        self.timestamp = time.strftime("%Y%m%d_%H%M%S")
        
        # Test configurations
        self.test_cases = [
            {"file": "input_36_7_Re100.nml", "name": "Re100", "description": "Low Reynolds Number"},
            {"file": "input_36_7_Re1000.nml", "name": "Re1000", "description": "High Reynolds Number"}
        ]
        
        # Executables to test
        self.executables = {
            "baseline": {"file": "SEM_2D_Baseline", "name": "Baseline", "color": "red"},
            "blas": {"file": "SEM_2D_BLAS", "name": "BLAS-Optimized", "color": "blue"},
            "openmp": {"file": "SEM_2D_OpenMP", "name": "OpenMP", "color": "green"},
            "hybrid": {"file": "SEM_2D_Hybrid", "name": "BLAS+OpenMP", "color": "purple"}
        }
        
        self.results = {}
    
    def run_single_test(self, executable: str, input_file: str, num_runs: int = 5) -> Dict:
        """Run a single benchmark test multiple times and collect statistics."""
        if not Path(executable).exists():
            raise FileNotFoundError(f"Executable {executable} not found")
        
        if not Path(input_file).exists():
            raise FileNotFoundError(f"Input file {input_file} not found")
        
        wall_times = []
        user_times = []
        sys_times = []
        memory_usage = []
        
        print(f"  Running {num_runs} iterations...", end="", flush=True)
        
        for run in range(num_runs):
            # Run with time measurement
            cmd = ["/usr/bin/time", "-f", "%e %U %S %M", f"./{executable}"]
            
            with open(input_file, 'r') as input_f:
                try:
                    result = subprocess.run(
                        cmd, 
                        stdin=input_f,
                        capture_output=True,
                        text=True,
                        timeout=300  # 5 minute timeout
                    )
                    
                    if result.returncode == 0:
                        # Parse timing output
                        timing_line = result.stderr.strip().split('\n')[-1]
                        times = timing_line.split()
                        
                        if len(times) >= 4:
                            wall_times.append(float(times[0]))
                            user_times.append(float(times[1]))
                            sys_times.append(float(times[2]))
                            memory_usage.append(float(times[3]))
                            print(".", end="", flush=True)
                        else:
                            print("E", end="", flush=True)  # Error parsing
                    else:
                        print("F", end="", flush=True)  # Failed run
                        
                except subprocess.TimeoutExpired:
                    print("T", end="", flush=True)  # Timeout
                except Exception as e:
                    print("X", end="", flush=True)  # Exception
        
        print(" done")
        
        if not wall_times:
            return None
        
        return {
            'wall_time': {
                'mean': statistics.mean(wall_times),
                'median': statistics.median(wall_times),
                'stdev': statistics.stdev(wall_times) if len(wall_times) > 1 else 0,
                'min': min(wall_times),
                'max': max(wall_times),
                'values': wall_times
            },
            'user_time': {
                'mean': statistics.mean(user_times),
                'median': statistics.median(user_times),
                'stdev': statistics.stdev(user_times) if len(user_times) > 1 else 0
            },
            'sys_time': {
                'mean': statistics.mean(sys_times),
                'median': statistics.median(sys_times),
                'stdev': statistics.stdev(sys_times) if len(sys_times) > 1 else 0
            },
            'memory': {
                'mean': statistics.mean(memory_usage),
                'median': statistics.median(memory_usage),
                'stdev': statistics.stdev(memory_usage) if len(memory_usage) > 1 else 0
            },
            'num_runs': len(wall_times)
        }
    
    def run_all_benchmarks(self, num_runs: int = 5):
        """Run all benchmark combinations."""
        print("Starting comprehensive benchmark suite...")
        
        for exec_key, exec_info in self.executables.items():
            executable = exec_info["file"]
            
            if not Path(executable).exists():
                print(f"Skipping {exec_info['name']} - executable not found")
                continue
            
            print(f"\nBenchmarking {exec_info['name']}:")
            self.results[exec_key] = {}
            
            for test_case in self.test_cases:
                input_file = test_case["file"]
                case_name = test_case["name"]
                
                if not Path(input_file).exists():
                    print(f"  Skipping {case_name} - input file not found")
                    continue
                
                print(f"  Test case: {case_name}")
                
                try:
                    result = self.run_single_test(executable, input_file, num_runs)
                    if result:
                        self.results[exec_key][case_name] = result
                        print(f"    Wall time: {result['wall_time']['mean']:.3f}s ± {result['wall_time']['stdev']:.3f}s")
                    else:
                        print(f"    Failed to collect timing data")
                except Exception as e:
                    print(f"    Error: {e}")
    
    def save_results(self):
        """Save results to JSON and CSV files."""
        # Save JSON
        json_file = self.results_dir / f"benchmark_results_{self.timestamp}.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        # Save CSV
        csv_file = self.results_dir / f"benchmark_summary_{self.timestamp}.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Executable', 'TestCase', 'WallTime_Mean', 'WallTime_StDev',
                'UserTime_Mean', 'SysTime_Mean', 'Memory_Mean', 'NumRuns'
            ])
            
            for exec_key, exec_results in self.results.items():
                exec_name = self.executables[exec_key]["name"]
                for case_name, case_results in exec_results.items():
                    writer.writerow([
                        exec_name, case_name,
                        case_results['wall_time']['mean'],
                        case_results['wall_time']['stdev'],
                        case_results['user_time']['mean'],
                        case_results['sys_time']['mean'],
                        case_results['memory']['mean'],
                        case_results['num_runs']
                    ])
        
        print(f"\nResults saved:")
        print(f"  JSON: {json_file}")
        print(f"  CSV:  {csv_file}")
    
    def generate_performance_report(self):
        """Generate detailed performance analysis report."""
        report_file = self.results_dir / f"performance_report_{self.timestamp}.txt"
        
        with open(report_file, 'w') as f:
            f.write("SEM Performance Analysis Report\n")
            f.write("=" * 50 + "\n")
            f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # System information
            f.write("System Information:\n")
            f.write("-" * 20 + "\n")
            try:
                import platform
                f.write(f"OS: {platform.system()} {platform.release()}\n")
                f.write(f"CPU: {platform.processor()}\n")
                f.write(f"Python: {platform.python_version()}\n")
            except:
                f.write("System info not available\n")
            f.write("\n")
            
            # Performance comparison
            if 'baseline' in self.results:
                baseline_results = self.results['baseline']
                
                for case_name in baseline_results:
                    f.write(f"Test Case: {case_name}\n")
                    f.write("-" * 30 + "\n")
                    
                    baseline_time = baseline_results[case_name]['wall_time']['mean']
                    baseline_memory = baseline_results[case_name]['memory']['mean']
                    
                    f.write(f"{'Method':<15} {'Time (s)':<10} {'Speedup':<8} {'Memory (KB)':<12} {'Efficiency':<10}\n")
                    f.write("-" * 65 + "\n")
                    
                    for exec_key, exec_results in self.results.items():
                        if case_name in exec_results:
                            exec_name = self.executables[exec_key]["name"]
                            wall_time = exec_results[case_name]['wall_time']['mean']
                            memory = exec_results[case_name]['memory']['mean']
                            
                            speedup = baseline_time / wall_time if wall_time > 0 else 0
                            efficiency = speedup  # Could be adjusted for number of cores
                            
                            f.write(f"{exec_name:<15} {wall_time:<10.3f} {speedup:<8.2f} {memory:<12.0f} {efficiency:<10.2f}\n")
                    
                    f.write("\n")
            
            # Statistical summary
            f.write("Statistical Summary:\n")
            f.write("-" * 20 + "\n")
            for exec_key, exec_results in self.results.items():
                exec_name = self.executables[exec_key]["name"]
                f.write(f"\n{exec_name}:\n")
                
                for case_name, case_results in exec_results.items():
                    wall_stats = case_results['wall_time']
                    f.write(f"  {case_name}: {wall_stats['mean']:.3f}s ± {wall_stats['stdev']:.3f}s "
                           f"(min: {wall_stats['min']:.3f}s, max: {wall_stats['max']:.3f}s)\n")
        
        print(f"Performance report: {report_file}")
        return report_file
    
    def plot_results(self):
        """Generate performance visualization plots."""
        if not HAS_PLOTTING:
            print("Plotting libraries not available. Skipping visualization.")
            return
        
        if not self.results:
            print("No results to plot.")
            return
        
        # Prepare data for plotting
        plot_data = {}
        for exec_key, exec_results in self.results.items():
            exec_name = self.executables[exec_key]["name"]
            plot_data[exec_name] = {}
            for case_name, case_results in exec_results.items():
                plot_data[exec_name][case_name] = case_results['wall_time']['mean']
        
        # Create plots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Bar plot comparison
        test_cases = list(next(iter(plot_data.values())).keys())
        x = np.arange(len(test_cases))
        width = 0.2
        
        for i, (exec_name, data) in enumerate(plot_data.items()):
            times = [data.get(case, 0) for case in test_cases]
            color = self.executables.get(
                next((k for k, v in self.executables.items() if v["name"] == exec_name), ""),
                {"color": f"C{i}"}
            )["color"]
            ax1.bar(x + i * width, times, width, label=exec_name, color=color)
        
        ax1.set_xlabel('Test Cases')
        ax1.set_ylabel('Wall Time (seconds)')
        ax1.set_title('Performance Comparison by Test Case')
        ax1.set_xticks(x + width * (len(plot_data) - 1) / 2)
        ax1.set_xticklabels(test_cases)
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Speedup plot (relative to baseline)
        if 'Baseline' in plot_data:
            baseline_data = plot_data['Baseline']
            
            for exec_name, data in plot_data.items():
                if exec_name != 'Baseline':
                    speedups = []
                    for case in test_cases:
                        if case in baseline_data and case in data and data[case] > 0:
                            speedup = baseline_data[case] / data[case]
                            speedups.append(speedup)
                        else:
                            speedups.append(0)
                    
                    color = self.executables.get(
                        next((k for k, v in self.executables.items() if v["name"] == exec_name), ""),
                        {"color": "gray"}
                    )["color"]
                    ax2.plot(test_cases, speedups, 'o-', label=exec_name, color=color, linewidth=2, markersize=8)
            
            ax2.set_xlabel('Test Cases')
            ax2.set_ylabel('Speedup Factor')
            ax2.set_title('Speedup Relative to Baseline')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            ax2.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        plot_file = self.results_dir / f"performance_comparison_{self.timestamp}.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Performance plot: {plot_file}")
        
        # Show plot if in interactive mode
        if hasattr(sys, 'ps1'):  # Interactive mode
            plt.show()
        
        plt.close()
    
    def analyze_scaling(self):
        """Analyze performance scaling characteristics."""
        if not self.results:
            return
        
        print("\nPerformance Scaling Analysis:")
        print("=" * 40)
        
        if 'baseline' in self.results:
            baseline_results = self.results['baseline']
            
            for case_name in baseline_results:
                print(f"\nTest Case: {case_name}")
                print("-" * 20)
                
                baseline_time = baseline_results[case_name]['wall_time']['mean']
                
                for exec_key, exec_results in self.results.items():
                    if exec_key != 'baseline' and case_name in exec_results:
                        exec_name = self.executables[exec_key]["name"]
                        wall_time = exec_results[case_name]['wall_time']['mean']
                        speedup = baseline_time / wall_time if wall_time > 0 else 0
                        efficiency = speedup  # Could be normalized by core count
                        
                        print(f"{exec_name:15}: {speedup:.2f}x speedup, {efficiency:.2f} efficiency")


def main():
    parser = argparse.ArgumentParser(description="SEM Performance Benchmark Tool")
    parser.add_argument("--runs", type=int, default=5, help="Number of benchmark runs per test")
    parser.add_argument("--results-dir", default="benchmark_results", help="Results directory")
    parser.add_argument("--plot", action="store_true", help="Generate performance plots")
    parser.add_argument("--report-only", action="store_true", help="Generate report from existing results")
    parser.add_argument("--load-results", help="Load results from JSON file")
    
    args = parser.parse_args()
    
    benchmark = SEMBenchmark(args.results_dir)
    
    if args.load_results:
        # Load existing results
        with open(args.load_results, 'r') as f:
            benchmark.results = json.load(f)
        print(f"Loaded results from {args.load_results}")
    elif not args.report_only:
        # Run benchmarks
        benchmark.run_all_benchmarks(args.runs)
        benchmark.save_results()
    
    # Generate analysis
    if benchmark.results:
        benchmark.generate_performance_report()
        benchmark.analyze_scaling()
        
        if args.plot:
            benchmark.plot_results()
    else:
        print("No results available for analysis.")

if __name__ == "__main__":
    main()
