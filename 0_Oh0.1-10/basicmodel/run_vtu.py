import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import sys
import shutil

# Parameters are defined within the script
tmax = 50
tsnap = 0.01
output_dir = "paraviewResults"
num_threads = 10

def compile_program():
    """Compiles the C program."""
    executable_name = "dumpToVTU"
    # Remove the executable if it exists
    if os.path.exists(executable_name):
        os.remove(executable_name)
        print(f"Removed existing executable: {executable_name}")

    # Compile the C program
    compile_command = "qcc -O2 -w -fopenmp -Wall dumpToVTU.c -lm -o dumpToVTU"
    print("Compiling the C program...")
    compile_result = subprocess.run(compile_command, shell=True, text=True)
    if compile_result.returncode != 0:
        print("Compilation failed.")
        return False
    print("Compilation successful.")
    return True

def run_process(t_real):
    """Runs the compiled C program with given parameters."""
    command = f"./dumpToVTU {t_real} {tsnap} {output_dir}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

def main():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(current_dir)   
    if not compile_program():
        return 1

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)  # Removes the directory and all contents
        print(f"Deleted existing directory: {output_dir}")
    os.makedirs(output_dir)  # Recreate the directory

    # Set environment variable for OpenMP
    os.environ['OMP_NUM_THREADS'] = str(num_threads)
    print(f"Set OMP_NUM_THREADS to {num_threads}.")

    # Create a list of time steps
    t_reals = [i * tsnap for i in range(int(tmax / tsnap) + 1)]

    # Run the C program in parallel for each time step
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        executor.map(run_process, t_reals)

    print("All tasks completed successfully.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
