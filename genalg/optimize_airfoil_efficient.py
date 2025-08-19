import pygad
import os
import shutil
import subprocess
import re
import time
from dotenv import set_key, find_dotenv, load_dotenv
import csv

# --- Import functions from our refactored case runner ---
import case_runner

# --- DEFINE THE SCRIPT'S ABSOLUTE DIRECTORY ---
# This is the single source of truth for our case path.
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# --- Load .env and find its path ---
load_dotenv()
ENV_FILE_PATH = find_dotenv()
MESH_CACHE_DIR = os.getenv("MESH_CACHE_DIR")

def wait_for_file(filepath, timeout=10):
    """
    Waits for a file to exist for a maximum number of seconds (timeout).
    Returns True if the file is found in time, False otherwise.
    """
    start_time = time.time()
    while not os.path.exists(filepath):
        if time.time() - start_time > timeout:
            print(f"  [Error] Timeout: Waited {timeout}s for file '{filepath}' but it was not found.")
            return False
        time.sleep(0.5) # Check every half-second
    print(f"  File '{filepath}' found. Proceeding.")
    return True

def parse_coeffs():
    """Parses the coefficient.dat file using an absolute path."""
    try:
        # Use the absolute path to ensure we're looking in the right place.
        coeffs_file = os.path.join(SCRIPT_DIR, "postProcessing", "forceCoeffs", "0", "coefficient.dat")
        with open(coeffs_file, "r") as f:
            lines = f.readlines()
            last_line = ""
            for line in reversed(lines):
                if not line.startswith("#"):
                    last_line = line
                    break
            if not last_line: return None, None
            parts = re.split(r'\s+', last_line.strip())
            cd = float(parts[2])
            cl = float(parts[3])
            return cl, cd
    except (FileNotFoundError, IndexError, ValueError):
        return None, None

def fitness_func(ga_instance, solution, solution_idx):
    """Fitness function with explicit directory control."""
    aoa, scale = solution
    print(f"\n--- Evaluating Solution {solution_idx}: AoA = {aoa:.2f} deg, Scale = {scale:.4f} ---")

    # 1. Update .env and module variables
    set_key(ENV_FILE_PATH, "ANGLE_OF_ATTACK", str(aoa))
    set_key(ENV_FILE_PATH, "AIRFOIL_SCALE", str(scale))
    load_dotenv(override=True)
    case_runner.ANGLE_OF_ATTACK = float(os.getenv("ANGLE_OF_ATTACK"))
    case_runner.AIRFOIL_SCALE = float(os.getenv("AIRFOIL_SCALE"))
    case_runner.INLET_VELOCITY = float(os.getenv("INLET_VELOCITY"))
    print(f"  Confirmed: AoA = {case_runner.ANGLE_OF_ATTACK:.2f}°, Scale = {case_runner.AIRFOIL_SCALE:.4f}")

    # 2. MESHING STAGE (with caching)
    scale_id = f"scale_{scale:.1f}"
    cached_mesh_path = os.path.join(MESH_CACHE_DIR, scale_id)
    if os.path.exists(cached_mesh_path):
        # CACHE HIT
        print(f"  CACHE HIT: Found existing mesh for scale {scale:.4f}.")
        # Run cleaning in the correct directory
        subprocess.run(["foamCleanTutorials"], cwd=SCRIPT_DIR, capture_output=True)
        case_runner.prepare_directories()
        
        # --- FIX: Add dirs_exist_ok=True to prevent error if directory exists ---
        shutil.copytree(os.path.join(cached_mesh_path, "polyMesh"), os.path.join(SCRIPT_DIR, "constant/polyMesh"), dirs_exist_ok=True)
        shutil.copytree(os.path.join(cached_mesh_path, "triSurface"), os.path.join(SCRIPT_DIR, "constant/triSurface"), dirs_exist_ok=True)
        # -------------------------------------------------------------------------

        print("  Writing solver configuration files for cache hit...")
        case_runner.write_controlDict()
        case_runner.write_fvSchemes()
        case_runner.write_fvSolution()
        case_runner.write_transportProperties()
        case_runner.write_turbulenceProperties()
        case_runner.write_decomposeParDict()
    else:
        # CACHE MISS
        print(f"  CACHE MISS: No mesh found for scale {scale:.4f}. Generating new mesh...")
        # Run cleaning in the correct directory
        subprocess.run(["foamCleanTutorials"], cwd=SCRIPT_DIR, capture_output=True)
        case_runner.prepare_directories()
        coords = case_runner.read_and_scale_airfoil(case_runner.AIRFOIL_FILE, scale)
        case_runner.create_stl(coords, case_runner.STL_OUTPUT)
        # Write all config files
        case_runner.write_controlDict()
        case_runner.write_blockMeshDict()
        case_runner.write_surfaceFeatureExtractDict()
        case_runner.write_decomposeParDict()
        case_runner.write_snappyHexMeshDict()
        case_runner.write_fvSchemes()
        case_runner.write_fvSolution()
        case_runner.write_transportProperties()
        case_runner.write_turbulenceProperties()
        case_runner.write_boundary_conditions()
        # Run meshing in the correct directory
        mesh_success = case_runner.run_meshing_workflow(cwd=SCRIPT_DIR)
        if not mesh_success:
            print("  !!! Meshing failed. Assigning poor fitness. !!!")
            return -1000.0
        print(f"  Saving new mesh to cache: {cached_mesh_path}")
        os.makedirs(cached_mesh_path, exist_ok=True)
        shutil.copytree(os.path.join(SCRIPT_DIR, "constant/polyMesh"), os.path.join(cached_mesh_path, "polyMesh"), dirs_exist_ok=True)
        shutil.copytree(os.path.join(SCRIPT_DIR, "constant/triSurface"), os.path.join(cached_mesh_path, "triSurface"), dirs_exist_ok=True)

    # 3. SOLVER STAGE
    case_runner.write_boundary_conditions()
    # Run solver in the correct directory
    solver_success = case_runner.run_solver_workflow(cwd=SCRIPT_DIR)
    if not solver_success:
        print("  !!! Solver failed. Assigning poor fitness. !!!")
        return -1000.0

    # 4. RESULTS & FITNESS CALCULATION
    print("  Solver finished. Waiting for results file to be written...")
    # Use the absolute path to wait for the file
    results_filepath = os.path.join(SCRIPT_DIR, "postProcessing", "forceCoeffs", "0", "coefficient.dat")
    if not wait_for_file(results_filepath):
        print("  !!! Results file not found after waiting. Assigning poor fitness. !!!")
        return -1000.0

    cl, cd = parse_coeffs()
    if cl is None or cd is None:
        print("  !!! Failed to parse results. Assigning poor fitness. !!!")
        return -1000.0

    fitness = cl / cd if cd > 1e-6 else 0.0
    set_key(ENV_FILE_PATH, "LAST_CL", f"{cl:.5f}")
    set_key(ENV_FILE_PATH, "LAST_CD", f"{cd:.5f}")
    set_key(ENV_FILE_PATH, "LAST_CL_CD_RATIO", f"{fitness:.5f}")
    print(f"--- Result: Cl={cl:.4f}, Cd={cd:.4f}, Fitness (L/D) = {fitness:.4f} ---")
    return fitness

def on_gen_callback(ga_instance):
    """PyGAD callback function to check for simulation failures and stop the GA."""
    last_fitness = ga_instance.last_generation_fitness
    if -1000.0 in last_fitness:
        print("\n!!! A simulation failed. Stopping the genetic algorithm. !!!")
        ga_instance.stop = True  # Set our custom flag for the final check
        return "stop"            # This is the correct way to terminate the run


if __name__ == '__main__':
    os.makedirs(MESH_CACHE_DIR, exist_ok=True)
    gene_space = [
        {'low': -10.0, 'high': 10.0, 'step': 0.5},
        {'low': 0.5, 'high': 1.5, 'step': 0.1}
    ]
    ga_instance = pygad.GA(
        num_generations=20,
        num_parents_mating=4,
        sol_per_pop=8,
        num_genes=2,
        fitness_func=fitness_func,
        gene_space=gene_space,
        crossover_type="single_point",
        mutation_type="random",
        mutation_percent_genes=40,
        on_generation=on_gen_callback
    )

    # --- FIX 1: Initialize the custom attribute before running ---
    ga_instance.stop = False

    ga_instance.run()
    
    # This check will now work correctly in all scenarios
    if ga_instance.stop:
        print("\n================= OPTIMIZATION TERMINATED =================")
        print("Optimization was stopped due to a simulation failure.")
        print("=======================================================")
    else:
        solution, solution_fitness, solution_idx = ga_instance.best_solution()
        print("\n================= OPTIMIZATION COMPLETE =================")
        print(f"Best solution: AoA = {solution[0]:.4f}, Scale = {solution[1]:.4f}")
        print(f"Highest L/D Ratio (Fitness): {solution_fitness:.4f}")
        print("=======================================================")
        
        # --- NEW CODE FOR PLOTTING AND SAVING ---

        # 1. Save the plot to a file instead of showing it
        ga_instance.plot_fitness(title="GA Fitness Progression (L/D Ratio vs. Generation)",
                                 save_dir="fitness_plot.png")
        print("\nFitness progression plot saved to 'fitness_plot.png'")

        # 2. Save the full GA object for later analysis
        ga_instance.save(filename="ga_instance.pkl")
        print("Full GA instance saved to 'ga_instance.pkl'")

        # 3. Save the generation-by-generation data to a CSV file
        best_solutions_data = []
        for generation, sol, fit in zip(range(ga_instance.generations_completed),
                                        ga_instance.best_solutions,
                                        ga_instance.best_solutions_fitness):
            best_solutions_data.append([generation + 1, fit, sol[0], sol[1]])

        csv_filename = "ga_progression_data.csv"
        with open(csv_filename, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Generation", "Fitness_LD_Ratio", "AoA", "Scale"]) # Header
            writer.writerows(best_solutions_data)

        print(f"Progression data saved to '{csv_filename}'")
