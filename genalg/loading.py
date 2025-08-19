import pygad
import matplotlib.pyplot as plt # Often useful for custom plotting

# 1. Define the filename
filename = "ga_instance"

# 2. Load the saved GA instance from the file
print(f"Loading GA instance from '{filename}'...")
loaded_ga_instance = pygad.load(filename=filename)
print("Load complete.")

# 3. Now you can access all the data and methods from the completed run
print("\n--- Best Solution ---")
solution, solution_fitness, solution_idx = loaded_ga_instance.best_solution()
print(f"Parameters (AoA, Scale): {solution}")
print(f"Fitness (L/D Ratio): {solution_fitness}")

# You can also re-plot the fitness graph
print("\nRegenerating fitness plot...")
loaded_ga_instance.plot_fitness(title="Saved GA Fitness Progression",
                                save_dir="fitness_plot_from_load.png")
print("Plot saved to 'fitness_plot_from_load.png'")

# You can inspect other properties, like the fitness value of the best
# solution at each generation
print("\n--- Fitness Progression per Generation ---")
print(loaded_ga_instance.best_solutions_fitness)