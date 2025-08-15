import numpy as np
import os
import re

def extract_cl_cd_from_openfoam(case_directory='.'):
    """
    Extract Cl and Cd from OpenFOAM forceCoeffs function object
    
    Returns:
    dict: Contains time series data and final converged values
    """
    
    # Try forceCoeffs first (direct coefficients)
    coeff_file = os.path.join(case_directory, 'postProcessing/forceCoeffs/0/coefficient.dat')
    
    if os.path.exists(coeff_file):
        return extract_from_forceCoeffs(coeff_file)
    else:
        # Fallback to forces file and manual calculation
        force_file = os.path.join(case_directory, 'postProcessing/forces/0/forces.dat')
        if os.path.exists(force_file):
            print("forceCoeffs file not found, using forces.dat for manual calculation")
            return extract_from_forces(force_file, case_directory)
        else:
            print("Error: Neither forceCoeffs nor forces files found!")
            print("Make sure you have run the simulation with the updated controlDict")
            return None


def extract_from_forceCoeffs(coeff_file):
    """
    Extract coefficients directly from forceCoeffs function object output
    FIXED: Correct column indices for OpenFOAM forceCoeffs output format
    """
    try:
        data = []
        with open(coeff_file, 'r') as f:
            for line in f:
                if not line.startswith('#') and line.strip():
                    values = line.split()
                    if len(values) >= 12:  # Updated minimum columns needed
                        data.append([float(v) for v in values])
        
        if not data:
            print("Error: No valid data found in forceCoeffs file")
            return None
        
        data = np.array(data)
        
        # FIXED COLUMN INDICES based on OpenFOAM forceCoeffs format:
        # Time  Cd  Cd(f)  Cd(r)  Cl  Cl(f)  Cl(r)  CmPitch  CmRoll  CmYaw  Cs  Cs(f)  Cs(r)
        #   0   1    2      3     4    5      6       7        8       9     10   11     12
        
        time = data[:, 0]
        Cd = data[:, 1]     # Total drag coefficient
        Cl = data[:, 4]     # Total lift coefficient (FIXED!)
        Cs = data[:, 10]    # Total side force coefficient
        CmPitch = data[:, 7] # Pitching moment coefficient
        
        # Also extract sub-components for analysis
        Cd_friction = data[:, 2]  # Cd(f) - friction drag
        Cd_pressure = data[:, 3]  # Cd(r) - pressure drag  
        Cl_friction = data[:, 5]  # Cl(f) - friction lift
        Cl_pressure = data[:, 6]  # Cl(r) - pressure lift
        
        # Get final converged values (last 100 time steps)
        n_avg = min(100, len(time))
        Cd_final = np.mean(Cd[-n_avg:])
        Cl_final = np.mean(Cl[-n_avg:])
        Cd_std = np.std(Cd[-n_avg:])
        Cl_std = np.std(Cl[-n_avg:])
        CmPitch_final = np.mean(CmPitch[-n_avg:])
        
        # Additional analysis
        Cd_friction_final = np.mean(Cd_friction[-n_avg:])
        Cd_pressure_final = np.mean(Cd_pressure[-n_avg:])
        Cl_friction_final = np.mean(Cl_friction[-n_avg:])
        Cl_pressure_final = np.mean(Cl_pressure[-n_avg:])
        
        print("=== FORCE COEFFICIENTS RESULTS (CORRECTED) ===")
        print(f"Data points: {len(time)}")
        print(f"Time range: {time[0]:.1f} to {time[-1]:.1f}")
        print(f"")
        print(f"Final Converged Values (averaged over last {n_avg} time steps):")
        print(f"  Cl (Lift Coefficient):    {Cl_final:.6f} ± {Cl_std:.6f}")
        print(f"  Cd (Drag Coefficient):    {Cd_final:.6f} ± {Cd_std:.6f}")
        print(f"  CmPitch (Moment Coeff):   {CmPitch_final:.6f}")
        print(f"")
        print(f"Component Breakdown:")
        print(f"  Cl_friction:  {Cl_friction_final:.6f}")
        print(f"  Cl_pressure:  {Cl_pressure_final:.6f}")
        print(f"  Cd_friction:  {Cd_friction_final:.6f}")
        print(f"  Cd_pressure:  {Cd_pressure_final:.6f}")
        
        if Cd_final != 0:
            print(f"")
            print(f"  Cl/Cd Ratio:             {Cl_final/Cd_final:.3f}")
        else:
            print(f"  Cl/Cd Ratio:             ∞ (zero drag)")
        
        return {
            'time': time,
            'Cl': Cl,
            'Cd': Cd,
            'Cs': Cs,
            'CmPitch': CmPitch,
            'Cl_friction': Cl_friction,
            'Cl_pressure': Cl_pressure,
            'Cd_friction': Cd_friction,
            'Cd_pressure': Cd_pressure,
            'Cl_final': Cl_final,
            'Cd_final': Cd_final,
            'CmPitch_final': CmPitch_final,
            'Cl_std': Cl_std,
            'Cd_std': Cd_std,
            'converged': Cl_std < 0.001 and Cd_std < 0.001  # Convergence check
        }
        
    except Exception as e:
        print(f"Error reading forceCoeffs file: {e}")
        return None


def extract_from_forces(force_file, case_directory, inlet_velocity=14.6285, 
                       chord_length=1.0, span=0.01, rho=1.225):
    """
    Manual calculation from forces.dat file
    """
    try:
        data = []
        with open(force_file, 'r') as f:
            for line in f:
                if not line.startswith('#') and line.strip():
                    values = line.split()
                    if len(values) >= 7:
                        data.append([float(v) for v in values[:7]])
        
        if not data:
            print("Error: No valid data found in forces file")
            return None
        
        data = np.array(data)
        
        time = data[:, 0]
        # Total forces
        Fx = data[:, 1] + data[:, 4]  # Pressure + Viscous in X
        Fy = data[:, 2] + data[:, 5]  # Pressure + Viscous in Y
        
        # Calculate coefficients
        q_inf = 0.5 * rho * inlet_velocity**2
        S_ref = chord_length * span
        
        # Assuming small angle of attack: Drag ≈ Fx, Lift ≈ Fy
        Cd = Fx / (q_inf * S_ref)
        Cl = Fy / (q_inf * S_ref)
        
        # Final values
        n_avg = min(100, len(time))
        Cd_final = np.mean(Cd[-n_avg:])
        Cl_final = np.mean(Cl[-n_avg:])
        Cd_std = np.std(Cd[-n_avg:])
        Cl_std = np.std(Cl[-n_avg:])
        
        print("=== FORCE COEFFICIENTS (Manual Calculation) ===")
        print(f"Dynamic Pressure: {q_inf:.4f} Pa")
        print(f"Reference Area: {S_ref:.6f} m²")
        print(f"")
        print(f"Final Converged Values:")
        print(f"  Cl (Lift Coefficient):  {Cl_final:.6f} ± {Cl_std:.6f}")
        print(f"  Cd (Drag Coefficient):  {Cd_final:.6f} ± {Cd_std:.6f}")
        if Cd_final != 0:
            print(f"  Cl/Cd Ratio:           {Cl_final/Cd_final:.3f}")
        
        return {
            'time': time,
            'Cl': Cl,
            'Cd': Cd,
            'Fx': Fx,
            'Fy': Fy,
            'Cl_final': Cl_final,
            'Cd_final': Cd_final,
            'Cl_std': Cl_std,
            'Cd_std': Cd_std,
            'converged': Cl_std < 0.001 and Cd_std < 0.001
        }
        
    except Exception as e:
        print(f"Error reading forces file: {e}")
        return None


def plot_convergence(results, save_plot=True):
    """
    Plot convergence of Cl and Cd over time
    """
    try:
        import matplotlib.pyplot as plt
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # Lift coefficient
        ax1.plot(results['time'], results['Cl'], 'b-', linewidth=1.5, label='Cl (Total)')
        if 'Cl_friction' in results:
            ax1.plot(results['time'], results['Cl_friction'], 'g--', alpha=0.7, label='Cl (Friction)')
            ax1.plot(results['time'], results['Cl_pressure'], 'r--', alpha=0.7, label='Cl (Pressure)')
        ax1.axhline(y=results['Cl_final'], color='b', linestyle=':', alpha=0.7, 
                   label=f'Converged: {results["Cl_final"]:.4f}')
        ax1.set_ylabel('Lift Coefficient (Cl)')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        ax1.set_title('Lift Coefficient Convergence')
        
        # Drag coefficient
        ax2.plot(results['time'], results['Cd'], 'r-', linewidth=1.5, label='Cd (Total)')
        if 'Cd_friction' in results:
            ax2.plot(results['time'], results['Cd_friction'], 'g--', alpha=0.7, label='Cd (Friction)')
            ax2.plot(results['time'], results['Cd_pressure'], 'orange', linestyle='--', alpha=0.7, label='Cd (Pressure)')
        ax2.axhline(y=results['Cd_final'], color='r', linestyle=':', alpha=0.7,
                   label=f'Converged: {results["Cd_final"]:.4f}')
        ax2.set_ylabel('Drag Coefficient (Cd)')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        ax2.set_title('Drag Coefficient Convergence')
        
        # Cl/Cd ratio
        if results['Cd_final'] != 0:
            cl_cd_ratio = results['Cl'] / results['Cd']
            ax3.plot(results['time'], cl_cd_ratio, 'purple', linewidth=1.5)
            ax3.axhline(y=results['Cl_final']/results['Cd_final'], color='purple', 
                       linestyle=':', alpha=0.7, 
                       label=f'Final: {results["Cl_final"]/results["Cd_final"]:.2f}')
            ax3.set_ylabel('Lift-to-Drag Ratio (Cl/Cd)')
            ax3.grid(True, alpha=0.3)
            ax3.legend()
            ax3.set_title('Aerodynamic Efficiency')
        
        # Moment coefficient (if available)
        if 'CmPitch' in results:
            ax4.plot(results['time'], results['CmPitch'], 'orange', linewidth=1.5)
            ax4.axhline(y=results['CmPitch_final'], color='orange', linestyle=':', alpha=0.7,
                       label=f'Final: {results["CmPitch_final"]:.4f}')
            ax4.set_ylabel('Pitching Moment Coefficient (Cm)')
            ax4.grid(True, alpha=0.3)
            ax4.legend()
            ax4.set_title('Pitching Moment Convergence')
        
        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_xlabel('Time (iterations)')
        
        plt.tight_layout()
        
        if save_plot:
            plt.savefig('force_coefficients_convergence.png', dpi=300, bbox_inches='tight')
            print("Convergence plot saved as 'force_coefficients_convergence.png'")
        
        plt.show()
        
    except ImportError:
        print("Matplotlib not available. Install with: pip install matplotlib")


if __name__ == "__main__":
    # Extract Cl and Cd from simulation results
    results = extract_cl_cd_from_openfoam('.')

    if results is not None:
        print(f"\n=== SUMMARY ===")
        print(f"Final Cl: {results['Cl_final']:.6f}")
        print(f"Final Cd: {results['Cd_final']:.6f}")
        if results['Cd_final'] != 0:
            print(f"Final L/D: {results['Cl_final']/results['Cd_final']:.2f}")
        
        # Plot convergence
        plot_convergence(results)
    else:
        print("No simulation results found.")