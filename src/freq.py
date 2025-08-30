import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import pandas as pd

def read_output_file(filename="output.txt"):
    """
    Reads the simulation output file and returns the time steps and data.
    :param filename: Name of the output file (default: "output.txt").
    :return: time_steps (1D array), data (3D array where [vessel][variable][time])
    """
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            
            # Parse time steps
            time_steps = []
            vessels_data = []
            
            for line in lines:
                # Convert to a list of float values
                values = list(map(float, line.split()))
                time_steps.append(values[0])  # First column is time
                
                # Remaining columns are vessel data
                data = values[1:]
                num_vessels = len(data) // 3  # Each vessel has 3 variables
                
                if not vessels_data:
                    vessels_data = [[] for _ in range(num_vessels)]
                
                for i in range(num_vessels):
                    vessels_data[i].append(data[i * 3:(i + 1) * 3])
            
            # Convert to numpy arrays for easy manipulation
            time_steps = np.array(time_steps)
            vessels_data = [np.array(vessel) for vessel in vessels_data]
            return time_steps, vessels_data
    
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None, None
    except Exception as e:
        print(f"Error: {e}")
        return None, None


def calculate_frequency(time_steps, displacement_data):
    """
    Calculate the frequency of oscillation based on displacement data.
    :param time_steps: Time steps of the simulation.
    :param displacement_data: Displacement data of a vessel.
    :return: Frequency in Hz.
    """
    # Find peaks in displacement data and compute periods
    peaks = np.where((displacement_data[1:] > displacement_data[:-1]) & 
                     (displacement_data[1:] > displacement_data[2:]))[0] + 1
    
    if len(peaks) < 2:
        return 0
    
    # Calculate periods between peaks
    periods = np.diff(time_steps[peaks])
    
    # Calculate frequency (1 / period)
    frequency = 1 / np.mean(periods)
    return frequency


def investigate_frequency_change():
    """
    Investigates how frequency changes with different simulation parameters.
    """
    # Define a range of parameters to investigate (e.g., vessel mass, resistance)
    # For simplicity, assume we are changing one parameter at a time
    parameter_values = np.linspace(1.0, 10.0, 10)  # Example range for a parameter
    frequencies = []

    for param_value in parameter_values:
        # Simulate with the current parameter value
        output_file = f"output_{param_value}.txt"
        
        # Read data
        time_steps, vessels_data = read_output_file(output_file)
        
        if time_steps is None or len(vessels_data) == 0:
            continue
        
        # Calculate frequency for the first vessel (change based on your analysis)
        frequency = calculate_frequency(time_steps, vessels_data[0][:, 0])  # Displacement (y)
        frequencies.append(frequency)
    
    # Plot frequency change with respect to the parameter
    plt.figure(figsize=(10, 6))
    plt.plot(parameter_values, frequencies, marker='o', linestyle='-', color='b')
    plt.title("Frequency Change with Parameter Variation")
    plt.xlabel("Parameter Value")
    plt.ylabel("Frequency (Hz)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("frequency_change.png")
    plt.show()


def save_final_data_table(time_steps, vessels_data, output_filename="final_data.csv"):
    """
    Saves the final displacement, velocity, and water level of vessels to a CSV file.
    :param time_steps: Array of time steps.
    :param vessels_data: List of vessel data arrays (one per vessel).
    :param output_filename: Name of the output file.
    """
    if not vessels_data or not time_steps.any():
        print("No data to save.")
        return

    # Create a DataFrame and save the final data
    final_data = []
    for vessel_idx, vessel_data in enumerate(vessels_data):
        final_row = vessel_data[-1]
        final_data.append([f"Vessel {vessel_idx + 1}", time_steps[-1], *final_row])

    df = pd.DataFrame(final_data, columns=["Vessel", "Final Time", "Final Displacement (y)", "Final Velocity (v)", "Final Water Level (h)"])
    df.to_csv(output_filename, index=False)
    print(f"Final data saved to '{output_filename}'.")


def main():
    # Run the investigation on how frequency changes with a parameter
    investigate_frequency_change()


if __name__ == "__main__":
    main()
