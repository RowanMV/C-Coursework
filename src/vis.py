import matplotlib.pyplot as plt
import numpy as np
import os
import csv

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

    with open(output_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Vessel", "Final Time", "Final Displacement (y)", "Final Velocity (v)", "Final Water Level (h)"])
        
        for vessel_idx, vessel_data in enumerate(vessels_data):
            final_row = vessel_data[-1]
            writer.writerow([f"Vessel {vessel_idx + 1}", time_steps[-1], *final_row])
    
    print(f"Final data saved to '{output_filename}'.")


def plot_vessel_data(time_steps, vessels_data):
    """
    Plots the displacement for all vessels.
    :param time_steps: Array of time steps.
    :param vessels_data: List of vessel data arrays (one per vessel).
    """
    if not vessels_data or not time_steps.any():
        print("No data to plot.")
        return

    # Only plot displacement (y)
    variable = "Displacement (y)"
    
    plt.figure(figsize=(12, 8))

    for vessel_idx, vessel_data in enumerate(vessels_data):
        plt.plot(time_steps, vessel_data[:, 0], label=f"Vessel {vessel_idx + 1}")  # Only the first column for displacement
    
    plt.title(variable)
    plt.xlabel("Time (s)")
    plt.ylabel(variable)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig("vessel_displacement.png")
    plt.show()



def main():
    output_file = "output.txt"
    if not os.path.exists(output_file):
        print(f"Error: '{output_file}' not found in the current working directory.")
    else:
        # Read data
        time_steps, vessels_data = read_output_file(output_file)
        
        # Save final data to table
        save_final_data_table(time_steps, vessels_data)
        
        # Plot data
        plot_vessel_data(time_steps, vessels_data)


if __name__ == "__main__":
    main()
