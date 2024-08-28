import time
import pandas as pd
import concurrent.futures

# Function to append data to CSV in a separate thread
def append_to_csv(data, filename):
    df = pd.DataFrame(data)
    df.to_csv(filename, mode='a', header=False, index=False)

# Initialize the CSV file with headers
filename = 'simulation_data.csv'
columns = ['Column1', 'Column2', 'Column3']
pd.DataFrame(columns=columns).to_csv(filename, index=False)

# Example simulation loop
data_batch = []
batch_size = 10  # Adjust as needed

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    for i in range(1000):  # Assuming 1000 iterations
        
        # Step 1: Simulation - Replace this with your actual simulation code
        new_data = {'Column1': i, 'Column2': i * 2, 'Column3': i * 3}
        data_batch.append(new_data)

        # Step 2: Append to CSV in the background when batch is ready
        if len(data_batch) >= batch_size:
            # Offload CSV writing to a separate thread
            executor.submit(append_to_csv, data_batch, filename)
            data_batch = []  # Clear batch for next iteration

        # Step 3: Proceed to next iteration - The loop continues immediately

    # Final step: Write any remaining data in the last batch
    if data_batch:
        executor.submit(append_to_csv, data_batch, filename)
