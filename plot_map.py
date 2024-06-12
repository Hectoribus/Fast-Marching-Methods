import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the data from the Excel file
file_path = 'C:/Users/cheto/Desktop/UC3M/4th_Course/TFG/Maps/RESULT_MAPS/FM2_result.xlsx'  # Replace with your file path
df = pd.read_excel(file_path, header=None)

# Convert the DataFrame to a NumPy array and handle non-numeric data
data = df.apply(pd.to_numeric, errors='coerce').fillna(np.inf).to_numpy()

# Plotting the 2D graph
plt.figure()
plt.imshow(data, cmap='viridis', aspect='auto')
plt.colorbar(label='Value')
plt.show()

# Plotting the 3D graph
#fig = plt.figure(1)
#ax = fig.add_subplot(projection='3d')

# Create grid for 3D plot
#X, Y = np.meshgrid(np.arange(data.shape[0]), np.arange(data.shape[1]))

# Plot the surface
#ax.plot_surface(X, Y, data, rstride=2, cstride=2, cmap='viridis')
#plt.show()
