import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import sobel

# Create a 2D potential field with a local minimum
x, y = np.meshgrid(np.arange(0, 10, 0.1), np.arange(0, 10, 0.1))
potential_field = (np.sin(x) - 5)**2 + (y - 5)**2

# Set the initial position
initial_position = np.array([2, 2])

# Set learning rate and number of iterations
learning_rate = 0.05
num_iterations = 200

# Move towards the minimum
current_position = initial_position.copy().astype(float)  # Convert to floats
trajectory = [current_position]

for _ in range(num_iterations):
    gradient = -sobel(potential_field, mode='constant')
    
    # Correct the indexing of the gradient array and round to the nearest integer
    update = learning_rate * gradient[current_position[1].astype(int), current_position[0].astype(int)]
    current_position += np.round(update)
    
    current_position = np.clip(current_position, 0, np.array(potential_field.shape) - 1)
    trajectory.append(current_position.copy())


# Visualize the potential field and trajectory
plt.figure(figsize=(8, 8))
plt.contourf(x, y, potential_field, levels=20, cmap='viridis')
plt.colorbar(label='Potential Field')
plt.plot(*zip(*trajectory), marker='o', color='red', label='Trajectory')
plt.scatter(*initial_position, color='blue', label='Initial Position')
plt.scatter(*trajectory[-1], color='green', label='Final Position')
plt.title('Movement towards Local Minimum')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.legend()
plt.show()

