import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import sobel

def potential_field(x, y):
    # Define your 2D potential field function here
    return x**2 + y**2  # Example: a simple quadratic field

def gradient_sobel(Z):
    # Calculate the gradient of the potential field using the Sobel operator
    grad_x = sobel(Z, axis=1, mode='reflect')
    grad_y = sobel(Z, axis=0, mode='reflect')
    return grad_x, grad_y

# Create a 2D grid
x_values = np.linspace(-5, 5, 100)
y_values = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x_values, y_values)

# Calculate potential field values for each point on the grid
Z = potential_field(X, Y)

# Calculate the gradient of the potential field using the Sobel operator
grad_x, grad_y = gradient_sobel(Z)

# Set parameters for gradient descent
starting_point = (4, 4)
learning_rate = 0.1
num_iterations = 50

# Perform gradient descent using Sobel gradients
current_point = np.array(starting_point)
points = [current_point]

for _ in range(num_iterations):
    current_point[0] -= learning_rate * grad_x[int((current_point[1] + 5) / 10 * 100), int((current_point[0] + 5) / 10 * 100)]
    current_point[1] -= learning_rate * grad_y[int((current_point[1] + 5) / 10 * 100), int((current_point[0] + 5) / 10 * 100)]
    points.append(current_point.copy())

points = np.array(points)

# Plot the 2D potential field and the gradient descent path
plt.contourf(X, Y, Z, cmap='viridis', levels=20)
plt.colorbar(label='Potential Field Value')
plt.scatter(points[:, 0], points[:, 1], color='red', label='Sobel Gradient Descent Path')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('2D Potential Field Map with Sobel Gradient Descent')
plt.legend()
plt.show()


