import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import heapq



def is_any_zero(matrix):
    for row in matrix:
        for element in row:
            if element != 0:
                return True
    return False

def minTDim(x, Time, dim):
    i, j = x
    if dim == 1:
        # En la primera dimensión, definimos los dos vecinos
        nx_left = i - 1
        nx_right = i + 1
        # Si ninguno se sale del mapa, devolvemos el qe tenga menor valor de tiempo
        if 0 <= nx_left < grid.shape[0] and 0 <= nx_right < grid.shape[0]:
            return min(Time[(i + 1, j)], Time[(i - 1, j)])
        elif 0 > nx_left  or  nx_left >= grid.shape[0]:
            return Time[(i + 1, j)]
        elif 0 > nx_right  or  nx_right >= grid.shape[0]:
            return Time[(i - 1, j)]
            
    elif dim == 2:
        # En la segunda simensión lo mismo
        ny_top = j + 1
        ny_bottom = j - 1
        if 0 <= ny_top < grid.shape[1] and 0 <= ny_bottom < grid.shape[1]:
            return min(Time[(i, j + 1)], Time[(i, j - 1)])
        elif 0 > ny_top  or  ny_top >= grid.shape[1]:
            return Time[(i, j - 1)]
        elif 0 > ny_bottom  or  ny_bottom >= grid.shape[1]:
            return Time[(i, j + 1)]


def solveNDims(x, dim, T_values, V_field):
    if dim == 1:
        filtered_queue = []
    
        for element in T_values:
            time, dimension = element
            if dimension == 1:
                heapq.heappush(filtered_queue, (time, dimension))
                s_time, dimension = heapq.heappop(filtered_queue)

                return s_time + h/V_field[x]
        

    sum_T = 0
    sum_T2 = 0

    for element in T_values:
        time, dimension = element
        sum_T = sum_T + time
        sum_T2 = sum_T2 + time*time

    a = dim
    b = -2*sum_T
    c = sum_T2 - h**2/V_field[x]
    q = b**2 - 4*a*c

    if q < 0:
        return np.inf
    else:      
        return (-b + np.sqrt(q))/(2*a)



def solveEikonal(x, Time, V_field):
    
    a = N
    T_values = []

    for dim in range(1, N+1):
        min_T = minTDim(x, Time, dim)
        #print("min_T: ", min_T)
        
        if min_T != np.inf and min_T < Time[x]:
            heapq.heappush(T_values, (min_T, dim)) 
            #print("time x: ", Time[x])
            #print("enter");  
        else:
            a = a - 1
        
    if a == 0:
        return np.inf
    
    for dim in range(1, a+1):
        new_Time = solveNDims(x, dim, T_values, V_field)
        if dim == a:  # Need to change for higher dimensions (Thesis shows)
            break

    return new_Time


def get_neighbors(point, shape):
    x = point[0]
    y = point[1]
    neighbors = []
    for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < shape[0] and 0 <= ny < shape[1]: # Check if the neighors belong to the grid
            neighbors.append((nx, ny))
    return neighbors



def FMM(grid, Time, V_field, sources):
    Unknown = np.ones(grid.shape)
    Narrow = np.zeros(grid.shape)
    Frozen = np.zeros(grid.shape)

    Time = np.inf * np.ones(grid.shape)

    for x in sources:
        Time[x] = 0
        Unknown[x] = 0
        Narrow[x] = 1

    counter = 0
    
    while is_any_zero(Narrow): # Mientras que no valga todo 0 en Narrow
        #print("NARROW: ",Narrow,"\n\n")
        min_t = None
        x_min = None

        rows, cols = grid.shape
        for x in range(0, rows):
            for y in range(0, cols):
                if Narrow[x][y] == 1:
                    if min_t is None or Time[x][y] < min_t:
                        min_t = Time[x][y]
                        x_min = (x, y)
        if x_min is not None:
            neighbors = get_neighbors(x_min, grid.shape)
        #print(x_min,": siguiente")#--------------------------------------------------------------
        for x in neighbors:
            if Frozen[x] == 0:   # check every neighbor that is not Frozen
                if grid[x] == 0: # check that it's not an obstacle
                    new_Time = solveEikonal(x, Time, V_field)
                    #print(x," and time ",new_Time)#--------------------------------------------------
                    if new_Time < Time[x]:
                        #print("new min")#------------------------------------------------------------
                        Time[x] = new_Time
                    
                    if Unknown[x] == 1:
                        #print("discover")#-----------------------------------------------------------
                        Narrow[x] = 1
                        Unknown[x] = 0
                
        Narrow[x_min] = 0
        Frozen[x_min] = 1

        counter = counter + 1
        #if counter%50 == 0:
            #visualize(Time)
            #visualize(Narrow)
    return Time, Narrow



def visualize(map):
    plt.imshow(map, cmap='viridis', origin='lower')
    plt.colorbar()
    plt.title('Map')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

def visualize_grid(grid):
    plt.imshow(grid, cmap='viridis', origin='lower')
    plt.title('Grid with Obstacles')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

def update_obstacles(event):
    global grid
    if event.button == 1:  # Left mouse button
        x, y = int(event.ydata), int(event.xdata)
        grid[x, y] = 1  # Set the clicked cell as an obstacle
        visualize_grid(grid)

N = 2
h = 1

grid_shape = (10, 10)
grid = np.zeros(grid_shape)
Times = np.zeros(grid_shape)
Velocities = np.ones(grid_shape)
source = [(1, 1), (9, 9)]

fig, ax = plt.subplots()
ax.imshow(np.zeros(grid_shape), cmap='viridis', origin='lower')
ax.set_title('Left-click to Add Obstacles')
ax.set_xlabel('X')
ax.set_ylabel('Y')
fig.canvas.mpl_connect('button_press_event', update_obstacles)
plt.show()

# Wait for obstacle selection
# input("Press Enter to calculate the distance map...")

# Set obstacle cells
obstacle_cells = set()

for cell in obstacle_cells:
    grid[cell] = 1  # 1 represents an obstacle

print ("hello")

new_Times, Narrow = FMM(grid, Times, Velocities, source)
#print(new_Times)
visualize(new_Times)

#X, Y = np.meshgrid(np.arange(grid_shape[0]), np.arange(grid_shape[1]))

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#surf = ax.plot_surface(X, Y, new_Times, cmap='viridis')
#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.show()