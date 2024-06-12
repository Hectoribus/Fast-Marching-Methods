import numpy as np
import heapq
import matplotlib.pyplot as plt


def velocity_field(grid, obstacles):
    # Initialize the distance field, assigning every cell a value of inf (to ensure that they are not considered in the wave propagation)
    velocity_F = np.inf * np.ones(grid.shape)
    status = np.zeros(grid.shape)  # 0 for unprocessed, 1 for accepted

    # Priority queue for active cells
    pq = []

    # Initialize the distance values and status for the sources (obstacles)
    for obstacle in obstacles:
        velocity_F[obstacle] = 0
        status[obstacle] = 1
        pq.append((0, obstacle))

    while pq:
        # Obtain the first element of the queue
        h, current = heapq.heappop(pq)

        # Update neighbors
        neighbors = get_neighbors(current, grid.shape)

        for neighbor in neighbors:
            if status[neighbor] == 0:  # Check if neighbor is unprocessed
                if grid[neighbor] == 0:  # Check if neighbor is not an obstacle

                    # Calculate a new possible distance, based on the value of the current cell + the
                    # travel cost to the neighbor (Euclidean distance).
                    new_distance = calculate_distance(velocity_F, current, neighbor)

                    # Update new distance
                    if new_distance < velocity_F[neighbor]:
                        velocity_F[neighbor] = new_distance
                        status[neighbor] = 1

                        # Insert the neighbor to the queue
                        heapq.heappush(pq, (new_distance, neighbor))

    return velocity_F

def calculate_distance(distance, current, neighbor):
    return distance[current] + euclidean_distance(current, neighbor)

def euclidean_distance(p1, p2):
    return np.linalg.norm(np.array(p1) - np.array(p2))


def visualize(map):
    plt.imshow(map, cmap='viridis', origin='lower')
    plt.colorbar()
    plt.title('Map')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()



def get_neighbors(point, shape):
    x, y = point
    neighbors = []
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            new_x, new_y = x + dx, y + dy
            if 0 <= new_x < shape[0] and 0 <= new_y < shape[1]:
                neighbors.append((new_x, new_y))
    return neighbors

def fast_marching_with_obstacles(source_point, target_point, velocity_map, shape, obstacles):
    # Initialize the grid and distance field
    distance_field = np.full(shape, np.inf)
    source_processed = np.zeros(shape, dtype=bool) 
    target_processed = np.zeros(shape, dtype=bool)

    # Create a priority queue for wavefront propagation
    source_wavefront = []
    target_wavefront = []

    distance_field[source_point] = 0
    distance_field[target_point] = 0
    heapq.heappush(source_wavefront, (0, source_point))
    heapq.heappush(target_wavefront, (0, target_point))
    
    current_source_point = source_point
    current_target_point = target_point

    break_outer = False

    while source_processed[current_target_point] != 1   and   target_processed[current_source_point] != 1:
        # Pop the point with the smallest distance
        source_distance, current_source_point = heapq.heappop(source_wavefront)
        target_distance, current_target_point = heapq.heappop(target_wavefront)

        if source_processed[current_source_point] or target_processed[current_target_point]:
            continue

        source_processed[current_source_point] = True
        target_processed[current_target_point] = True

        for neighbor_point in get_neighbors(current_source_point, shape):
            if target_processed[neighbor_point] != 1:
                if not source_processed[neighbor_point]:
                    # Check if the neighbor is an obstacle
                    if neighbor_point in obstacles:
                        continue
                    
                    # Calculate the time it takes to reach the neighbor
                    x, y = neighbor_point
                    dx = x - current_source_point[0]
                    dy = y - current_source_point[1]
                    distance_to_neighbor = np.sqrt(dx*dx + dy*dy)
                    dx = 1.0 / velocity_map[x, y]
                    tentative_distance = source_distance + dx * distance_to_neighbor

                    if tentative_distance < distance_field[x, y]:
                        distance_field[x, y] = tentative_distance
                        heapq.heappush(source_wavefront, (tentative_distance, neighbor_point))
            else:
                break_outer = True
                break

        if break_outer:
            break

        for neighbor_point in get_neighbors(current_target_point, shape):
            if source_processed[neighbor_point] != 1:
                if not target_processed[neighbor_point]:
                    # Check if the neighbor is an obstacle
                    if neighbor_point in obstacles:
                        continue
                    
                    # Calculate the time it takes to reach the neighbor
                    x, y = neighbor_point
                    dx = x - current_target_point[0]
                    dy = y - current_target_point[1]
                    distance_to_neighbor = np.sqrt(dx*dx + dy*dy)
                    dx = 1.0 / velocity_map[x, y]
                    tentative_distance = target_distance + dx * distance_to_neighbor

                    if tentative_distance < distance_field[x, y]:
                        distance_field[x, y] = tentative_distance
                        heapq.heappush(target_wavefront, (tentative_distance, neighbor_point))
            else:
                break_outer = True
                break
        if break_outer:
            break

        # visualize(distance_field)
        


    return distance_field


def find_shortest_path_2dir(arrival_time_Map, source, target):
    # Initiate the path from the target 

    for i in range(len(arrival_time_Map)):
        for j in range(len(arrival_time_Map[i])):
            if arrival_time_Map[i][j] == np.inf:
                arrival_time_Map[i][j] = 0


    path_source = [source]
    path_target = [target]
    current_source = source
    current_target = target
    
    while current_source != current_target:
        neighbors_source = get_neighbors(current_source, arrival_time_Map.shape)
        neighbors_target = get_neighbors(current_target, arrival_time_Map.shape)
        
        max_neighbor_source = None
        max_distance_source = arrival_time_Map[current_source]
        max_neighbor_target = None
        max_distance_target = arrival_time_Map[current_target]
        
        for neighbor in neighbors_source:
            if arrival_time_Map[neighbor] > max_distance_source:
                max_distance_source = arrival_time_Map[neighbor]
                max_neighbor_source = neighbor
        
        for neighbor in neighbors_target:
            if arrival_time_Map[neighbor] > max_distance_target:
                max_distance_target = arrival_time_Map[neighbor]
                max_neighbor_target = neighbor


        if max_neighbor_source is not None:
            current_source = max_neighbor_source
            path_source.append(current_source)
        else:
            break

        if max_neighbor_target is not None:
            current_target = max_neighbor_target
            path_target.append(current_target)
        else:
            break
    

    for i in range(len(arrival_time_Map)):
        for j in range(len(arrival_time_Map[i])):
            if arrival_time_Map[i][j] == 0   and  (i, j) != source   and   (i, j) != target:
                arrival_time_Map[i][j] = np.inf


    path_target.reverse()

    path = path_source + path_target
    return np.array(path)


def visualize(map):
    plt.imshow(map, cmap='viridis', origin='lower')
    plt.colorbar()
    plt.title('Map')
    plt.xlabel('X')
    plt.ylabel('Y')
    #plt.show()
    plt.draw()

def visualize_path(grid, path):
    plt.imshow(grid, cmap='viridis', origin='lower')
    plt.colorbar()
    plt.title('Path Visualization')
    plt.xlabel('X')
    plt.ylabel('Y')
    
    x, y = zip(*path)
    plt.plot(y, x, color='red', linewidth=2)
    
    plt.show()


grid_shape = (100, 100)
grid = np.zeros(grid_shape)

obstacles = [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8),
             (5, 1), (5, 9), (5, 2), (5, 3), (5, 5), (5, 6), (5, 7), (5, 8),
             (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 6), (8, 7), (8, 8),
             (0, 10), (1, 10), (2, 10), (3, 10), (4, 10), (5, 10)]


seed_point = (0, 0)
target_point = (50, 50)


velocity_map = np.ones(grid_shape)

velocity_map = velocity_field(grid, obstacles) #NORMALIZE

visualize(velocity_map)


arrival_time_map = fast_marching_with_obstacles(seed_point, target_point, velocity_map, grid_shape, obstacles)

visualize(arrival_time_map)


shortest_path_2dir = find_shortest_path_2dir(arrival_time_map, (0, 0), (50, 50))

visualize_path(arrival_time_map, shortest_path_2dir)
