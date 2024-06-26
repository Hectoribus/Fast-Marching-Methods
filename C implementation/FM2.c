#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define N 2
#define H 1
#define GRID_SIZE 500
#define PI 3.14159265359

typedef struct {
    int x;
    int y;
} Point;

typedef struct {
    double x;
    double y;
} pathPoint;

typedef struct {
    float min_time;
    int dimension;
} TValue;


//-----------------------------------------------------------------------------------------------------------------------------
//-------------------------------------                     H·E·A·P                --------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------

typedef struct {
    int x;
    int y;
    float time;
} Narrow_Point;

// Declare a heap structure
typedef struct{
    Narrow_Point* arr;
    int size;
    int capacity;
} heap;

// forward declarations
heap* createHeap(int capacity);
void insertHelper(heap* h, int index);
void heapify(heap* h, int index);
Narrow_Point extractMin(heap* h);
void insert(heap* h, Narrow_Point data);
void resizeHeap(heap* h, int newCapacity);

// createHeap function
heap* createHeap(int capacity)
{
    // Allocating memory to heap h
    heap* h = (heap*)malloc(sizeof(heap));

    // Checking if memory is allocated to h or not
    if (h == NULL) {
        printf("Memory error");
        return NULL;
    }
    // set the values to size and capacity
    h->size = 0;
    h->capacity = capacity;

    // Allocating memory to array
    if (capacity > 0) {
        h->arr = (Narrow_Point*)malloc(capacity * sizeof(Narrow_Point));
        // Checking if memory is allocated to h->arr or not
        if (h->arr == NULL) {
            printf("Memory error");
            return NULL;
        }
    } else {
        h->arr = NULL; // Initialize as NULL if capacity is zero
    }

    return h;
}

// InsertHelper function
void insertHelper(heap* h, int index)
{
    // Store parent of element at index
    // in parent variable
    int parent = (index - 1) / 2;

    if (h->arr[parent].time > h->arr[index].time) {
        // Swapping when child is smaller
        // than parent element
        Narrow_Point temp = h->arr[parent];
        h->arr[parent] = h->arr[index];
        h->arr[index] = temp;

        // Recursively calling insertHelper
        insertHelper(h, parent);
    }
}

void heapify(heap* h, int index)
{
    int left = index * 2 + 1;
    int right = index * 2 + 2;
    int min = index;

    // Checking whether our left or child element
    // is at right index or not to avoid index error
    if (left >= h->size || left < 0)
        left = -1;
    if (right >= h->size || right < 0)
        right = -1;

    // store left or right element in min if
    // any of these is smaller than its parent
    if (left != -1 && h->arr[left].time < h->arr[index].time)
        min = left;
    if (right != -1 && h->arr[right].time < h->arr[min].time)
        min = right;

    // Swapping the nodes
    if (min != index) {
        Narrow_Point temp = h->arr[min];
        h->arr[min] = h->arr[index];
        h->arr[index] = temp;

        // recursively calling for their child elements
        // to maintain min heap
        heapify(h, min);
    }
}

Narrow_Point extractMin(heap* h)
{
    Narrow_Point deleteItem;

    // Checking if the heap is empty or not
    if (h->size == 0) {
        printf("\nHeap is empty.");
        return (Narrow_Point){-1, -1, -1}; // Return a dummy Point in case of empty heap
    }

    // Store the node in deleteItem that
    // is to be deleted.
    deleteItem = h->arr[0];

    // Replace the deleted node with the last node
    h->arr[0] = h->arr[h->size - 1];
    // Decrement the size of heap
    h->size--;

    // Call heapify for 0th index
    // to maintain the heap property
    heapify(h, 0);
    return deleteItem;
}

Narrow_Point getMin(heap* h)
{
    // Checking if the heap is empty or not
    if (h->size == 0) {
        printf("\nHeap is empty.");
        return (Narrow_Point){-1, -1, -1.0}; // Return a dummy Point in case of empty heap
    }

    // Return the root element, which is the minimum
    return h->arr[0];
}

// Insert function
void insert(heap* h, Narrow_Point data)
{
    // Checking if heap is full or not
    if (h->size >= h->capacity) {
        // Resize the heap if full
        int newCapacity = (h->capacity == 0) ? 1 : h->capacity * 2;
        resizeHeap(h, newCapacity);
    }

    // Inserting data into an array
    h->arr[h->size] = data;
    // Calling insertHelper function
    insertHelper(h, h->size);
    // Incrementing size of array
    h->size++;
}

void resizeHeap(heap* h, int newCapacity)
{
    h->arr = (Narrow_Point*)realloc(h->arr, newCapacity * sizeof(Narrow_Point));
    if (h->arr == NULL) {
        printf("Memory error");
        exit(1); // Exit if memory allocation fails
    }
    h->capacity = newCapacity;
}

void printHeap(heap* h)
{
    for (int i = 0; i < h->size; i++) {
        printf("(%d, %d, %.2f) ", h->arr[i].x, h->arr[i].y, h->arr[i].time);
    }
    printf("\n");
}

void freeHeap(heap* h)
{
    if (h->arr != NULL) {
        free(h->arr);
    }
    free(h);
}


//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------



void read_matrix(const char *filename, double **matrix, int rows, int cols) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (fscanf(file, "%lf,", &matrix[i][j]) != 1) {
                perror("Error reading file");
                fclose(file);
                exit(EXIT_FAILURE);
            }
            matrix[i][j] = (matrix[i][j] == 5.0) ? 1 : 0;
        }
    }

    fclose(file);
}


float minTDim(Point point, int dim, float *Time_MA) { // DEVOLVEMOS EL VALOR MÁS PEQUEÑO DE TIEMPO EN EL ARRAY TIME PARA CADA DIMENSIÓN (UP-DOWN, LEFT-RIGHT)
    int i = point.x, j = point.y;
    if (dim == 1) {
        int nx_left = i - 1;
        int nx_right = i + 1;
        if (0 <= nx_left && nx_left < GRID_SIZE && 0 <= nx_right && nx_right < GRID_SIZE) {
            return fmin(Time_MA[(i+1)*GRID_SIZE + j], Time_MA[(i-1)*GRID_SIZE + j]);
        } else if (0 > nx_left || nx_left >= GRID_SIZE) {
            return Time_MA[(i+1)*GRID_SIZE + j];
        } else if (0 > nx_right || nx_right >= GRID_SIZE) {
            return Time_MA[(i-1)*GRID_SIZE + j];
        }
    } else if (dim == 2) {
        int ny_top = j + 1;
        int ny_bottom = j - 1;
        if (0 <= ny_top && ny_top < GRID_SIZE && 0 <= ny_bottom && ny_bottom < GRID_SIZE) {
            return fmin(Time_MA[i*GRID_SIZE + (j+1)], Time_MA[i*GRID_SIZE + (j-1)]);
        } else if (0 > ny_top || ny_top >= GRID_SIZE) {
            return Time_MA[i*GRID_SIZE + (j-1)];
        } else if (0 > ny_bottom || ny_bottom >= GRID_SIZE) {
            return Time_MA[i*GRID_SIZE + (j+1)];
        }
    }
    return 0.0;
}

float solveNDims(Point point, int dim, TValue T_Values[2], float *V_field_MA) {
    int num_dim = 0;
    int i;

    for(i = 0; i < 2; i++){ //DESCUBRIR QUÉ T_VALUE ES EL QUE NO VALE, PARA ANALIZAR EL OTRO.
        if(T_Values[i].dimension == 0){
            num_dim = 1;
            break;
        }
    }

    if (num_dim == 1) { // CASO PARA CUANDO SÓLO HAYA VALORES EN UNA DE LAS DIMENSIONES.
        if(i == 0){
            return T_Values[1].min_time + H/V_field_MA[point.x*GRID_SIZE + point.y];
        }
        else if(i == 1){
            return T_Values[0].min_time + H/V_field_MA[point.x*GRID_SIZE + point.y];
        }
    }

    // CASO PARA CUANDO HAY VALORES VÁLIDOS EN AMBAS DIMENSIONES Y HAY QUE HACER LOS SIGUIENTES CÁLCULOS PARA EL NUEVO TIEMPO
    float sum_T = 0.0, sum_T2 = 0.0;

    for(int i = 0; i < 2; i++) {
        sum_T = sum_T + T_Values[i].min_time;
        sum_T2 = sum_T2 + (T_Values[i].min_time * T_Values[i].min_time);
    }

    float a = dim;
    float b = -2 * sum_T;
    float c = sum_T2 - pow(H, 2) / V_field_MA[point.x*GRID_SIZE + point.y];
    float q = b * b - 4 * a * c;

    if (q < 0) {
        return INFINITY;
    } else {
        return (-b + sqrt(q)) / (2 * a);
    }

    return 0.0;
}

float solveEikonal(Point point, float *V_field_MA, float *Time_MA) {
    int a = N;
    float new_Time = 0.0;

    TValue T_Values[2] = {{INFINITY, 0}, {INFINITY, 0}}; //INICIALIZAMOS T_VALUES

    for (int dim = 1; dim <= N; dim++) {
        float min_T = minTDim(point, dim, Time_MA);

        if (min_T != INFINITY && min_T < Time_MA[point.x*GRID_SIZE + point.y]) { // METER EN T_VALUES LOS TIEMPOS QUE SEAN MENORES QUE LOS ANTERIORES
            for(int i = 0; i < 2; i++){
                if(T_Values[i].min_time == INFINITY){
                    T_Values[i].min_time = min_T;
                    T_Values[i].dimension = dim;

                    break;
                }
            }
        } else {
            a--;
        }
    }
    
    if (a == 0) {
        return INFINITY;
    }
    
    for (int dim = 1; dim <= a; dim++) {
        new_Time = solveNDims(point, dim, T_Values, V_field_MA);
        if(dim == a){
            break;
        }
    }
    
    return new_Time;
}

Point* appendToArray(Point *arr, int *size, Point element) { 
    // Increment the size of the array 
    (*size)++; 
     
    // Allocate memory for the new array 
    arr = (Point*)realloc(arr, (*size) * sizeof(Point));
     
    if (arr == NULL) { 
        // Handle memory allocation failure 
        printf("Memory allocation failed.\n"); 
        exit(1); 
    } 
    
    // Add the new struct at the end of the array
    arr[(*size) - 1] = element;
 
    // Return the updated array 
    return arr; 
}

pathPoint* appendToPath(pathPoint *arr, int *size, pathPoint element) { 
    // Increment the size of the array 
    (*size)++; 
     
    // Allocate memory for the new array 
    arr = (pathPoint*)realloc(arr, (*size) * sizeof(pathPoint));
     
    if (arr == NULL) { 
        // Handle memory allocation failure 
        printf("Memory allocation failed.\n"); 
        exit(1); 
    } 
    
    // Add the new struct at the end of the array
    arr[(*size) - 1] = element;
 
    // Return the updated array 
    return arr; 
} 

Point* removeFromArray(Point *arr, int *size, Point elementToRemove) {
    int i, j;
    int found = 0;

    // Search for the element in the array
    for (i = 0; i < *size; i++) {
        if (arr[i].x == elementToRemove.x && arr[i].y == elementToRemove.y) {
            found = 1;
            break;
        }
    }

    // If the element is found, remove it
    if (found) {
        // Shift elements to fill the gap
        for (j = i; j < *size - 1; j++) {
            arr[j] = arr[j + 1];
        }
        // Decrease the size of the array
        (*size)--;

        // Resize the array or return NULL if size becomes 0
        if (*size > 0) {
            arr = (Point*)realloc(arr, (*size) * sizeof(Point));

            if (arr == NULL) {
                // Handle memory allocation failure
                printf("Memory allocation failed.\n");
                exit(1);
            }
        } else {
            free(arr);
            arr = NULL;
        }
    } else {
        printf("Point not found in the array.\n");
    }

    // Return the updated array
    return arr;
}


void FMM(double **Grid_MA, float *V_field_MA, float *Time_MA, Point *sources, int num_sources, Point target) {
    
    bool *Unknown_ma = (bool *)malloc(GRID_SIZE * GRID_SIZE * sizeof(bool));
    bool *Frozen_ma = (bool *)malloc(GRID_SIZE * GRID_SIZE * sizeof(bool));

    heap* HEAP_Narrow_band = createHeap(0);

    printf("Start!\n");

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            Unknown_ma[i*GRID_SIZE + j] = 1;
            Frozen_ma[i*GRID_SIZE + j] = 0;
        }
    }

    Point source = {0, 0};

    for (int k = 0; k < num_sources; k++) {
        source.x = sources[k].x;
        source.y = sources[k].y;
        Time_MA[source.x*GRID_SIZE + source.y] = 0;
        Unknown_ma[source.x*GRID_SIZE + source.y] = 0;

        insert(HEAP_Narrow_band, (Narrow_Point){source.x, source.y, Time_MA[source.x*GRID_SIZE + source.y]});
        //printHeap(HEAP_Narrow_band);
    }

    long counter_points = 0;

    Narrow_Point minPoint;

    while (HEAP_Narrow_band->size != 0) { // MIENTRAS QUE EL ARRAY NARROW NO ESTÉ VACÍO (0 EN TODAS LAS POSICIONES)
        
        // BUSCAMOS EL PUNTO DEL FRENTE DE ONDA (NARROW) CON MENOR TIEMPO PARA ANALIZARLO (EXPANDIR LA ONDA)
        minPoint = getMin(HEAP_Narrow_band);

        counter_points = counter_points + 1;

        if((minPoint.y == target.y) && (minPoint.x == target.x)){
            break;
        }
        
        Point neighbors[4];
        
        int n_num = 0;
        int nx, ny;
        for(int dx = -1; dx <= 1; dx = dx+2){
            nx = minPoint.x + dx;
            ny = minPoint.y;
            if((0 <= nx) && (nx < GRID_SIZE) && (0 <= ny) && (ny < GRID_SIZE)){
                neighbors[n_num].x = nx;
                neighbors[n_num].y = ny;
                n_num++;
            }
        }
        for(int dy = -1; dy <= 1; dy = dy+2){
            nx = minPoint.x;
            ny = minPoint.y + dy;
            if((0 <= nx) && (nx < GRID_SIZE) && (0 <= ny) && (ny < GRID_SIZE)){
                neighbors[n_num].x = nx;
                neighbors[n_num].y = ny;
                n_num++;
            }
        }

        for(int i = 0; i < n_num; i++) {
            Point neighbor = neighbors[i];
            //ANALIZAMOS CADA VECINO SIEMPRE QUE NO SEA UN OBSTÁCULO O YA LO HAYAMOS ANALIZADO ANTES
            if (Frozen_ma[neighbor.x*GRID_SIZE + neighbor.y] == 0 && Grid_MA[neighbor.x][neighbor.y] == 0) { 

                float new_Time = solveEikonal(neighbor, V_field_MA, Time_MA); //OBTENEMOS NUEVO TIEMPO PARA EL VECINO

                if (new_Time < Time_MA[neighbor.x*GRID_SIZE + neighbor.y]) {
                    Time_MA[neighbor.x*GRID_SIZE + neighbor.y] = new_Time;
                }
                if (Unknown_ma[neighbor.x*GRID_SIZE + neighbor.y] == 1) {
                    insert(HEAP_Narrow_band, (Narrow_Point){neighbor.x, neighbor.y, new_Time});
                    Unknown_ma[neighbor.x*GRID_SIZE + neighbor.y] = 0;
                }
            }
        }

        extractMin(HEAP_Narrow_band);
        Frozen_ma[minPoint.x*GRID_SIZE + minPoint.y] = 1;
    }

    free(Unknown_ma);
    free(Frozen_ma);
    freeHeap(HEAP_Narrow_band);
    printf("Counter_Points: %li\n", counter_points);
}

void findPath(float *Time_MA, pathPoint source, pathPoint target, double step){
    
    int counter = 0;

    int size = 0;
    pathPoint* path = NULL;

    path = appendToPath(path, &size, target);

    
    pathPoint p_i = {target.x, target.y};
    pathPoint p_next_i;
    int x = target.x;
    int y = target.y;

    double right, left, up, down;
    double grad;
    double gx, gy;

    //(floor(p_i.x) != source.x) && (floor(p_i.y) != source.y) ||  && (counter < 20)

    while(((ceil(p_i.x) != source.x) || (ceil(p_i.y) != source.y)) && (counter < 200)){
        /*
                        _________>  Y
                        |
                        |
                        |
                    X  \/
        */

        int x = floor(p_i.x);
        int y = floor(p_i.y);

        //printf("Current point: [%i, %i];   Time: %lf\n", x, y, Time_MA[x*GRID_SIZE + y]);
        
        
        right = (p_i.y + 1 < GRID_SIZE) ? Time_MA[x*GRID_SIZE + y + 1] : 0;
        left = (p_i.y - 1 >= 0) ? Time_MA[x*GRID_SIZE + y - 1] : 0;
        up = (p_i.x - 1 >= 0) ? Time_MA[(x-1)*GRID_SIZE + y] : 0;
        down = (p_i.x + 1 < GRID_SIZE) ? Time_MA[(x+1)*GRID_SIZE + y] : 0;

        //printf("up: %lf, down: %lf, right: %lf, left: %lf\n", up, down, right, left);
        
        
        if(right == INFINITY || left == INFINITY || up == INFINITY || down == INFINITY){
            if (right == INFINITY){
                right = Time_MA[x*GRID_SIZE + y];
                //printf("right inf\n");
                gy = (up-down);
                gx = (right-left);
                grad = fabs(atan(gy/gx));
                //printf("grad: %lf\n", grad*180/PI);
                if(up < down){
                    p_next_i.x = p_i.x - step*sin(grad);
                    p_next_i.y = p_i.y - step*cos(grad);
                }
                else{
                    p_next_i.x = p_i.x + step*sin(grad);
                    p_next_i.y = p_i.y - step*cos(grad);
                }
            }
            else if(left == INFINITY){
                left = Time_MA[x*GRID_SIZE + y];
                //printf("left inf\n");
                gy = (up-down);
                gx = (right-left);
                grad = fabs(atan(gy/gx));
                //printf("grad: %lf\n", grad*180/PI);
                if(up < down){
                    p_next_i.x = p_i.x - step*sin(grad);
                    p_next_i.y = p_i.y + step*cos(grad);
                }
                else{
                    p_next_i.x = p_i.x + step*sin(grad);
                    p_next_i.y = p_i.y + step*cos(grad);
                }
            }
            if(up == INFINITY){
                up = Time_MA[x*GRID_SIZE + y];
                //printf("up inf\n");
                gy = (up-down);
                gx = (right-left);
                grad = fabs(atan(gy/gx));
                //printf("grad: %lf\n", grad*180/PI);
                if(right < left){
                    p_next_i.x = p_i.x + step*sin(grad);
                    p_next_i.y = p_i.y + step*cos(grad);
                }
                else{
                    p_next_i.x = p_i.x + step*sin(grad);
                    p_next_i.y = p_i.y - step*cos(grad);
                }
            }
            else if(down == INFINITY){
                down = Time_MA[x*GRID_SIZE + y];
                //printf("down inf\n");
                gy = (up-down);
                gx = (right-left);
                grad = fabs(atan(gy/gx));
                //printf("grad: %lf\n", grad*180/PI);
                if(right < left){
                    p_next_i.x = p_i.x - step*sin(grad);
                    p_next_i.y = p_i.y + step*cos(grad);
                }
                else{
                    p_next_i.x = p_i.x - step*sin(grad);
                    p_next_i.y = p_i.y - step*cos(grad);
                }
            }
        }
        else{
            
            gy = (up-down);
            gx = (right-left);

            grad = fabs(atan(gy/gx));

            //printf("grad: %lf\n", grad*180/PI);

            if(up > down){
                p_next_i.x = p_i.x + step*sin(grad);
            }
            else{
                p_next_i.x = p_i.x - step*sin(grad);
            }
            if(right > left){
                p_next_i.y = p_i.y - step*cos(grad);
            }
            else{
                p_next_i.y = p_i.y + step*cos(grad);
            }
        }
        //printf("up: %lf, down: %lf, right: %lf, left: %lf\n", up, down, right, left);


        path = appendToPath(path, &size, p_next_i);
        //printf("[%lf, %lf]\n", p_next_i.x, p_next_i.y);

        p_i.x = p_next_i.x;
        p_i.y = p_next_i.y;

        counter++;
        //printf("%i\n", counter);

        /*
        if (right == INFINITY){
            right = arrival_time_map[x][y];
            printf("right inf\n");
        }
        if(left == INFINITY){
            left = arrival_time_map[x][y];
            printf("left inf\n");
        }
        if(up == INFINITY){
            up = arrival_time_map[x][y];
            printf("up inf\n");
        }
        if(down == INFINITY){
            down = arrival_time_map[x][y];
            printf("down inf\n");
        }
        */
    }

    
    
    for (int i = 0; i < size; i++) {
        printf("%lf, %lf\n", path[i].x, path[i].y);
    }

    free(path);
    
}

int main(){

    double **obs_grid = malloc(GRID_SIZE * sizeof(double *));
    for (int i = 0; i < GRID_SIZE; i++) {
        obs_grid[i] = malloc(GRID_SIZE * sizeof(double));
    }

    read_matrix("C:/Users/cheto/Desktop/UC3M/4th_Course/TFG/Maps/mapas_fer/TESTS/FM2_500x500.csv", obs_grid, GRID_SIZE, GRID_SIZE);

    float *times_ma = (float *)malloc(GRID_SIZE * GRID_SIZE * sizeof(float));
    float *times_ma_2 = (float *)malloc(GRID_SIZE * GRID_SIZE * sizeof(float));
    float *velocities_ma = (float *)malloc(GRID_SIZE * GRID_SIZE * sizeof(float));

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            times_ma[i*GRID_SIZE + j] = INFINITY;
            times_ma_2[i*GRID_SIZE + j] = INFINITY;
            velocities_ma[i*GRID_SIZE + j] = 1;
        }
    }

    clock_t start, end;
    double time_elapsed;

    int size_s = 0;
    int obs_size = 0;
    Point* sources = NULL;
    //Point source1 = {1, 1};
    //sources = appendToArray(sources, &size_s, source1);

    Point target = {598, 598};

    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            if(obs_grid[i][j] != 0.0){
                Point source1 = {i, j};
                sources = appendToArray(sources, &size_s, source1);
                obs_size = obs_size + 1;
            }
        }
    }
    

    double **obstacles = malloc(GRID_SIZE * sizeof(double *));
    for (int i = 0; i < GRID_SIZE; i++) {
        obstacles[i] = malloc(GRID_SIZE * sizeof(double));
    }

    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            obstacles[i][j] = 0;
        }
    }


    

    start = clock();

    FMM(obstacles, velocities_ma, times_ma, sources, obs_size, target);

    end = clock();
    time_elapsed = (double) (end - start)/CLOCKS_PER_SEC; 

    double max = 0;

    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            //printf("%.2f ", times_ma[i*GRID_SIZE + j]);
            if(times_ma[i*GRID_SIZE + j] > max){
                max = times_ma[i*GRID_SIZE + j];
            }
        }
        //printf("\n");
    }
    /*
    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            times_ma[i*GRID_SIZE + j] = times_ma[i*GRID_SIZE + j]/max;
        }
    }
    
    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            printf("%.2f ", times_ma[i*GRID_SIZE + j]);
        }
        printf("\n");
    }*/

    
    int size_s2 = 0;
    Point* sources2 = NULL;
    Point source2 = {1, 1};
    sources2 = appendToArray(sources2, &size_s2, source2);


    FMM(obs_grid, times_ma, times_ma_2, sources2, 1, target);


    
    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            printf("%.2f ", times_ma_2[i*GRID_SIZE + j]);
        }
        printf("\n");
    }
    
    

    pathPoint pathSource = {50, 50};
    pathPoint pathTarget = {50, 97};

    //findPath(times_ma, pathSource, pathTarget, 0.8);
    
    free(times_ma);
    free(velocities_ma);
    free(obs_grid);

    printf("Finished!\nFMM Execution time: %lf\n", time_elapsed);

    return 0;
}
