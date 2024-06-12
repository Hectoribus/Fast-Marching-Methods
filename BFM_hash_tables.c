#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define N 2
#define H 1
#define GRID_SIZE 1004

typedef struct {
    int x;
    int y;
    bool type;
} Point;

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
    bool type;
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

// Define a createHeap function
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

// Defining insertHelper function
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

// Define an insert function
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

bool containsPoint(heap* h, int x, int y, int type)
{
    for (int i = 0; i < h->size; i++) {
        if ((abs(h->arr[i].x - x) <= 1) && (abs(h->arr[i].y - y) <= 1) && (h->arr[i].type == type)) {
            printf("Connection points: (%i, %i)     (%i, %i)\n", x, y, h->arr[i].x, h->arr[i].y);
            return true; // Point found
        }
    }
    return false; // Point not found
}

/*
int arePointsClose(int x0, int y0, int x1, int y1) {
    int dx = abs(x0 - x1);
    int dy = abs(y0 - y1);
    return (dx <= 1 && dy <= 1);
}
*/

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







//-----------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------          H·A·S·H    T·A·B·L·E·S       -------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------

// Hash table key structure (for coordinates)
typedef struct {
    int x;
    int y;
} PointKey;


// Hash table node structure
typedef struct HashNode {
    PointKey key;
    bool value;
    struct HashNode* next;
} HashNode;


// Hash table structure
typedef struct {
    HashNode** buckets;
    size_t size;
} HashTable;




// Hash function to compute the hash value of a PointKey
size_t hash_function(PointKey key, size_t size) {
    return (size_t)((key.x * 73856093) ^ (key.y * 19349663)) % size;
}


// Function to create a new hash table
HashTable* create_table(size_t size) {
    HashTable* table = (HashTable*)malloc(sizeof(HashTable));
    table->buckets = (HashNode**)calloc(size, sizeof(HashNode*));
    table->size = size;
    return table;
}


// Function to insert a point into the hash table
void insert_point(HashTable* table, int x, int y) {
    PointKey key = {x, y};
    size_t hash_index = hash_function(key, table->size);
    HashNode* new_node = (HashNode*)malloc(sizeof(HashNode));
    new_node->key = key;
    new_node->value = true;
    new_node->next = table->buckets[hash_index];
    table->buckets[hash_index] = new_node;
}


// Function to delete a point from the hash table
void delete_point(HashTable* table, int x, int y) {
    PointKey key = {x, y};
    size_t hash_index = hash_function(key, table->size);
    HashNode* current = table->buckets[hash_index];
    HashNode* prev = NULL;
    while (current) {
        if (current->key.x == x && current->key.y == y) {
            if (prev) {
                prev->next = current->next;
            } else {
                table->buckets[hash_index] = current->next;
            }
            free(current);
            return;
        }
        prev = current;
        current = current->next;
    }
}


// Function to check if a point exists in the hash table
bool contains_point_hash(HashTable* table, int x, int y) {
    PointKey key = {x, y};
    size_t hash_index = hash_function(key, table->size);
    HashNode* current = table->buckets[hash_index];
    while (current) {
        if (current->key.x == x && current->key.y == y) {
            return true;
        }
        current = current->next;
    }
    return false;
}


// Function to free the hash table
void free_table(HashTable* table) {
    for (size_t i = 0; i < table->size; ++i) {
        HashNode* current = table->buckets[i];
        while (current) {
            HashNode* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(table->buckets);
    free(table);
}


void process_point(Point x_min, HashTable* narrow_band_wave0, HashTable* narrow_band_wave1) {
    // Check for neighboring points in the other wavefront's hash table
    for (int dx = -1; dx <= 1; dx=dx+2) {
        for (int dy = -1; dy <= 1; dy=dy+2) {
            int nx = x_min.x + dx;
            int ny = x_min.y + dy;
            if (x_min.type == 0) {
                if (contains_point_hash(narrow_band_wave1, nx, ny)) {
                    printf("Connection found at (%d, %d)\n", nx, ny);
                    // Handle connection (e.g., stop the algorithm)
                    exit(0);
                }
            } else {
                if (contains_point_hash(narrow_band_wave0, nx, ny)) {
                    printf("Connection found at (%d, %d)\n", nx, ny);
                    // Handle connection (e.g., stop the algorithm)
                    exit(0);
                }
            }
        }
    }
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

/*
bool insideNarrowBand(heap* HEAP_Narrow_Band, int size, int nx, int ny, int narrow_type){
    for(int point=0; point < size; point++){
        if((Narrow_Band[point].x == nx) && (Narrow_Band[point].y == ny) && (Narrow_Band[point].type == narrow_type)){
            return true;
        }
    }
    return false;
}
*/

bool get_neighbors(Point x_min, Point neighbors[4], int *n_num, HashTable* narrow_band_wave0, HashTable* narrow_band_wave1){
    *n_num = 0;
    int nx, ny;

    int x = x_min.x;
    int y = x_min.y;
    bool type = x_min.type;

    //printf("neighbors:  ");
    
    for(int dx = -1; dx <= 1; dx = dx+2){
        nx = x + dx;
        ny = y;
        if((0 <= nx) && (nx < GRID_SIZE) && (0 <= ny) && (ny < GRID_SIZE)){
            neighbors[*n_num].x = nx;
            neighbors[*n_num].y = ny;
            neighbors[*n_num].type = type;
            //printf("(%i, %i)   ", nx, ny);
            (*n_num)++;
            if(type == 0 && contains_point_hash(narrow_band_wave1, nx, ny)){
                printf("Connection points: (%i, %i)     (%i, %i)\n", x_min.x, x_min.y, nx, ny);
                return true;
            }
            else if(type == 1 && contains_point_hash(narrow_band_wave0, nx, ny)){
                printf("Connection points: (%i, %i)     (%i, %i)\n", x_min.x, x_min.y, nx, ny);
                return true;
            }
        }
    }
    for(int dy = -1; dy <= 1; dy = dy+2){
        nx = x;
        ny = y + dy;
        if((0 <= nx) && (nx < GRID_SIZE) && (0 <= ny) && (ny < GRID_SIZE)){
            neighbors[*n_num].x = nx;
            neighbors[*n_num].y = ny;
            neighbors[*n_num].type = type;
            //printf("(%i, %i)   ", nx, ny);
            (*n_num)++;
            if(type == 0 && contains_point_hash(narrow_band_wave1, nx, ny)){
                printf("Connection points: (%i, %i)     (%i, %i)\n", x_min.x, x_min.y, nx, ny);
                return true;
            }
            else if(type == 1 && contains_point_hash(narrow_band_wave0, nx, ny)){
                printf("Connection points: (%i, %i)     (%i, %i)\n", x_min.x, x_min.y, nx, ny);
                return true;
            }
        }
    }
    //printf("\n");

    return false;
}


float calculateDistance(int x0, int y0, int x1, int y1){
    float a = abs(x0 - x1)^2;
    float b = abs(y0 - y1)^2;

    float c = sqrt(a + b);

    return c;
}


void FMM(double **Grid_MA, float *V_field_MA, Point *sources, int num_sources, float *Time_MA, Point connection_points[2]) {
    
    bool *Unknown_ma = (bool *)malloc(GRID_SIZE * GRID_SIZE * sizeof(bool));
    bool *Frozen_ma = (bool *)malloc(GRID_SIZE * GRID_SIZE * sizeof(bool));

    int size = 0;

    heap* HEAP_Narrow_band = createHeap(0);

    // Initialize hash tables for both wavefronts
    size_t initial_table_size = 1024;  // Starting size, can be adjusted based on expected usage
    HashTable* narrow_band_wave0 = create_table(initial_table_size);
    HashTable* narrow_band_wave1 = create_table(initial_table_size);

    printf("Start!\n");

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            Unknown_ma[i*GRID_SIZE + j] = 1;
            Frozen_ma[i*GRID_SIZE + j] = 0;
            Time_MA[i*GRID_SIZE + j] = INFINITY;
        }
    }

    Point source = {0, 0, 0};

    for (int k = 0; k < num_sources; k++) {
        source.x = sources[k].x;
        source.y = sources[k].y;
        source.type = sources[k].type;
        Time_MA[source.x*GRID_SIZE + source.y] = 0;
        Unknown_ma[source.x*GRID_SIZE + source.y] = 0;
        insert(HEAP_Narrow_band, (Narrow_Point){source.x, source.y, source.type, Time_MA[source.x*GRID_SIZE + source.y]});
        if (source.type == 0) {
            insert_point(narrow_band_wave0, source.x, source.y);
        } else {
            insert_point(narrow_band_wave1, source.x, source.y);
        }
    }

    float min_t = -1;
    Point x_min = {-1, -1, 0};
    long counter_points = 0;
    bool connection = 0;

    Narrow_Point minPoint;

    while ((HEAP_Narrow_band->size != 0) && (connection == 0) && (counter_points < 600000)) { // MIENTRAS QUE EL ARRAY NARROW NO ESTÉ VACÍO (0 EN TODAS LAS POSICIONES)

        minPoint = getMin(HEAP_Narrow_band);

        x_min.x = minPoint.x;
        x_min.y = minPoint.y;
        x_min.type = minPoint.type;

        //printf("X_MIN: (%i, %i)   ", x_min.x, x_min.y);


        // Add the point to the corresponding hash table
        if (x_min.type == 0) {
            insert_point(narrow_band_wave0, x_min.x, x_min.y);
        } else {
            insert_point(narrow_band_wave1, x_min.x, x_min.y);
        }

        
        //Point c_points[2];

        /*
        if(x_min.type == 0 && containsPoint(HEAP_Narrow_band, x_min.x, x_min.y, 1)){
            connection = 1;
        }
        */

        /*
        connection_points[0].x = c_points[0].x;
        connection_points[0].y = c_points[0].y;

        connection_points[1].x = c_points[1].x;
        connection_points[1].y = c_points[1].y;
        */

        counter_points = counter_points + 1;
        
        Point neighbors[4];
        int n_num = 0;

        connection = get_neighbors(x_min, neighbors, &n_num, narrow_band_wave0, narrow_band_wave1);


        for(int i = 0; i < n_num; i++) {
            Point neighbor = neighbors[i];
            if (Frozen_ma[neighbor.x*GRID_SIZE + neighbor.y] == 0 && Grid_MA[neighbor.x][neighbor.y] == 0) { //ANALIZAMOS CADA VECINO SIEMPRE QUE NO SEA UN OBSTÁCULO O YA LO HAYAMOS ANALIZADO ANTES

                float new_Time = solveEikonal(neighbor, V_field_MA, Time_MA); //OBTENEMOS NUEVO TIEMPO PARA EL VECINO

                if (new_Time < Time_MA[neighbor.x*GRID_SIZE + neighbor.y]) {
                    Time_MA[neighbor.x*GRID_SIZE + neighbor.y] = new_Time;
                }
                if (Unknown_ma[neighbor.x*GRID_SIZE + neighbor.y] == 1) {
                    insert(HEAP_Narrow_band, (Narrow_Point){neighbor.x, neighbor.y, neighbor.type, new_Time});
                    if (neighbor.type == 0) {
                        insert_point(narrow_band_wave0, neighbor.x, neighbor.y);
                    } else if(neighbor.type == 1) {
                        insert_point(narrow_band_wave1, neighbor.x, neighbor.y);
                    }
                    Unknown_ma[neighbor.x*GRID_SIZE + neighbor.y] = 0;
                }
            }
        }

        extractMin(HEAP_Narrow_band);
        if (x_min.type == 0) {
            delete_point(narrow_band_wave0, x_min.x, x_min.y);
        } else if(x_min.type == 1) {
            delete_point(narrow_band_wave1, x_min.x, x_min.y);
        }

        Frozen_ma[x_min.x*GRID_SIZE + x_min.y] = 1;
    }

    // Free the hash tables
    free_table(narrow_band_wave0);
    free_table(narrow_band_wave1);
    
    free(Unknown_ma);
    free(Frozen_ma);
    freeHeap(HEAP_Narrow_band);
    printf("Counter_Points: %li\n", counter_points);
}

int main(){

    double **obs_grid = malloc(GRID_SIZE * sizeof(double *));
    for (int i = 0; i < GRID_SIZE; i++) {
        obs_grid[i] = malloc(GRID_SIZE * sizeof(double));
    }

    read_matrix("C:/Users/cheto/Desktop/UC3M/4th_Course/TFG/Maps/mapas_fer/TESTS/city_map_1004x1004.csv", obs_grid, GRID_SIZE, GRID_SIZE);


    clock_t start, end;
    double time_elapsed;

    float *times_ma = (float *)malloc(GRID_SIZE * GRID_SIZE * sizeof(float));
    float *velocities_ma = (float *)malloc(GRID_SIZE * GRID_SIZE * sizeof(float));

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            times_ma[i*GRID_SIZE + j] = 0;
            velocities_ma[i*GRID_SIZE + j] = 1;
        }
    }


    int size_s = 0;
    Point* sources = NULL;
    Point source1 = {10, 10, 0};
    sources = appendToArray(sources, &size_s, source1);
    Point source2 = {994, 994, 1};
    sources = appendToArray(sources, &size_s, source2);

    Point connection_points[2] = {{-1, -1, 0}, {-1, -1, 0}};

    start = clock();

    FMM(obs_grid, velocities_ma, sources, size_s, times_ma, connection_points);

    end = clock();
    time_elapsed = (double) (end - start)/CLOCKS_PER_SEC; 

    
    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            printf("%.2f ", times_ma[i*GRID_SIZE + j]);
        }
        printf("\n");
    }
    

    //printf("\nConnection points: [%i, %i]   [%i, %i]\n", connection_points[0].x, connection_points[0].y, connection_points[1].x, connection_points[1].y);
    
    free(times_ma);
    free(velocities_ma);
    free(obs_grid);

    printf("Finished!\nFMM Execution time: %lf\n", time_elapsed);
    return 0;
}