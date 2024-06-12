#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define N 2
#define H 1
#define GRID_SIZE 20

typedef struct {
    int x;
    int y;
} Point;

typedef struct {
    float min_time;
    int dimension;
} TValue;


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


void FMM(bool *Grid_MA, float *V_field_MA, Point *sources, int num_sources, float *Time_MA) {
    
    bool *Unknown_ma = (bool *)malloc(GRID_SIZE * GRID_SIZE * sizeof(bool));
    bool *Frozen_ma = (bool *)malloc(GRID_SIZE * GRID_SIZE * sizeof(bool));

    int size = 0;
    Point* Narrow_band = NULL;
    Point Narrow_point = {0, 0};

    printf("Start!\n");

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            Unknown_ma[i*GRID_SIZE + j] = 1;
            Frozen_ma[i*GRID_SIZE + j] = 0;
            Time_MA[i*GRID_SIZE + j] = INFINITY;
        }
    }

    Point source = {0, 0};
    for (int k = 0; k < num_sources; k++) {
        source.x = sources[k].x;
        source.y = sources[k].y;
        Time_MA[source.x*GRID_SIZE + source.y] = 0;
        Unknown_ma[source.x*GRID_SIZE + source.y] = 0;
        Narrow_band = appendToArray(Narrow_band, &size, source);
    }

    while (Narrow_band != NULL) { // MIENTRAS QUE EL ARRAY NARROW NO ESTÉ VACÍO (0 EN TODAS LAS POSICIONES)

        float min_t = -1;
        Point x_min = {-1, -1};

        // BUSCAMOS EL PUNTO DEL FRENTE DE ONDA (NARROW) CON MENOR TIEMPO PARA ANALIZARLO (EXPANDIR LA ONDA)
        
        for(int x=0; x < size; x++){
            if (min_t == -1 || Time_MA[Narrow_band[x].x*GRID_SIZE + Narrow_band[x].y] < min_t) {
                min_t = Time_MA[Narrow_band[x].x*GRID_SIZE + Narrow_band[x].y];
                x_min.x = Narrow_band[x].x;
                x_min.y = Narrow_band[x].y;      
            }
        }
        
        Point neighbors[4];
        
        int n_num = 0;
        int nx, ny;
        for(int dx = -1; dx <= 1; dx = dx+2){
            nx = x_min.x + dx;
            ny = x_min.y;
            if((0 <= nx) && (nx < GRID_SIZE) && (0 <= ny) && (ny < GRID_SIZE)){
                neighbors[n_num].x = nx;
                neighbors[n_num].y = ny;
                n_num++;
            }
        }
        for(int dy = -1; dy <= 1; dy = dy+2){
            nx = x_min.x;
            ny = x_min.y + dy;
            if((0 <= nx) && (nx < GRID_SIZE) && (0 <= ny) && (ny < GRID_SIZE)){
                neighbors[n_num].x = nx;
                neighbors[n_num].y = ny;
                n_num++;
            }
        }

        for(int i = 0; i < n_num; i++) {
            Point neighbor = neighbors[i];
            if (Frozen_ma[neighbor.x*GRID_SIZE + neighbor.y] == 0 && Grid_MA[neighbor.x*GRID_SIZE + neighbor.y] == 0) { //ANALIZAMOS CADA VECINO SIEMPRE QUE NO SEA UN OBSTÁCULO O YA LO HAYAMOS ANALIZADO ANTES

                float new_Time = solveEikonal(neighbor, V_field_MA, Time_MA); //OBTENEMOS NUEVO TIEMPO PARA EL VECINO

                if (new_Time < Time_MA[neighbor.x*GRID_SIZE + neighbor.y]) {
                    Time_MA[neighbor.x*GRID_SIZE + neighbor.y] = new_Time;
                }
                if (Unknown_ma[neighbor.x*GRID_SIZE + neighbor.y] == 1) {
                    Narrow_band = appendToArray(Narrow_band, &size, neighbor);
                    Unknown_ma[neighbor.x*GRID_SIZE + neighbor.y] = 0;
                }
            }
        }

        Narrow_band = removeFromArray(Narrow_band, &size, x_min);
        Frozen_ma[x_min.x*GRID_SIZE + x_min.y] = 1;
    }

    free(Narrow_band);
    free(Unknown_ma);
    free(Frozen_ma);
}

int main(){

    clock_t start, end;
    double time_elapsed;

    start = clock();
    
    // ----------------------------------------- SET-UP -----------------------------------------

    float *times_ma = (float *)malloc(GRID_SIZE * GRID_SIZE * sizeof(float));
    float *times_ma_2 = (float *)malloc(GRID_SIZE * GRID_SIZE * sizeof(float));
    float *velocities_ma = (float *)malloc(GRID_SIZE * GRID_SIZE * sizeof(float));
    bool *grid_ma = (bool *)malloc(GRID_SIZE * GRID_SIZE * sizeof(bool));

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            grid_ma[i*GRID_SIZE + j] = 0;
            times_ma[i*GRID_SIZE + j] = 0;
            times_ma_2[i*GRID_SIZE + j] = 0;
            velocities_ma[i*GRID_SIZE + j] = 1;
        }
    }


    //--------------------------------------- OBSTACLES ----------------------------------------------

    int size = 0;
    Point* obstacles = NULL;
    Point obstacle = {0, 0};

    for(int i=0; i<GRID_SIZE; i++){// OBSTACLES
        for(int j=0; j<GRID_SIZE; j++){
            if ((i>5) && (i<15) && (j>5) && (j<15)){
                if(i=j){ // || ((15-i) == j) || (15-j == i)
                    obstacle.x = i;
                    obstacle.y = j;
                    grid_ma[i*GRID_SIZE + j] = 1;
                    obstacles = appendToArray(obstacles, &size, obstacle);
                }
            }
        }
    }

    //------------------------------------ FIRST FAST MARCHING ---------------------------------

    FMM(grid_ma, velocities_ma, obstacles, size, times_ma); 
    
    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            printf("%.2f ", times_ma[i*GRID_SIZE + j]);
        }
        printf("\n");
    }


    //------------------------------------- SATURATION -----------------------------------------

    double exp = 2;
    float max = 0;
    float min = INFINITY;

    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            times_ma[i*GRID_SIZE + j] = pow(times_ma[i*GRID_SIZE + j], exp);
            if(times_ma[i*GRID_SIZE + j] < min){
                min = times_ma[i*GRID_SIZE + j];
            }
            if(times_ma[i*GRID_SIZE + j] > max){
                max = times_ma[i*GRID_SIZE + j];
            }
        }
    }


    //------------------------------------- NORMALIZATION -----------------------------------------

    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            times_ma[i*GRID_SIZE + j] = (times_ma[i*GRID_SIZE + j] - min) / (max - min);
        }
    }

    printf("\n Min: %lf,   Max: %lf\n\n", min, max);
    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            printf("%.4f ", times_ma[i*GRID_SIZE + j]);
        }
        printf("\n");
    }


    //----------------------------------- SECOND FAST MARCHING -------------------------------------

    int size_s = 0;
    Point* sources = NULL;
    Point source = {1, 1};
    sources = appendToArray(sources, &size_s, source);


    FMM(grid_ma, times_ma, sources, size_s, times_ma_2);
    
    for(int i=0; i<GRID_SIZE; i++){
        for(int j=0; j<GRID_SIZE; j++){
            printf("%.3f ", times_ma_2[i*GRID_SIZE + j]);
        }
        printf("\n");
    }

    // ----------------------------------------------------------------------------------------------
    
   
    free(obstacles);
    free(times_ma);
    free(times_ma_2);
    free(velocities_ma);
    free(grid_ma);

    end = clock();
    time_elapsed = (double) (end - start)/CLOCKS_PER_SEC;

    printf("Finished!\n");

    printf("FMM Execution time: %lf\n", time_elapsed);

     

    return 0;
}
