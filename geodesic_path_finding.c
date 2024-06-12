#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define N 2
#define H 1
#define GRID_SIZE 50
#define PI 3.14159265359

double grid[GRID_SIZE][GRID_SIZE];
double times[GRID_SIZE][GRID_SIZE];
double velocities[GRID_SIZE][GRID_SIZE];

typedef struct {
    int x;
    int y;
} Point;

typedef struct {
    double x;
    double y;
} pathPoint;

typedef struct {
    double min_time;
    int dimension;
} TValue;


int is_any_zero(int map[GRID_SIZE][GRID_SIZE]) {
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            if (map[i][j] != 0) {
                return 1;
            }
        }
    }
    return 0;
}

double minTDim(Point point, double Time[GRID_SIZE][GRID_SIZE], int dim) { // DEVOLVEMOS EL VALOR MÁS PEQUEÑO DE TIEMPO EN EL ARRAY TIME PARA CADA DIMENSIÓN (UP-DOWN, LEFT-RIGHT)
    int i = point.x, j = point.y;
    if (dim == 1) {
        int nx_left = i - 1;
        int nx_right = i + 1;
        if (0 <= nx_left && nx_left < GRID_SIZE && 0 <= nx_right && nx_right < GRID_SIZE) {
            return fmin(Time[i + 1][j], Time[i - 1][j]);
        } else if (0 > nx_left || nx_left >= GRID_SIZE) {
            return Time[i + 1][j];
        } else if (0 > nx_right || nx_right >= GRID_SIZE) {
            return Time[i - 1][j];
        }
    } else if (dim == 2) {
        int ny_top = j + 1;
        int ny_bottom = j - 1;
        if (0 <= ny_top && ny_top < GRID_SIZE && 0 <= ny_bottom && ny_bottom < GRID_SIZE) {
            return fmin(Time[i][j + 1], Time[i][j - 1]);
        } else if (0 > ny_top || ny_top >= GRID_SIZE) {
            return Time[i][j - 1];
        } else if (0 > ny_bottom || ny_bottom >= GRID_SIZE) {
            return Time[i][j + 1];
        }
    }
    return 0.0;
}

double solveNDims(Point point, int dim, TValue T_Values[2], double V_field[GRID_SIZE][GRID_SIZE]) {
    //printf("Dentro de solveNDims\n"); //-----------------------------------------------------------------
    
    int num_dim = 0;
    int i;

    for(i = 0; i < 2; i++){ //DESCUBRIR QUÉ T_VALUE ES EL QUE NO VALE, PARA ANALIZAR EL OTRO.
        if(T_Values[i].dimension == 0){
            num_dim = 1;
            break;
        }
    }

    if (num_dim == 1) { // CASO PARA CUANDO SÓLO HAYA VALORES EN UNA DE LAS DIMENSIONES.
        //printf("Dentro del caso num_dim = 1\n");//-----------------------------------------------------------------
        //printf("%i\n", i);//-----------------------------------------------------------------
        if(i == 0){
            //printf("1. min_time: %f \n", T_Values[1].min_time + H/V_field[point.x][point.y]);//-----------------------------------------------------------------
            return T_Values[1].min_time + H/V_field[point.x][point.y];
        }
        else if(i == 1){
            //printf("2. min_time: %f \n", T_Values[0].min_time + H/V_field[point.x][point.y]);//-----------------------------------------------------------------
            return T_Values[0].min_time + H/V_field[point.x][point.y];
        }
    }

    // CASO PARA CUANDO HAY VALORES VÁLIDOS EN AMBAS DIMENSIONES Y HAY QUE HACER LOS SIGUIENTES CÁLCULOS PARA EL NUEVO TIEMPO
    double sum_T = 0.0, sum_T2 = 0.0;

    for(int i = 0; i < 2; i++) {
        sum_T = sum_T + T_Values[i].min_time;
        sum_T2 = sum_T2 + (T_Values[i].min_time * T_Values[i].min_time);
    }

    double a = dim;
    double b = -2 * sum_T;
    double c = sum_T2 - pow(H, 2) / V_field[point.x][point.y];
    double q = b * b - 4 * a * c;

    if (q < 0) {
        return INFINITY;
    } else {
        return (-b + sqrt(q)) / (2 * a);
    }

    return 0.0;
}

double solveEikonal(Point point, double Time[GRID_SIZE][GRID_SIZE], double V_field[GRID_SIZE][GRID_SIZE]) {
    //printf("Dentro de solveEikonal\n");//-----------------------------------------------------------------
    
    int a = N;
    double new_Time = 0.0;

    TValue T_Values[2] = {{INFINITY, 0}, {INFINITY, 0}}; //INICIALIZAMOS T_VALUES

    for (int dim = 1; dim <= N; dim++) {
        double min_T = minTDim(point, Time, dim);

        //printf("Dimension %i, min time: %f\n", dim, min_T);//-----------------------------------------------------------------

        if (min_T != INFINITY && min_T < Time[point.x][point.y]) { // METER EN T_VALUES LOS TIEMPOS QUE SEAN MENORES QUE LOS ANTERIORES
            //printf("New minimum time found! \n");//-----------------------------------------------------------------
            for(int i = 0; i < 2; i++){
                if(T_Values[i].min_time == INFINITY){
                    T_Values[i].min_time = min_T;
                    T_Values[i].dimension = dim;

                    break;
                }
            }
            //printf("T_Values: [val: %f, dim: %i] [val: %f, dim: %i]\n", T_Values[0].min_time, T_Values[0].dimension, T_Values[1].min_time, T_Values[1].dimension);//-----------------------------------------------------------------
            //printf("time x: %f \n", Time[point.x][point.y]);
            //printf("enter\n");
            
        } else {
            a--;
        }
    }
    
    if (a == 0) {
        //printf("a = 0\n");//-----------------------------------------------------------------
        return INFINITY;
    }
    
    for (int dim = 1; dim <= a; dim++) {
        //printf("Solving for new time...\n");//-----------------------------------------------------------------
        new_Time = solveNDims(point, dim, T_Values, V_field);
        if(dim == a){
            break;
        }
    }
    
    return new_Time;
}

void visualize(double map[GRID_SIZE][GRID_SIZE]) {
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            printf("%.2f ", map[i][j]);
        }
        printf("\n");
    }
}

void visualize_narrow(int map[GRID_SIZE][GRID_SIZE]) {
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            printf("%i ", map[i][j]);
        }
        printf("\n");
    }
}

void FMM(double grid[GRID_SIZE][GRID_SIZE], double Time[GRID_SIZE][GRID_SIZE], double V_field[GRID_SIZE][GRID_SIZE], Point sources, int num_sources) {
    int Unknown[GRID_SIZE][GRID_SIZE];
    int Narrow[GRID_SIZE][GRID_SIZE];
    int Frozen[GRID_SIZE][GRID_SIZE];

    printf("Start!\n");

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            Unknown[i][j] = 1;
            Narrow[i][j] = 0;
            Frozen[i][j] = 0;
            Time[i][j] = INFINITY;
        }
    }

    for (int k = 0; k < num_sources; k++) {
        Point source = sources;
        Time[source.x][source.y] = 0;
        Unknown[source.x][source.y] = 0;
        Narrow[source.x][source.y] = 1;
    }

    int counter = 0;

    while (is_any_zero(Narrow)) { // MIENTRAS QUE EL ARRAY NARROW NO ESTÉ VACÍO (0 EN TODAS LAS POSICIONES)
        //printf("inside the loop\n");//-----------------------------------------------------------------
        double min_t = -1;
        Point x_min = {-1, -1};


        // BUSCAMOS EL PUNTO DEL FRENTE DE ONDA (NARROW) CON MENOR TIEMPO PARA ANALIZARLO (EXPANDIR LA ONDA)
        for (int x = 0; x < GRID_SIZE; x++) {
            for (int y = 0; y < GRID_SIZE; y++) {
                if (Narrow[x][y] == 1) {
                    if (min_t == -1 || Time[x][y] < min_t) {
                        min_t = Time[x][y];
                        x_min.x = x;
                        x_min.y = y;
                        
                    }
                }
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

        //printf("Point: [%i, %i]\n", x_min.x, x_min.y);//-----------------------------------------------------------------
        //printf("Neighbor 1: [%i, %i]\n", neighbors[0].x, neighbors[0].y);//-----------------------------------------------------------------
        //printf("Neighbor 2: [%i, %i]\n", neighbors[1].x, neighbors[1].y);//-----------------------------------------------------------------
        //printf("Neighbor 3: [%i, %i]\n", neighbors[2].x, neighbors[2].y);//-----------------------------------------------------------------
        //printf("Neighbor 4: [%i, %i]\n", neighbors[3].x, neighbors[3].y);//-----------------------------------------------------------------

        for(int i = 0; i < n_num; i++) {
            Point neighbor = neighbors[i];
            if (Frozen[neighbor.x][neighbor.y] == 0 && grid[neighbor.x][neighbor.y] == 0) { //ANALIZAMOS CADA VECINO SIEMPRE QUE NO SEA UN OBSTÁCULO O YA LO HAYAMOS ANALIZADO ANTES
                
                //printf("Calculando new time para vecino %i... [%i, %i] \n", i+1, neighbor.x, neighbor.y);//-----------------------------------------------------------------
                
                double new_Time = solveEikonal(neighbor, Time, V_field); //OBTENEMOS NUEVO TIEMPO PARA EL VECINO

                //printf("new_Time: %0.2f \n", new_Time);//-----------------------------------------------------------------

                if (new_Time < Time[neighbor.x][neighbor.y]) {
                    Time[neighbor.x][neighbor.y] = new_Time;
                }
                if (Unknown[neighbor.x][neighbor.y] == 1) {
                    Narrow[neighbor.x][neighbor.y] = 1;
                    Unknown[neighbor.x][neighbor.y] = 0;
                }
            }
        }

        Narrow[x_min.x][x_min.y] = 0;
        Frozen[x_min.x][x_min.y] = 1;

        
        counter++;
        if(counter % 1 == 0){ //VER LA EXPANSIÓN PROGRESIVAMENTE
            //visualize(Time);
            //printf("\n");
        }

        
    }

}

pathPoint* appendToArray(pathPoint *arr, int *size, pathPoint element) { 
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


void findPath(double arrival_time_map[GRID_SIZE][GRID_SIZE], pathPoint source, pathPoint target, double step){
    
    int counter = 0;

    int size = 0;
    pathPoint* path = NULL;

    path = appendToArray(path, &size, target);

    
    pathPoint p_i = {target.x, target.y};
    pathPoint p_next_i;
    int x = target.x;
    int y = target.y;

    double right, left, up, down;
    double grad;
    double gx, gy;

    //(floor(p_i.x) != source.x) && (floor(p_i.y) != source.y) ||  && (counter < 20)

    while(((ceil(p_i.x) != source.x) || (ceil(p_i.y) != source.y)) && (counter < 500)){
        /*
                        _________>  Y
                        |
                        |
                        |
                    X  \/
        */

        int x = floor(p_i.x);
        int y = floor(p_i.y);

        printf("Current point: [%i, %i];   Time: %lf\n", x, y, arrival_time_map[x][y]);
        
        
        right = (p_i.y + 1 < GRID_SIZE) ? arrival_time_map[x][y + 1] : 0;
        left = (p_i.y - 1 >= 0) ? arrival_time_map[x][y - 1] : 0;
        up = (p_i.x - 1 >= 0) ? arrival_time_map[x - 1][y] : 0;
        down = (p_i.x + 1 < GRID_SIZE) ? arrival_time_map[x + 1][y] : 0;

        printf("up: %lf, down: %lf, right: %lf, left: %lf\n", up, down, right, left);
        
        
        if(right == INFINITY || left == INFINITY || up == INFINITY || down == INFINITY){
            if (right == INFINITY){
                right = arrival_time_map[x][y];
                printf("right inf\n");
                gy = (up-down);
                gx = (right-left);
                grad = fabs(atan(gy/gx));
                printf("grad: %lf\n", grad*180/PI);
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
                left = arrival_time_map[x][y];
                printf("left inf\n");
                gy = (up-down);
                gx = (right-left);
                grad = fabs(atan(gy/gx));
                printf("grad: %lf\n", grad*180/PI);
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
                up = arrival_time_map[x][y];
                printf("up inf\n");
                gy = (up-down);
                gx = (right-left);
                grad = fabs(atan(gy/gx));
                printf("grad: %lf\n", grad*180/PI);
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
                down = arrival_time_map[x][y];
                printf("down inf\n");
                gy = (up-down);
                gx = (right-left);
                grad = fabs(atan(gy/gx));
                printf("grad: %lf\n", grad*180/PI);
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

            printf("grad: %lf\n", grad*180/PI);

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
        printf("up: %lf, down: %lf, right: %lf, left: %lf\n", up, down, right, left);


        path = appendToArray(path, &size, p_next_i);
        printf("[%lf, %lf]\n", p_next_i.x, p_next_i.y);

        p_i.x = p_next_i.x;
        p_i.y = p_next_i.y;

        counter++;
        printf("%i\n", counter);

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


// ----------- findPath but i return the final path to the main, instead of printing it inside the function itself. -------------
/*

Point* findPath(double arrival_time_map[GRID_SIZE][GRID_SIZE], Point source, Point target, double step, int *final_size){
    
    //Point path[] = {{source.x, source.y}};
    int size = 0;
    double *path_x = NULL;
    double *path_y = NULL;

    path_x = appendToArray(path_x, &size, source.x);
    path_y = appendToArray(path_y, &size, source.y);


    Point p_i = {target.x, target.y};
    Point p_next_i;
    int x = target.x;
    int y = target.y;

    double right, left, up, down;
    double grad;
    double gx, gy;

    while((floor(p_i.x) != source.x) && (floor(p_i.y) != source.y)){
        
                        _________>  Y
                        |
                        |
                        |
                    X  \/

        
        right = (p_i.y + 1 < GRID_SIZE) ? arrival_time_map[p_i.x][p_i.y + 1] : 0;
        left = (p_i.y - 1 >= 0) ? arrival_time_map[p_i.x][p_i.y - 1] : 0;
        up = (p_i.x - 1 >= 0) ? arrival_time_map[p_i.x - 1][p_i.y] : 0;
        down = (p_i.x + 1 >= GRID_SIZE) ? arrival_time_map[p_i.x  + 1][p_i.y] : 0;

        if (right == INFINITY || left == INFINITY){
            if(up < down){
                grad = PI/2;
                p_next_i.x = p_i.x - step*sin(grad);
                p_next_i.y = p_i.y;
            }
            else{
                grad = -PI/2;
                p_next_i.x = p_i.x - step*sin(grad);
                p_next_i.y = p_i.y;
            }
        }
        else if(up == INFINITY || down == INFINITY){
            if(right < left){
                grad = PI;
                p_next_i.x = p_i.x;
                p_next_i.y = p_i.y - step*cos(grad);
            }
            else{
                grad = 0;
                p_next_i.x = p_i.x;
                p_next_i.y = p_i.y - step*cos(grad);
            }
        }
        else{
            gy = (up-down);
            gx = (right-left);

            grad = fabs(atan(gy/gx));

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

        path_x = appendToArray(path_x, &size, p_next_i.x);
        path_y = appendToArray(path_y, &size, p_next_i.y);

        p_i.x = p_next_i.x;
        p_i.y = p_next_i.y;


    }

    *final_size = size;

    Point* path = (Point*)malloc(size * sizeof(Point));
    
    if (path == NULL) {
        // Handle memory allocation failure
        printf("Memory allocation failed.\n");
        exit(1);
    }
    
    // Populate the array of structs with some data
    for (int i = 0; i < size; i++) {
        path[i].x = path_x[i];
        path[i].y = path_y[i];
    }

    free(path_x);
    free(path_y);
    
    // Return a pointer to the array of structs
    return path;
}

*/






int main(){

    Point source = {48, 48};

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            grid[i][j] = 0;
            times[i][j] = 0;
            velocities[i][j] = 1;
        }
    }

    
    for(int i=0; i<35; i++){ // OBSTACLES
        for(int j=0; j<10; j++){
            grid[i][15] = 1;
            grid[i+1][16] = 1;
            grid[i+2][17] = 1;
            grid[i+3][18] = 1;
            grid[i+4][19] = 1;
            grid[i+5][20] = 1;
            grid[i+6][20+j] = 1;

            grid[i+15][37] = 1;
        }
    }

    FMM(grid, times, velocities, source, 1);

    visualize(times);

    pathPoint pathSource = {48, 48};
    pathPoint pathTarget = {1, 1};

    findPath(times, pathSource, pathTarget, 1); 
 
    printf("Finished!");

    return 0;
}
