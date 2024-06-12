#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

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

int main() {
    int rows = 508, cols = 508;
    double **matrix = malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        matrix[i] = malloc(cols * sizeof(double));
    }

    read_matrix("C:/Users/cheto/Desktop/UC3M/4th_Course/TFG/Maps/mapas_fer/TESTS/city_map.csv", matrix, rows, cols);

    // Print the matrix to verify
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.f ", matrix[i][j]);
        }
        printf("\n");
    }

    // Free allocated memory
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}
