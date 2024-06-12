#include <stdio.h>
#include "libxl.h" // Include the libxl library header

#define MAX_ROWS 30
#define MAX_COLS 30

int main() {
    BookHandle book = xlCreateBook(); // Create a book handle
    if (book) {
        const char* filename = "C:/Users/cheto/Desktop/UC3M/4th_Course/TFG/Maps/Book1.xlsx";
        if (xlBookLoad(book, filename)) { // Load Excel file
            SheetHandle sheet = xlBookGetSheet(book, 0); // Get the first sheet
            if (sheet) {
                int numRows = xlSheetLastRow(sheet);
                int numCols = xlSheetLastCol(sheet);
                
                // Ensure we don't exceed the array dimensions
                numRows = numRows > MAX_ROWS ? MAX_ROWS : numRows;
                numCols = numCols > MAX_COLS ? MAX_COLS : numCols;
                
                // Define a 2D array to store the data
                int data[MAX_ROWS][MAX_COLS];
                
                // Read Excel data into the array
                for (int i = 0; i < numRows; i++) {
                    for (int j = 0; j < numCols; j++) {
                        data[i][j] = xlSheetReadNum(sheet, i + 1, j + 1, NULL); // Read numeric value
                    }
                }
                
                // Print the values stored in the 2D array
                printf("Values read from Excel:\n");
                for (int i = 0; i < numRows; i++) {
                    for (int j = 0; j < numCols; j++) {
                        printf("%d\t", data[i][j]);
                    }
                    printf("\n");
                }
            } else {
                printf("Failed to open sheet.\n");
            }
        } else {
            printf("Failed to load Excel file.\n");
        }
        
        xlBookRelease(book); // Release the book handle
    } else {
        printf("Failed to create book.\n");
    }
    
    return 0;
}
