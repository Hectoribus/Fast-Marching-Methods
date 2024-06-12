/*
#include <stdio.h>
#include <stdlib.h>

// Define the Point struct
typedef struct {
    int x;
    int y;
    float time;
} Point;

// Declare a heap structure
typedef struct{
    Point* arr;
    int size;
    int capacity;
} heap;



// forward declarations
heap* createHeap(int capacity, Point* points);
void insertHelper(heap* h, int index);
void heapify(heap* h, int index);
Point extractMin(heap* h);
void insert(heap* h, Point data);

// Define a createHeap function
heap* createHeap(int capacity, Point* points)
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
    h->arr = (Point*)malloc(capacity * sizeof(Point));

    // Checking if memory is allocated to h or not
    if (h->arr == NULL) {
        printf("Memory error");
        return NULL;
    }
    int i;
    for (i = 0; i < capacity; i++) {
        h->arr[i] = points[i];
    }

    h->size = i;
    i = (h->size - 2) / 2;
    while (i >= 0) {
        heapify(h, i);
        i--;
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
        Point temp = h->arr[parent];
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
    // any of these is smaller that its parent
    if (left != -1 && h->arr[left].time < h->arr[index].time)
        min = left;
    if (right != -1 && h->arr[right].time < h->arr[min].time)
        min = right;

    // Swapping the nodes
    if (min != index) {
        Point temp = h->arr[min];
        h->arr[min] = h->arr[index];
        h->arr[index] = temp;

        // recursively calling for their child elements
        // to maintain min heap
        heapify(h, min);
    }
}

Point extractMin(heap* h)
{
    Point deleteItem;

    // Checking if the heap is empty or not
    if (h->size == 0) {
        printf("\nHeap is empty.");
        return (Point){-1, -1, -1}; // Return a dummy Point in case of empty heap
    }

    // Store the node in deleteItem that
    // is to be deleted.
    deleteItem = h->arr[0];

    // Replace the deleted node with the last node
    h->arr[0] = h->arr[h->size - 1];
    // Decrement the size of heap
    h->size--;

    // Call minheapify_top_down for 0th index
    // to maintain the heap property
    heapify(h, 0);
    return deleteItem;
}

// Define a insert function
void insert(heap* h, Point data)
{
    // Checking if heap is full or not
    if (h->size < h->capacity) {
        // Inserting data into an array
        h->arr[h->size] = data;
        // Calling insertHelper function
        insertHelper(h, h->size);
        // Incrementing size of array
        h->size++;
    }
}

void printHeap(heap* h)
{
    for (int i = 0; i < h->size; i++) {
        printf("%.2f    ",h->arr[i].time);
    }
    printf("\n");
}

int main()
{
    Point points[9] = { {9, 0, 9.0}, {8, 0, 8.0}, {7, 0, 7.0}, {6, 0, 6.0}, {5, 0, 5.0}, {4, 0, 4.0}, {3, 0, 3.0}, {2, 0, 2.0}, {1, 0, 1.0} };
    heap* hp = createHeap(15, points);

    printHeap(hp);
    extractMin(hp);
    printHeap(hp);
    insert(hp, (Point){10, 0, 6.5});
    printHeap(hp);
    insert(hp, (Point){10, 0, 6.5});
    printHeap(hp);
    insert(hp, (Point){10, 0, 6.5});
    printHeap(hp);
    insert(hp, (Point){10, 0, 6.5});
    printHeap(hp);
    insert(hp, (Point){10, 0, 6.5});
    printHeap(hp);
    insert(hp, (Point){10, 0, 6.5});
    printHeap(hp);


    return 0;
}




#include <stdio.h>
#include <stdlib.h>

// Define the Point struct
typedef struct {
    int x;
    int y;
    float time;
} Point;

// Declare a heap structure
struct Heap {
    Point* arr;
    int size;
    int capacity;
};

// define the struct Heap name
typedef struct Heap heap;

// forward declarations
heap* createHeap(int capacity);
void insertHelper(heap* h, int index);
void heapify(heap* h, int index);
Point extractMin(heap* h);
void insert(heap* h, Point data);
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
        h->arr = (Point*)malloc(capacity * sizeof(Point));
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
        Point temp = h->arr[parent];
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
    // any of these is smaller that its parent
    if (left != -1 && h->arr[left].time < h->arr[index].time)
        min = left;
    if (right != -1 && h->arr[right].time < h->arr[min].time)
        min = right;

    // Swapping the nodes
    if (min != index) {
        Point temp = h->arr[min];
        h->arr[min] = h->arr[index];
        h->arr[index] = temp;

        // recursively calling for their child elements
        // to maintain min heap
        heapify(h, min);
    }
}

Point extractMin(heap* h)
{
    Point deleteItem;

    // Checking if the heap is empty or not
    if (h->size == 0) {
        printf("\nHeap is empty.");
        return (Point){-1, -1, -1}; // Return a dummy Point in case of empty heap
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

// Define a insert function
void insert(heap* h, Point data)
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
    h->arr = (Point*)realloc(h->arr, newCapacity * sizeof(Point));
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

int main()
{
    heap* hp = createHeap(0);

    insert(hp, (Point){9, 0, 9.0});
    printHeap(hp);
    insert(hp, (Point){8, 0, 8.0});
    printHeap(hp);
    insert(hp, (Point){7, 0, 7.0});
    printHeap(hp);
    insert(hp, (Point){6, 0, 6.0});
    printHeap(hp);
    insert(hp, (Point){5, 0, 5.0});
    insert(hp, (Point){4, 0, 4.0});
    insert(hp, (Point){3, 0, 3.0});
    insert(hp, (Point){2, 0, 2.0});
    printHeap(hp);
    insert(hp, (Point){1, 0, 1.0});

    printHeap(hp);
    extractMin(hp);
    printHeap(hp);
    insert(hp, (Point){10, 0, 6.5});
    printHeap(hp);

    return 0;
}




#include <stdio.h>
#include <stdlib.h>

// Define the Point struct
typedef struct {
    int x;
    int y;
    float time;
} Point;

// Declare a heap structure
struct Heap {
    Point* arr;
    int size;
    int capacity;
};

// define the struct Heap name
typedef struct Heap heap;

// forward declarations
heap* createHeap(int capacity);
void insertHelper(heap* h, int index);
void heapify(heap* h, int index);
Point extractMin(heap* h);
void insert(heap* h, Point data);
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
        h->arr = (Point*)malloc(capacity * sizeof(Point));
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
        Point temp = h->arr[parent];
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
        Point temp = h->arr[min];
        h->arr[min] = h->arr[index];
        h->arr[index] = temp;

        // recursively calling for their child elements
        // to maintain min heap
        heapify(h, min);
    }
}

Point extractMin(heap* h)
{
    Point deleteItem;

    // Checking if the heap is empty or not
    if (h->size == 0) {
        printf("\nHeap is empty.");
        return (Point){-1, -1, -1.0}; // Return a dummy Point in case of empty heap
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

// Define an insert function
void insert(heap* h, Point data)
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
    h->arr = (Point*)realloc(h->arr, newCapacity * sizeof(Point));
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

int main()
{
    heap* hp = createHeap(0);

    insert(hp, (Point){9, 0, 9.0});
    insert(hp, (Point){8, 0, 8.0});
    insert(hp, (Point){7, 0, 7.0});
    insert(hp, (Point){6, 0, 6.0});
    insert(hp, (Point){5, 0, 5.0});
    insert(hp, (Point){4, 0, 4.0});
    insert(hp, (Point){3, 0, 3.0});
    insert(hp, (Point){2, 0, 2.0});
    insert(hp, (Point){1, 0, 1.0});

    printHeap(hp);

    Point minPoint = extractMin(hp);
    printf("Extracted: (%d, %d, %.2f)\n", minPoint.x, minPoint.y, minPoint.time);

    printHeap(hp);

    insert(hp, (Point){10, 0, 6.5});
    printHeap(hp);

    return 0;
}

*/


#include <stdio.h>
#include <stdlib.h>

// Define the Point struct
typedef struct {
    int x;
    int y;
    float time;
} Point;

// Declare a heap structure
struct Heap {
    Point* arr;
    int size;
    int capacity;
};

// define the struct Heap name
typedef struct Heap heap;

// forward declarations
heap* createHeap(int capacity);
void insertHelper(heap* h, int index);
void heapify(heap* h, int index);
Point extractMin(heap* h);
Point getMin(heap* h);
void insert(heap* h, Point data);
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
        h->arr = (Point*)malloc(capacity * sizeof(Point));
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
        Point temp = h->arr[parent];
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
        Point temp = h->arr[min];
        h->arr[min] = h->arr[index];
        h->arr[index] = temp;

        // recursively calling for their child elements
        // to maintain min heap
        heapify(h, min);
    }
}

Point extractMin(heap* h)
{
    Point deleteItem;

    // Checking if the heap is empty or not
    if (h->size == 0) {
        printf("\nHeap is empty.");
        return (Point){-1, -1, -1}; // Return a dummy Point in case of empty heap
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

// Define the getMin function
Point getMin(heap* h)
{
    // Checking if the heap is empty or not
    if (h->size == 0) {
        printf("\nHeap is empty.");
        return (Point){-1, -1, -1.0}; // Return a dummy Point in case of empty heap
    }

    // Return the root element, which is the minimum
    return h->arr[0];
}

// Define an insert function
void insert(heap* h, Point data)
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
    h->arr = (Point*)realloc(h->arr, newCapacity * sizeof(Point));
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

int main()
{
    heap* hp = createHeap(0);

    insert(hp, (Point){9, 0, 9.0});

    printHeap(hp);

    if(hp->size != 0){
        printf("All correct\n");
    }
    else{
        printf("Wrong\n");
    }

    extractMin(hp);

    if(hp->size != 0){
        printf("All correct\n");
    }
    else{
        printf("Wrong\n");
    }

    return 0;
}
