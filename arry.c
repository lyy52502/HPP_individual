#include <stdio.h>
#include <stdlib.h>
// learn from lab9


static inline void swap(int* restrict a, int* restrict b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Loop unrooling, learn from lab7_task1
void sort_array(int* restrict arr, int size) {
    if (size < 2) return;

    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(&arr[j], &arr[j + 1]);
            }
        }
    }
}

void print_array(int arr[], int N) {
    for (int i = 0; i < N; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

int main(int argc, char* argv[]) {
	int size = atoi(argv[1]);
	int* arr = (int*) malloc(size * sizeof(int));
    	for (int i = 0; i < size; i++) {
        	arr[i] = rand() % 1000;
    	}

	printf("Unsorted array:\n");
	print_array(arr, size);
	sort_array(arr,size);
	printf("Sorted array:\n");
    	print_array(arr, size);
	for (int i = 0; i < size - 1; i++) {
        	if (arr[i] > arr[i + 1]) {
            	printf("Error! array not sorted!\n");
            	free(arr);
            	return -1;
        	}
    	}
	printf("OK, array is sorted!\n");

	free(arr);
	return 0;
}
