#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <string.h>

/* Reference: QuicSort.pdf, lab8, lab9, lab11, Assignment4, https://iq.opengenus.org/parallel-quicksort/, chatgpt, deepseek */
typedef struct {
    int* data;
    int size;
} ThreadData;

// learn from lab10_task5
static inline double get_wall_seconds() {
    return omp_get_wtime();
}

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

// Reference: QuickSort.pdf
// Here I use strategy 3:Sort the medians and select the mean value of the two middlemost medians in each processor set and step.
int select_pivot(ThreadData* threadData,int num_threads, int num_groups) {

    /*thread ID, 1, 2, 3, 4;
    the array size is 20
    1. every thread has 5 data, here we have a local median as pivot and split into two group;
    2. group 0, group 1;
    3. every group has a median, so we need to split the data into group 2, group 3, group4 group 5  */
    int myid = omp_get_thread_num();// the thread ID
    int group_size = num_threads / num_groups; // the number of threads within a group


    // Local ID within the group and group ID
    int locid = myid % group_size; //every group has 2 threads now, and we need to get the id of thread within group, 1%2=1, 2%2=0
    int group = myid / group_size;

    // 1. Sort the data in each thread (local array)
    int* data = threadData[myid].data; // the data within one thread
    int data_size = threadData[myid].size;// the number of data within one thread
    sort_array(data, data_size);

    // 2. Compute the median of each thread
    int median;
    if (data_size % 2 == 0) {
        median = (data[data_size / 2 - 1] + data[data_size / 2]) / 2;
    } else {
        median = data[data_size / 2];
    }


    // Ensures only one thread executes this block, the use of single learn  from chatgpt
    int** medians = (int**) malloc(num_groups * sizeof(int*));
    int* pivots = (int*) malloc(num_groups * sizeof(int));

    for (int g = 0; g < num_groups; g++) {
        medians[g] = (int*) malloc(group_size * sizeof(int));
    }

    #pragma omp critical
    {
        medians[group][locid] = median;
    }

    #pragma omp barrier 

    // 3. Sorting the medians array and get the median of each group as pivot
    // Ensuring only one thread per group does sorting, avoding race conditions, I learn this method from chatgpt
    if (locid == 0) {
        sort_array(medians[group], group_size);
        if (group_size % 2 == 0) {
            pivots[group] = (medians[group][group_size / 2 - 1] + medians[group][group_size / 2]) / 2;
        } else {
            pivots[group] = medians[group][group_size / 2];
        }
    }
    #pragma omp barrier
   
    #pragma omp single
    {
        for (int g = 0; g < num_groups; g++) {
            free(medians[g]);
        }
        free(medians);
    }

    return pivots[group];
}

// Here I split the data into different group based on the comparison with pivot, I made modification based on the advice of chatgpt.
int findsplit(int* restrict data, int data_size, int pivot) {
    int splitpoint = 0;

    // Phase 1: Count elements <= pivot
    #pragma omp parallel for reduction(+:splitpoint)
    for (int i = 0; i < data_size; i++) {
        if (data[i] <= pivot) splitpoint++;
        printf(" the data before spliting:%d\n", data[i]);
    }

    // Phase 2: Partition
    int* temp = (int*) malloc(data_size * sizeof(int));
    int left = 0, right = splitpoint;

    #pragma omp parallel
    {
        int local_left = 0, local_right = 0;
        int* local_temp = (int*) malloc(data_size * sizeof(int));

        #pragma omp for nowait
        for (int i = 0; i < data_size; i++) {
            if (data[i] <= pivot) 
                local_temp[local_left++] = data[i];
            else 
                local_temp[splitpoint + local_right++] = data[i];
        }

        // Copy local partition back to temp
        #pragma omp critical
        {
            memcpy(temp + left, local_temp, local_left * sizeof(int));
            left += local_left;
            memcpy(temp + right, local_temp + splitpoint, local_right * sizeof(int));
            right += local_right;
        }
        for (int i = 0; i < data_size; i++) {
            printf(" the data after spliting:%d\n", data[i]);
        }

        free(local_temp);
    }

    // Phase 3: Copy back
    memcpy(data, temp, data_size * sizeof(int));
    free(temp);
    
    return splitpoint;
}
// Here is a basic merge function, I learn from https://www.geeksforgeeks.org/merge-sort/
void merge( int* restrict merged_data,  int* restrict data1, int size1, int* restrict data2, int size2, int lowerparts) {
    int* temp = (int*) malloc((size1 + size2) * sizeof(int));
    int i = 0, j = 0, k = 0;

 
    printf("Merging: size1=%d, size2=%d\n", size1, size2);
    printf("data1: ");
    for (int i = 0; i < size1; i++) printf("%d ", data1[i]);
    printf("\ndata2: ");
    for (int i = 0; i < size2; i++) printf("%d ", data2[i]);
    printf("\n");
            

    // Merge based on lowerparts flag
    if (lowerparts) {
        // Merge in ascending order
        while (i < size1 && j < size2) {
            if (data1[i] <= data2[j]) {
                temp[k++] = data1[i++];
            } else {
                temp[k++] = data2[j++];
            }
        }
    } else {
        // Merge in descending order
        while (i < size1 && j < size2) {
            if (data1[i] >= data2[j]) {
                temp[k++] = data1[i++];
            } else {
                temp[k++] = data2[j++];
            }
        }
    }
    


    // Copy remaining elements from data1
    while (i < size1) {
        temp[k++] = data1[i++];
    }

    // Copy remaining elements from data2
    while (j < size2) {
        temp[k++] = data2[j++];
    }

    // Copy merged data back to original arrays
    #pragma omp parallel for
    for (i = 0; i < size1 + size2; i++) {
        merged_data[i] = temp[i]; 
    }

    free(temp);

}


void global_sort(ThreadData* threadData, int myid, int num_threads, int num_groups) {
    int* data = threadData[myid].data;
    int data_size = threadData[myid].size;
    if (num_groups <= 1) {
        sort_array(data, data_size);
        return;
    }
    int group_size = num_threads / num_groups;
    if (group_size < 2) group_size = 1;
    int locid = myid % group_size;
    int group = myid / group_size;
    int pivot=select_pivot(threadData, num_threads,num_groups);
    findsplit(data, data_size, pivot);
    printf("group_size:%d\n",group_size);
    printf("split finished.....\n");
    
    
    #pragma omp barrier
    
    // Merge the data based on the split
    if (group_size>=2){
        printf("............\n");
        int merge_part=-1;
        int new_size = threadData[myid].size + threadData[merge_part].size;
        int* merged_data = (int*)malloc(new_size * sizeof(int));
        if (locid < group_size / 2) {
            merge_part=myid+group_size/2;
            printf("Thread %d (group %d) merging with %d (upperpart) \n", myid, group, merge_part);
            merge(merged_data, data, data_size, threadData[myid + group_size / 2].data, threadData[myid + group_size / 2].size, 1);
        } else {
            merge_part=myid-group_size/2;
            printf("Thread %d (group %d) merging with %d (lowerpart) \n", myid, group, merge_part);
            merge(merged_data, data, data_size,  threadData[myid - group_size / 2].data, threadData[myid - group_size / 2].size, 0);
        }
       

    }
    
    #pragma omp barrier

    // Recursively sort
    if (num_groups > 1) {
        int next_num_groups = num_groups / 2;  // Reduce the number of groups
        if (next_num_groups > 0) {  // Ensure recursion is valid
            global_sort(threadData, myid, num_threads, next_num_groups);
        }
    }
}

void print_array(int arr[], int N) {
    for (int i = 0; i < N; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}
// count the times of a number occur, learn from lab10_task9
static int count_values(const int* arr, int n, int x) {
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (arr[i] == x)
            count++;
    }
    return count;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: %s <array_size> <num_threads>\n", argv[0]);
        return 1;
    }

    int size = atoi(argv[1]);
    int num_threads = atoi(argv[2]);
    omp_set_num_threads(num_threads);


    int* arr = (int*) malloc(size * sizeof(int));
    for (int i = 0; i < size; i++) {
        arr[i] = rand() % 1000;
    }

    int count7 = count_values(arr, size, 7);
    printf("Before sort: the number 7 occurs %d times in the list.\n", count7);

    printf("Unsorted array:\n");
    print_array(arr, size);

    ThreadData* threadData = (ThreadData*) malloc(num_threads * sizeof(ThreadData));
    // Dividing data into different threads
    int chunk_size = size / num_threads;
    int remainder = size % num_threads;
    int offset = 0;
    for (int i = 0; i < num_threads; i++) {
        threadData[i].size = (i < remainder) ? chunk_size + 1 : chunk_size;
        threadData[i].data = arr + offset; //pionter
        offset += threadData[i].size;
    }

    double startTime = get_wall_seconds();

    // int num_groups = 1;
    // while (num_groups * 2 <= num_threads) num_groups *= 2;
    // printf("num_threads: %d, num_groups: %d\n", num_threads, num_groups);

    #pragma omp parallel num_threads(num_threads)
    {
        int myid = omp_get_thread_num();
        global_sort(threadData, myid, num_threads, num_threads);
    }

    double endTime = get_wall_seconds();

    int count7_again = count_values(arr, size, 7);
    printf("After sort : the number 7 occurs %d times in the list.\n", count7_again);

    printf("Sorted array:\n");
    print_array(arr, size);

    // Verify that the list is sorted, learn from lab10_task9
    for (int i = 0; i < size - 1; i++) {
        if (arr[i] > arr[i + 1]) {
            printf("Error! array not sorted!\n");
            free(arr);      
            return -1;
        }
    }
    printf("OK, array is sorted!\n");

    printf("Execution Time: %f seconds\n", endTime - startTime);

    free(threadData);
    free(arr);

    return 0;
}