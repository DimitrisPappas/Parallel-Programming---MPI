#include "mpi.h"
#include "omp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define FILE_NAME "binaryFile.bin"
#define K 5                             //  K nearest neighbors
#define MAX_NUM_OF_ELEMENTS 100
#define NUM_OF_DIMENSIONS 6
#define ROWS 60000             // file contains 60000 rows
#define COLS 30             // file contains 30 columns
#define MASTER 0            // Master task
#define P 4                 // number of processes


typedef struct {
    double distance;
    int element;
}knn_struct;


int Asc (const void * a, const void * b);
void Initialization(knn_struct knn[ROWS][K]);

// =============================================================================================================

int main(int argc, char *argv[]){

    double coordinates_1[COLS];         // i row of file
    double coordinates_2[COLS];         // j row of file
    knn_struct knn[ROWS][K];            // Each row has the K nearest elements in Ascended order(element) - File

    
    Initialization(knn);

    int i=0, j=0, c=0, k=0;
    double dis;                 // distance between elements    
    FILE *fpi;                  // pointer for i element
    FILE *fpj;                  // pointer for j element
    double file_value_1;
    double file_value_2;
    int pnumber, pid, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Status myStatus; 
    double start, end, time;          // start & end of mpi communication time, time = overall time data passing
    double pstart,pend,ptime;
    MPI_Request send_request, recv_request;
    
    MPI_Init(&argc, &argv);                     // Initializes the MPI execution environment
    MPI_Comm_size(MPI_COMM_WORLD, &pnumber);    // Returns the total number of MPI processes
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);        // Returns the rank (id) of the calling MPI process 
    MPI_Barrier(MPI_COMM_WORLD);                // Wait until all processes arrive
    pstart = MPI_Wtime();                        // Time when mpi started
    
    // Struct Derived Data Type
    MPI_Datatype knn_type, oldtypes[2];     // required variables for the send/recv
    int blockcounts[2];                     // 2 different types included in the struct
    MPI_Aint offsets[2], extent;            // MPI_Aint type used to be consistent with syntax of MPI_Type_extent routine
    
    // Setup description of the 1 MPI_DOUBLE field distance
    offsets[0] = 0;
    oldtypes[0] = MPI_DOUBLE;       
    blockcounts[0] = 1;
    
    // Setup description of the 1 MPI_INT field element
    // Need to first figure offset by getting size of MPI_DOUBLE
    MPI_Type_extent(MPI_DOUBLE, &extent);
    offsets[1] = 1 * extent;
    oldtypes[1] = MPI_INT;       
    blockcounts[1] = 1;
    
    // Define structured type and commit it
    MPI_Type_struct(2, blockcounts, offsets, oldtypes, &knn_type);      // count, blocklens[], offsets[], old_type, &newtype
    MPI_Type_commit(&knn_type);                                         // Commits new datatype to the system
    
    printf ("MPI task %d has started...\n", pid);
    
    
    // Open first file pointer (for the i element)
    fpi = fopen(FILE_NAME, "rb");

    if (!fpi){
    printf("Unable to open file! \nfpi error! \n");
    exit(1);
    }

    fseek(fpi, (ROWS/P)*COLS*pid*sizeof(double), SEEK_SET);         // setting fpi according to the prosses's id


    // Open second file pointer (for the j element)
    fpj = fopen(FILE_NAME, "rb");

    if (!fpj){
        printf("Unable to open file! \nfpj error! \n");
        exit(1);
    }

    #pragma omp parallel for
    for (i=0 ; i<ROWS ; i++){       // First element to compare
        k = 0;

        // Store the coordinates of the i element
        for (c=0 ; c<COLS ; c++){
            fread(&file_value_1, sizeof(file_value_1), 1, fpi);     // coordinates of i row (i element)
            coordinates_1[c] = file_value_1;
            //printf("(i) %d element: coordinates[%d] = %lf\n", i, c, coordinates_1[c]);
        }

        fseek(fpj, 0, SEEK_SET);         // setting fpj in the start of file

        for (j=0 ; j<ROWS ; j++){          // Second element to compare

            dis = 0;                    // distance

            // Store the coordinates of j element and calculate distance between i, j (dis = |i-j|^2)
            for (c=0 ; c<COLS ; c++){                      // c = number of dimensions

                fread(&file_value_2, sizeof(file_value_2), 1, fpj);
                coordinates_2[c] = file_value_2;
                //printf("(i=%d) %d element: coordinates[%d] = %lf\n", i, j, c, coordinates_2[c]);

                // no need for square root (makes my code slower)
                dis = dis + (coordinates_1[c] - coordinates_2[c])*(coordinates_1[c] - coordinates_2[c]);   
            }
            
            if (dis == 0){      // don't want to calculate the distance of each element with itself (dis = 0)
                continue;
            }

            if (k<=K-1){
                knn[i][k].distance = dis;
                knn[i][k].element = j;
                qsort(&knn[i][0], K, sizeof(knn_struct), Asc);
                //printf("i=%d, j=%d, k=%d, knn[%d][%d].distance = %lf\n", i, j,k,i,k, knn[i][k].distance);
                k++;          // column counter of knn
            }else{
                if (dis < knn[i][K-1].distance){
                    knn[i][K-1].distance = dis;
                    knn[i][K-1].element = j;
                    qsort(&knn[i][0], K, sizeof(knn_struct), Asc);
                }
            }
        }
        start = MPI_Wtime();
        if (pid != MASTER){
            MPI_Irecv(&knn[i][0], K, knn_type, pid-1, 0, MPI_COMM_WORLD, &recv_request);                // pid receives i row from pid - 1
            MPI_Wait(&recv_request, &myStatus);
            printf("Process %d received from process %d\n", pid, pid - 1);
        }
                
        MPI_Isend(&knn[i][0], K, knn_type, (pid + 1) % pnumber, 0, MPI_COMM_WORLD, &send_request);      // pid sends to pid + 1 (last process sends to MASTER=0)
        MPI_Wait(&send_request, &myStatus);
                
        // Now MASTER process can receive from the last process
        if (pid == MASTER){
            MPI_Irecv(&knn[i][0], K, knn_type, pnumber-1, 0, MPI_COMM_WORLD, &recv_request);            // MASTER receives from last process 
            MPI_Wait(&recv_request, &myStatus);
            printf(" MASTER Process %d received from process %d\n", pid, pnumber - 1);
        }
        end = MPI_Wtime(); 
        time += end - start;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    pend = MPI_Wtime();              // Time when mpi ended
    ptime = pend - pstart;
    
    fclose(fpj);
    fclose(fpi);
    
    MPI_Type_free(&knn_type);       // Free datatype when done using it
    MPI_Finalize();                 // Done with MPI  

    
    // Results of knn 
    if (pid == MASTER){
        printf("knn array:\n\n");
        for (i=0 ; i<ROWS ; i++){
            printf("%d element's %d nearest neighbors are: ", i, K);
            for (k=0 ; k<K ; k++){
                printf(" %d + %lf   ", knn[i][k].element, knn[i][k].distance);
            }
            printf("\n\n");
        }

        printf("Clock time of MASTER's data passing (pid=0) = %lf\n", time);
        printf("Clock time of program = %lf", ptime);
        printf("\n\nEnd\n\n==========================================================================\n\n");
    }
    
     return 0;
}

// =============================================================================================================


// Sorting in Ascended order a struct
int Asc(const void * a, const void * b)
{

  knn_struct *knnA = (knn_struct *)a;
  knn_struct *knnB = (knn_struct *)b;

  if (knnA->distance < knnB->distance){
    return -1;
  }else{
    return 1;
  }
}


// ------------------------------------------------------------------------------------------------------------


void Initialization(knn_struct knn[ROWS][K]){
    int i=0, k=0;

    // Setting knn.distance = infinite and knn.element = -1
    for (i=0 ; i<ROWS ; i++){
        for (k=0 ; k<K ; k++){
            knn[i][k].distance = INFINITY;
            knn[i][k].element = -1;
        }
    }
}





