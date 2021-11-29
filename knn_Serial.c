//  Δημήτρης Παππάς 8391
//
//  Εργασία 2
//
//  Serial

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define FILE_NAME "binaryFile.bin"
#define K 5                             //  K nearest neighbors
#define MAX_NUM_OF_ELEMENTS 100
#define NUM_OF_DIMENSIONS 6
#define ROWS 60000             // file contains 60000 rows
#define COLS 30             // file contains 30 columns


struct timeval startwtime, endwtime;
double seq_time;


typedef struct {
    double distance;
    int element;
}knn_struct;



int Asc (const void * a, const void * b);
void knnFinder(double coordinates_1[COLS], double coordinates_2[COLS], knn_struct knn[ROWS][K]);
void Initialization(knn_struct knn[ROWS][K]);



// =============================================================================================================

int main(int argc, char *argv[]){


    //      --- Serial knn code ---

    double coordinates_1[COLS];         // i row of file
    double coordinates_2[COLS];         // j row of file
    knn_struct knn[ROWS][K];            // Each row has the K nearest elements in Ascended order(element) - File

    // Calculating the distances and finding the K nearest neighbors of each element of file
    gettimeofday (&startwtime, NULL);
    knnFinder(coordinates_1, coordinates_2, knn);
    gettimeofday (&endwtime, NULL);

    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

    printf("Serial clock time = %f\n", seq_time);


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

  //return ( knnA->distance - knnB->distance );
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


// ------------------------------------------------------------------------------------------------------------


void knnFinder(double coordinates_1[COLS], double coordinates_2[COLS], knn_struct knn[ROWS][K]){

    Initialization(knn);

    int i=0, j=0, c=0, k=0;
    double dis;                 // distance between elements

    FILE *fpi;
    FILE *fpj;
    double file_value_1;
    double file_value_2;

    // Open first file pointer (for the i element)
    fpi = fopen(FILE_NAME, "rb");

    if (!fpi){
    printf("Unable to open file! \nfpi error! \n");
    exit(1);
    }

    fseek(fpi, 0, SEEK_SET);         // setting fpi in the start of file


    // Open second file pointer (for the j element)
    fpj = fopen(FILE_NAME, "rb");

    if (!fpj){
        printf("Unable to open file! \nfpj error! \n");
        exit(1);
    }


    for (i=0 ; i<ROWS ; i++){       // First element to compare
        k = 0;

        // Store the coordinates of the i element
        for (c=0 ; c<COLS ; c++){
            fread(&file_value_1, sizeof(file_value_1), 1, fpi);     // coordinates of i row (i element)
            coordinates_1[c] = file_value_1;
            //printf("(i) %d element: coordinates[%d] = %lf\n", i, c, coordinates_1[c]);
        }

        fseek(fpj, 0, SEEK_SET);         // setting fpj in the start of file

        for (j=0 ; j<ROWS ; j++){       // Second element to compare

            dis = 0;

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

        /*printf("\n\nAscending order for %d element: \n", i);
        for (k=0 ; k<K ; k++){
            printf(" knn[%d][%d].distance=%lf (knn[%d][%d].element=%d) \n", i, k, knn[i][k].distance, i, k, knn[i][k].element);
        }
        printf("\n");
        */
    }

    fclose(fpj);
    fclose(fpi);

    // Results of knn 
    printf("knn program:\n\n");
    for (i=0 ; i<ROWS ; i++){
        printf("%d element's %d nearest neighbors are: ", i, K);
        for (k=0 ; k<K ; k++){
            printf(" %d + %lf   ", knn[i][k].element, knn[i][k].distance);
        }
        printf("\n\n");
    }

    printf("\n\nEnd\n\n==========================================================================\n\n");

}





