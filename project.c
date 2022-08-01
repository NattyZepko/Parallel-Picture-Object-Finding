// Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "mpi.h"

typedef struct Picture
{
    int ID;
    int size;
    int **matrix;
} Picture;

/*
# return 1 when found
# return 0 if nothing found
# return -1 if error
*/
int findSubMatrix(Picture *picture, Picture *object, double matchingValue, int *x, int *y)
{
    if (object->size > picture->size)
        return -1;

    int repreats = (picture->size) - (object->size) + 1;
    double val = 0;
    for (int i = 0; i < repreats; i++)
    {
        for (int j = 0; j < repreats; j++)
        {
            val = 0;
            /*See if object fits*/

            for (int k = i; k < i + object->size; k++)
            {
                for (int l = j; l < j + object->size; l++)
                {
                    val += (abs(picture->matrix[k][l] - object->matrix[k - i][l - j])) / (double)(picture->matrix[k][l]);
                }
                if (val > matchingValue) // if we already surpassed the thereshold, we can stop trying this position, and go to the next position.
                    break;
            }
            if (matchingValue > val) // found
            {
                *x = i;
                *y = j;
                return 1;
            }
        }
    }
    *x = -1;
    *y = -1;
    return 0;
}

void readFromFile(double *matchingValue, int *numberOfPictures,
                  int *numberOfObjects, Picture ***objects, Picture ***pictures,
                  char *fileName)
{
    FILE *fp;
    if ((fp = fopen(fileName, "r")) == 0)
    {
        printf("cannot open file %s for reading\n", fileName);
        exit(0);
    }
    // matching value
    fscanf(fp, "%lf", matchingValue);
    // printf("%lf\n", *matchingValue); // print the matching value
    //  number of pictures
    fscanf(fp, "%d", numberOfPictures);
    // printf("%d\n", *numberOfPictures); // print the number of pictures
    pictures[0] = (Picture **)malloc(sizeof(Picture *) * (*numberOfPictures)); // assign space for pointers to all pictures
    if (!pictures[0])
    {
        printf("error in malloc readFromFile all pictures\n");
        exit(0);
    }
    for (int i = 0; i < (*numberOfPictures); i++) // Assign pictures in a loop
    {
        pictures[0][i] = (Picture *)malloc(sizeof(Picture)); // pointer to the i-th picture
        if (!pictures[0][i])
        {
            printf("error in malloc readFromFile all pictures[i]\n");
            exit(0);
        }
        fscanf(fp, "%d", &pictures[0][i]->ID);                                         // Read ID
        fscanf(fp, "%d", &pictures[0][i]->size);                                       // Read Size
        pictures[0][i]->matrix = (int **)malloc(sizeof(int *) * pictures[0][i]->size); // init matrix
        if (!pictures[0][i]->matrix)
        {
            printf("error in malloc readFromFile picture[i]->matrix\n");
            exit(0);
        }
        for (int k = 0; k < pictures[0][i]->size; k++)
        {
            pictures[0][i]->matrix[k] = (int *)malloc(
                sizeof(int) * pictures[0][i]->size);
            if (!pictures[0][i]->matrix[k])
            {
                printf("error in malloc readFromFile picture[i][k]\n");
                exit(0);
            }
            for (int t = 0; t < pictures[0][i]->size; t++)
            {
                fscanf(fp, "%d", &pictures[0][i]->matrix[k][t]);
            }
        }
    }
    // read objects
    fscanf(fp, "%d", numberOfObjects); // number of objects
    objects[0] = (Picture **)malloc(sizeof(Picture *) * (*numberOfObjects));
    for (int i = 0; i < (*numberOfObjects); i++)
    {
        objects[0][i] = (Picture *)malloc(sizeof(Picture));
        if (!objects[0][i])
        {
            printf("error in malloc readFromFile all objects[i]\n");
            exit(0);
        }
        fscanf(fp, "%d", &(objects[0][i]->ID));                                        // read ID
        fscanf(fp, "%d", &(objects[0][i]->size));                                      // read Size
        objects[0][i]->matrix = (int **)malloc(sizeof(int *) * (objects[0][i]->size)); // init matrix
        if (!objects[0][i]->matrix)
        {
            printf("error in malloc readFromFile obj[i]->matrix\n");
            exit(0);
        }
        for (int k = 0; k < (objects[0][i]->size); k++)
        {
            objects[0][i]->matrix[k] = (int *)malloc(
                sizeof(int) * (objects[0][i]->size));
            if (!objects[0][i]->matrix[k])
            {
                printf("error in malloc readFromFile obj[i][k]\n");
                exit(0);
            }
            for (int t = 0; t < (objects[0][i]->size); t++)
            {
                fscanf(fp, "%d", &objects[0][i]->matrix[k][t]);
            }
        }
    }
    fclose(fp);
}

int main(int argc, char *argv[])
{
    // General MPI variables
    int rank, processesCount;
    MPI_Status status;

    // All processes' variables
    double tStart, tEnd;
    double tRead, tSending, tWork;
    double matchingValue = 0;
    int numberOfPictures = 0;
    int numberOfObjects = 0;
    char fileName[] = "input.txt";

    Picture **pictures;
    Picture **objects;

    // MPI Init
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processesCount);

    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(8); // Use exactly 8 threads for all consecutive parallel regions. 8 turns out to be the optimal amount by trial, however, both lines are dismisable, and the program runs close to optimally without it as well.

    if (rank == 0)
    {
        // Validations
        if (processesCount != 2)
        {
            printf("Configure the program to run with exactly 2 processes!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(0);
        }
        else
            printf("--------------PROGRAM START--------------\n");

        /* STEP 1: READ FROM FILE */

        printf("Step 1: Starting to read from file.\n");
        tStart = MPI_Wtime();
        readFromFile(&matchingValue, &numberOfPictures, &numberOfObjects, &objects, &pictures, fileName);
        tEnd = MPI_Wtime();
        tRead = tEnd - tStart;
        printf("Finished reading. Total reading time: %f\n", tRead);

        /* STEP 2: SEND INFO TO THE OTHER PC */

        printf("Step 2: Sending info to proccess 1.\n");
        tStart = MPI_Wtime();
        MPI_Send(&matchingValue, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD); // send matching value to #1
        MPI_Send(&numberOfPictures, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&numberOfObjects, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        // Send half the pictures. (rounded down in case of odd number)
        for (int i = 0; i < numberOfPictures / 2; i++)
        {
            MPI_Send(&(pictures[i]->ID), 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            MPI_Send(&(pictures[i]->size), 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            for (int j = 0; j < pictures[i]->size; j++)
                MPI_Send(pictures[i]->matrix[j], pictures[i]->size, MPI_INT, 1, 0, MPI_COMM_WORLD);
        }
        // Send all objects.
        for (int i = 0; i < numberOfObjects; i++)
        {
            MPI_Send(&(objects[i]->ID), 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            MPI_Send(&(objects[i]->size), 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            for (int j = 0; j < objects[i]->size; j++)
                MPI_Send((objects[i]->matrix[j]), objects[i]->size, MPI_INT, 1, 0, MPI_COMM_WORLD);
        }
        tEnd = MPI_Wtime();
        tSending = tEnd - tStart;
        printf("Finished Sending. Total sending time: %f\n", tSending);
    }
    else // rank 1:
    {
        // recieve info from other pc goes here
        MPI_Recv(&matchingValue, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&numberOfPictures, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&numberOfObjects, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        pictures = (Picture **)malloc(sizeof(Picture *) * (numberOfPictures / 2));
        for (int i = 0; i < numberOfPictures / 2; i++)
        {
            pictures[i] = (Picture *)malloc(sizeof(Picture));
        }
        objects = (Picture **)malloc(sizeof(Picture *) * (numberOfObjects));
        for (int i = 0; i < numberOfObjects; i++)
        {
            objects[i] = (Picture *)malloc(sizeof(Picture));
        }
        // recieve half the pictures (again, rounded down)
        for (int i = 0; i < numberOfPictures / 2; i++)
        {
            MPI_Recv(&(pictures[i]->ID), 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&(pictures[i]->size), 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            pictures[i]->matrix = (int **)malloc(sizeof(int *) * (pictures[i]->size));
            for (int j = 0; j < pictures[i]->size; j++)
            {
                pictures[i]->matrix[j] = (int *)malloc(sizeof(int *) * (pictures[i]->size));
                MPI_Recv(pictures[i]->matrix[j], pictures[i]->size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            }
        }
        // recieve all the objects
        for (int i = 0; i < numberOfObjects; i++)
        {
            MPI_Recv(&(objects[i]->ID), 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&(objects[i]->size), 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            objects[i]->matrix = (int **)malloc(sizeof(int *) * ((objects[i]->size)));
            for (int j = 0; j < objects[i]->size; j++)
            {
                objects[i]->matrix[j] = (int *)malloc(sizeof(int *) * (objects[i]->size));
                MPI_Recv(objects[i]->matrix[j], objects[i]->size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            }
        }
    }

    /* STEP 3: WORK ON FINDING OBJECTS IN PICTURES */

    if (rank == 0)
    {
        printf("Step 3: Find Objects in Pictures.\n");
        tStart = MPI_Wtime();
    }

    // open a file for output.
    char outputFile[] = "output.txt";
    FILE *fp;
    if ((fp = fopen(outputFile, "w")) == 0)
    {
        printf("cannot open file %s for reading\n", outputFile);
        exit(0);
    }

    // Set boundry parameters, and reception parameters.
    int firstPicture = 0;
    int lastPicture = numberOfPictures;
    int result, x, y, tempPicID, tempObjID;
    if (rank == 0)
        firstPicture = numberOfPictures / 2; // Rank 0 starts at the middle
    else
        lastPicture = numberOfPictures / 2; // Rank 1 end at the middle (EXCLUDING!!)

    // loop for each picture
    #pragma omp parallel for private(result, x, y, tempPicID, tempObjID)
    for (int i = firstPicture; i < lastPicture; i++)
    {
        result = 0;                               // set as not found
        for (int j = 0; j < numberOfObjects; j++) // For each picture, match each object
        {
            result = findSubMatrix(pictures[i], objects[j], matchingValue, &x, &y); // 1 = found, other value otherwise
            // printf("I am proccess %d, and I finished working on picture %d, with object %d\n", rank, pictures[i]->ID, objects[j]->ID); //Enable line to show progress
            if (result) // FOUND!
            {
                if (rank == 0) // Just print
                    fprintf(fp, "Picture %d found Object %d in Position (%d , %d)\n", pictures[i]->ID, objects[j]->ID, x, y);
                else // I send it to 0 to print
                {
                    MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(&(pictures[i]->ID), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(&(objects[j]->ID), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(&x, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(&y, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
                // After I found the object, I am done with this picture, therefore:
                j = numberOfObjects; // exit the loop of object, go to next picture
            }
        }
        // after the loop of all objects didn't find, print that none were found:
        if (!result)
        {
            if (rank == 0) // proccess 0 just prints.
                fprintf(fp, "Picture %d no objects were found\n", pictures[i]->ID);
            else // proccess 1 signals a 0 for not found.
            {
                MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                MPI_Send(&(pictures[i]->ID), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
    }
    // printf("Proccess %d finished working.\n", rank); enable line for progress following
    // Proccess 0 will recieve data from proccess 1 about their picture, and print the result.
    if (rank == 0)
    {
        for (int i = 0; i < numberOfPictures / 2; i++)
        {
            MPI_Recv(&result, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status); // get result, if 0 - none were found, if 1 - an object was found.
            if (result)
            {
                MPI_Recv(&tempPicID, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);                               // get the picture ID
                MPI_Recv(&tempObjID, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);                               // get the object ID
                MPI_Recv(&x, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);                                       // x coordinate
                MPI_Recv(&y, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);                                       // y coordinate
                fprintf(fp, "Picture %d found Object %d in Position (%d , %d)\n", tempPicID, tempObjID, x, y); // print where it was found
            }
            else
            {
                MPI_Recv(&tempPicID, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status); // get the picture ID
                fprintf(fp, "Picture %d no objects were found\n", tempPicID);    // pring that none were found
            }
        }
    }

    /* STEP 4: End the program */
    if (rank == 0)
    {
        tEnd = MPI_Wtime();
        tWork = tEnd - tStart;
        printf("Total calculation time: %f\n", tWork);
        printf("Total time: %f\n", (tRead + tSending + tWork));
    }

    MPI_Finalize();
    return 0;
}
