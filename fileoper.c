#include <stdlib.h>
#include <stdio.h>

/*************************************************************************/

float **ReadInputFile(filename,dimension,numpoints)

/*************************************************************************/

char *filename;
int *dimension,*numpoints;

{
  int j,k;
  FILE *f;
  float **A;
  float *center;

  f=fopen(filename,"r"); 

  if (f==NULL) {
    *dimension = -1;
    printf("file not found\n");
    return(NULL);
  }
  
  fscanf(f,"%d\n",numpoints);  /* read number of points */
  fscanf(f,"%d\n",dimension);  /* read dimension */

  center = (float *)malloc((*dimension)*sizeof(float));
  for (j=0; j < *dimension; j++) {
    center[j]=0.0;
  } /* used for center of gravity calculations */

  A = (float **)malloc((*numpoints)*sizeof(float *));

  for (k=0; k < *numpoints; k++) {
    A[k] = (float *)malloc((*dimension)*sizeof(float));
  }

  for (k=0;k<*numpoints;k++) {
    for (j=0; j < *dimension; j++) {
      fscanf(f,"%f\n",&A[k][j]);
      center[j]+=A[k][j];
    }
  }

  for (j=0; j < *dimension; j++) {
    center[j]=center[j]/(*numpoints);
  } 
  /* I have no idea why they were doing this -SCL */
  /* for (k=0;k<*numpoints;k++) { */
  /*   for (j=0; j < *dimension; j++) { */
  /*     A[k][j]=A[k][j]-center[j]; */
  /*   } */
  /* } /\* make center of gravity the origin *\/ */
  
  fclose(f);
  return(A);
}


