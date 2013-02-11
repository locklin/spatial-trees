#include <stdlib.h>
#include "nn.h"

/* from metric.c import */
extern float EuclidDist2();
extern float ManhattDist();
extern float LInfinityDist();
extern float LGeneralDist();

/***************************************************************************
*                                                                          *
* This module contains naive implementations of nearest neighbor           *
* and rectangular queries.                                                 *
*                                                                          *
****************************************************************************/


/***************************************************************************/

int *NaiveRectQuery(Points,NumPoints,Dimension,RectQuery)

/***************************************************************************
 *                                                                         *
 *  Returns a pointer to an array of array indices of points in the array  *
 *  Points falling within the hyper-rectangle defined by RectQuery.        *
 *  The first element of the array contains the number of points falling   *
 *  inside the query.                                                      * 
 *                                                                         *
 ***************************************************************************/

float **Points,**RectQuery;
int NumPoints,Dimension;

{
  int count,j,k, InsideRange,
      *found;            /* start at 50, and realloc if necessary */ 
 
  count = 0;

  found = (int *) malloc(NumPoints*sizeof(int)); 
  for (j=0; j < NumPoints; j++) {
    InsideRange = 1;
    for (k=0; (k < Dimension); k++) {
      if (Points[j][k] < RectQuery[k][0] || Points[j][k] > RectQuery[k][1]) {
        InsideRange = 0;
        break;
      }
    }
    if (InsideRange) {
      count++;
      found[count] = j;
    }
  }
  found[0] = count;
  return(found);
}



/***************************************************************************/

void selection(a, perm, N, k)

/***************************************************************************/

/* Makes the perm partition the array Values along the element k.          */
/* Adapted from Sedgewick's Algorithms in C (p. 128)                       */

float *a;
int *perm, N,k;

{
float v;
int t,i,j,l,r;

l=0;r=N-1;

while(r>l) {
  v=a[perm[r]]; i=l-1; j=r;
  for (;;) {
    while (a[perm[++i]] < v);
    while (a[perm[--j]] > v && j>l) ;
    if (i >= j) break;
    t=perm[i]; perm[i] = perm[j]; perm[j]=t;
  }
  t=perm[i]; perm[i] = perm[r]; perm[r]=t;
  if (i>=k) r=i-1;
  if (i<=k) l=i+1;
}
}

/***************************************************************************/

int *NaiveNNQuery(Points,NumPoints,Dimension,NNQPoint,NumNN,Metric,MinkP)

/***************************************************************************
 *                                                                         *
 *  Returns a pointer to an array of the NumNN indices of the array points *
 *  closest to the query point q.                                          *
 *  The first element of the array contains the number of points asked for,*
 *  NumNN.                                                                 *
 *                                                                         *
 ***************************************************************************/

float **Points,*NNQPoint;
int NumPoints,Dimension,NumNN,Metric,MinkP;

{
  int j,*perm,*found; 

  float *Dist;

  float (*Distance)();


  switch(Metric) {
    case EUCLIDEAN  : Distance = EuclidDist2;
                      break;
    case MANHATTAN  : Distance = ManhattDist;
                      break;
    case L_INFINITY : Distance = LInfinityDist;
                      break;
    case L_P        : Distance = LGeneralDist;
                      break;
    default         : Distance = EuclidDist2;
        
  }

  Dist = (float *) malloc((NumPoints)*sizeof(float)); 
  perm = (int *) malloc((NumPoints)*sizeof(int));   

  for (j=0; j < NumPoints; j++) {
    Dist[j]=Distance(Points,j,NNQPoint,Dimension,MinkP);
    perm[j] = j;
  }

  selection(Dist,perm,NumPoints,NumNN);

  free(Dist);

  found = (int *) malloc((NumNN+1)*sizeof(int)); 
  found[0] = NumNN;
  for (j=1; j <=NumNN; j++) {
    found[j] = perm[j-1];
  }
   free(perm);
  return(found);
}



