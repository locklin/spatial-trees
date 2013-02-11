/**************************************************************************/
/* Program to perform orthogonal range searches and nearest neighbor      */
/* querys in a more sophisticated k-d tree.  In this implementation the,  */
/* nodes on any given level of the tree do not have the same              */
/* discriminating dimension as the discrimiator is chosen based on the    */
/* dimension with   the "maxspead."                                       */
/*                                                                        */
/* References:  J.H. Friedman, J.L. Bentley, R.A. Finkel  "An Algorithm   */
/* for Finding Best Matches in Logarithmic Expected Time."                */
/* ACM Transactions on Mathematical Software, Vol 3 No. 3 Sept. 1977      */
/* pp. 209-226.                                                           */
/**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "optkd.h"
#include "nn.h"

/* Used to create a new tree in the k-d tree */
#define TESTTREE(PP)  ((PP) = (optkdNode *)malloc(sizeof(optkdNode)))
#define NEWTREE(PP)  if (TESTTREE(PP)==NULL) \
                         {printf("memory error\n");return;}

/*from pqueue.c */
extern void PQInsert();
extern void PQreplace();

/* from metric.c */
extern float ManhattDist();
extern float LInfinityDist();
extern float LGeneralDist();

static int *perm;  /* permutation array */
int *optfound;
static int Count=0;
float *nndist;

float (*Distance)();
extern double fabs();

extern float pradius;


/***************************************************************************/

void Selection(a, l,N, k,discrim)

/***************************************************************************/
/* Makes the perm partition the array Values along the element k.          */
/* Adapted from Sedgewick's Algorithms in C (p. 128)                       */
/***************************************************************************/

float **a;
int l,N,k,discrim;

{
float v;
int t,i,j,r;

r=N;

while(r>l) {
  v=a[perm[r]][discrim]; i=l-1; j=r;
  for (;;) {
    while (a[perm[++i]][discrim] < v);
    while (a[perm[--j]][discrim] > v && j>l); 
    if (i >= j) break;
    t=perm[i]; perm[i] = perm[j]; perm[j]=t;
  }
  t=perm[i]; perm[i] = perm[r]; perm[r]=t;
  if (i>=k) r=i-1;
  if (i<=k) l=i+1;
}
}

/****************************************************************************/

int findmaxspread(l,u,dimension,points)

/****************************************************************************/

int l,u,dimension;
float **points;

{
int i,j,maxdim;
float max       =-999999999.0,
      min       = 999999999.0,
      maxspread =-999999999.0;

for (i=0; i < dimension; i++) {
  max =-999999999.0;
  min = 999999999.0;
  for (j=l; j <= u; j++) {
    if (max < points[perm[j]][i]) { 
      max = points[perm[j]][i];
    }
    if (min > points[perm[j]][i]) { 
      min = points[perm[j]][i];
    }
    if (maxspread < fabs(max-min)) {
      maxspread = fabs(max-min);
      maxdim = i;
    }
  }
}
return(maxdim);
}

/*******************************************************************************/

optkdNode *BuildkdTree(points,l,u,dimension)

/*******************************************************************************/

int l,u;
float **points;

{
  optkdNode *p;
 int m;

/*  NEWTREE(p); again with the weird global thing -SCL */
 p=(optkdNode *)malloc(sizeof(optkdNode));
   if (u-l+1 <= BUCKETSIZE) {
     p->bucket = 1;
     p->lopt = l;
     p->hipt = u;
     p->loson = NULL;
     p->hison = NULL;
   } else {
     p->bucket =0;
     p->discrim = findmaxspread(l,u,dimension,points);
     m=(l+u)/2;
     Selection(points,l,u,m,p->discrim);
     p->cutval = points[perm[m]][p->discrim];
     p->loson = BuildkdTree(points,l,m,dimension);
     p->hison = BuildkdTree(points,m+1,u,dimension);
   }
 return(p);
}

/*******************************************************************************/

optkdNode *BuildOptTree(points,numPoints,dimension)

/*******************************************************************************/
int dimension,numPoints;
float **points;
{

int j;

  /* initialize perm array */
  perm = (int *) malloc(numPoints*sizeof(int));
  for (j=0; j < numPoints; j++) {
    perm[j]=j;
  }
  return(BuildkdTree(points,0,numPoints-1,dimension));
}


/*******************************************************************************/

void rnnEuclidean(p,querpoint,points,dimension,numpoints)

/*******************************************************************************/

/* special searching algorithm to take advantage of the fact that square roots
   do not need to be evaulated */

optkdNode *p;
float *querpoint;
float **points;
int dimension,numpoints;

{
  int i,j;
  float d,thisdist,val;

  if (p->bucket) {
    for (i=p->lopt; i <= p->hipt; i++) {
      thisdist=0.0;
      for (j=0; j<dimension; j++) {
        d=(querpoint[j]-points[perm[i]][j]);
        thisdist=thisdist+d*d;
      }        

      if (optfound[0] < numpoints) {
        PQInsert(thisdist,perm[i],nndist,optfound);
      } else {
        PQreplace(thisdist,nndist,optfound,perm[i]);
      }
    }
  } else {
    val = querpoint[p->discrim] - p->cutval;
    if (val < 0) {
      rnnEuclidean(p->loson,querpoint,points,dimension,numpoints);
      if (nndist[1] >= val*val) {
        rnnEuclidean(p->hison,querpoint,points,dimension,numpoints);
      }
    } else {
      rnnEuclidean(p->hison,querpoint,points,dimension,numpoints);
      if (nndist[1] >= val*val) {
        rnnEuclidean(p->loson,querpoint,points,dimension,numpoints);
      }
    }
  }
}

/*******************************************************************************/

void rnnGeneral(p,querpoint,points,dimension,numpoints,MinkP)

/*******************************************************************************/


optkdNode *p;
float *querpoint;
float **points;
int dimension,numpoints,MinkP;

{
  int i;
  float thisdist,val,thisx;

  if (p->bucket) {
    for (i=p->lopt; i <= p->hipt; i++) {
      thisdist=Distance(points,perm[i],querpoint,dimension,MinkP);

      if (optfound[0] < numpoints) {
        PQInsert(thisdist,perm[i],nndist,optfound);
      } else {
        PQreplace(thisdist,nndist,optfound,perm[i]);
      }
    }
  } else {
    val = p->cutval;
    thisx=querpoint[p->discrim];
    if (thisx < val) {
      rnnGeneral(p->loson,querpoint,points,dimension,numpoints,MinkP);
      if (thisx + nndist[1] > val) {
        rnnGeneral(p->hison,querpoint,points,dimension,numpoints,MinkP);
      }
    } else {
      rnnGeneral(p->hison,querpoint,points,dimension,numpoints,MinkP);
      if (thisx - nndist[1] < val) {
        rnnGeneral(p->loson,querpoint,points,dimension,numpoints,MinkP);
      }
    }
  }
}


/*******************************************************************************/

int *kdOptNNQuery(points,dimension, querpoint,numNN,Metric,root,MinkP)

/*******************************************************************************/

optkdNode *root;
float *querpoint, **points;
int dimension,numNN,MinkP;

{
  int j;

  /* set up found array */
  optfound = (int *) malloc((numNN+1)*sizeof(int));
  optfound[0]=1;  /* for now */

  /* nndist is a priority queue of the distances of the nearest neighbors found */
  nndist = (float *)malloc((numNN+1)*sizeof(float));
  for (j=0; j < numNN+1; j++) {
    nndist[j] = 99999999999.0;
  }

  switch(Metric) {
    case EUCLIDEAN : rnnEuclidean(root,querpoint,points,dimension,numNN);
                      break;
    case MANHATTAN : Distance=ManhattDist;
                      rnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
    case L_INFINITY: Distance=LInfinityDist;
                      rnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
    case L_P       : Distance=LGeneralDist;
                      rnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
  }
  free(nndist);
  return(optfound);
}

/*******************************************************************************/

int *kdOptNNQueryDist(points,dimension, querpoint,nndist,numNN,Metric,root,MinkP)

/*******************************************************************************/

optkdNode *root;
float *querpoint, *nndist, **points;
int dimension,numNN,MinkP;

{
  int j;

  /* set up found array */
  optfound = (int *) malloc((numNN+1)*sizeof(int));
  optfound[0]=1;  /* for now */

  /* nndist is a priority queue of the distances of the nearest neighbors found */

  for (j=0; j < numNN+1; j++) {
    nndist[j] = 99999999999.0;
  }

  switch(Metric) {
    case EUCLIDEAN : rnnEuclidean(root,querpoint,points,dimension,numNN);
                      break;
    case MANHATTAN : Distance=ManhattDist;
                      rnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
    case L_INFINITY: Distance=LInfinityDist;
                      rnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
    case L_P       : Distance=LGeneralDist;
                      rnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
  }
  return(optfound);
}


/***************************************************************************/

void KillOptTree(P)

/***************************************************************************/

/*  Kills a kd-tree to avoid memory holes.   */


optkdNode *P;

{
  if (perm != NULL) {
    free(perm);
  }  /* free permutation array */

  if (P==NULL) {
    return;
  } /* just to be sure */
  if (P->loson != NULL) {
    KillOptTree(P->loson);
  }

  if (P->hison != NULL) {
    KillOptTree(P->hison);
  }

  free(P);

}


/***************************************************************************/

void optInRegion(P,Dimension,Points,RectQuery)

/***************************************************************************/
/* Determines if the treenode P falls inside the rectangular query         */
/* RectQuery.  If so, adds the array index of the point to the found       */
/* array.                                                                  */
/***************************************************************************/

optkdNode *P;
int Dimension;
float **Points, **RectQuery;

{
  int index,dc,InsideRange;

  
  for (index=P->lopt;index<=P->hipt;index++) {
    InsideRange=1;
/*    circle(Points[perm[index]][0],Points[perm[index]][1],pradius); */

    for (dc=0; dc < Dimension; dc ++) {
      if ((Points[perm[index]][dc] < RectQuery[dc][0] || Points[perm[index]][dc] > 
            RectQuery[dc][1])) {   /* P is in the region */
         InsideRange=0;
         break;
      }
    }
    if (InsideRange) {
      Count++;
      if ((Count % 4000) == 0) {
         realloc(optfound,(4000+Count)*sizeof(int));
      }
      if (optfound == NULL) {
        printf("we have a memory problem\n");
      }
      optfound[Count] = perm[index];
    }
  }
}

/***************************************************************************/

void optAddRegion(P,Dimension,Points,RectQuery)

/***************************************************************************/
/* Adds the array index of each point in the bucket the point to the found */
/* array.  There is no need to check if the points are in it because we    */
/* have proven so already.                                                 */
/***************************************************************************/

optkdNode *P;
int Dimension;
float **Points, **RectQuery;

{
  int index;

  for (index=P->lopt;index<=P->hipt;index++) {
    Count++;
    if ((Count % 4000) == 0) {
       realloc(optfound,(4000+Count)*sizeof(int));
    }
    if (optfound == NULL) {
      printf("we have a memory problem\n");
    }
    optfound[Count] = perm[index];
  }
}

/***************************************************************************/

int optBoundsIntersectRegion(B,RectQuery,Dimension)

/***************************************************************************/
/* Returns true iff the hyper-rectangle defined by bounds array B          */
/* intersects the rectangular query RectQuery.                             */
/***************************************************************************/

float *B,**RectQuery;
int Dimension;

{
  int dc;

  for (dc=0; dc < Dimension; dc++) {
    if (B[2*dc] > RectQuery[dc][1] || B[2*dc+1] < RectQuery[dc][0]) {
      return(0);
    }
  }
  return(1);
}

/***************************************************************************/

int optBoundsContainsRegion(B,RectQuery,Dimension)

/***************************************************************************/
/* Returns true iff the hyper-rectangle defined by bounds array B          */
/* is completely contained inside the rectangular query RectQuery.         */
/***************************************************************************/

float *B,**RectQuery;
int Dimension;

{
  int dc;

  for (dc=0; dc < Dimension; dc++) {
    if (!(B[2*dc]   >= RectQuery[dc][0] &&
        B[2*dc+1] <= RectQuery[dc][1])) {
      return(0);
    }
  }
  return(1);
}


/***************************************************************************/

void optRangeSearch(P,Points,Dimension,RectQuery,B)

/***************************************************************************/

optkdNode *P;
float **RectQuery, **Points, *B;
int Dimension;


{
  int dc, disc;
  float *BHigh,*BLow;
  

  if (P==NULL) {printf("somehow a null pointer got sent here\n");}

  if (P->bucket) {   
    if (optBoundsContainsRegion(B,RectQuery,Dimension)) {
      optAddRegion(P,Dimension,Points,RectQuery);
    } else {
      optInRegion(P,Dimension,Points,RectQuery);
    }
    return;
  }

  /* Claim: P is not a bucket node */ 
  disc=P->discrim;
  BLow =  (float *)(malloc(2*Dimension*sizeof(float)));
  BHigh = (float *)(malloc(2*Dimension*sizeof(float)));

  if (BLow == NULL || BHigh == NULL) {
   printf("we have a memory error\n");
  }

  /* copy the region B into BLow, BHigh */
  for (dc=0; dc < 2*Dimension; dc++) {
    BLow[dc]  = B[dc];
    BHigh[dc] = B[dc];
  }

  /* Improve the Bounds for the subtrees */
  BLow[2*disc+1] = P->cutval;
  BHigh[2*disc] = P->cutval;

  if (optBoundsIntersectRegion(BLow,RectQuery,Dimension)) {
    optRangeSearch(P->loson,Points,Dimension,RectQuery,BLow);
  }
  free(BLow);
  if (optBoundsIntersectRegion(BHigh,RectQuery,Dimension)) {
    optRangeSearch(P->hison,Points,Dimension,RectQuery,BHigh);
  }
  free(BHigh);
}


/***************************************************************************/

int *kdOptRectQuery(root,Points,Dimension,RectQuery)

/***************************************************************************/

optkdNode *root;
float **RectQuery, **Points;
int Dimension;

{
float *B;
int dc;

  B =  (float *)(malloc(2*Dimension*sizeof(float)));
  if (B == NULL) {
    printf("We have a memory problem\n");
  }

  for (dc =0; dc < Dimension; dc++) {
    B[2*dc]   = RectQuery[dc][0];
    B[2*dc+1] = RectQuery[dc][1];
  }

  Count=0;
  optfound = (int *)(malloc(4000*sizeof(int)));
  if (optfound == NULL) {
    printf("We have a memory problem\n");
  }

  optRangeSearch(root,Points,Dimension,RectQuery,B);
  free(B);
  optfound[0] = Count;
  return(optfound);
}

