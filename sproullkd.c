/**************************************************************************/
/* Program to perform orthogonal range searches and nearest neighbor      */
/* querys in a Sproull k-d tree.  In this implementation,                 */
/* the nodes on any given level of the tree do not have orthogonal        */
/* discrimintating dimensions, rather there is an arbitary partition plane*/
/* chosen by the principal eigenvector of the covarience matrix.          */
/*                                                                        */
/* References:  R.F. Sproull "Refinements to Nearest-Neighbor Searching   */
/* in k-Dimensional Trees.  J. Algorithmica 1990.                         */
/* pp. 579-589.                                                           */
/**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "sproullkd.h"
#include "nn.h"

/* Used to create a new tree in the k-d tree */
#define TESTTREE(PP)  ((PP) = (sproullkdNode *)malloc(sizeof(sproullkdNode)))
#define NEWTREE(PP)  if (TESTTREE(PP)==NULL) \
                         {printf("memory error\n");return;}
#define EPS 1.0e-6

float *covar;  /* contains covarience matrix in lower triangular form */

/*from pqueue.c */
extern void PQInsert();
extern void PQreplace();

/* from metric.c */
extern float ManhattDist();
extern float LInfinityDist();
extern float LGeneralDist();

/* from eigens.c */
extern void eigens();

static int *perm;  /* permutation array */
int *sproullfound;
float *SPnndist;

float (*SPDistance)();

extern double fabs();

/***************************************************************************/

void SPSelection(points, l,N, k,plane,dimension)

/***************************************************************************/
/* Makes the perm partition the array Values along the plane, plane        */
/* Adapted from Sedgewick's Algorithms in C (p. 128)                       */
/***************************************************************************/

float **points,*plane;
int l,N,k;

{
float v,w,x;
int d,t,i,j,r;

r=N;

while(r>l) {
  v=0.0;
  for (d=0; d<dimension; d++) {
    v+=points[perm[r]][d]*plane[d];
  }
  i=l-1; j=r;
  for (;;) {
    do {
      i++; w=0.0;
      for (d=0; d<dimension; d++) {
        w+=points[perm[i]][d]*plane[d];
      }
    } while (w < v);

    do {
      j--; x=0.0;
      for (d=0; d<dimension; d++) {
        x+=points[perm[j]][d]*plane[d];
      }
    } while (x>v && j>l);

    if (i >= j) break;
    t=perm[i]; perm[i] = perm[j]; perm[j]=t;
  }
  t=perm[i]; perm[i] = perm[r]; perm[r]=t;
  if (i>=k) r=i-1;
  if (i<=k) l=i+1;
}
}

/****************************************************************************/

void ComputeCovarience(l,u,dimension,points)

/****************************************************************************/

/* compute the covariance matrix of the subsets of points defined by:       */
/* Points[perm[l]]..Points[perm[u]]                                         */

int l,u,dimension;
float **points;

{
int i,j,k,n,o;
float oneovern;

n=(u-l+1);
oneovern=1.0/(float)n;

for (i=0; i < dimension; i++) {
  for (j=0; j <= i; j++) {
    o=(i*i+i)/2 + j;
    covar[o]=0.0;
    for (k=l; k <= u; k++) {
      covar[o] += points[perm[k]][i]*points[perm[k]][j];
    }
    covar[o] = covar[o]*oneovern;
  }
}
}

/*******************************************************************************/

sproullkdNode *BuildSproullTree(points,l,u,dimension)

/*******************************************************************************/


int l,u;
float **points;

{
sproullkdNode *p;
int m,n, eigenindex;
float *eigenvects, *eigenvals, maxev;
/* NEWTREE(p); global weirds -SCL */
  p =(sproullkdNode *)malloc(sizeof(sproullkdNode));
  if (u-l+1 <= SPBUCKETSIZE) {
    p->bucket = 1;
    p->lopt = l;
    p->hipt = u;
    p->loson = NULL;
    p->hison = NULL;
  } else {
    p->bucket =0;
    ComputeCovarience(l,u,dimension,points);  /* compute covarience matrix */
    eigenvects = (float *)malloc(dimension*dimension*sizeof(float));
    eigenvals = (float *)malloc(dimension*sizeof(float));

    eigens(covar,eigenvects,eigenvals,dimension); /*compute eigenvectors using canned routine */
    maxev=-99999999.0;
    for (n=0; n < dimension; n++) {
      if (maxev < eigenvals[n]) {
        maxev=eigenvals[n];eigenindex=n;
      }
    }

    p->plane = (float *) malloc(dimension*sizeof(float));
    for (n=0; n < dimension; n++) {
      p->plane[n] = eigenvects[dimension*eigenindex + n];
    } /* data structure adjustment */

    free(eigenvects);free(eigenvals);

    m=(l+u)/2;
    SPSelection(points,l,u,m,p->plane,dimension);
    /* find median partition element */    

    p->cutval = 0.0;
    for (n=0; n<dimension; n++) {
     p->cutval += p->plane[n]*points[perm[m]][n];
    }

    p->loson = BuildSproullTree(points,l,m,dimension);
    p->hison = BuildSproullTree(points,m+1,u,dimension);
  }
  return(p);
}

/*******************************************************************************/

sproullkdNode *BuildSproullkdTree(points,numPoints,dimension)
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
  covar = (float *)malloc((dimension*(dimension+1)/2)*sizeof(float));

  return(BuildSproullTree(points,0,numPoints-1,dimension));
}


/*******************************************************************************/

void SPrnnEuclidean(p,querpoint,points,dimension,numpoints)

/*******************************************************************************/

/* special searching algorithm to take advantage of the fact that square roots
   do not need to be evaulated */

sproullkdNode *p;
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

      if (sproullfound[0] < numpoints) {
        PQInsert(thisdist,perm[i],SPnndist,sproullfound);
      } else {
        PQreplace(thisdist,SPnndist,sproullfound,perm[i]);
      }
    }
  } else {
    val = -p->cutval;
    for (i=0; i < dimension; i++) {
      val+=p->plane[i]*querpoint[i];
    }
    if (val < 0) {
      SPrnnEuclidean(p->loson,querpoint,points,dimension,numpoints);
      if (SPnndist[1] >= val*val) {
        SPrnnEuclidean(p->hison,querpoint,points,dimension,numpoints);
      }
    } else {
      SPrnnEuclidean(p->hison,querpoint,points,dimension,numpoints);
      if (SPnndist[1] >= val*val) {
        SPrnnEuclidean(p->loson,querpoint,points,dimension,numpoints);
      }
    }
  }
}

/*******************************************************************************/

void SPrnnGeneral(p,querpoint,points,dimension,numpoints,MinkP)

/*******************************************************************************/


sproullkdNode *p;
float *querpoint;
float **points;
int dimension,numpoints,MinkP;

{
  int i;
  float thisdist,val,thisx;

  if (p->bucket) {
    for (i=p->lopt; i <= p->hipt; i++) {
      thisdist=SPDistance(points,perm[i],querpoint,dimension,MinkP);

      if (sproullfound[0] < numpoints) {
        PQInsert(thisdist,perm[i],SPnndist,sproullfound);
      } else {
        PQreplace(thisdist,SPnndist,sproullfound,perm[i]);
      }
    }
  } else {
    val = -p->cutval; /* should this be negative ? -SCL */
    thisx=querpoint[p->discrim]; /* added this line, based mostly on optkd.o -SCL */
    for (i=0; i < dimension; i++) {
      val+=p->plane[i]*querpoint[i];
    }
    if (val < 0) {
      SPrnnGeneral(p->loson,querpoint,points,dimension,numpoints,MinkP);
      if (thisx + SPnndist[1] > val) {
        SPrnnGeneral(p->hison,querpoint,points,dimension,numpoints,MinkP);
      }
    } else {
      SPrnnGeneral(p->hison,querpoint,points,dimension,numpoints,MinkP);
      if (thisx - SPnndist[1] < val) {
        SPrnnGeneral(p->loson,querpoint,points,dimension,numpoints,MinkP);
      }
    }
  }
}


/*******************************************************************************/

int *kdSproullNNQuery(points,dimension, querpoint,numNN,Metric,root,MinkP)

/*******************************************************************************/

sproullkdNode *root;
float *querpoint, **points;
int dimension,numNN,MinkP;

{
  int j;

  /* set up found array */
  sproullfound = (int *) malloc((numNN+1)*sizeof(int));
  sproullfound[0]=1;  /* for now */

  /* SPnndist is a ordered list of the SPDistances of the nearest neighbors found */
  SPnndist = (float *)malloc((numNN+1)*sizeof(float));
  for (j=0; j < numNN+1; j++) {
    SPnndist[j] = 99999999999.0;
  }

  switch(Metric) {
    case EUCLIDEAN  : SPrnnEuclidean(root,querpoint,points,dimension,numNN);
                      break;
    case MANHATTAN  : SPDistance=ManhattDist;
                     SPrnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
    case L_INFINITY : SPDistance=LInfinityDist;
                     SPrnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
    case L_P        : SPDistance=LGeneralDist;
                     SPrnnGeneral(root,querpoint,points,dimension,numNN,MinkP);
                      break;
  }
  free(SPnndist);
  return(sproullfound);
}


/***************************************************************************/

void KillSproullTree(P)

/***************************************************************************/

/*  Kills a kd-tree to avoid memory holes.   */


sproullkdNode *P;

{
  if (perm != NULL) {
    free(perm);
  }  /* free permutation array */

  if (P==NULL) {
    return;
  } /* just to be sure */
  if (P->loson != NULL) {
    KillSproullTree(P->loson);
  }

  if (P->hison != NULL) {
    KillSproullTree(P->hison);
  }

  free(P);

}



