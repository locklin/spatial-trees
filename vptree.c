/**************************************************************************/
/* Program to perform nearest neighbor queries in a VPtree.               */
/* In this implementation, we use the simplest VPtree (Algorithm 1) with  */
/* four values retained indicating lower/uperbounds of each subspace as   */
/* seen by the vantage point.                                             */
/*                                                                        */
/* References:  P.N. Yianilos "Data Structures and Algorithms for Nearest */
/* Search in General Metric Spaces.  SODA '93.                            */
/*                                                                        */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "vptree.h"
#include "nn.h"

/* from metric.c */
extern float EuclidDist2();
extern double sqrt();
/*from pqueue.c */
extern void PQInsert();
extern void PQreplace();

/* Used to create a new tree in the k-d tree */
#define TESTTREE(PP)  ((PP) = (VPTreeNode*)malloc(sizeof(VPTreeNode)))
#define NEWTREE(PP)  if (TESTTREE(PP)==NULL) \
                         {printf("memory error\n");return;}

#define TESTINDEX(PP) ((PP)=(IndexRec *)malloc(sizeof(IndexRec)))
#define NEWINDEX(PP) if (TESTINDEX(PP)==NULL) \
                         {printf("memory error\n");return;}

/* Global Variables */
int *VPfound;
float *tau;

/***************************************************************************/

void MedianDistance(a,N,k)

/***************************************************************************/
/* Makes the VPPerm partition the array Values along the element k.        */
/* Adapted from Sedgewick's Algorithms in C (p. 128)                       */
/***************************************************************************/

float *a;
int k,N;

{
float v,t;
int i,j,r,l;

l=1;r=N;

while(r>l) {
  v=a[r]; i=l-1; j=r;
  for (;;) {
    while (a[++i] < v);
    while (a[--j] > v && j>l); 
    if (i >= j) break;
    t=a[i]; a[i] = a[j]; a[j]=t;
  }
  t=a[i]; a[i] = a[r]; a[r]=t;
  if (a[i]>=(float)k) r=i-1;
  if (a[i]<=(float)k) l=i+1;
}

}


/***************************************************************************/

IndexRec *Select_vp(Points,S,numPts,dimension,median)

/***************************************************************************/

float **Points, *median;
int dimension,numPts;
IndexRec *S;

{
  int k,m;
  float spread,best, bestmu,mu, *dist,ex,ex2;
  IndexRec *P,*D,*index;  /* P,D are supposed to be random sample of S, but
                      I've chosen to choose them all */


  if ((dist = (float *) malloc((numPts)*sizeof(float)))==NULL) {
    printf("select vp memory error\n");
  } 
  
  best=-1.0; P=S->Next;
  while (P!=NULL) {
    mu=0.0;m=-1; D=S->Next;
    while (D!=NULL) {
     if (P->Index!=D->Index) {
        m++;
        mu = (float)sqrt(EuclidDist2(Points,P->Index,Points[D->Index],dimension));
        dist[m] = mu/(1.0+mu); 
      }
     D=D->Next;
    } 
    MedianDistance(dist,m,m/2);
    mu = dist[m/2];
    ex=0.0;ex2=0.0;
    for (k=0; k<=m; k++) {
      ex+=(dist[k]-mu);
      ex2+=(dist[k]-mu)*(dist[k]-mu);
    }
    spread=ex2-ex*ex;  /* compute varience */
    if (spread > best) {
      best=spread;
      bestmu=mu;
      index=P;
    }
    P=P->Next;
  }
  *median=bestmu;
  free(dist);
  return(index);
}

/*************************************************************************/

VPTreeNode *MakeVPTree(Points,S,numPts,dimension)

/*************************************************************************/

float **Points;
int dimension,numPts;
IndexRec *S;
{
VPTreeNode *p;
int nl,nr;
float mu, dist;
IndexRec *Index,*prev,*next,*L,*R,*LPrev,*RPrev;

  if (S->Next==NULL) {
    return(NULL);
  }

  if (numPts == 1) {
    p = (VPTreeNode*)malloc(sizeof(VPTreeNode));
    /* NEWTREE(p); global weirds -SCL */
    p->Index=S->Next->Index;
    p->Mu = 0.0; 
    p->ll=5.0; p->lu=-5.0;
    p->rl=5.0; p->ru=-5.0;  /* useful to reduce search time */
    p->Left=NULL;p->Right=NULL;
    return(p);
  }  /* unroll recursion by one step */

   NEWTREE(p); /*global weirds -SCL */
  /*  p = (VPTreeNode*)malloc(sizeof(VPTreeNode)); */
  Index=Select_vp(Points,S,numPts,dimension,&mu);

  p->Index=Index->Index;
  p->Mu=mu;

  prev=S; next=S->Next;
    NEWINDEX(L);NEWINDEX(R); /*global and ugly weirds! */
  /* L = (IndexRec *)malloc(sizeof(IndexRec));*/
  /*R = (IndexRec *)malloc(sizeof(IndexRec));*/
  LPrev=L;RPrev=R; nl=0;nr=0;
  while (next!=NULL) {
    if (next->Index!=Index->Index) {
      dist=sqrt(EuclidDist2(Points,next->Index,Points[Index->Index],dimension));
      dist=dist/(1.0 + dist);
      if (dist < mu) { /* add to left list */
        if (p->ll > dist) {
          p->ll = dist;
        }
        if (p->lu < dist) {
          p->lu=dist;
        }
        LPrev->Next=next;
        LPrev=next;
        nl++;
      } 
      else {
        if (p->rl > dist) {
          p->rl=dist;
        }
        if (p->ru < dist) {
          p->ru = dist;
        }
        RPrev->Next=next;            
        RPrev=next;
        nr++;
      }
    } else {
      prev=next;
      free(prev); /* kill current point from list*/    
    }
    next=next->Next;
  }
  LPrev->Next=NULL;RPrev->Next=NULL;
  free(S);
  p->Left=MakeVPTree(Points,L,nl,dimension);
  p->Right=MakeVPTree(Points,R,nr,dimension);
  return(p);
}


/*************************************************************************/

VPTreeNode *BuildVPTree(Points,numPoints,dimension)

/*************************************************************************/

int dimension,numPoints;
float **Points;

{
int j;
VPTreeNode *root;
IndexRec *FirstIndex,*PrevIndex,*NextIndex;

/* NEWINDEX(PrevIndex); global weirds -SCL */
  PrevIndex = (IndexRec *)malloc(sizeof(IndexRec));
  FirstIndex=PrevIndex;  /* dummy header */
  for (j=0; j<numPoints; j++) {
    /* NEWINDEX(NextIndex); global weirds -SCL */
    NextIndex = (IndexRec *)malloc(sizeof(IndexRec));
    NextIndex->Index=j;
    PrevIndex->Next=NextIndex;
    PrevIndex=NextIndex;
  }
  NextIndex->Next=NULL;

  root=MakeVPTree(Points,FirstIndex,numPoints,dimension);
  return(root);
}


/*******************************************************************************/

void VPnnSearch(p,querpoint,points,dimension,numpoints,MinkP)

/*******************************************************************************/


VPTreeNode *p;
float *querpoint;
float **points;
int dimension,numpoints,MinkP;

{
  float x,middle;

  if (p==NULL) {
    return;
  }
  x=sqrt(EuclidDist2(points,p->Index,querpoint,dimension,MinkP));
  x=x/(1.0+x);

  if (VPfound[0] < numpoints) {
    PQInsert(x,p->Index,tau,VPfound);
  } else {
    PQreplace(x,tau,VPfound,p->Index);
  }

  middle=((p->lu)+(p->rl))/2.0;

  if (x < middle) {
    if ((x > (p->ll - tau[1])) && (x < (p->lu + tau[1]))) {
      VPnnSearch(p->Left,querpoint,points,dimension,numpoints,MinkP);
    }
     
    if ((x > (p->rl - tau[1])) && (x < (p->ru + tau[1]))) {
      VPnnSearch(p->Right,querpoint,points,dimension,numpoints,MinkP);
    }
  }
  else {
    if ((x > (p->rl - tau[1])) && (x < (p->ru + tau[1]))) {
      VPnnSearch(p->Right,querpoint,points,dimension,numpoints,MinkP);
    }
    if ((x > (p->ll - tau[1])) && (x < (p->lu + tau[1]))) {
      VPnnSearch(p->Left,querpoint,points,dimension,numpoints,MinkP);
    }
  }
}


/*****************************************************************************/

int *vpNNQuery(points,dimension, querpoint,numNN,Metric,root,MinkP)

/*****************************************************************************/

VPTreeNode *root;
float *querpoint, **points;
int dimension,numNN,MinkP;

{
  int j;

  /* set up found array */
  VPfound = (int *) malloc((numNN+1)*sizeof(int));
  VPfound[0]=1;  /* for now */
   
  /* tau is a priority queue of the distances of the nearest neighbors found */
  tau = (float *)malloc((numNN+1)*sizeof(float));
  for (j=0; j < numNN+1; j++) {
    tau[j] = 99999999999.0;
  }

  switch(Metric) {
    case EUCLIDEAN : 
    case MANHATTAN : 
    case L_INFINITY: 
    case L_P       : 
                      VPnnSearch(root,querpoint,points,dimension,numNN,MinkP);
                      break;
  }
  free(tau);
  return(VPfound);
}


/***************************************************************************/

void KillVPTree(P)

/***************************************************************************/

/*  Kills a kd-tree to avoid memory holes.   */


VPTreeNode *P;

{
  if (P==NULL) {
    return;
  } /* just to be sure */
  if (P->Left != NULL) {
    KillVPTree(P->Left);
  }

  if (P->Right != NULL) {
    KillVPTree(P->Right);
  }

  free(P);

}
