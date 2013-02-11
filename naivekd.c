/**************************************************************************/
/* Program to perform orthogonal range searches in a very naive k-d tree. */
/* In this implementation the nodes on any given level of the tree all    */
/* have the same discriminating dimension and the discriminator is chosen */
/* as NextDisc(i) = i+1 mod k.                                            */
/*                                                                        */
/* Later refinements employ much more sophisticated methods for the       */
/* selection of the discriminator.                                        */
/*                                                                        */
/* References:  J.L. Bentley "Multidimensional Binary Search Trees used   */
/* for Associative Searching.  ACM Sept. 1975 Vol. 18 No. 9.              */
/*                                                                        */
/**************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "naivekd.h"
/* Used to create a new tree in the k-d tree */
#define TESTTREE(PP)  ((PP) = (TreePtr)malloc(sizeof(TreeRec)))
#define NEWTREE(PP)  if (TESTTREE(PP)==NULL) \
                         {printf("memory error\n");return;}


/* Global Variables */
int *found;
int count=0;

TreePtr Root;

/*contains a pointer to an integer array of indices of points that were
  found in the rectangular query */

/***************************************************************************/

void InRegion(P,Dimension,Points,RectQuery)

/***************************************************************************/
/* Determines if the treenode P falls inside the rectangular query         */
/* RectQuery.  If so, adds the array index of the point to the found       */
/* array.                                                                  */
/***************************************************************************/

TreePtr P;
int Dimension;
float **Points, **RectQuery;

{
  int index, dc;

  index = P->Index;

  for (dc=0; dc < Dimension; dc ++) {
    if (Points[index][dc] < RectQuery[dc][0] || 
        Points[index][dc] > RectQuery[dc][1]) {
         return;
     }
  }
  /* P is in the region */
  count++;
  if ((count % 4000) == 0) {
    realloc(found,(4000+count)*sizeof(int));
  }
  if (found == NULL) {
    printf("we have a memory problem\n");
  }

  found[count] = index;
}

/***************************************************************************/

int BoundsIntersectRegion(B,RectQuery,Dimension)

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

void RangeSearch(P,Points,Dimension,RectQuery,B)

/***************************************************************************/

TreePtr P;
float **RectQuery, **Points, *B;
int Dimension;


{
  int dc, disc;
  float *BHigh,*BLow;
  

  if (P==NULL) {printf("somehow a null pointer got sent here\n");}
  disc=P->Discrim;

  /* Check to see if the P is in the region */
  InRegion(P,Dimension,Points,RectQuery);

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
  BLow[2*disc+1] = Points[P->Index][disc];
  BHigh[2*disc] = Points[P->Index][disc];

  if (P->Left != NULL && BoundsIntersectRegion(BLow,RectQuery,Dimension)) {
    RangeSearch(P->Left,Points,Dimension,RectQuery,BLow);
  }
  free(BLow);
  if (P->Right != NULL && BoundsIntersectRegion(BHigh,RectQuery,Dimension)) {
    RangeSearch(P->Right,Points,Dimension,RectQuery,BHigh);
  }
  free(BHigh);
}


/***************************************************************************/

int *kdNaiveRectQuery(root,Points,Dimension,RectQuery)

/***************************************************************************/

TreePtr root;
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

  count=0;
  found = (int *)(malloc(4000*sizeof(int)));
  if (found == NULL) {
    printf("We have a memory problem\n");
  }

  RangeSearch(root,Points,Dimension,RectQuery,B);
  free(B);
  found[0] = count;
  return(found);
}


/**************************************************************************/

int Successor(P,Q,Dimension,Points)

/**************************************************************************/

TreePtr P,Q;
int Dimension;
float **Points;

{
int dc,disc;
  disc = Q -> Discrim;
  for (dc=0; dc < Dimension; dc++) {
    if (Points[P->Index][(disc+dc) % Dimension] <
        Points[Q->Index][(disc+dc) % Dimension]) 
      {
        return(-1);
      }
    if (Points[P->Index][(disc+dc) % Dimension] >
        Points[Q->Index][(disc+dc) % Dimension]) 
      {
        return(1);
      }
  }
  return(0);  /* they must be equal */

}


/***************************************************************************/

void InsertNode(P,Points,Dimension)

/***************************************************************************/

/* Inserts a node into the naive k-d tree */
float **Points;
TreePtr P;
{
  int succ;
  TreePtr Son,Q;

  /* check for root */
  if (Root == NULL) {
    Root = P;
    P->Discrim = 0;
    P->Left=NULL;
    P->Right=NULL;
    return;
  }

  Son=Root;

  do {
    Q=Son;
    succ = Successor(P,Q,Dimension,Points);
    switch(succ) {
      case -1: Son=Q->Left;
               break;
      case  1: Son=Q->Right;
               break;
      case  0: return;  /* don't insert the point */
    }
  } while (Son != NULL);

  /* ASSERT: Q points to the leaf of the tree that needs to be added */
  if (succ==-1) {
    Q->Left = P;
  }
  else {
    Q->Right =P;
  }

  P->Discrim = (Q->Discrim + 1) % Dimension;
  P->Left = NULL;
  P->Right=NULL;  
}

/**************************************************************************/

TreePtr kdNaiveBuildTree(Points,NumPoints,Dimension)

/**************************************************************************/
float **Points;
int NumPoints,Dimension;

{
int k;
TreePtr Node;

Root=NULL;
for (k=0; k < NumPoints; k++) {
  /* NEWTREE(Node); this is some weird kind of global thing; fuckit */
  Node = (TreePtr)malloc(sizeof(TreeRec)); 
  Node->Index = k;
  InsertNode(Node,Points,Dimension);
}

if (Root == NULL) 
{
printf("returning null tree\n");
}
return(Root);

}

/***************************************************************************/

void KillNaiveTree(P)

/***************************************************************************/

/*  Kills a kd-tree to avoid memory holes.   */


TreePtr P;

{

  if (P==NULL) {
    return;
  } /* just to be sure */
  if (P->Left != NULL) {
    KillNaiveTree(P->Left);
  }

  if (P->Right != NULL) {
    KillNaiveTree(P->Right);
  }

  free(P);

}
