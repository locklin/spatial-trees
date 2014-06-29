/*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "nn.h"
#include "optkd.h"
#include "naivekd.h"
#include "sproullkd.h"
#include "vptree.h"

/* file reading */
extern float **ReadInputFile();


/* from naive.c */
extern int *NaiveRectQuery();
extern int *NaiveNNQuery();

/* from naivekd.c */
extern TreePtr kdNaiveBuildTree();
extern int *kdNaiveRectQuery();
extern void KillNaiveTree();

/* from optkd.c */
extern optkdNode *BuildOptTree();
extern int *kdOptNNQuery();
extern int *DrawkdOptNNQuery2d();
extern int *DrawkdOptNNQuery3d();
extern void KillOptTree();
extern int *kdOptRectQuery();

/* from sproullkd.c */
extern sproullkdNode *BuildSproullkdTree();
extern int *kdSproullNNQuery();
extern void KillSproullTree();

/* from vptree.c */
extern VPTreeNode *BuildVPTree();
extern int *vpNNQuery();
extern void KillVPtree();
float **(*DistribFun)();

float **Create_Points_Array(n,d)

/*******************************************************************************/
/*  									       */
/*  Creates a dynamically allocated two dimensional array consisting of the    */
/*  number of points, n, each with a corresponding d-vector.                   */
/*									       */
/*******************************************************************************/

int 
     n,   /* number of points */
     d;   /* dimension of the number of points */

{
  int i;
  float **x;

  x=(float **)malloc(n * sizeof(float *));
  for (i =0; i < n; i++) {
    x[i] = (float *)malloc(d*sizeof(float));
  }
  return(x);
}
  
int main(int argc, char* argv[]) {
  /*  float **points,*point; */
  float **fpoints,*point,*dists,**qpts;
  int dim,num,nn,select;
  int *retdx;
  char *filename;
  /*  TreePtr *tree; */
  optkdNode *tree;  
  /* sproullkdNode *tree; */
  /* VPTreeNode *tree; */

  if((argc==1)|(argc>4)) {
    printf("usage: ./test filename\n./test nn filename\n./test selected nn filename\ndefaults to selected=10,nn=2");
    return(-1);
  } else {
    select=10;
    nn =2;
    if(argc==2) {
      filename=argv[1];
    } else if (argc==3) {
      filename=argv[2];
      nn=atoi(argv[1]);
    } else if (argc==4) {
      select=atoi(argv[1]);
      nn=atoi(argv[2]);
      filename=argv[3];
    }

    fpoints = ReadInputFile(filename,&dim,&num);
    printf("%d dim    %d num\n",dim,num); 
    point=fpoints[select];
    printf("%f     %f\n",point[0],point[1]);  
    /* retdx = NaiveNNQuery(fpoints,num,dim,fpoints[select],nn,EUCLIDEAN,1); */
    /* tree = kdNaiveBuildTree(fpoints,num,dim); */
    /* printf("built tree\n"); */
    /* qpts= (float **)malloc((nn+1)*sizeof(float)); */
    /* qpts[0] = fpoints[select]; */
    /* qpts[1] = fpoints[select]; */
    /* retdx = kdNaiveRectQuery(tree,fpoints,dim,fpoints); */
    /*  KillNaiveTree(tree); */
    /*  printf("killed tree\n");         */
  
    tree = BuildOptTree(fpoints,num,dim);
    retdx=kdOptNNQuery(fpoints,dim,fpoints[select],nn,EUCLIDEAN,tree,1);
   /* dists = (float *)malloc((nn+1)*sizeof(float)); */
    /* retdx=kdOptNNQueryDist(fpoints,dim,fpoints[select],dists,nn,EUCLIDEAN,tree,1);  */

    /* seems to return index in reverse order */
     /* tree = BuildSproullTree(fpoints,num,dim);  */
     /* retdx = kdSproullNNQuery(fpoints,dim,fpoints[select],nn,EUCLIDEAN,tree,1);  */
    /*the way this is defined, l..u covariance rather than num*/
    /* seems to not work -no surprise there */
    /* tree = BuildVPTree(fpoints,num,dim); */
    /* printf("built tree\n");       */
    /* retdx=vpNNQuery(fpoints,dim, fpoints[select],nn,EUCLIDEAN,tree,1); */
    printf("%d number returned %d index returned\n",retdx[0],retdx[1]);
    if(nn>1) {
      printf("next nearest index is %d\n",retdx[2]);
    }
    /* now a kdtree search ... */

  return(1);
  }
}
