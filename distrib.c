/*************************************************************************/
/*  Code to generate the non-uniform data sets as described by:          */
/*								         */
/*      Bentley, J.L. "k-D Trees for Semidynamic Point Sets.             */
/*							        	 */
/*  The distributions generated are as follows:			         */
/*     								  	 */
/*       uniform  - uniform on the unit hypercube 		         */
/*       annulus  - uniform on the edge of a circle		         */
/*       arith    - x0 = 0,1,4,...n^2    x1=x2=...xd-1 = 0               */
/*       ball     - uniform inside a circle			         */
/*       clusnorm - Normal(0.05) at 10 points in the unit hypercube      */
/*       cubediam - x0=x1=...=xd-1  where x0 comes from U[0,1]	         */
/*       cubeedge - x0=U[0,1], x1=...=xk-1 =0			         */
/*       corners  - x0 and x1 are taken Uniformily from rectangular      */
/*                  regions around (0,0),(2,0), (0,2), and (2,2).        */
/*                  x2..xd-1 are from the unit interval [0,1].           */
/*       grid     - N points chosen from a grid hypercube with 1.3N      */
/*                  partitions.                                          */
/*       normal   - each dimension from Normal(1)                        */
/*       spokes   - N/D at (U[0,1],1/2,1/2,...,1/2)                      */
/*                  N/D at (1/2,U[0,1],1/2,...,1/2)                      */
/*							   	         */
/*************************************************************************/
#include <stdlib.h>
#define MAXINT 2147483647.0   /* (2^31)-1 */
#define PI 3.1415926535

extern float **Create_Points_Array();
extern double sin(),cos(),sqrt(),floor(),log(),pow();

/***************************************************************************/

float Normal(std)

/***************************************************************************/
/*                                                                         */
/* returns a normally distributed random number with standard deviation    */
/* std.  Based on Numerical Recipies in C.                                 *
/*                                                                         */
/***************************************************************************/

float
  std;

{
static int iset=0;
static float gset;
float fac,r,v1,v2;
float ran1();

if (iset==0) {
  do {
    v1=2.0*(random()/MAXINT) - 1.0;
    v2=2.0*(random()/MAXINT) - 1.0;
    r=v1*v1+v2*v2;
  } while (r>=1.0 || r == 0.0);

  fac = sqrt(-2.0*log(r)/r);
  gset=v1*fac;
  iset=1;
  return(v2*fac*std);
}
else {
  iset=0;
  return(gset*std);
}
}


/*****************************************************************************/

float **Create_Points_Array(n,d)

/*****************************************************************************/
/*   								             */
/*  Creates a dynamically allocated two dimensional array consisting of the  */
/*  number of points, n, each with a corresponding d-vector.                 */
/*									     */
/*****************************************************************************/

int
  n,   /* number of points */
  d;   /* dimension of the number of points */

{
  int i;
  float **x;
  /*  char *malloc();  this fucks things up -SCL*/

  x=(float **)malloc(n * sizeof(float *));
  for (i =0; i < n; i++) {
    x[i] = (float *)malloc(d*sizeof(float));
  }
  return(x);
}


/*************************************************************************/

float **GenerateUniform(n,d)

/*************************************************************************/
/*							                 */
/*  Returns a pointer to a two dimensional array with points distributed */
/*  uniformly in the unit hypercube of dimension d.  The first dimension */
/*  of the array is the 						 */
/*  number of points, n.  The second dimension is of the array is the    */
/*  d-vectors (each point) in R^d.			                 */
/*									 */
/*************************************************************************/

int n,d;

{
float **arr;
int i,j;
  arr=Create_Points_Array(n,d);
  
  for (i=0; i < n; i++) {
    for(j=0; j < d; j++) {
      arr[i][j]=(random()/MAXINT) - 0.5;
    }
  }
 
  return(arr);

}

/**************************************************************************/

float **GenerateAnnulus(n,d)

/**************************************************************************/
/*							                  */
/*  Returns a pointer to a two dimensional array containing n d-dim points*/
/*  A point is generated as follows:                                      */
/*    x0,x1 are uniformly distributed on a circle.                        */
/*    x2,..xd-1 are uniformly random on the unit interval.                */
/*							                  */
/**************************************************************************/


int n,d;

{
float **arr;
float theta;

int i,j;
  arr=Create_Points_Array(n,d);
  
  for (i=0; i < n; i++) {
    theta=(random()/MAXINT)*2.0*PI;

    arr[i][0]=sin(theta); 
    arr[i][1]=cos(theta);   /* x0,x1 are on a circle */


    for (j=2; j < d; j++) {
      arr[i][j]=(random()/MAXINT) - 0.5;
    }  /* the rest are uniformly random */
  }
  return(arr);
}


/**************************************************************************/

float **GenerateArith(n,d)

/**************************************************************************/
/*							                  */
/*  Returns a pointer to a two dimensional array containing n d-dim points*/
/*  A point is generated as follows:                                      */
/*    x0 = 0,1,4,16,...n^2                                                */
/*    x1=x2=x3..xd-1=0                                                    */
/*							                  */
/**************************************************************************/


int n,d;

{
float **arr;

int i,j;
  arr=Create_Points_Array(n,d);
  
  for (i=0; i < n; i++) {
    arr[i][0]=i*i; /* the first coord is the point number squared */

    for (j=1; j < d; j++) {
      arr[i][j]=0;
    }  /* the rest are zero */
  }
  return(arr);
}

/**************************************************************************/

float **GenerateBall(n,d)

/**************************************************************************/
/*							                  */
/*  Returns a pointer to a two dimensional array containing n d-dim points*/
/*  A point is uniformly chosen from inside a d-dimensional sphere.       */
/*                                                                        */
/**************************************************************************/

int n,d;

{
float **arr;
float magnitude,distance;
float *directionvector;

int i,j;
  arr=Create_Points_Array(n,d);


  directionvector=(float *)malloc(d*sizeof(float));

  for (i=0; i < n; i++) {
    distance = pow((random()/MAXINT),1.0/(float)d);
    magnitude=0.0;

    for (j=0; j < d; j++) {
      directionvector[j]=Normal(1.0); /* (2.0*(random()/MAXINT) - 1.0); */
      magnitude=magnitude+directionvector[j]*directionvector[j];
    }  

    magnitude=sqrt(magnitude);
    
    for (j=0; j < d; j++) {
      arr[i][j] = (directionvector[j]/magnitude)*distance;
    }
  }

/*
 i=0;
 while (i < n) {
   magnitude=0;
   for (j=0; j<d; j++) {
     arr[i][j] = 2.0*(random()/MAXINT)-1.0;
     magnitude = magnitude + arr[i][j]*arr[i][j];
   }
   if (magnitude < 1) { i++;}
 }
*/
  return(arr);
}

/**************************************************************************/

float **GenerateClusnorm(n,d)

/**************************************************************************/
/*							                  */
/*  Returns a pointer to a two dimensional array containing n d-dim points*/
/*  A point is normally chosen to be a distance of 0.05 from a one of ten */
/*  uniformly chosen point inside the unit hypercube.                     */
/*                                                                        */
/**************************************************************************/

int n,d;

{
float **arr;
float **cluspoints;
int i,j;

  arr=Create_Points_Array(n,d);
  cluspoints=Create_Points_Array(10,d);
  for (i=0; i < 10; i++) {
    for (j=0; j < d; j++) {
      cluspoints[i][j] = (random()/MAXINT) - 0.5;
    }
  } /* create random points */

  for (i=0; i < n; i++) {
    for (j=0; j < d; j++) {
      arr[i][j] = cluspoints[(i % 10)][j] + Normal(0.05);
    }
  }
  return(arr);
}


/**************************************************************************/

float **GenerateCubediam(n,d)

/**************************************************************************/
/*							                  */
/*  Returns a pointer to a two dimensional array containing n d-dim points*/
/*  A point is uniformly chosen from the line x0=x1=x2=x3=...xd-1         */
/*							                  */
/**************************************************************************/


int n,d;

{
float **arr;

int i,j;
  arr=Create_Points_Array(n,d);
  
  for (i=0; i < n; i++) {
    arr[i][0] = random()/MAXINT;
    for (j=1; j < d; j++) {
      arr[i][j]=arr[i][0];
    }  
  }
  return(arr);
}


/**************************************************************************/

float **GenerateCubeedge(n,d)

/**************************************************************************/
/*							                  */
/*  Returns a pointer to a two dimensional array containing n d-dim points*/
/*  A point is generated as follows:                                      */
/*    							                  */
/*    x0 is selected from U[0,1].                                         */
/*    x1,x2,...xd-1 are all zero                                          */
/**************************************************************************/


int n,d;

{
float **arr;

int i,j;
  arr=Create_Points_Array(n,d);
  
  for (i=0; i < n; i++) {
    arr[i][0]=random()/MAXINT; /* the first coord is in U[0,1] */

    for (j=1; j < d; j++) {
      arr[i][j]=0;
    }  /* the rest are zero */
  }
  return(arr);
}


/**************************************************************************/

float **GenerateCorners(n,d)

/**************************************************************************/
/*							                  */
/*  Returns a pointer to a two dimensional array containing n d-dim points*/
/*  The x0, and x1 coordinates of a point are chosen uniformly from the   */
/*  rectangle of length one centered at one of the following points:      */
/*  (-1,-1), (-1,1), (1,-1), (1,1). The rest of the dimensions are chosen */
/*  uniformly  from inside the unit interval [0,1].                       */
/*                                                                        */
/**************************************************************************/

int n,d;

{
float **arr,*center;
int i,j;


  arr=Create_Points_Array(n,d);
  center = (float *)malloc(d*sizeof(float));
  
  for (j=0; j<d; j++) {
    center[j] = 0.0;
  } /* used to compute the center of gravity */

  for (i=0; i < n; i++) {
    for (j=0; j < 2; j++) {
      if (random()/MAXINT < 0.5) {
        arr[i][j] = random()/MAXINT - 1.0;
        center[j]+=arr[i][j];
      }
      else {
        arr[i][j] = random()/MAXINT + 1.0;
        center[j]+=arr[i][j];
      }
    }
    for (j=2; j < d; j++) {
      arr[i][j] = random()/MAXINT;
      center[j]+=arr[i][j];
    }
  }
  
  for (j=0; j<d; j++) {
    center[j] = center[j]/(float)n;
  }

  for (i=0; i < n; i++) {
    for (j=0; j < d; j++) {
      arr[i][j]=arr[i][j]-center[j];
    }
  }

  return(arr);
}

/**************************************************************************/

float **GenerateGrid(n,d)

/**************************************************************************/
/*							                  */
/*  Returns a pointer to a two dimensional array containing n d-dim points*/
/*  The points are chosen from a square grid of 1.3N points.              */
/**************************************************************************/

int n,d;

{
float **arr,gridsize;
int i,j;

  arr=Create_Points_Array(n,d);
  gridsize = pow((1.3*(float)n),1.0/((float)d));
  for (i=0; i < n; i++) {
    for (j=0; j < d; j++) {
      arr[i][j]=floor(gridsize*(random()/MAXINT))-(gridsize/2.0);
    }
  }
  return(arr);
}


/*************************************************************************/

float **GenerateNormal(n,d)

/*************************************************************************/
/*							                 */
/*  Returns a pointer to a two dimensional array with points distributed */
/*  normally in the unit hypercube of dimension d.  The first dimension  */
/*  of the array is the 						 */
/*  number of points, n.  The second dimension is of the array is the    */
/*  d-vectors (each point) in R^d.			                 */
/*									 */
/*************************************************************************/

int n,d;

{
float **arr;
int i,j;

  arr=Create_Points_Array(n,d);
  
  for (i=0; i < n; i++) {
    for(j=0; j < d; j++) {
      arr[i][j]=Normal(1.0);
    }
  }
  return(arr);
}

/**************************************************************************/

float **GenerateSpokes(n,d)

/**************************************************************************/
/*							                  */
/* Returns a pointer to a two dimensional array containing n d-dim points */
/*    N/K points are at (U[0,1],1/2,...,1/2),                             */
/*    N/K points are at (1/2,U[0,1],...,1/2),                             */
/*    ...                                                                 */
/**************************************************************************/

int n,d;

{
float **arr;
int i,j;

  arr=Create_Points_Array(n,d);

  for (i=0; i < n; i++) {
    for (j=0; j < d; j++) {
      arr[i][j]=0.0;
    }
    arr[i][i % d] = (random()/MAXINT)-0.5;
  }
  return(arr);
}

/**************************************************************************/

float **GenerateHenon(n,d)

/**************************************************************************/
/*							                  */
/* Returns a pointer to a two dimensional array containing n d-dim points */
/*    Each points come from the henon attractor.                          */
/**************************************************************************/

int n,d;

{
float **arr,*center;
register int i,j;
double x1, y1, x0, y0;
int init=1000;

  arr=Create_Points_Array(n,d);
  center = (float *)malloc(d*sizeof(float));

  x0 = 0.2; y0 = 0.2;
  for (i=0;i < init; ++i)
  {
    x1 = 1.0 - 1.4*x0*x0 + y0;
    y1 = 0.3*x0;
    x0 = x1; y0 = y1;
  }

  for (i=0; i< n; i++) {
    for (j=0; j<d; j++) {
      x1 = 1.0 - 1.4*x0*x0 + y0;
      y1 = 0.3*x0;
      x0 = x1; y0 = y1; 
      arr[i][j]=(float)x1;
      center[j]+=(float)x1;
    }
  }
  
  for (j=0; j<d; j++) {
    center[j] = center[j]/(float)n;
  }

  for (i=0; i < n; i++) {
    for (j=0; j < d; j++) {
      arr[i][j]=arr[i][j]-center[j];
    }
  } /* center of gravity computations */

  return(arr);
}




