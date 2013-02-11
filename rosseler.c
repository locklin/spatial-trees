/*************************************************************************/
/*	                                                                 */
/* PROGRAM 								 */
/*   naive.c								 */
/*									 */
/* PURPOSE								 */
/*   To return all points in a given hyper-rectangle in a d dimensional  */
/*   space using a naive brute force approach.                           */
/*                                                                       */
/*************************************************************************/

/*************************************************************************/
/*									 */
/*                        Preprocessing Directives                       */
/* 									 */
/*************************************************************************/


#include <stdio.h>
#include <math.h>

#define N        4000

#define TESTINDEX(PP)  ((PP) = (ARRPTR)malloc(sizeof(ARREC)))
#define NEWINDEX(PP)   if (TESTINDEX(PP)==NULL) \
                              {printf("memory error\n");return;}
/* Creates new array index record */


#define Dimension 2
#define Eps       5
#define P         10



/*************************************************************************/
/*									 */
/*                          Type Declarations                            */
/* 									 */
/*************************************************************************/


typedef float
   POINT[Dimension];         /* A point in R^Dimension */

typedef float
   REALREGION[Dimension*2];  /* Region definition */

typedef struct ArrIndexLst   /* Array Index List */
{
  int                   Index;          /* Index to array in time series */
  struct ArrIndexLst    *Next;          /* Next point in the linked list */
} ARREC, *ARRPTR;



/*************************************************************************/
/*									 */
/*                            Global Variables                           */
/* 									 */
/*************************************************************************/

int
  DC,                          /* Dimension Counter */
  CellLength,
  K;

float
  A[N],                        /* Time series array */
  J,                           /* Counter */

  Low  =  9999999999.0,        /* Smallest number in time series */
  High = -999999999.0;         /* Largest number in time series */
			       /* Low and high are initially +- infinity */

POINT
  EP;

REALREGION
  RecDef;                       /* Range to search */

ARRPTR
  AListF,
  AListL,
  FoundF,                       /* First Array Index in range */
  FoundL,                       /* Last Array Index in range */
  Index;


/***************************************************************************/

void ReadTimeSeries()

/***************************************************************************/
/*                                                                         */
/* Insert code here to compute or read from a file the time series to be   */
/* analyzed.  								   */
/*								           */
/* Make sure High and Low are initialized to the highest and lowest        */
/* numbers in the time series.  				           */
/*                                                                         */
/***************************************************************************/

/* Computes Rossler sprial using 4th order Runge Kutta.  */
/* Not implemented intelligently. */

{
  float
    a=0.2,
    b=0.2,
    c=5.7,   /* a,b,c parameters */
    h=0.05,
    w[3],
    k[4][3],
    x,y,z;

  int
    j;


  w[0]=0.2;
  w[1]=0.1;
  w[2]=0;    /* Set initial conditions */

  for(K=0; K<1000; K++)   /* Discard first 1000 points */
  {

    x = w[0];
    y = w[1];
    z = w[2];
    k[0][0] = h*(z-y);
    k[0][1] = h*(x+a*y);
    k[0][2] = h*(b*x - c*z + x*z);

    x = w[0] + 0.5*(k[0][0]);
    y = w[1] + 0.5*(k[0][1]);
    z = w[2] + 0.5*(k[0][2]);
    k[1][0] = h*(z-y);
    k[1][1] = h*(x+a*y);
    k[1][2] = h*(b*x - c*z + x*z);

    x = w[0] + 0.5*(k[1][0]);
    y = w[1] + 0.5*(k[1][1]);
    z = w[2] + 0.5*(k[1][2]);
    k[2][0] = h*(z-y);
    k[2][1] = h*(x+a*y);
    k[2][2] = h*(b*x - c*z + x*z);

    x = w[0] + (k[2][0]);
    y = w[1] + (k[2][1]);
    z = w[2] + (k[2][2]);
    k[3][0] = h*(z-y);
    k[3][1] = h*(x+a*y);
    k[3][2] = h*(b*x - c*z + x*z);

    for (j=0; j<3; j++)
    {
      w[j] = w[j] + (k[0][j] + 2*k[1][j] + 2*k[2][j] + k[3][j])/6;
    }
  }


  for(K=0; K<N; K++)
  {

    x = w[0];
    y = w[1];
    z = w[2];
    k[0][0] = h*(z-y);
    k[0][1] = h*(x+a*y);
    k[0][2] = h*(b*x - c*z + x*z);

    x = w[0] + 0.5*(k[0][0]);
    y = w[1] + 0.5*(k[0][1]);
    z = w[2] + 0.5*(k[0][2]);
    k[1][0] = h*(z-y);
    k[1][1] = h*(x+a*y);
    k[1][2] = h*(b*x - c*z + x*z);

    x = w[0] + 0.5*(k[1][0]);
    y = w[1] + 0.5*(k[1][1]);
    z = w[2] + 0.5*(k[1][2]);
    k[2][0] = h*(z-y);
    k[2][1] = h*(x+a*y);
    k[2][2] = h*(b*x - c*z + x*z);

    x = w[0] + (k[2][0]);
    y = w[1] + (k[2][1]);
    z = w[2] + (k[2][2]);
    k[3][0] = h*(z-y);
    k[3][1] = h*(x+a*y);
    k[3][2] = h*(b*x - c*z + x*z);

    for (j=0; j<3; j++)
    {
      w[j] = w[j] + (k[0][j] + 2*k[1][j] + 2*k[2][j] + k[3][j])/6;
    }

    A[K] = w[1];

    if (A[K] > High)
      High=A[K];

    if (A[K] < Low)
      Low=A[K];
  }

}


/**************************************************************************/

int main()
{
  int k, j;
  FILE *f;

  k=0;
  ReadTimeSeries();
  f=fopen("rosseler2d.rg","w");

  if (f==NULL) {
    return(k);
  }

  fprintf(f,"%d\n",N-(k+(Dimension-1)*P));  /* write number of points */
  fprintf(f,"%d\n",Dimension);  /* write dimension */

  while (k+(Dimension-1)*P < N) {
    for (j=0; j < Dimension; j++) {
      fprintf(f,"%f\n",A[(j*P)+k]);
    }
    k++;     
  }
  fclose(f);
}


