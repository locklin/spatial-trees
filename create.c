/*******************************************************************************/
#include <stdlib.h>
double **Create_Points_Array(n,d)

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
  double **x;

  x=(double **)malloc(n * sizeof(double *));
  for (i =0; i < n; i++) {
    x[i] = (double *)malloc(d*sizeof(double));
  }
  return(x);
}
  
void main() {
  double **points,*point;

  points=Create_Points_Array(4,2);
  points[2][0]=1.0;points[2][1]=2.0;

  point=points[2];
  printf("%f     %f\n",point[0],point[1]);
}
