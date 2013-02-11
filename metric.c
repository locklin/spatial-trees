/* module contains implementations of the distance function for the 
   Euclidean, Manhattan, and L_Infinity metric */

extern double fabs(),pow();

/****************************************************************************/

float EuclidDist2(Points,Index,NNQPoint,dimension)

/****************************************************************************/

/* returns the square of the Euclidean distance of Points[Index] and
   NNQPoint */

float **Points,*NNQPoint;
int Index,dimension;

{
float dist,d;
int j;

dist=0;
for (j=0; j<dimension; j++) {
  d=Points[Index][j]-NNQPoint[j];
  dist += d*d;
}
return(dist);

}
/****************************************************************************/

float ManhattDist(Points,Index,NNQPoint,dimension)

/****************************************************************************/

/* returns the manhattan distance between Points[Index] and
   NNQPoint */

float **Points,*NNQPoint;
int Index,dimension;

{
float dist,d;
int j;

dist=0;
for (j=0; j<dimension; j++) {
  d=Points[Index][j]-NNQPoint[j];
  dist = dist + fabs(d);
}

return(dist);
}

/****************************************************************************/

float LInfinityDist(Points,Index,NNQPoint,dimension)

/****************************************************************************/

/* returns the square of the Euclidean distance of Points[Index] and
   NNQPoint */

float **Points,*NNQPoint;
int Index,dimension;

{
float dist,d;
int j;

dist=-999999999.0;
for (j=0; j<dimension; j++) {
  d=fabs(Points[Index][j]-NNQPoint[j]);
  if (dist < d) {
    dist = d;
  }
}
return(dist);
}


/****************************************************************************/

float LGeneralDist(Points,Index,NNQPoint,dimension,MinkP)

/****************************************************************************/

/* returns the square of the Euclidean distance of Points[Index] and
   NNQPoint */

float **Points,*NNQPoint;
int Index,dimension,MinkP;

{
float dist;
int j;

dist=0;
for (j=0; j<dimension; j++) {
  dist +=fabs(pow(Points[Index][j]-NNQPoint[j],(float)MinkP));
}
return(pow(dist,1.0/(float)MinkP));
}
















