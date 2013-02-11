#include <stdio.h>
#include <stdlib.h>

/*****************************************************************************
*  This module contains code to parse the query expressions given by the
*  user
******************************************************************************/

/*****************************************************************************/

char *StripBlanks(as)

/*****************************************************************************/

char *as;

{
int i=0,j=0;
char *bs;

  bs = (char *) malloc((strlen(as)+1)*sizeof(char));
  while (as[i] != '\0') {
    if (as[i] != ' ') {
      bs[j++]=as[i];
    }
    i++;
  }
  free(as);
  bs[j]+='\0';
  return(bs);
}

/*****************************************************************************/

float **ParseRectangular(rectstring, dimension)

/*****************************************************************************/

char *rectstring;
int dimension;

{
  int stat,k;
  char brac[1],lparen[1],rparen[1],comma[1],rest[500];
  float **rectquery;

  rectquery =(float **)malloc(dimension *sizeof(float *));
  for (k=0; k < dimension; k++) {
    rectquery[k] = (float *)malloc(2*sizeof(float));
  }  /* allocate array that holds rectangular query */

  rectstring = StripBlanks(rectstring);
  sscanf(rectstring,"%c%s",brac,rest); /* get first bracket */
  rectstring=rest;
  if (strncmp(brac,"[",1) != 0) {
    DisplayError("Rectangular Query is invalid");
    return(NULL);
  }
  for (k=0; k <dimension;k++) {
    stat=sscanf(rectstring,"%c%f%c%f%c%s",lparen,&rectquery[k][0],comma,
                &rectquery[k][1],rparen,rest);
    rectstring=rest;
    if ((strncmp(lparen,"(",1) !=0) || (strncmp(comma,",",1)) !=0 ||
        (strncmp(rparen,")",1) != 0)) {
      DisplayError("Rectangular Query is invalid");
      return(NULL);
    }
    if (k < dimension -1 ) {
      sscanf(rectstring,"%c%s",comma,rest);
      rectstring=rest;
      if (strncmp(comma,",",1) != 0) {
        DisplayError("Rectangular Query is invalid");
        return(NULL);
      }
    }
  }
  sscanf(rectstring,"%c",brac); /* get last bracket */
  if (strncmp(brac,"]",1) !=0) {
    DisplayError("Rectangular Query is invalid");
    return(NULL);
  }
  return(rectquery);
}  


/*****************************************************************************/

float *ParseNearestNeighbor(NNstring, dimension)

/*****************************************************************************/

char *NNstring;
int dimension;

{
int j;
char lparen[1],rparen[1],comma[1],rest[500];
float *point;

point = (float *) malloc((dimension)*sizeof(float));


NNstring = StripBlanks(NNstring);

if (sscanf(NNstring,"(%s",rest) != 1) {
  return(NULL);
}

NNstring=rest;

for (j=0; j < dimension-1; j++) {
  if (sscanf(NNstring,"%f%s",&point[j],rest) != 2) {
    return(NULL);
  }
  NNstring = rest;  

  sscanf(NNstring,"%c%s",comma,rest);

  if (strncmp(comma,",",1) != 0) {
    return(NULL);
  }
  NNstring = rest;  
}

if (sscanf(NNstring,"%f,",&point[dimension-1]) !=1) {
    DisplayError("Nearest Neighbor Query is invalid   3\n");
    return(NULL);
  }
  for (j=0;j < dimension; j++) {
    printf("Dimension %d, %f\n",j,point[j]);
  }
  return(point);
}
