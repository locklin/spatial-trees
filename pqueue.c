/*************************************************************************************/
/*  Code to implement the abstract data type priority queue for use in j nearest     */
/*  neighbor searching.  Actual implementation is done using heaps.                  */
/*                                                                                   */
/*  Adapted from Sedgewick's: Algorithms in C p. 148-160.                            */
/*************************************************************************************/

/* 
   The heap data structure consists of two priority queues.  One for the j-smallest
   distances encountered, one to keep the indexes into the points array of the  
   points corresponding to the j-smallest distances. 
*/



/*************************************************************************************/

void PQupheap(DistArr,FoundArr,k)

/*************************************************************************************/

float *DistArr;  /* j-smallest distances encountered */

int *FoundArr,k;

{
float v;
int j;

v=DistArr[k]; DistArr[0] = 999999999999999.0;
j=FoundArr[k];

while(DistArr[k/2] <= v) {
  DistArr[k] = DistArr[k/2];
  FoundArr[k] = FoundArr[k/2];
  k=k/2;
}
DistArr[k] = v;
FoundArr[k] = j;
}

/*************************************************************************************/

void PQInsert(distance,index,DistArr,FoundArr)

/*************************************************************************************/

float distance,*DistArr;
int   index, *FoundArr;

{
  FoundArr[0]=FoundArr[0]+1;
  DistArr[FoundArr[0]] = distance;
  FoundArr[FoundArr[0]] = index;
  PQupheap(DistArr,FoundArr,FoundArr[0]);
}



/*************************************************************************************/

void PQdownheap(DistArr,FoundArr,k,index)

/*************************************************************************************/

float *DistArr;  /* j-smallest distances encountered */

int *FoundArr,k,index;

{

int j,N;
float v;

v=DistArr[k];

N = FoundArr[0];  /* tricky patch to maintain the data structure */
FoundArr[0]=index;

while (k <= N/2) {
  j=k+k;
  if (j < N && DistArr[j] <DistArr[j+1]) j++;
  if (v>=DistArr[j]) break;
  DistArr[k]=DistArr[j]; 
  FoundArr[k]=FoundArr[j];
  k=j;
}

DistArr[k] = v;
FoundArr[k]= index;
FoundArr[0]=N;  /* restore data struct */


}

/*************************************************************************************/

void PQreplace(distance,DistArr,FoundArr,index)

/*************************************************************************************/

float *DistArr,distance;
int *FoundArr;

{
  DistArr[0]=distance;
  PQdownheap(DistArr,FoundArr,0,index);
}
