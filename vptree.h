#define SAMPLE_SIZE 100


typedef struct VPTreeNode   /* Tree Record */
{
  int     Index;   /* Stores index of point in the points array */
  float   Mu;      /* Stores median distance */
  float   ll,lu,rl,ru; /* Stores left and right child's subbounds */
  struct  VPTreeNode   *Left, *Right; /*Pointers to sons of tree node */
} VPTREEREC, VPTreeNode, *VPTreePtr;

typedef struct IndexRec
{
  int Index; /* index to points array */
  struct IndexRec *Next;
} IndexRec;
