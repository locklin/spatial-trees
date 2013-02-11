#define SPBUCKETSIZE 50

typedef struct sproullkdNode {
  int bucket;  /* true if the node is a bucket node */
  int discrim; /* discriminator of node */
  float cutval; 
  float *plane;     /* normal vector along principal axis */
  struct sproullkdNode *loson, *hison;
  int lopt,hipt;
} sproullkdNode;

