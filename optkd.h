#define BUCKETSIZE 50

typedef struct optkdNode {
  int bucket;  /* true if the node is a bucket node */
  int discrim; /* discriminator of node */
  float cutval; 
  struct optkdnode *loson, *hison;
  int lopt,hipt;
} optkdNode;

