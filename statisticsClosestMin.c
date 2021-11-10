/*
  Test local sensitivity of the assigning closest smallest center heuristic.
  Input should be processed from hl{k}-{r}-{kd}.hash-w2 by the following two
  steps:
  1. remove -1 lines (unassigned k-mers)
  2. remove the third field of each line (dist to center)
  so that each line of the file contains a string representation of a k-mer
  and the index of the center assigned to it, separated by a space.

  The r parameter in the hash file is the largest possible distance from a
  k-mer and its assigned center. By triangle inequality, if two k-mers s and 
  t are assigned to the same center, dist(s, t) <= 2r.

  For d = 1, ..., 2r, a million pairs of k-mers (s, t) with dist(s, t) = d
  are randomly generated. The hash values h(s) and h(t) are then compared.
  We want to see what percentage of the pairs have the same hash value.
  Higher percentage indicates a better local sensitivity. For d=0, the value
  is 100%; for d>2r, the value is 0%.

  By: Ke@PSU
  Last edited: 11/10/2021
*/
#include "util.h"
#include <time.h>

#define TRIALS 1000000lu

int main(int argc, char* argv[]){
    if(argc != 4){
	printf("usage: statisticsClosestMin.out k r hash_file\n");
	return 1;
    }

    int k = atoi(argv[1]);
    int r = atoi(argv[2]);
    char* hash_file = argv[3];

    size_t NUM_KMERS = (1<<(k<<1));
    int* h = malloc_harder(sizeof *h *NUM_KMERS);
    readKMerHashFromFile(hash_file, k, h);
    
    srand(time(0));

    r <<= 1;
    if(r>k) r = k;
    int d;
    long unsigned match;
    long unsigned i;
    kmer s, t;
    for(d=1; d<=r; d+=1){
	match = 0lu;
	for(i=0lu; i<TRIALS; i+=1){
	    s = randomKMer(k);
	    t = randomEdit(s, k, d);
	    if(h[s] == h[t]) match += 1;
	}
	printf("d=%d %.2f\n", d, match*100.0/TRIALS);
    }
    
    free(h);
    return 0;
}
