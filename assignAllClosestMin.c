/*
  **************************************************************************
  Different from v1 (assignAll.c), we only keep the closest center with the
  smallest index for each kmer. v1 produces more information but requires
  too much memory for large k and r.
  **************************************************************************
  Given a list of centers C={c_1, c_2, ..., c_m} in the k-mer space S, and 
  a fixed radius r, we assign to each k-mer s the center c_{s} within 
  radius r that has the smallest dist(s, c_i) and smallest index (if tie in 
  distance) among all centers.

  If C is a maximal independent set of the k-mer space with parameter d 
  (i.e., dist(c_i, c_j) > d and for each k-mer s in S, there exists a 
  center c_i such that dist(c_i, s) <= d), and r >= d, then each k-mer s 
  is guaranteed to be assigned a center.
 
  By: Ke@PSU
  Last edited: 11/09/2021
*/

#include "util.h"
#include "ArrayList.h"
#include "HashTable.h"
//#include <time.h>
//#include <unistd.h>

typedef int bool;
#define TRUE 1
#define FALSE 0

//mask for (k-1)-mers
#define luMSB 0x8000000000000000lu

typedef struct{
    size_t center_idx;
    int dist;
}NeighborCenter;

static inline void assignNeighborCenter(NeighborCenter* nc, size_t idx, int d){
    // not assiged or this is a closer center
    if((nc->center_idx & luMSB) || nc->dist > d){
	nc->center_idx = idx;
	nc->dist = d;
    }
}

/*
  Do bfs for r layers from the center c, for each visited kmer s, if s is
  unassigned or dist(s, c_s) > dist(s, c), then assign c to s (if dist is
  the same, c must have a higher index then c_s).
*/
void bfsNeighborsRadius(kmer c, int c_idx, int k, int r, NeighborCenter** h){
    HashTable visited;
    HTableInit(&visited);

    ArrayList cur_layer, next_layer;
    AListInit(&cur_layer);
    AListInit(&next_layer);

    HTableInsert(&visited, c);
    AListInsert(&cur_layer, (void*)c);

    //assign the center itself
    assignNeighborCenter(h[c], c_idx, 0);

    size_t i, j;
    int depth = 1;
    kmer s, head, body, tail, x, m;

    for(depth=1; depth<=r; depth+=1){
	for(i=0; i<cur_layer.used; i+=1){
	    s = (kmer) cur_layer.arr[i];
	    //(k-1)-mer, no need to ^luMSB as the head will shift MSB out
	    if(s>=luMSB){
		//insertion
		for(j=0; j<k; j+=1){
		    head = (s>>(j<<1))<<((j+1)<<1);
		    tail = ((1lu<<(j<<1))-1) & s;
		    for(m=0; m<4; m+=1){
			body = m<<(j<<1);
			x = head|body|tail;
			//x is a k-mer
			//add to next_layer if not visited
			if(!HTableSearch(&visited, x)){
			    AListInsert(&next_layer, (void*) x);
			    HTableInsert(&visited, x);
			    //try to assign c to x
			    assignNeighborCenter(h[x], c_idx, depth);
			}
		    }
		}
	    }
	    //k-mer
	    else{
		//deletion
		for(j=0; j<k; j+=1){
		    head = (s>>((j+1)<<1))<<(j<<1);
		    tail = ((1lu<<(j<<1))-1) & s;
		    x = head|tail|luMSB;
		    //x is a (k-1)-mer
		    //add to next_layer if not visited
		    if(!HTableSearch(&visited, x)){
			AListInsert(&next_layer, (void*)x);
			HTableInsert(&visited, x);
		    }
		    
		}
		//substitution
		for(j=1; j<=k; j+=1){
		    head = (s>>(j<<1))<<(j<<1);
		    tail = ((1lu<<((j-1)<<1))-1) & s;
		    for(m=0; m<4; m+=1){
			body = m<<((j-1)<<1);
			x = head|body|tail;
			//x is a k-mer
			//add to next_layer if not visited
			if(!HTableSearch(&visited, x)){
			    AListInsert(&next_layer, (void*)x);
			    HTableInsert(&visited, x);
			    //try to assign c to x
			    assignNeighborCenter(h[x], c_idx, depth);
			}
		    }
		}
	    }//end k-mer
	}//end for each in cur_layer

	AListClear(&cur_layer, NULL);
	AListSwap(&cur_layer, &next_layer);
    }//end for depth from 1 to r

    HTableFree(&visited);
    AListFree(&cur_layer, NULL);
    AListFree(&next_layer, NULL);
    
}//end bfsNeighborsRadius


int main(int argc, char* argv[]){
    if(argc != 4){
	printf("usage: assignAllClosestMin.out k r centers_file\n");
	return 1;
    }

    int k = atoi(argv[1]);
    int r = atoi(argv[2]);
    char* centers_file = argv[3];

    size_t NUM_CENTERS;
    kmer* centers = readCentersFromFile(centers_file, k, &NUM_CENTERS);

    size_t NUM_KMERS = (1<<(k<<1));
    NeighborCenter** h = malloc_harder(sizeof *h *NUM_KMERS);

    size_t i;
    for(i=0; i<NUM_KMERS; i+=1){
	h[i] = malloc_harder(sizeof *h[i]);
	h[i] -> center_idx = luMSB;
    }

    for(i=0; i<NUM_CENTERS; i+=1){
	bfsNeighborsRadius(centers[i], i, k, r, h);
    }
    free(centers);
    
    char output_filename[50];
    sprintf(output_filename, "hl%d-%d-%.*s.hash-w2", k, r, 4, centers_file);
    FILE* fout = fopen(output_filename, "w");

    for(i=0; i<NUM_KMERS; i+=1){
	fprintf(fout, "%.*s ", k, decode(i, k, output_filename));
	if(h[i]->center_idx & luMSB){
	    fprintf(fout, "%d\n", -1);
	}else{
	    fprintf(fout, "%zu %d\n", h[i]->center_idx, h[i]->dist);
	}
	free(h[i]);
    }

    fclose(fout);
    free(h);
    return 0;
}
