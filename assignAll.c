/*
  Given a list of centers C={c_1, c_2, ..., c_m} in the k-mer space S, and 
  a fixed radius r, 
  
  ********Deprecated****************************************************
  we assign to each k-mer s a list of weighted centers 
  L_s={(c_{si}, w_{si})} such that for each i, dist(s, c_{si}) <= r and
  sum(w_{si) = 1.
  **********************************************************************

  we assign to each k-mer s a list of centers L_s={c_{si}} within radius
  r of s. L_s is sorted in ascending order of dist(s, c_{si}) and c_{si}
  (for centers with the same distance).

  If C is a maximal independent set of the k-mer space with parameter d 
  (i.e., dist(c_i, c_j) > d and for each k-mer s in S, there exists a 
  center c_i such that dist(c_i, s) <= d), and r >= d, then each k-mer s 
  is guaranteed to be assigned a non-empty list.
 
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

static inline void setNeighborCenter(NeighborCenter* nc, size_t idx, int d){
    nc->center_idx = idx;
    nc->dist = d;
}

/*
  Do bfs for r layers from the center c, add c to the center list of each
  visited kmer.
*/
void bfsNeighborsRadius(kmer c, int c_idx, int k, int r, ArrayList* h){
    HashTable visited;
    HTableInit(&visited);

    ArrayList cur_layer, next_layer;
    AListInit(&cur_layer);
    AListInit(&next_layer);

    HTableInsert(&visited, c);
    AListInsert(&cur_layer, (void*)c);

    //assign hash list for the center itself
    NeighborCenter *entry = malloc_harder(sizeof *entry);
    setNeighborCenter(entry, c_idx, 0);
    AListInsert(h+c, entry);

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
			    //add c to the hash list for x
			    entry = malloc_harder(sizeof *entry);
			    setNeighborCenter(entry, c_idx, depth);
			    AListInsert(h+x, entry);
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
			    //add c to the hash list for x
			    entry = malloc_harder(sizeof *entry);
			    setNeighborCenter(entry, c_idx, depth);
			    AListInsert(h+x, entry);
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

int cmpNeighborCenterByDistAsc(const void* c1, const void* c2){
    NeighborCenter* a = *(NeighborCenter**) c1;
    NeighborCenter* b = *(NeighborCenter**) c2;
    int dist = a->dist - b->dist;
    if(dist != 0) return dist;
    else return (a->center_idx > b->center_idx ? 1 : -1);
}

//this version does not work well for Jaccard Similarity, will output raw dist instead
//weight can then be calculated as needed
/*
void fprintNeighborCenters(FILE* fout, ArrayList* list){
    int ct = list->used;
    if(ct == 0){
	fprintf(fout, "\n");
	return;
    }
    if(ct == 1){
	fprintf(fout, " %zu %.05f\n",
		((NeighborCenter*)list->arr[0])->center_idx, 1.0);
	return;
    }

    qsort(list->arr, ct, sizeof(NeighborCenter*), cmpNeighborCenterByDistAsc);
    
    double normalizer = 0.0;
    double weights[ct];
    size_t i;
    for(i=0; i<ct; i+=1){
	weights[i] = 1.0 / ((NeighborCenter*)list->arr[i])->dist;
        normalizer += weights[i];
    }
    normalizer = 1.0 / normalizer;

    double sum_weights = 0.0;
    for(i=0; i<ct; i+=1){
	weights[i] *= normalizer;
	weights[i] = ((int)(weights[i] * 100000 + .5)) / 100000.0;
	sum_weights += weights[i];
    }

    //slightly favor the first closest center
    //or penalize the last furthest center
    //if there is any rounding error
    weights[sum_weights < 1.0 ? 0 : ct-1] += 1.0 - sum_weights;
    
    for(i=0; i<ct; i+=1){
	fprintf(fout, " %zu %.05f",
		((NeighborCenter*)list->arr[i])->center_idx,
		weights[i]);
    }
    fprintf(fout, "\n");

}
*/
void fprintNeighborCenters(FILE* fout, ArrayList* list){
    int ct = list->used;
    if(ct == 0){
	fprintf(fout, "\n");
	return;
    }

    qsort(list->arr, ct, sizeof(NeighborCenter*), cmpNeighborCenterByDistAsc);

    int i;
    NeighborCenter* cur;
    for(i=0; i<ct; i+=1){
	cur = (NeighborCenter*)list->arr[i];
	fprintf(fout, " %zu %d",
	        cur->center_idx, cur->dist);
    }
    fprintf(fout, "\n");

}


int main(int argc, char* argv[]){
    if(argc != 4){
	printf("usage: assignAll.out k r centers_file\n");
	return 1;
    }

    int k = atoi(argv[1]);
    int r = atoi(argv[2]);
    char* centers_file = argv[3];

    size_t NUM_CENTERS;
    kmer* centers = readCentersFromFile(centers_file, k, &NUM_CENTERS);

    size_t NUM_KMERS = (1<<(k<<1));
    ArrayList* h = malloc_harder(sizeof *h *NUM_KMERS);

    size_t i;
    for(i=0; i<NUM_KMERS; i+=1){
	AListInitSize(h+i, r); //xxx: Is r a good choice for init size?
    }

    for(i=0; i<NUM_CENTERS; i+=1){
	bfsNeighborsRadius(centers[i], i, k, r, h);
    }
    free(centers);
    
    char output_filename[50];
    sprintf(output_filename, "hl%d-%d-%s.hash", k, r, centers_file);
    FILE* fout = fopen(output_filename, "w");

    for(i=0; i<NUM_KMERS; i+=1){
	fprintf(fout, "%.*s %zu", k, decode(i, k, output_filename), h[i].used);
	fprintNeighborCenters(fout, h+i);
	AListFree(h+i, free);
    }

    fclose(fout);
    free(h);
    return 0;
}
