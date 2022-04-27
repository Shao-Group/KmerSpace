/*
  We try to identify a largest clique among all r-neighbors of a given kmer s,
  that is, a largest subset of r-neighbors of s with mutual distance at most r.
  The following ILP is used.

  For i=1...|N_d(s)|, binary variable x_i=1 iff the i-th neighbor n_i is in 
  the clique.
  Maximize: sum(x_i)
  Subject to: x_i+x_j<=1 if edit(n_i, n_j)>r.

  By: Ke@PSU
  Last edited: 04/26/2022
*/

#include "util.h"
#include "ArrayList.h"
#include "HashTable.h"
#include "gurobi_c.h"

//mask for (k-1)-mers
#define luMSB 0x8000000000000000lu

typedef struct{
    kmer x;
    int dist_to_s;
} Neighbor;

static inline void setNeighbor(Neighbor* n, kmer s, int d){
    n->x = s;
    n->dist_to_s = d;
}

/*
  Do bfs for r layers from the given kmer c, return a list of neightbors
*/
void bfsNeighborsRadius(kmer c, int k, int r, ArrayList* neighbors){
    HashTable visited;
    HTableInit(&visited);

    ArrayList cur_layer, next_layer;
    AListInit(&cur_layer);
    AListInit(&next_layer);

    HTableInsert(&visited, c);
    AListInsert(&cur_layer, (void*)c);

    size_t i, j;
    int depth = 1;
    kmer s, head, body, tail, x, m;
    Neighbor *entry;

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
			    //add x to the neighbor list
			    entry = malloc_harder(sizeof *entry);
			    setNeighbor(entry, x, depth);
			    AListInsert(neighbors, entry);
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
			    //add x to the neighbor list
			    entry = malloc_harder(sizeof *entry);
			    setNeighbor(entry, x, depth);
			    AListInsert(neighbors, entry);
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


int main(int argc, char *argv[]){
    if(argc != 4){
	printf("usage: maxClique_ILP k r s\n");
	return 1;
    }

    int k = atoi(argv[1]);
    int r = atoi(argv[2]);
    kmer s = encode(argv[3], k);

    size_t i, j;

    ArrayList neighbors;
    AListInit(&neighbors);

    bfsNeighborsRadius(s, k, r, &neighbors);
    
    GRBenv *env = NULL;
    GRBmodel *model = NULL;
    int error = 0;
    double *obj = NULL;
    char* vtype = NULL;

    /* Create environment */
    error = GRBemptyenv(&env);
    if (error) goto QUIT;

    char filename[100];
    sprintf(filename, "maxClique-%d-%d-%.*s.log", k, r, k, argv[3]);
    error = GRBsetstrparam(env, "LogFile", filename);
    if (error) goto QUIT;

    //number of threads to use
    error = GRBsetintparam(env, "Threads", 64);
    if (error) goto QUIT;
    //disable screen output
    error = GRBsetintparam(env, "LogToConsole", 0);
    if (error) goto QUIT;

    error = GRBstartenv(env);
    if (error) goto QUIT;

    /* Create an empty model */
    error = GRBnewmodel(env, &model, "maxClique", 0,
			NULL, NULL, NULL, NULL, NULL);
    if (error) goto QUIT;

    /* Add variables */
    obj = malloc_harder(sizeof *obj *(neighbors.used));
    vtype = malloc_harder(sizeof *vtype *(neighbors.used));
    for(i=0; i<(neighbors.used); i+=1){
	obj[i] = 1;
	vtype[i] = GRB_BINARY;
    }

    error = GRBaddvars(model, neighbors.used, 0, NULL, NULL, NULL,
		       obj, NULL, NULL, vtype, NULL);
    if (error) goto QUIT;

    /* Change objective sense to maximization */
    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
    if (error) goto QUIT;


    /* Add constraints */
    //for edit(n_i, n_j) > r, add x_i + x_j <= 1
    double val[2] = {1, 1};
    int ind[2];
    int edit;
    Neighbor *ni, *nj;
    
    for(i=0; i<neighbors.used; i+=1){
	ni = neighbors.arr[i];
	for(j=i+1; j<neighbors.used; j+=1){
	    nj = neighbors.arr[j];
	    //safe by triangle inequality
	    if(ni->dist_to_s + nj->dist_to_s <= r) continue;
	    edit = editDist(ni->x, nj->x, k, r+1);
	    if(edit > r){
		ind[0] = i;
		ind[1] = j;
		error = GRBaddconstr(model, 2, ind, val,
				     GRB_LESS_EQUAL, 1.0, NULL);
		if (error) goto QUIT;
	    }
	}
    }

    /* Write model to '.lp' */
    sprintf(filename, "maxClique-%d-%d-%.*s.lp", k, r, k, argv[3]);
    error = GRBwrite(model, filename);
    if (error) goto QUIT;
    
    /* Optimize model */
    error = GRBoptimize(model);
    if (error) goto QUIT;

    /* Capture solution information */
    int optimstatus;
    double objval;
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
    if (error) goto QUIT;

    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if (error) goto QUIT;

    //reuse obj array as sol
    double *sol = obj;
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, neighbors.used, sol);
    if (error) goto QUIT;

    printf("\nOptimization complete\n");

    sprintf(filename, "maxClique-%d-%d-%.*s.result", k, r, k, argv[3]);
    FILE* fout = fopen(filename, "w");
    if (optimstatus == GRB_OPTIMAL) {
	fprintf(fout, "Total neighbors %zu, Optimal objective: %lu\n", neighbors.used, (long unsigned)objval);
	for(i=0; i<neighbors.used; i+=1){
	    if(sol[i] > 0){
		ni = neighbors.arr[i];
		fprintf(fout, "%.*s\n", k,
			decode(ni->x, k, filename));
	    }
	}
    } else if (optimstatus == GRB_INF_OR_UNBD) {
	fprintf(fout, "Model is infeasible or unbounded\n");
    } else {
	fprintf(fout, "Optimization was stopped early\n");
    }
    fclose(fout);

QUIT:

    AListFree(&neighbors, free);
    
    free(obj);
    free(vtype);
    
    /* Error reporting */
    if (error) {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(1);
    }

    /* Free model */
    GRBfreemodel(model);

    /* Free environment */
    GRBfreeenv(env);

    return 0;
}
