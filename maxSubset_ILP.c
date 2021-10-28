/*
  We try to identify a largest possible subset A of diameter d 
  in the k-mer space with the following ILP:
  For i=1...4^k, binary variable x_i=1 iff k-mer i is in A.
  Maximize: sum(x_i)
  Subject to: (x_i+x_j-1)edit(i,j)<=d.

  By: Ke@PSU
  Last edited: 10/28/2021
*/

#include "util.h"
//#include "ArrayList.h"
//#include "HashTable.h"
#include "gurobi_c.h"

int main(int argc, char *argv[]){
    if(argc != 3){
	printf("usage: maxSubset k d\n");
	return 1;
    }

    int k = atoi(argv[1]);
    int d = atoi(argv[2]);

    size_t NUM_KMERS = 1<<(k<<1);
    size_t i, j;
    
    GRBenv *env = NULL;
    GRBmodel *model = NULL;
    int error = 0;
    double *obj = NULL;
    char* vtype = NULL;

    /* Create environment */
    error = GRBemptyenv(&env);
    if (error) goto QUIT;

    char filename[50];
    sprintf(filename, "maxSubset-%d-%d.log", k, d);
    error = GRBsetstrparam(env, "LogFile", filename);
    if (error) goto QUIT;

    //number of threads to use
    error = GRBsetintparam(env, "Threads", 64);
    if (error) goto QUIT;

    error = GRBstartenv(env);
    if (error) goto QUIT;

    /* Create an empty model */
    error = GRBnewmodel(env, &model, "maxSubset", 0,
			NULL, NULL, NULL, NULL, NULL);
    if (error) goto QUIT;

    /* Add variables */
    obj = malloc_harder(sizeof *obj *NUM_KMERS);
    vtype = malloc_harder(sizeof *vtype *NUM_KMERS);
    for(i=0; i<NUM_KMERS; i+=1){
	obj[i] = 1;
	vtype[i] = GRB_BINARY;
    }

    error = GRBaddvars(model, NUM_KMERS, 0, NULL, NULL, NULL,
		       obj, NULL, NULL, vtype, NULL);
    if (error) goto QUIT;

    /* Change objective sense to maximization */
    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
    if (error) goto QUIT;


    /* Add constraints */
    //for edit(i,j) > d, add x_i + x_j <= 1
    double val[2] = {1, 1};
    int ind[2];
    int edit;
    
    for(i=0; i<NUM_KMERS; i+=1){
	for(j=i+1; j<NUM_KMERS; j+=1){
	    edit = editDist(i, j, k, d+1);
	    if(edit > d){
		ind[0] = i;
		ind[1] = j;
		error = GRBaddconstr(model, 2, ind, val,
				     GRB_LESS_EQUAL, 1.0, NULL);
		if (error) goto QUIT;
	    }
	}
    }

    /* Write model to '.lp' */
    sprintf(filename, "maxSubset-%d-%d.lp", k, d);
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
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, NUM_KMERS, sol);
    if (error) goto QUIT;

    printf("\nOptimization complete\n");

    sprintf(filename, "maxSubset-%d-%d.result", k, d);
    FILE* fout = fopen(filename, "w");
    if (optimstatus == GRB_OPTIMAL) {
	fprintf(fout, "Optimal objective: %lu\n", (long unsigned)objval);
	for(i=0; i<NUM_KMERS; i+=1){
	    if(sol[i] > 0){
		fprintf(fout, "%.*s\n", k, decode(i, k, filename));
	    }
	}
    } else if (optimstatus == GRB_INF_OR_UNBD) {
	fprintf(fout, "Model is infeasible or unbounded\n");
    } else {
	fprintf(fout, "Optimization was stopped early\n");
    }
    fclose(fout);

QUIT:

    free(obj);
    free(vtype);

    /* Free model */
    GRBfreemodel(model);

    /* Free environment */
    GRBfreeenv(env);
    
    /* Error reporting */
    if (error) {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(1);
    }

    return 0;
}
