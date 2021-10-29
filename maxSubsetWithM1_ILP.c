/*
  We try to identify a largest possible subset A of diameter d 
  in the k-mer + (k-1)-mer space with the following ILP:
  For i=1...4^k+4^{k-1}, binary variable x_i=1 iff k-mer i (i<=4^k)
  or (k-1)-mer i-4^k (i>4^k) is in A.
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
	printf("usage: maxSubsetWithM1_ILP k d\n");
	return 1;
    }

    int k = atoi(argv[1]);
    int d = atoi(argv[2]);

    size_t NUM_KMERS = 1<<(k<<1);
    size_t NUM_KM1MERS = NUM_KMERS >> 2;
    size_t total = NUM_KMERS + NUM_KM1MERS;
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
    sprintf(filename, "maxSubsetWithM1-%d-%d.log", k, d);
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
    error = GRBnewmodel(env, &model, "maxSubset", 0,
			NULL, NULL, NULL, NULL, NULL);
    if (error) goto QUIT;

    /* Add variables */
    obj = malloc_harder(sizeof *obj *total);
    vtype = malloc_harder(sizeof *vtype *total);
    for(i=0; i<total; i+=1){
	obj[i] = 1;
	vtype[i] = GRB_BINARY;
    }

    error = GRBaddvars(model, total, 0, NULL, NULL, NULL,
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

    kmer s1, s2;
    int l1, l2;
    
    for(i=0; i<total; i+=1){
	if(i>=NUM_KMERS){
	    s1 = i ^ NUM_KMERS;
	    l1 = k-1;
	}else{
	    s1 = i;
	    l1 = k;
	}
	for(j=i+1; j<total; j+=1){
	    if(j>=NUM_KMERS){
		s2 = j ^ NUM_KMERS;
		l2 = k-1;
	    }else{
		s2 = j;
		l2 = k;
	    }
	    edit = editDist2(s1, l1, s2, l2, d+1);
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
    sprintf(filename, "maxSubsetWithM1-%d-%d.lp", k, d);
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
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, total, sol);
    if (error) goto QUIT;

    printf("\nOptimization complete\n");

    sprintf(filename, "maxSubsetWithM1-%d-%d.result", k, d);
    FILE* fout = fopen(filename, "w");
    if (optimstatus == GRB_OPTIMAL) {
	fprintf(fout, "Optimal objective: %lu\n", (long unsigned)objval);
	for(i=0; i<NUM_KMERS; i+=1){
	    if(sol[i] > 0){
		fprintf(fout, "%.*s\n", k, decode(i, k, filename));
	    }
	}
	double* solM1 = sol + NUM_KMERS;
	for(i=0; i<NUM_KM1MERS; i+=1){
	    if(solM1[i] > 0){
		fprintf(fout, "%.*s\n", k-1, decode(i, k-1, filename));
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
