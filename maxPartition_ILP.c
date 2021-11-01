/*
  A partition of the k-mer space S is a collection of subsets C={A_1, A_2, ..., 
  A_n, U}, where union(A_1, ..., A_n, U) = S and intersection(X,Y)=empty for
  all X!=Y in C. We say k-mers in A_i are assigned to the i-th island; k-mers
  in U are unassigned.

  A partition is called locality sensitive if it satisfies the following
  two conditions for the given parameters p and q:
  1. the diameter of each A_i is at most q, i.e., edit(x, y) <= q for all 
  x, y in A_i;
  2. the A_i's are at least p away from each other, i.e., edit(x, y) >= p
  for all x in A_i, y in A_j, i!=j.
  
  We try to identify a maximum size locality sensitive partition of the 
  k-mer space with the following ILP:
  Let N = 4^k. For r=1...N, binary variable x_r=1 iff k-mer r is assigned to
  one of the A_i's. For 1<=r<s<=N, binary variable x_{r,s}=1 iff k-mers r 
  and s are assigned to the same island.

  Maximize: sum(x_r)
  Subject to: 
  1. x_r + x_s -2 x_{r,s}>=0, for 1<=r<s<=N
  2. x_{r,s}+x_{s,t}-x_{r,t}<=1, for 1<=r<s<t<=N
  3. x_{r,s}-x_{s,t}+x_{r,t}<=1, for 1<=r<s<t<=N
  4. -x_{r,s}+x_{s,t}+x_{r,t}<=1, for 1<=r<s<t<=N
  5. edit(r, s)x_{r,s} <= q, for 1<=r<s<=N
  6. edit(r, s)(p(x_{r,s}-x_r-x_s+2)+1)>=p, for 1<=r<s<=N

  By: Ke@PSU
  Last edited: 10/29/2021
*/

#include "util.h"
//#include "ArrayList.h"
//#include "HashTable.h"
#include "gurobi_c.h"

int main(int argc, char *argv[]){
    if(argc != 4){
	printf("usage: maxPartition_ILP k p q\n");
	return 1;
    }

    int k = atoi(argv[1]);
    int p = atoi(argv[2]);
    int q = atoi(argv[3]);
    
    size_t NUM_KMERS = 1<<(k<<1);
    size_t r, s, t;
    
    GRBenv *env = NULL;
    GRBmodel *model = NULL;
    int error = 0;
    double *obj = NULL;
    char* vtype = NULL;
    double *sol = NULL;

    /* Create environment */
    error = GRBemptyenv(&env);
    if (error) goto QUIT;

    char filename[50];
    sprintf(filename, "maxPartition-%d-%d-%d.log", k, p, q);
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

    /* Add variables x_r */
    obj = malloc_harder(sizeof *obj *NUM_KMERS);
    vtype = malloc_harder(sizeof *vtype *NUM_KMERS);
    for(r=0; r<NUM_KMERS; r+=1){
	obj[r] = 1;
	vtype[r] = GRB_BINARY;
    }

    error = GRBaddvars(model, NUM_KMERS, 0, NULL, NULL, NULL,
		       obj, NULL, NULL, vtype, NULL);
    if (error) goto QUIT;

    free(obj);
    obj = NULL;
    free(vtype);
    vtype = NULL;
    
    /* Add variables x_{r,s} */
    size_t NUM_PAIRS = (1<<((k<<2)-1))-(1<<((k<<1)-1));
    vtype = malloc_harder(sizeof *vtype * NUM_PAIRS);
    for(r=0; r<NUM_PAIRS; r+=1){
	vtype[r] = GRB_BINARY;
    }

    error = GRBaddvars(model, NUM_PAIRS, 0, NULL, NULL, NULL,
		       NULL, NULL, NULL, vtype, NULL);

    if (error) goto QUIT;
    
    free(vtype);
    vtype = NULL;

    /* Change objective sense to maximization */
    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
    if (error) goto QUIT;


    /* Add constraints */
    double val[3];
    int ind[3];
    int edit;
    size_t i1, i2;
    
    for(r=0; r<NUM_KMERS; r+=1){
	i1 = NUM_KMERS*(r+1)-((r*(r+1))>>1); //index of x_{r, r+1}
	for(s=r+1; s<NUM_KMERS; s+=1){
	    //constraint 1: x_r + x_s - 2x_{r,s} >= 0
	    ind[0] = r;
	    ind[1] = s;
	    ind[2] = i1+s-r-1;
	    val[0] = 1;
	    val[1] = 1;
	    val[2] = -2;
	    error = GRBaddconstr(model, 3, ind, val,
				 GRB_GREATER_EQUAL, 0.0, NULL);
	    if (error) goto QUIT;

	    edit = editDist(r, s, k, -1);
	    if(edit < p){
		//constraint 6: if edit(r,s) < p, x_r + x_s - x_{r,s}<=1
		//reuse from above
		//ind[0] = r;
		//ind[1] = s;
		//ind[2] = i1+s-r-1;
		//val[0] = 1;
		//val[1] = 1;
		val[2] = -1;
		error = GRBaddconstr(model, 3, ind, val,
				     GRB_LESS_EQUAL, 1.0, NULL);
		if (error) goto QUIT;
	    }
	    if(edit > q){
		//constraint 5: if edit(r,s) > q, x_{r,s} = 0
		ind[0] = i1+s-r-1;
		val[0] = 1;
		error = GRBaddconstr(model, 1, ind, val,
				     GRB_EQUAL, 0.0, NULL);
		if (error) goto QUIT;
	    }
	    
	    
	    i2 = NUM_KMERS*(s+1)-((s*(s+1))>>1); //index of x_{s, s+1}
	    
	    for(t=s+1; t<NUM_KMERS; t+=1){
		//constraints 2-4
		ind[0] = i1+s-r-1; //x_{r,s}
		ind[1] = i2+t-s-1; //x_{s,t}
		ind[2] = i1+t-r-1; //x_{r,t}
		val[0] = 1;
		val[1] = 1;
		val[2] = -1;
		error = GRBaddconstr(model, 3, ind, val,
				     GRB_LESS_EQUAL, 1.0, NULL);
		if (error) goto QUIT;

		val[1] = -1;
		val[2] = 1;
		error = GRBaddconstr(model, 3, ind, val,
				     GRB_LESS_EQUAL, 1.0, NULL);
		if (error) goto QUIT;

		val[0] = -1;
		val[1] = 1;
		error = GRBaddconstr(model, 3, ind, val,
				     GRB_LESS_EQUAL, 1.0, NULL);
		if (error) goto QUIT;
	    }
	}
    }

    /* Write model to '.lp' */
    sprintf(filename, "maxPartition-%d-%d-%d.lp", k, p, q);
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

    t = NUM_KMERS+NUM_PAIRS;
    sol = malloc_harder(sizeof *sol *t);
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X,
			       0, t, sol);
    if (error) goto QUIT;

    printf("\nOptimization complete\n");

    sprintf(filename, "maxPartition-%d-%d-%d.result", k, p, q);
    FILE* fout = fopen(filename, "w");
    if (optimstatus == GRB_OPTIMAL) {
	fprintf(fout, "Optimal objective: %lu\n", (long unsigned)objval);
	size_t ct = 0;
	for(r=0; r<NUM_KMERS; r+=1){
	    if(sol[r] > 0){
                // r is assigned, this is a new island
		fprintf(fout, "island %zu\n%.*s\n",
			ct, k, decode(r, k, filename)); 
		i1 = NUM_KMERS*(r+1)-((r*(r+1))>>1); //index of x_{r, r+1}
		for(s=r+1; s<NUM_KMERS; s+=1, i1+=1){
		    if(sol[i1] > 0){ //r, s in the same island
			sol[s] = 0; //won't output s again
			fprintf(fout, "%.*s\n", k, decode(s, k, filename));
		    }
		}
		ct += 1;
	    }
	}
    } else if (optimstatus == GRB_INF_OR_UNBD) {
	fprintf(fout, "Model is infeasible or unbounded\n");
    } else {
	fprintf(fout, "Optimization was stopped early\n");
    }
    fclose(fout);

QUIT:

    if(obj) free(obj);
    if(vtype) free(vtype);
    if(sol) free(sol);

    /* Error reporting */
    if (error) {
	printf("ERROR: %d - %s\n", error, GRBgeterrormsg(env));
	exit(1);
    }

    /* Free model */
    GRBfreemodel(model);

    /* Free environment */
    GRBfreeenv(env);
    

    return 0;
}
