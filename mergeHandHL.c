/*
  Merge the island based PLSHashing (h{k}-{p}-{q}-{kd}.hash[-v*]) 
  with the closest min center hashing (hl{k}-{r}-{kd}.hash-w2c).
  When r=d, the CMCHashing guarantees a value for each k-mer.
  If both assigned values to a k-mer, they should agree (remove -1 entries
  from PLSHashing first). 
  
  So this is essentially using PLSHashing to mark in the CMCHashing 
  which ones are within a proper island as defined by the PLSHashing.
  1 - this kmer is in an island
  0 - otherwise
 
  By: Ke@PSU
  Last edited: 11/21/2021
*/

#include "util.h"
#include <string.h>


int main(int argc, char* argv[]){
    if(argc != 4){
	printf("usage: mergeHandHL.out k PLSHashing.hash[-v*] CMCHashing.hash-w2c\n");
	return 1;
    }

    int k = atoi(argv[1]);
    char* pls_file = argv[2];
    char* cmc_file = argv[3];

    FILE* pls_in = fopen(pls_file, "r");
    FILE* cmc_in = fopen(cmc_file, "r");

    char merge_file[20];
    sprintf(merge_file, "merge%d.hash", k);
    FILE* fout = fopen(merge_file, "w");


    char format[20];
    sprintf(format, "%%%ds %%d\n", k);
    
    size_t N = 1<<(k<<1);
    size_t i;
    char cmc_cur[k+1];
    int cmc_hash;
    char pls_cur[k+1];
    int pls_hash;
    int stat;
    int matched=1;
    for(i=0; i<N; i+=1){
	stat = fscanf(cmc_in, format, &cmc_cur, &cmc_hash);
	if(matched){//if not match, keep scanning cmc until a match is found
	    stat = fscanf(pls_in, format, &pls_cur, &pls_hash);
	    if(stat == EOF){//no more hash in pls
		fprintf(fout, "%*s %d 0\n", k, cmc_cur, cmc_hash);
		i += 1;
		break;
	    }
	    matched = 0;
	}
	
	if(strcmp(cmc_cur, pls_cur) == 0){
	    matched = 1;
	    fprintf(fout, "%*s %d 1\n", k, cmc_cur, cmc_hash);
	}else{
	    fprintf(fout, "%*s %d 0\n", k, cmc_cur, cmc_hash);
	}
    }

    for(; i<N; i+=1){
	stat = fscanf(cmc_in, format, &cmc_cur, &cmc_hash);
	fprintf(fout, "%*s %d 0\n", k, cmc_cur, cmc_hash);
    }

    fclose(cmc_in);
    fclose(pls_in);
    fclose(fout);
    return 0;
}
