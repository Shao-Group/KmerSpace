#include "util.h"

void printIntArray(const int* x, const int len){
  int i;
  for(i=0; i<len; i+=1){
    printf("%2d ", x[i]);
  }
  printf("\n");
}

int editDist(const kmer s1, const kmer s2, const int k, const int max_d){
    int* row1 = malloc(sizeof *row1 *(k+1));
    int* row2 = malloc(sizeof *row2 *(k+1));
    int* tmptr;
    int i, j, tmp;
    for(i=0; i<k+1; i+=1){
	row1[i] = i;
    }

    kmer s1_copy, s2_copy;
    for(i=1, s1_copy=s1; i<k+1; i+=1, s1_copy>>=2){
	row2[0] = i;
	
	for(j=1, s2_copy=s2; j<k+1; j+=1, s2_copy>>=2){
	    //substitution
	    row2[j] = row1[j-1] +
		((s1_copy & 3) == (s2_copy & 3) ? 0 : 1);
	    //deletion
	    tmp = row1[j] + 1;
	    row2[j] = row2[j] > tmp ? tmp : row2[j];
	    //insertion
	    tmp = row2[j-1] +1;
	    row2[j] = row2[j] > tmp ? tmp : row2[j];	
	}

	if(max_d >= 0 && row2[i] >= max_d){
	    tmp = row2[i];
	    free(row1);
	    free(row2);
	    return tmp;
	}
	tmptr = row1;
	row1 = row2;
	row2 = tmptr;
    }

    tmp = row1[k];
    free(row1);
    free(row2);
    return tmp;
}

kmer encode(const char* str, const int k){
    kmer enc = 0;
    int i, x;
    for(i=0; i<k; i+=1){
	switch(str[i]){
	case 'A': x=0; break;
	case 'C': x=1; break;
	case 'G': x=2; break;
	case 'T': x=3; break;
	}
	enc = (enc << 2)|x;
    }
    return enc;
}

char* decode(const kmer enc, const int k, char* str){
    if(str == NULL){
	str = malloc(sizeof *str *k);
    }
    kmer enc_copy = enc;
    char base[] = {'A', 'C', 'G', 'T'};
    int i;
    for(i=k-1; i>=0; i-=1){
	str[i] = base[enc_copy & 3];
	enc_copy >>= 2;
    }
    return str;
}

kmer* readCentersFromFile(const char* filename, const int k, size_t* numOfCenters){
  FILE* fin = fopen(filename, "r");

  fscanf(fin, "%zu\n", numOfCenters);
  kmer* centers = malloc(sizeof *centers *(*numOfCenters));

  int buffersize = k+2;//for \n and \0
  char* kmer_str = malloc(sizeof *kmer_str *buffersize);

  int i;
  for(i=0; i<*numOfCenters; i+=1){
      fgets(kmer_str, buffersize, fin);
      centers[i] = encode(kmer_str, k);
  }

  free(kmer_str);
  
  fclose(fin);
  return centers;
}
