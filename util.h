/*
  Utility functions used by the k-mer partition project.
  By: Ke@PSU
  Last edited: 11/10/2021
*/

#ifndef _UTIL_H
#define _UTIL_H 1

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> //for sleep

/*
  Leran@PSU provides the centers of the partitions. 
  Following her, each k-mer is represented by a long unsigned int 
  with the encoding A-00, C-01, G-10, T-11.
*/
typedef long unsigned kmer;

/*
  Calculate Levenshtein distance between two strings using Wagner-Fischer algorithm.
  If max_d is nonnegative, the calculation may stop earlier if a diagonal entry
  reaches max_d.
*/
int editDist3(const char* s1, const int l1, const char* s2, const int l2, const int max_d);

/*
  Calculate Levenshtein distance between two x-mers using Wagner-Fischer algorithm.
  If max_d is nonnegative, the calculation may stop earlier if a diagonal entry
  reaches max_d.
*/
int editDist2(const kmer s1, const int k1, const kmer s2, const int k2, const int max_d);
/*
  Special case where |s1| = |s2|, kept for backwards compatibility
*/
int editDist(const kmer s1, const kmer s2, const int k, const int max_d);

/*
  Encode the string representation of a k-mer.
*/
kmer encode(const char* c, const int k);

/*
  Decode an k-mer into its string representation.
  If str is not null, it is used to store the resulting string; 
  otherwise a new char array is allocated.
*/
char* decode(const kmer enc, const int k, char* str);

/*
  Read in a list of centers (strings representing k-mers) from a file.
  The first line of the file is the number of centers.
  The following lines each contains (a string representation of) a k-mer.
  Use "make filename.MaxID" to convert Leran's output into this format.
*/
kmer* readCentersFromFile(const char* filename, const int k, size_t* numOfCenters);

/*
  Read in a list of centers (k-mers cliques) from a file.
  The first line of the file is the number of centers.
  The following lines each contains a list of (string representations of) 
  k-mers and (k-1)-mers.
  Use "make clique{k}{d}.MaxID" to convert Leran's output into this format.
  A (k-1)-mer x is stored as (x | km1Mask) to distinguish with k-mers.
  Return an int array for the sizes of the cliques.
*/
int* readCliquesFromFile(const char* filename, const int k, kmer km1Mask,
			  kmer*** centers, size_t* numOfCenters);

/*
  Read in a list of hash value for k-mers from a file and update
  the provided hash array.
  Each line of the file contains a k-mer and its hash value, separated by space.
  The output generated by partitionByLayers.c and 
  partitionByLayersCheckByNeighbors.c can be used.
*/
void readKMerHashFromFile(const char* filename, const int k, int* h);

/*
  If memory allocation failed, sleep for 30s and try again until success.
*/
void* malloc_harder(size_t size);
void* calloc_harder(size_t num, size_t size);
void* realloc_harder(void* ptr, size_t new_size);

/*
  Generate a random k-mer. User is responsible for seeding with srand().
*/
kmer randomKMer(int k);

/*
  Given a k-mer s, randomly generate a k-mer t with dist(s,t) = d. 
  User is responsible for seeding with srand().
*/
kmer randomEdit(kmer s, int k, int d);

#endif // util.h
