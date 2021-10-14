/*
  Utility functions used by the k-mer partition project.
  By: Ke@PSU
  Last edited: 10/11/2021
*/

#ifndef _UTIL_H
#define _UTIL_H 1

#include <stdlib.h>
#include <stdio.h>

/*
  Leran@PSU provides the centers of the partitions. 
  Following her, each k-mer is represented by an long unsigned int 
  with the encoding A-00, C-01, G-10, T-11.
*/
typedef long unsigned kmer;

/*
  Calculate Levenshtein distance between two k-mers using Wagner-Fischer algorithm.
  If max_d is nonnegative, the calculation may stop earlier if a diagonal entry
  reaches max_d.
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
  Generate all neighbors of distance 1 from the given k-mer.
 */
void getNeighbors(const kmer enc, const int k);

#endif // util.h
