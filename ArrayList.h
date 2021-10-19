/*
  A dynamic array based list for k-mers (long unsigned int).
  By Ke@PSU
  Last modified: 10/11/2021
*/

#ifndef _ARRAYLIST_H
#define _ARRAYLIST_H 1

#include "util.h"
//#include <stdlib.h>

typedef struct {
    long unsigned* arr;
    size_t size;
    size_t used;
} ArrayList;

/*
  Initialize an array to size 16.
*/
void AListInit(ArrayList* list);

/*
  Initialize an array to the given size.
*/
void AListInitSize(ArrayList* list, size_t size);

/*
  Free the entire array.
*/
void AListFree(ArrayList* list);

/*
  Trim the list to the minimum size needed.
*/
void AListTrim(ArrayList* list);

/*
  Insert a new k-mer into the array.
  If already full, the size of the array will be doubled before the insertion.
*/
void AListInsert(ArrayList* list, long unsigned enc);

/*
  Clear used count so the array can be used as if a new one.
*/
void AListClear(ArrayList* list);

/*
  Swap two lists.
*/
void AListSwap(ArrayList* la, ArrayList* lb);

#endif // ArrayList.h
