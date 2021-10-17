#include "ArrayList.h"

void AListInit(ArrayList* list){
    AListInitSize(list, 16);
}

void AListInitSize(ArrayList* list, size_t size){
    list->arr = malloc(sizeof *(list->arr) *size);
    list->size = size;
    list->used = 0;
}

void AListFree(ArrayList* list){
    free(list->arr);
    list->arr = NULL;
    list->used = 0;
    list->size = 0;
}

static inline void AListResize(ArrayList* list, size_t size){
    list->arr = realloc(list->arr, sizeof *(list->arr) *size);
    list->size = size;
}

void AListTrim(ArrayList* list){
    AListResize(list, list->used);
}

void AListInsert(ArrayList* list, long unsigned enc){
    if(list->used == list->size){
        AListResize(list, list->size*2);
    }
    list->arr[list->used++] = enc;
}

void AListClear(ArrayList* list){
    list->used = 0;
}

void AListSwap(ArrayList* la, ArrayList* lb){
    size_t tmp = la->size;
    la->size = lb->size;
    lb->size = tmp;

    tmp = la->used;
    la->used = lb->used;
    lb->used = tmp;

    long unsigned* arr_ptr = la->arr;
    la->arr = lb->arr;
    lb->arr = arr_ptr;
}
