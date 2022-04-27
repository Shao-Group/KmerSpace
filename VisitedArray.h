#ifndef VISITEDARRAY_H
#define VISITEDARRAY_H

/*
 * A class for the visited array used in search algorithms. The array is 
 * chopped into subarrays to avoid failures of dynamic allocation.
 *
 * Author: Leran Ma (lkm5463@psu.edu)
 * Date:   9:47 PM, Tuesday, January 11, 2022
 */
class VisitedArray
{
private:
    char **aptrs;               // An array of pointers to subarrays
    unsigned long int size;     // The capacity of the whole array in elements
    unsigned long int num_subs; // Number of subarrays
    unsigned long int sub_size; // Size of a subarray in bytes

public:
    /*
     * Constructor
     *
     * s: total number of elements
     */
    VisitedArray( unsigned long int s );

    /*
     * Constructor
     *
     * s : total number of elements
     * ss: 2^ss is the size of a subarray in bytes
     */
    VisitedArray( unsigned long int s, int ss );

    /*
     * Destructor
     */
    ~VisitedArray();

    /*
     * Overload [] operator to return the visited value indexed by sub
     * 0 - unvisited
     * 1 - visited
     *
     * sub: The index of the element to be extract
     */
    bool operator[]( const unsigned long int &sub );
    
    /*
     * Set the element indexed by sub as visited
     *
     * sub: The index of the element
     */
    void setVisit( const unsigned long int sub );
};

#endif
