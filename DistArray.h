/*
 * A class for the dist array used in BFS. The array is chopped into subarrays to avoid
 * failures of dynamic allocation.
 *
 * Author: Leran Ma (lkm5463@psu.edu)
 * Date  : 1:58 PM, Wednesday, November 3, 2021
 */
class DistArray
{
private:
    char **aptrs;               // An array of pointers to subarrays
    unsigned long int size;     // The capacity of the whole array in elements
    unsigned long int num_subs; // Number of subarrays
    unsigned long int sub_size; // Size of a subarray in bytes
    unsigned int *values;       // Possible dist values for an element
    unsigned int **masks;       // Masks used for easier caulculation

public:
    /*
     * Constructor
     *
     * s: total number of elements
     * d: maximum edit distance allowed
     */
    DistArray( unsigned long int s, unsigned int d );

    /*
     * Constructor
     *
     * s : total number of elements
     * d : maximum edit distance allowed
     * ss: 2^ss is the size of a subarray in bytes
     */
    DistArray( unsigned long int s, unsigned int d, int ss );

    /*
     * Destructor
     */
    ~DistArray();

    /*
     * Overload [] operator to return the dist value indexed by sub
     *
     * sub: The index of the element to be extract
     */
    unsigned int operator[]( const unsigned long int &sub );

    /*
     * Set the dist value for an element indexed by sub
     *
     * sub : The index of the element
     * dist: The dist value of the element
     */
    void setDist( const unsigned long int sub, int dist );
};
