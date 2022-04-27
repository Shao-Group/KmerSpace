#include <cstdlib>
#include <cstring>
#include "DistArray.h"

/*
 * Constructor
 *
 * s: total number of elements
 * d: maximum edit distance allowed
 */
DistArray::DistArray( unsigned long int s, unsigned int d )
{
    sub_size = 1ul << 20; // Default sub_size is 2^20 bytes.
    num_subs = s / 4 / sub_size; // Each element occupies 2 bits (1/4 bytes).
    if ( (s / 4) % sub_size != 0 )
    {
        num_subs++;
    }
    aptrs = (char **) calloc( num_subs, sizeof(char *) );
    for (int i = 0; i < num_subs; ++i)
    {
        aptrs[i] = (char *) calloc( sub_size, 1 );
    }
    size = s;

    values = (unsigned int *) malloc( 4 * sizeof(unsigned int) );
    values[0] = d + 1;
    values[1] = values[0] / 2;
    values[2] = values[1] + 1;
    values[3] = values[2] + 1;

    // Generate masks for the 4 bit offets
    // Each offset has four possible value
    masks = (unsigned int **) malloc( 4 * sizeof(unsigned int *) );
    for (int i = 0; i < 4; ++i)
    {
        masks[i] = (unsigned int *) malloc( 4 * sizeof(unsigned int) );
        for (unsigned int j = 0; j < 4; ++j)
        {
            masks[i][j] = j << (6 - 2*i);
        }
    }
}

/*
 * Constructor
 *
 * s : total number of elements
 * d : maximum edit distance allowed
 * ss: 2^ss is the size of a subarray in bytes
 */
DistArray::DistArray( unsigned long int s, unsigned int d, int ss )
{
    sub_size = 1ul << ss;
    num_subs = s / 4 / sub_size; // Each element occupies 2 bits (1/4 bytes).
    if ( (s / 4) % sub_size != 0 )
    {
        num_subs++;
    }
    aptrs = (char **) calloc( num_subs, sizeof(char *) );
    for (int i = 0; i < num_subs; ++i)
    {
        aptrs[i] = (char *) calloc( sub_size, 1 );
    }
    size = s;

    values = (unsigned int *) malloc( 4 * sizeof(unsigned int) );
    values[0] = d + 1;
    values[1] = values[0] / 2;
    values[2] = values[1] + 1;
    values[3] = values[2] + 1;

    // Generate masks for the 4 bit offets
    // Each offset has four possible value
    masks = (unsigned int **) malloc( 4 * sizeof(unsigned int *) );
    for (int i = 0; i < 4; ++i)
    {
        masks[i] = (unsigned int *) malloc( 4 * sizeof(unsigned int) );
        for (unsigned int j = 0; j < 4; ++j)
        {
            masks[i][j] = j << (6 - 2*i);
        }
    }
}

/*
 * Destructor
 */
DistArray::~DistArray()
{
    for (int i = 0; i < num_subs; ++i)
    {
        free( aptrs[i] );
    }
    free( aptrs );
    free( values );
    for (int i = 0; i < 4; ++i)
    {
        free( masks[i] );
    }
    free( masks );
}

/*
 * Overload [] operator to return the dist value indexed by sub
 *
 * sub: The index of the element to be extract
 */
unsigned int DistArray::operator[]( const unsigned long int &sub )
{
    // Find the bit offset
    int offset = sub % 4;

    // Find the byte in which the required element is located
    unsigned long int bytePos = sub / 4;

    unsigned int temp;
    memcpy( &temp, &aptrs[bytePos/sub_size][bytePos%sub_size], 1 );

    temp = temp & masks[offset][3];
    for (int i = 0; i < 4; ++i)
    {
        if ( temp == masks[offset][i] )
        {
            return values[i];
        }
    }
    return values[0];
}

/*
 * Set the dist value for an element indexed by sub
 *
 * sub : The index of the element
 * dist: The dist value of the element
 */
void DistArray::setDist( const unsigned long int sub, int dist )
{
    for (int i = 1; i < 4; ++i)
    {
        if ( dist <= values[i] )
        {
            dist = i;
            break;
        }
    }
    if ( dist >= values[0] )
    {
        dist = 0;
    }

    // Find the bit offset
    int offset = sub % 4;

    // Find the byte in which the required element is located
    unsigned long int bytePos = sub / 4;

    unsigned int tempByte;
    memcpy( &tempByte, &aptrs[bytePos/sub_size][bytePos%sub_size], 1 );
    tempByte = tempByte & (~masks[offset][3]);
    tempByte = tempByte | masks[offset][dist];
    memcpy( &aptrs[bytePos/sub_size][bytePos%sub_size], &tempByte, 1 );
}
