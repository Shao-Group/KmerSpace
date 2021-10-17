/*
 * This program addresses the following problem from Prof. Mingfu Shao (mxs2589@psu.edu).
 *
 * --------------------------------------------------------------------------------------
 * Let K be the set of all kmers over alphabet {A, C, G, T}. |K| = 4^k. We build a graph
 * G = (K, E) with K be the set of vertices. Let x and y be two kmers in K. We add an
 * edge to E iff the edit distance between x and y is at most d. We want to calculate the
 * an maximal independent set (MIS) of G that covers all vertices.
 *
 * The input for this problem are two integers k and d; the output is an independent set
 * of G constructed above.
 *
 * We can try a greedy algorithm. We first generate K. We then iteratively pick a kmer x
 * from K, and remove x and all other kmers whose edit distance with x is at most d; we
 * repeat this procedure until K becomes empty. In order to save time, we do not do
 * simple pairwise comparison to find the edit distance. Instead, we build a graph
 * containing all possible k-mers, (k-1)-mers, and (k+1)-mers as vertices. If the edit
 * distance between any two vertices is 1, we add an edge of length 1 between the two
 * vertices. Then, we add the picked vertex to the MIS, explore all vertices at a
 * distance less than or equal to d from it using BFS, and remove them from K. We would
 * not do redundant exploration if we detected that a vertex was already explored via a
 * shorter path from another vertex in the MIS.
 * --------------------------------------------------------------------------------------
 *
 * In the program, the author uses a binary encoding {00, 01, 10, 11} for the alphabet
 * {A, C, G, T} in order to save space. Thus, any k-mer that is not longer than 32
 * characters can be represented by a 64-bit binary number stored in an 8-byte slot (an
 * unsigned long int). Besides, the author uses the last 2 bits of the 8-byte slot to
 * represent the length of the k-mer. 00 represents k - 1; 01 represents k; 11 represents
 * k + 1.
 *
 * Author: Leran Ma (lkm5463@psu.edu)
 * Date:   7:15 PM, Sunday, October 3, 2021
 */

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <unistd.h>
#include <ios>
#include <fstream>
#include <string>

using namespace std;

/*
 * A class for the dist array used in BFS. The array is chopped into
 * subarrays to avoid failures of dynamic allocation.
 */
class DistArray
{
private:
    char **aptrs;               // An array of pointers to subarrays
    unsigned long int size;     // The capacity of the whole array in elements
    unsigned long int num_subs; // Number of subarrays
    unsigned long int sub_size; // Size of a subarray in bytes
    unsigned int max_d;         // Maximum edit distance allowed

public:
    /*
     * Constructor
     *
     * s: total number of elements
     */
    DistArray( unsigned long int s, unsigned int d )
    {
        sub_size = 1;
        sub_size = sub_size << 30;
        num_subs = s / 4 / sub_size; // Each element occupies 1/4 bytes.
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
        max_d = d;
    }

    /*
     * Destructor
     */
    ~DistArray()
    {
        for (int i = 0; i < num_subs; ++i)
        {
            free( aptrs[i] );
        }
        free( aptrs );
    }

    /*
     * Overload [] operator to return the dist value indexed by sub
     *
     * sub: The index of the element to be extract
     */
    unsigned int operator[]( const unsigned long int &sub )
    {
        // Find the bit offset
        int offset = (sub % 4) * 2;

        // Find the byte in which the required element is located
        unsigned long int bytePos = sub / 4;

        unsigned int dist;
        memcpy( &dist, &aptrs[bytePos/sub_size][bytePos%sub_size], 1 );

        dist = (dist << (24 + offset)) >> 30;
        return dist + (max_d + 1)/2 - 1;
    }
    
    /*
     * Set the dist value for an element indexed by sub
     *
     * sub : The index of the element
     * dist: The dist value of the element
     */
    void setDist( const unsigned long int sub, int dist )
    {
        dist = dist + 1 - (max_d + 1)/2;
        if (dist < 2)
        {
            dist = 1;
        }

        // Find the byte in which the required element is located
        unsigned long int bytePos = sub / 4;

        // Find the bit offset
        int offset = (sub % 4) * 2;

        unsigned int tempByte;
        memcpy( &tempByte, &aptrs[bytePos/sub_size][bytePos%sub_size], 1 );
        unsigned int head = tempByte >> (8 - offset);
        head = head << (8 - offset);
        unsigned int tail = tempByte << 24 << (offset + 2);
        tail = tail >> 24 >> (offset + 2);
        tempByte = tempByte & (head | tail);
        unsigned int body = dist;
        body = body << (6 - offset);
        tempByte = tempByte | body;
        memcpy( &aptrs[bytePos/sub_size][bytePos%sub_size], &tempByte, 1 );
    }
};

/*
 * Prints a k-mer given its binary encoding
 *
 * enc: The binary encoding of the k-mer
 * k  : The length of the k-mer
 */
void printKmer( unsigned long int enc, int k )
{
    char base[4] = {'A', 'C', 'G', 'T'};
    char kmer[k + 1];
    kmer[k] = '\0';
    for (int i = k - 1; i >= 0; --i)
    {
        enc = enc >> 2;
        kmer[i] = base[enc & 3];
    }
    cerr << kmer;
}

/*
 * Gets neighbors of a vertex
 *
 * enc: The binary encoding of the k-mer
 * k  : The length of the k-mer
 * n  : A unordered_set to hold the neighbors
 */
void getNeighbor( unsigned long int enc, int k, unordered_set<unsigned long int> &n )
{
    unsigned long int i = enc >> 2;

    // Populate the adjacency list for kmers
    if ( (enc & 3) == 1 )
    {
        // Handle deletion
        for (int j = 1; j <= k; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * (j - 1));
            unsigned long int tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            n.emplace( (head + tail) << 2 );
        }

        // Handle insertion
        for (int j = 0; j <= k; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * (j + 1));
            unsigned long int tail = (i << 1 << (63 - 2 * j)) >> (63 - 2 * j) >> 1;
            for (unsigned long int l = 0; l < 4; ++l)
            {
                unsigned long int body = l << (2 * j);
                unsigned long int node = head + body + tail;
                n.emplace( (node << 2) | 2 );
            }
        }

        // Handle substitution
        for (int j = 1; j <= k; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * j);
            unsigned long int tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            for (unsigned long int l = 0; l < 4; ++l)
            {
                unsigned long int body = l << (2 * (j - 1));
                unsigned long int node = head + body + tail;
                if ( node != i )
                {
                    n.emplace( (node << 2) | 1 );
                }
            }
        }
    }

    // Populate the adjacency list for kMinus1mers
    else if ( (enc & 3) == 0 )
    {
        // Handle insertion
        for (int j = 0; j <= k - 1; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * (j + 1));
            unsigned long int tail = (i << 1 << (63 - 2 * j)) >> (63 - 2 * j) >> 1;
            for (unsigned long int l = 0; l < 4; ++l)
            {
                unsigned long int body = l << (2 * j);
                unsigned long int node = head + body + tail;
                n.emplace( (node << 2) | 1 );
            }
        }
    }

    // Populate the adjacency list for kPlus1mers
    else
    {
        // Handle deletion
        for (int j = 1; j <= k + 1; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * (j - 1));
            unsigned long int tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            unsigned long int node = head + tail;
            n.emplace( (node << 2) | 1 );
        }
    }
}

/*
 * Reports the time and space usage
 */
void reportPerformance()
{
    string temp;
    unsigned long int utime, stime, vmpeak, vmhwm;

    // Find time usage in "/proc/self/stat"
    ifstream time_stream("/proc/self/stat", ios_base::in);

    // Skip all irrelevant attributes
    for (int i = 0; i < 13; ++i)
    {
        time_stream >> temp;
    }

    time_stream >> utime >> stime;
    time_stream.close();

    // Find space usage in "/proc/self/status"
    ifstream space_stream("/proc/self/status", ios_base::in);
    while ( temp.compare("VmPeak:") != 0 )
    {
        space_stream >> temp;
    }
    space_stream >> vmpeak;
    while ( temp.compare("VmHWM:") != 0 )
    {
        space_stream >> temp;
    }
    space_stream >> vmhwm;
    space_stream.close();

    cerr << "Time in user mode:        " << utime / sysconf(_SC_CLK_TCK) << " sec\n"
         << "Time in kernel mode:      " << stime / sysconf(_SC_CLK_TCK) << " sec\n"
         << "Peak virtual memory size: " << vmpeak << " kB\n"
         << "Peak resident set size:   " << vmhwm << " kB\n\n";
}

int main()
{
    int k;
    int d;
    cerr << "Please enter k: ";
    cin >> k;
    cerr << k << endl;
    cerr << "Plesae enter d: ";
    cin >> d;
    cerr << d << endl;

    // Initialize dist arrays for BFS
    unsigned long int num_kmers = 1;
    num_kmers = num_kmers << (2 * k);
    DistArray dist_kmer(num_kmers, d);

    unsigned long int num_kMinus1mers = 1;
    num_kMinus1mers = num_kMinus1mers << (2 * (k-1));
    DistArray dist_kMinus1mer(num_kMinus1mers, d);

    unsigned long int num_kPlus1mers = 1;
    num_kPlus1mers = num_kPlus1mers << (2 * (k+1));
    DistArray dist_kPlus1mer(num_kPlus1mers, d);

    unsigned long int num_indep_nodes = 0;
    cerr << "\nList of independent nodes: " << endl;
    for ( unsigned long int i = 0; i < num_kmers; ++i )
    {
        if ( dist_kmer[i] != (d + 1)/2 - 1 )
        {
            continue;
        }
        printKmer(i << 2, k);
        cerr << ' ';
        num_indep_nodes++;

        // Do BFS
        vector<unsigned long int> Q; // Initialize an empty queue
        Q.push_back( (i << 2) | 1 );
        unordered_map<unsigned long int, unsigned int> hist; // Keep the search history
        hist.emplace( (i << 2) | 1, 0 );

        dist_kmer.setDist(i, 0);
        while ( !Q.empty() )
        {
            auto q0 = hist.find( Q[0] );
            if ( q0->second + 1 > d )
            {
                break;
            }
            unordered_set<unsigned long int> neighbors;
            getNeighbor( Q[0], k, neighbors );

            if ( (Q[0] & 3) == 1 )
            {
                for ( auto &j : neighbors )
                {
                    if ( hist.find(j) != hist.end() )
                    {
                        continue;
                    }
                    unsigned int targetDist;
                    if ( (j & 3) == 1 )
                    {
                        targetDist = dist_kmer[j >> 2];
                        if ( targetDist == (d + 1)/2 - 1 ||
                             (targetDist != (d + 1)/2 &&
                              targetDist > q0->second + 1) )
                        {
                            dist_kmer.setDist(j >> 2, q0->second + 1);
                            Q.push_back(j);
                            hist.emplace( j, q0->second + 1 );
                        }
                    }
                    else if ( (j & 3) == 0 )
                    {
                        targetDist = dist_kMinus1mer[j >> 2];
                        if ( targetDist == (d + 1)/2 - 1 ||
                             (targetDist != (d + 1)/2 &&
                              targetDist > q0->second + 1) )
                        {
                            dist_kMinus1mer.setDist(j >> 2, q0->second + 1);
                            Q.push_back(j);
                            hist.emplace( j, q0->second + 1 );
                        }
                    }
                    else
                    {
                        targetDist = dist_kPlus1mer[j >> 2];
                        if ( targetDist == (d + 1)/2 - 1 ||
                             (targetDist != (d + 1)/2 && 
                              targetDist > q0->second + 1) )
                        {
                            dist_kPlus1mer.setDist(j >> 2, q0->second + 1);
                            Q.push_back(j);
                            hist.emplace( j, q0->second + 1 );
                        }
                    }
                }
            }
            else
            {
                for ( auto &j : neighbors )
                {
                    if ( hist.find(j) != hist.end() )
                    {
                        continue;
                    }
                    unsigned int targetDist = dist_kmer[j >> 2];
                    if ( targetDist == (d + 1)/2 - 1 || 
                         (targetDist != (d + 1)/2 &&
                          targetDist > q0->second + 1) )
                    {
                        dist_kmer.setDist(j >> 2, q0->second + 1);
                        Q.push_back(j);
                        hist.emplace( j, q0->second + 1 );
                    }
                }
            }
            Q.erase( Q.begin() );
        }
    }

    cerr << "\nThe graph has an independent set of size " << num_indep_nodes << ".\n\n";
    reportPerformance();

    return 0;
}
