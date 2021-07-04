/*
 * This program addresses the following problem from Prof. Mingfu Shao (mxs2589@psu.edu).
 *
 * --------------------------------------------------------------------------------------
 * Let K be the set of all kmers over alphabet {A, C, G, T}. |K| = 4^k. We build a graph
 * G = (K, E) with K be the set of vertices. Let x and y be two kmers in K. We add an
 * edge to E iff the edit distance between x and y is at most d. We want to calculate the
 * an maximal independent set (MIS) of G that covers all vertices.
 *
 * The input for this problem is two integers k and d; the output is an independent set
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
 * Date:   7:15 PM, Saturday, July 3, 2021
 */

#include <iostream>
#include <vector>
#include <map>

using namespace std;

/*
 * Prints a k-mer given its binary encoding
 *
 * enc: The binary encoding of the k-mer
 * k  : The length of the k-mer
 */
void printKmer(unsigned long int enc, int k)
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
 * enc : The binary encoding of the k-mer
 * k   : The length of the k-mer
 * temp: A vector to hold the neighbors
 */
void getNeighbor( unsigned long int enc, int k, vector<unsigned long int> &temp )
{
    unsigned long int i = enc >> 2;

    // First store neighbors in a map to avoid duplication
    map<unsigned long int, bool> m;

    // Populate the adjacency list for kmers
    if ( (enc & 3) == 1 )
    {
        // Handle deletion
        for (int j = 1; j <= k; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * (j - 1));
            unsigned long int tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            m.emplace( (head + tail) << 2, true );
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
                m.emplace( (node << 2) | 2, true );
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
                    m.emplace( (node << 2) | 1, true );
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
                m.emplace( (node << 2) | 1, true );
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
            m.emplace( (node << 2) | 1, true );
        }
    }

    for ( auto &j : m )
    {
        temp.push_back( j.first );
    }
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

    // Generate the k-mer space
    unsigned long int num_kmers = 1;
    num_kmers = num_kmers << (2 * k);
    vector<bool> kmerSpace(num_kmers, false);

    // Initialize the dist array for BFS
    vector<int> dist_kmer(num_kmers, d + 1);

    unsigned long int num_kMinus1mers = 1;
    num_kMinus1mers = num_kMinus1mers << (2 * (k-1));
    vector<int> dist_kMinus1mer(num_kMinus1mers, d + 1);

    unsigned long int num_kPlus1mers = 1;
    num_kPlus1mers = num_kPlus1mers << (2 * (k+1));
    vector<int> dist_kPlus1mer(num_kPlus1mers, d + 1);

    int num_indep_nodes = 0;
    cerr << "\nList of independent nodes: " << endl;
    for (unsigned long int i = 0; i < num_kmers; ++i)
    {
        if ( kmerSpace[i] )
        {
            continue;
        }
        printKmer(i << 2, k);
        cerr << ' ';
        num_indep_nodes++;
        kmerSpace[i] = true;

        // Do BFS
        vector<unsigned long int> Q; // Initialize an empty queue
        Q.push_back( (i << 2) | 1 );
        dist_kmer[i] = 0;
        while ( !Q.empty() )
        {
            vector<unsigned long int> neighbors;
            getNeighbor( Q[0], k, neighbors );

            if ( (Q[0] & 3) == 1 )
            {
                for ( auto &j : neighbors )
                {
                    if ( (j & 3) == 1 )
                    {
                        if (dist_kmer[j >> 2] > dist_kmer[Q[0] >> 2] + 1)
                        {
                            kmerSpace[j >> 2] = true;
                            dist_kmer[j >> 2] = dist_kmer[Q[0] >> 2] + 1;
                            Q.push_back(j);
                        }
                    }
                    else if ((j & 3) == 0 )
                    {
                        if (dist_kMinus1mer[j >> 2] > dist_kmer[Q[0] >> 2] + 1)
                        {
                            dist_kMinus1mer[j >> 2] = dist_kmer[Q[0] >> 2] + 1;
                            Q.push_back(j);
                        }
                    }
                    else
                    {
                        if (dist_kPlus1mer[j >> 2] > dist_kmer[Q[0] >> 2] + 1)
                        {
                            dist_kPlus1mer[j >> 2] = dist_kmer[Q[0] >> 2] + 1;
                            Q.push_back(j);
                        }
                    }
                }
            }
            else if ( (Q[0] & 3) == 0 )
            {
                for ( auto &j : neighbors )
                {
                    if (dist_kmer[j >> 2] > dist_kMinus1mer[Q[0] >> 2] + 1)
                    {
                        kmerSpace[j >> 2] = true;
                        dist_kmer[j >> 2] = dist_kMinus1mer[Q[0] >> 2] + 1;
                        Q.push_back(j);
                    }
                }
            }
            else
            {
                for ( auto &j : neighbors )
                {
                    if (dist_kmer[j >> 2] > dist_kPlus1mer[Q[0] >> 2] + 1)
                    {
                        kmerSpace[j >> 2] = true;
                        dist_kmer[j >> 2] = dist_kPlus1mer[Q[0] >> 2] + 1;
                        Q.push_back(j);
                    }
                }
            }
            Q.erase(Q.begin());
        }
    }

    cerr << "\nThe graph has an independent set of size " << num_indep_nodes << ".\n";

    return 0;
}
