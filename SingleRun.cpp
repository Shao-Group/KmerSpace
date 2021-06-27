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
 * We can try a greedy algorithm. We first generate K. We then iteratively and randomly
 * pick a kmer x from K, and remove x and all other kmers whose edit distance with x is
 * at most d; we repeat this procedure until K becomes empty. In order to save time, we
 * do not do simple pairwise comparison to find the edit distance. Instead, we build a
 * graph containing all possible k-mers, (k-1)-mers, and (k+1)-mers as vertices. If the
 * edit distance between any two vertices is 1, we add an edge of length 1 between the
 * two vertices. Then, we add the picked vertex to the MIS, explore all vertices at a
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
 */

#include <iostream>
#include <vector>
#include <map>
#include <ctime>
#include <algorithm>

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
    cout << kmer;
}

/*
 * Remove duplicate elements within a vector
 */
void removeDuplicates( vector<unsigned long int> &temp )
{
    map<unsigned long int, bool> m;
    for ( unsigned long &i : temp )
    {
        m.emplace( i, true );
    }
    temp.clear();
    for ( auto &i : m )
    {
        temp.push_back( i.first );
    }
}

/*
 * Find and erase an element from the k-mer space. In order to be faster, we start
 * finding from the position of the last target of this function and iterate forward
 * and backward.
 *
 * kmerSpace: The k-mer space
 * enc      : The binary encoding of the k-mer to be erased
 * lastFound: The index of the last target of this function
 */
void findAndErase( vector<unsigned long int> &kmerSpace,
                   unsigned long int enc,
                   unsigned long int &lastFound )
{
    if ( lastFound >= kmerSpace.size() )
    {
        for (unsigned long int i = kmerSpace.size() - 1; i <= 0 ; --i)
        {
            if ( kmerSpace[i] == enc )
            {
                kmerSpace.erase( kmerSpace.begin() + i );
                lastFound = i;
                return;
            }
        }
    }

    unsigned long int d = lastFound;
    if ( kmerSpace.size() - lastFound - 1 < d )
    {
        d = kmerSpace.size() - lastFound - 1;
    }

    for (unsigned long int i = 0; i <= d; ++i)
    {
        if ( kmerSpace[lastFound + i] == enc )
        {
            kmerSpace.erase( kmerSpace.begin() + lastFound + i );
            lastFound = lastFound + i;
            return;
        }
        if ( kmerSpace[lastFound - i] == enc )
        {
            kmerSpace.erase( kmerSpace.begin() + lastFound - i );
            lastFound = lastFound - i;
            return;
        }
    }

    if ( d == lastFound )
    {
        for (unsigned long int i = lastFound + d + 1; i < kmerSpace.size(); ++i)
        {
            if ( kmerSpace[i] == enc )
            {
                kmerSpace.erase( kmerSpace.begin() + i );
                lastFound = i;
                return;
            }
        }
    }
    else
    {
        for (unsigned long int i = lastFound - d - 1; i >= 0; --i)
        {
            if ( kmerSpace[i] == enc )
            {
                kmerSpace.erase( kmerSpace.begin() + i );
                lastFound = i;
                return;
            }
        }
    }
}

int main()
{
    int k;
    int d;
    cout << "Please enter k: ";
    cin >> k;
    cout << "Plesae enter d: ";
    cin >> d;

    // Generate the k-mer space
    unsigned long int num_kmers = 1;
    num_kmers = num_kmers << (2 * k);
    vector<unsigned long int> kmerSpace(num_kmers);
    for (unsigned long int i = 0; i < num_kmers; ++i)
    {
        kmerSpace[i] = i;
    }

    // Generate the adjacency list
    vector<unsigned long int> temp;
    vector< vector<unsigned long int> > kmers( num_kmers, temp );

    unsigned long int num_kMinus1mers = 1;
    num_kMinus1mers = num_kMinus1mers << (2 * (k - 1));
    vector< vector<unsigned long int> > kMinus1mers( num_kMinus1mers, temp );

    unsigned long int num_kPlus1mers = 1;
    num_kPlus1mers = num_kPlus1mers << (2 * (k + 1));
    vector< vector<unsigned long int> > kPlus1mers( num_kPlus1mers, temp );

    // Populate the adjacency list for kmers
    for (unsigned long int i = 0; i < num_kmers; ++i)
    {
        // Handle deletion
        for (int j = 1; j <= k; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * (j - 1));
            unsigned long int tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            kmers[i].push_back( (head + tail) << 2 );
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
                kmers[i].push_back( (node << 2) | 2 );
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
                    kmers[i].push_back( (node << 2) | 1 );
                }
            }
        }

        removeDuplicates( kmers[i] );
    }

    // Populate the adjacency list for kMinus1mers
    for ( unsigned long int i = 0; i < num_kMinus1mers; ++i )
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
                kMinus1mers[i].push_back( (node << 2) | 1 );
            }
        }
        removeDuplicates( kMinus1mers[i] );
    }

    // Populate the adjacency list for kPlus1mers
    for ( unsigned long int i = 0; i < num_kPlus1mers; ++i )
    {
        // Handle deletion
        for (int j = 1; j <= k + 1; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * (j - 1));
            unsigned long int tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            unsigned long int node = head + tail;
            kPlus1mers[i].push_back( (node << 2) | 1 );
        }
        removeDuplicates( kPlus1mers[i] );
    }

    // Initialize the dist array for BFS
    vector<int> dist_kmer(num_kmers, d + 1);
    vector<int> dist_kMinus1mer(num_kMinus1mers, d + 1);
    vector<int> dist_kPlus1mer(num_kPlus1mers, d + 1);

    srand( time(nullptr) );
    int num_indep_nodes = 0;
    unsigned long int lastFound;
    cout << "List of independent nodes: " << endl;
    while ( !kmerSpace.empty() )
    {
        int picked = rand() % kmerSpace.size();
        unsigned long int picked_kmer = kmerSpace[picked];
        printKmer(picked_kmer << 2, k);
        cout << ' ';
        num_indep_nodes++;
        kmerSpace.erase(kmerSpace.begin() + picked);
        lastFound = picked;

        // Do BFS
        vector<unsigned long int> Q; // Initialize an empty queue
        Q.push_back( (picked_kmer << 2) | 1 );
        dist_kmer[picked_kmer] = 0;
        while ( !Q.empty() )
        {
            if ( (Q[0] & 3) == 1 )
            {
                for ( auto &i : kmers[Q[0] >> 2] )
                {
                    if ( (i & 3) == 1 )
                    {
                        if (dist_kmer[i >> 2] > dist_kmer[Q[0] >> 2] + 1)
                        {
                            dist_kmer[i >> 2] = dist_kmer[Q[0] >> 2] + 1;
                            Q.push_back(i);
                            findAndErase( kmerSpace, i >> 2, lastFound );
                        }
                    }
                    else if ( (i & 3) == 0 )
                    {
                        if (dist_kMinus1mer[i >> 2] > dist_kmer[Q[0] >> 2] + 1)
                        {
                            dist_kMinus1mer[i >> 2] = dist_kmer[Q[0] >> 2] + 1;
                            Q.push_back(i);
                        }
                    }
                    else
                    {
                        if (dist_kPlus1mer[i >> 2] > dist_kmer[Q[0] >> 2] + 1)
                        {
                            dist_kPlus1mer[i >> 2] = dist_kmer[Q[0] >> 2] + 1;
                            Q.push_back(i);
                        }
                    }
                }
            }
            else if ( (Q[0] & 3) == 0 )
            {
                for ( auto &i : kMinus1mers[Q[0] >> 2] )
                {
                    if (dist_kmer[i >> 2] > dist_kMinus1mer[Q[0] >> 2] + 1)
                    {
                        dist_kmer[i >> 2] = dist_kMinus1mer[Q[0] >> 2] + 1;
                        Q.push_back(i);
                        findAndErase( kmerSpace, i >> 2, lastFound );
                    }
                }
            }
            else
            {
                for ( auto &i : kPlus1mers[Q[0] >> 2] )
                {
                    if (dist_kmer[i >> 2] > dist_kPlus1mer[Q[0] >> 2] + 1)
                    {
                        dist_kmer[i >> 2] = dist_kPlus1mer[Q[0] >> 2] + 1;
                        Q.push_back(i);
                        findAndErase( kmerSpace, i >> 2, lastFound );
                    }
                }
            }
            Q.erase(Q.begin());
        }
    }

    cout << "\nThe graph has an independent set of size " << num_indep_nodes << ".\n";

    return 0;
}