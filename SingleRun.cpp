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
 * characters can be represented by a 64-bit binary number stored in an 8-byte slot.
 *
 * Author: Leran Ma (lkm5463@psu.edu)
 */

#include <iostream>
#include <vector>
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
        kmer[i] = base[enc & 3];
        enc = enc >> 2;
    }
    cout << kmer;
}

struct node
{
    unsigned long int read; // The encoding of a read of length k-1, k, or k+1
    int length; // The lenght of the read

    /*
     * A constructor
     *
     * r: The encoding of a read of length k-1, k, or k+1
     * l: The lenght of the read
     */
    node( unsigned long int r, int l )
    {
        read = r;
        length = l;
    }

    /*
     * Overload the operator ==
     */
    bool operator==(node const &rhs) const
    {
        return read == rhs.read && length == rhs.length;
    }
};

int main()
{
    int k;
    int d;
    cout << "Please enter k: ";
    cin >> k;
    cout << "Plesae enter d: ";
    cin >> d;

    // Generate the k-mer space
    unsigned long int kmerSpaceSize = 1;
    kmerSpaceSize = kmerSpaceSize << (2 * k);
    vector<unsigned long int> kmerSpace(kmerSpaceSize);
    for (unsigned long int i = 0; i < kmerSpaceSize; ++i)
    {
        kmerSpace[i] = i;
    }

    // Generate the adjacency list
    vector<node> temp;
    unsigned long int num_kMinus1mers = 1;
    num_kMinus1mers = num_kMinus1mers << (2 * (k - 1));
    vector< vector<node> > kMinus1mers( num_kMinus1mers, temp );

    unsigned long int num_kmers = 1;
    num_kmers = num_kmers << (2 * k);
    vector< vector<node> > kmers( num_kmers, temp );

    unsigned long int num_kPlus1mers = 1;
    num_kPlus1mers = num_kPlus1mers << (2 * (k + 1));
    vector< vector<node> > kPlus1mers( num_kPlus1mers, temp );

    for (unsigned long int i = 0; i < num_kmers; ++i)
    {
        // Handle deletion
        for (int j = 1; j <= k; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * (j - 1));
            unsigned long int tail = (i << (64 - 2 * (j - 1))) >> (64 - 2 * (j - 1));
            node n( head + tail, k - 1 );
            auto it = find(kmers[i].begin(), kmers[i].end(), n); // Check for duplication
            if ( it == kmers[i].end() )
            {
                kmers[i].push_back( n );
                node m( i, k );
                kMinus1mers[head + tail].push_back( m );
            }
        }

        // Handle insertion
        for (int j = 0; j <= k; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * (j + 1));
            unsigned long int tail = (i << (64 - 2 * j)) >> (64 - 2 * j);
            for (unsigned long int l = 0; l < 4; ++l)
            {
                unsigned long int body = l << (2 * j);
                node n( head + body + tail, k + 1 );
                auto it = find(kmers[i].begin(), kmers[i].end(), n);
                if ( it == kmers[i].end() )
                {
                    kmers[i].push_back( n );
                    node m( i, k );
                    kPlus1mers[head + body + tail].push_back( m );
                }
            }
        }

        // Handle substitution
        for (int j = 1; j <= k; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * j);
            unsigned long int tail = (i << (64 - 2 * (j - 1))) >> (64 - 2 * (j - 1));
            for (unsigned long int l = 0; l < 4; ++l)
            {
                unsigned long int body = l << (2 * (j - 1));
                node n( head + body + tail, k );
                auto it = find(kmers[i].begin(), kmers[i].end(), n);
                if ( it == kmers[i].end() )
                {
                    kmers[i].push_back( n );
                }
            }
        }
    }

    // Initialize the dist array for BFS
    vector<int> dist_kmer(num_kmers, d + 1);
    vector<int> dist_kMinus1mer(num_kMinus1mers, d + 1);
    vector<int> dist_kPlus1mer(num_kPlus1mers, d + 1);

    srand( time(nullptr) );
    int num_indep_nodes = 0;
    cout << "List of independent nodes: " << endl;
    while ( !kmerSpace.empty() )
    {
        int picked = rand() % kmerSpace.size();
        unsigned long int picked_kmer = kmerSpace[picked];
        printKmer(picked_kmer, k);
        cout << ' ';
        num_indep_nodes++;
        kmerSpace.erase(kmerSpace.begin() + picked);

        // Do BFS
        vector<node> Q; // Initialize an empty queue
        node n(picked_kmer, k);
        Q.push_back(n);
        dist_kmer[picked_kmer] = 0;
        while ( !Q.empty() )
        {
            if ( Q[0].length == k )
            {
                for ( auto &i : kmers[Q[0].read] )
                {
                    if (i.length == k)
                    {
                        if (dist_kmer[i.read] > dist_kmer[Q[0].read] + 1)
                        {
                            dist_kmer[i.read] = dist_kmer[Q[0].read] + 1;
                            Q.push_back(i);
                            auto it = find(kmerSpace.begin(), kmerSpace.end(), i.read);
                            if ( it != kmerSpace.end() )
                            {
                                kmerSpace.erase(it);
                            }
                        }
                    }
                    else if (i.length == k - 1)
                    {
                        if (dist_kMinus1mer[i.read] > dist_kmer[Q[0].read] + 1) {
                            dist_kMinus1mer[i.read] = dist_kmer[Q[0].read] + 1;
                            Q.push_back(i);
                        }
                    }
                    else
                    {
                        if (dist_kPlus1mer[i.read] > dist_kmer[Q[0].read] + 1)
                        {
                            dist_kPlus1mer[i.read] = dist_kmer[Q[0].read] + 1;
                            Q.push_back(i);
                        }
                    }
                }
            }
            else if ( Q[0].length == k - 1 )
            {
                for ( auto &i : kMinus1mers[Q[0].read] )
                {
                    if (dist_kmer[i.read] > dist_kMinus1mer[Q[0].read] + 1)
                    {
                        dist_kmer[i.read] = dist_kMinus1mer[Q[0].read] + 1;
                        Q.push_back(i);
                        auto it = find(kmerSpace.begin(), kmerSpace.end(), i.read);
                        if ( it != kmerSpace.end() )
                        {
                            kmerSpace.erase(it);
                        }
                    }
                }
            }
            else
            {
                for ( auto &i : kPlus1mers[Q[0].read] )
                {
                    if (dist_kmer[i.read] > dist_kPlus1mer[Q[0].read] + 1)
                    {
                        dist_kmer[i.read] = dist_kPlus1mer[Q[0].read] + 1;
                        Q.push_back(i);
                        auto it = find(kmerSpace.begin(), kmerSpace.end(), i.read);
                        if ( it != kmerSpace.end() )
                        {
                            kmerSpace.erase(it);
                        }
                    }
                }
            }
            Q.erase(Q.begin());
        }
    }

    cout << "\nThe graph has an independent set of size " << num_indep_nodes << ".\n";

    return 0;
}