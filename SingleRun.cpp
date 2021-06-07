/*
 * This program addresses the following problem from Prof. Mingfu Shao (mxs2589@psu.edu).
 *
 * --------------------------------------------------------------------------------------
 * Let K be the set of all kmers over alphabet {A, C, G, T}. |K| = 4^k. We build a graph
 * G = (K, E) with K be the set of vertices. Let x and y be two kmers in K. We add an
 * edge to E iff the edit distance between x and y is at most d. We want to calculate the
 * an maximal independent set of G that covers all vertices.
 *
 * The input for this problem is two integers k and d; the output is an independent set
 * of G constructed above.
 *
 * We can try a simple greedy algorithm. We first generate K. We then iteratively and
 * randomly pick a kmer x from K, and remove x and all other kmers whose edit distance
 * with x is at most d; we repeat this procedure until K becomes empty.
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
#include <cstdlib>

using namespace std;

/*
 * Calculates the edit distance between 2 k-mers
 *
 * s1: The encoding of the first k-mer
 * s2: The encoding of the second k-mer
 * k : The length of the two k-mers
 */
int editDist(const unsigned long int s1, const unsigned long int s2, const int k)
{
    int DPtable[k + 1][k + 1];
    for (int i = 0; i < k + 1; ++i)
    {
        DPtable[0][i] = i;
        DPtable[i][0] = i;
    }
    for (int i = 1; i < k + 1; ++i)
    {
        for (int j = 1; j < k + 1; ++j)
        {
            if ( ((s1 >> (2*(i-1))) & 3) != ((s2 >> (2*(j-1))) & 3) )
            {
                DPtable[i][j] = DPtable[i-1][j-1] + 1;
            }
            else
            {
                DPtable[i][j] = DPtable[i-1][j-1];
            }
            if (DPtable[i-1][j] + 1 < DPtable[i][j])
            {
                DPtable[i][j] = DPtable[i-1][j] + 1;
            }
            if (DPtable[i][j-1] + 1 < DPtable[i][j])
            {
                DPtable[i][j] = DPtable[i][j-1] + 1;
            }
        }
    }
    return DPtable[k][k];
}

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

int main()
{
    int k;
    int d;
    cout << "Please enter k: ";
    cin >> k;
    cout << "Plesae enter d: ";
    cin >> d;

    // Generate the k-mer space
    vector<unsigned long int> kmerSpace;
    unsigned long int kmerSpaceSize = 1;
    kmerSpaceSize = kmerSpaceSize << (2 * k);
    for (unsigned long int i = 0; i < kmerSpaceSize; ++i)
    {
        kmerSpace.push_back(i);
    }

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

        // Compare the picked k-mer with the rest of the k-mers
        for (unsigned long int i = 0; i < kmerSpace.size(); ++i)
        {
            if (editDist(picked_kmer, kmerSpace[i], k) <= d)
            {
                kmerSpace.erase(kmerSpace.begin() + i);
                i--;
            }
        }
    }

    cout << "\nThe graph has an independent set of size " << num_indep_nodes << ".\n";

    return 0;
}