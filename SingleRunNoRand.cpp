/*
 * This program addresses the following problem from Prof. Mingfu Shao (mxs2589@psu.edu).
 *
 * Let K be the set of all kmers over alphabet {A, C, G, T}. |K| = 4^k. We build a graph
 * G = (K, E) with K be the set of vertices. Let x and y be two kmers in K. We add an
 * edge to E iff the edit distance between x and y is at most d. We want to calculate the
 * maximum independent set of G.
 *
 * The input for this problem is two integers k and d; the output is a maximum
 * independent set of G constructed above.
 *
 * We can try a simple greedy algorithm. We first generate K. We then iteratively and
 * randomly pick a kmer x from K, and remove x and all other kmers whose edit distance
 * with x is at most d; we repeat this procedure until K becomes empty.
 *
 *
 * Author: Leran Ma (lkm5463@psu.edu)
 */

#include <iostream>
#include <vector>
#include <cstring>

using namespace std;

/*
 * Generates the k-mer space
 *
 * bases    : an array of the four bases 'A', 'C', 'G', 'T'
 * kmerSpace: the k-mer space
 * kmer     : a temporary string to hold the k-mer that is being editing
 * k        : the length of a k-mer
 * pos      : the position that is currently being editing
 */
void generateKmerSpace(const char *bases, vector<char *> &kmerSpace, char *kmer, const int k, int pos)
{
    if (pos == k)
    {
        char *temp = new char[k];
        strcpy(temp, kmer);
        kmerSpace.push_back(temp);
        return;
    }
    for (int i = 0; i < 4; ++i)
    {
        kmer[pos] = bases[i];
        generateKmerSpace(bases, kmerSpace, kmer, k, pos + 1);
    }
}

/*
 * Calculates the edit distance between two strings of length k
 *
 * s1: The first string
 * s2: The second string
 * k : The length of the two strings
 */
int editDist(const char *s1, const char *s2, const int k)
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
            if (s1[i-1] != s2[j-1])
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

int main()
{
    char bases[4] = {'A', 'C', 'G', 'T'};
    vector<char *> kmerSpace; // The k-mer space
    int k;
    int d;

    cout << "Please enter k:";
    cin >> k;
    cout << "Plesae enter d:";
    cin >> d;

    char temp[k];
    generateKmerSpace(bases, kmerSpace, temp, k, 0); // Generate the k-mer space


    cout << "List of independent nodes: " << endl;
    int num_of_indep_nodes = 0;
    while (!kmerSpace.empty())
    {
        char picked_kmer[k + 1];
        strcpy(picked_kmer, kmerSpace[0]);
        picked_kmer[k] = '\0';
        cout << picked_kmer << " ";
        delete [] kmerSpace[0];
        kmerSpace.erase(kmerSpace.begin());
        num_of_indep_nodes++;
        for (unsigned int i = 0; i < kmerSpace.size(); ++i) // Compare the picked k-mer with the rest of the k-mers
        {
            if (editDist(picked_kmer, kmerSpace[i], k) <= d)
            {
                delete [] kmerSpace[i];
                kmerSpace.erase(kmerSpace.begin() + i);
                i--;
            }
        }
    }
    cout << "\nThe graph has an independent set of size " << num_of_indep_nodes << ".\n";

    return 0;
}
