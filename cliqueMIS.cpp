/*
 * Author: Leran Ma (lkm5463@psu.edu)
 * Date:   6:10 PM, Wednesday, October 27, 2021
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
 * A class for the dist array used in BFS. The array is chopped into subarrays to avoid
 * failures of dynamic allocation.
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
     * s : total number of elements
     * d : maximum edit distance allowed
     */
    DistArray( unsigned long int s, unsigned int d )
    {
        sub_size = 1;
        sub_size = sub_size << 20;
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
 * A data structure to store the coveage information for kmers. The array is chopped into
 * subarrays to avoid failures of dynamic allocation.
 */
class CoverageArray
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
     * s: capacity of the array
     */
    CoverageArray( unsigned long int s )
    {
        sub_size = 1;
        sub_size = sub_size << 20;
        num_subs = s * 5 / sub_size; // Each element occupies 5 bytes.
        if ( (s * 5) % sub_size != 0 )
        {
            num_subs++;
        }
        aptrs = (char **) calloc( num_subs, sizeof(char *) );
        for (int i = 0; i < num_subs; ++i)
        {
            // Add a few extra bytes to the end of each subarray to avoid potential out
            // of range access
            aptrs[i] = (char *) calloc( sub_size + 5, 1 );
        }
        size = s;
    }

    /*
     * Destructor
     */
    ~CoverageArray()
    {
        for (int i = 0; i < num_subs; ++i)
        {
            free( aptrs[i] );
        }
        free( aptrs );
    }

    /*
     * Overload [] operator to return the coverage of the kmer indexed by sub
     *
     * sub: The index of the element to be extract
     */
    unsigned long int operator[]( const unsigned long int &sub )
    {
        // Find the byte position where the required element is located
        // Use modulo to avoid index out of range
        unsigned long int bytePos = (sub % size) * 5;

        unsigned long int cov = 0;
        memcpy( &cov, &aptrs[bytePos/sub_size][bytePos%sub_size], 5);
        return cov;
    }
    
    /*
     * Set the coverage for an element indexed by sub
     *
     * sub: The index of the element
     * cov: The coverage of the element
     */
    void setCov( const unsigned long int sub, unsigned long int cov )
    {
        unsigned long int bytePos = (sub % size) * 5;
        memcpy( &aptrs[bytePos/sub_size][bytePos%sub_size], &cov, 5);
    }
};

/*
 * A class for the visited array. The array is chopped into subarrays to avoid failures
 * of dynamic allocation.
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
    VisitedArray( unsigned long int s )
    {
        sub_size = 1;
        sub_size = sub_size << 20;
        num_subs = s / 8 / sub_size; // Each element occupies 1/8 bytes.
        if ( (s / 8) % sub_size != 0 )
        {
            num_subs++;
        }
        aptrs = (char **) calloc( num_subs, sizeof(char *) );
        for (int i = 0; i < num_subs; ++i)
        {
            aptrs[i] = (char *) calloc( sub_size, 1 );
        }
        size = s;
    }

    /*
     * Destructor
     */
    ~VisitedArray()
    {
        for (int i = 0; i < num_subs; ++i)
        {
            free( aptrs[i] );
        }
        free( aptrs );
    }

    /*
     * Overload [] operator to return the visited value indexed by sub
     *
     * sub: The index of the element to be extract
     */
    bool operator[]( const unsigned long int &sub )
    {
        // Find the bit offset
        int offset = sub % 8;

        // Find the byte in which the required element is located
        unsigned long int bytePos = sub / 8;

        unsigned int visited;
        memcpy( &visited, &aptrs[bytePos/sub_size][bytePos%sub_size], 1 );

        visited = (visited << (24 + offset)) >> 31;
        return (visited == 1);
    }
    
    /*
     * Set the element indexed by sub as visited
     *
     * sub : The index of the element
     */
    void setVisit( const unsigned long int sub )
    {
        // Find the byte in which the required element is located
        unsigned long int bytePos = sub / 8;

        // Find the bit offset
        int offset = sub % 8;

        unsigned int tempByte;
        memcpy( &tempByte, &aptrs[bytePos/sub_size][bytePos%sub_size], 1 );
        unsigned int head = tempByte >> (8 - offset);
        head = head << (8 - offset);
        unsigned int tail = tempByte << 24 << (offset + 1);
        tail = tail >> 24 >> (offset + 1);
        unsigned int body = 1;
        body = body << (8 - offset - 1);
        tempByte = tempByte | body;
        memcpy( &aptrs[bytePos/sub_size][bytePos%sub_size], &tempByte, 1 );
    }
};

/*
 * Calculates the edit distance between 2 k-mers
 *
 * s1: The encoding of the first k-mer
 * s2: The encoding of the second k-mer
 * k1: The length of the first k-mer
 * k2: The length of the second k-mer
 * d : The maximum edit distance allowed
 */
int editDist( const unsigned long int s1, 
              const unsigned long int s2, 
              const int k1,
              const int k2,
              const int d )
{
    int DPtable[k1 + 1][k2 + 1];
    for (int i = 0; i < k1 + 1; ++i)
    {
        DPtable[i][0] = i;
    }
    for (int i = 0; i < k2 + 1; ++i)
    {
        DPtable[0][i] = i;
    }
    for (int i = 1; i < k1 + 1; ++i)
    {
        for (int j = 1; j < k2 + 1; ++j)
        {
            if ( ((s1 >> (2*(i-1))) & 3) != ((s2 >> (2*(j-1))) & 3) )
            {
                DPtable[i][j] = DPtable[i-1][j-1] + 1;
            }
            else
            {
                DPtable[i][j] = DPtable[i-1][j-1];
            }
            if ( DPtable[i-1][j] + 1 < DPtable[i][j] )
            {
                DPtable[i][j] = DPtable[i-1][j] + 1;
            }
            if ( DPtable[i][j-1] + 1 < DPtable[i][j] )
            {
                DPtable[i][j] = DPtable[i][j-1] + 1;
            }
        }
        if ( DPtable[i][i] > d && k1 == k2 )
        {
            return DPtable[i][i];
        }
    }
    return DPtable[k1][k2];
}

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
    else
    {
        // Handle substitution
        for (int j = 1; j <= k - 1; ++j)
        {
            unsigned long int head = (i >> (2 * j)) << (2 * j);
            unsigned long int tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            for (unsigned long int l = 0; l < 4; ++l)
            {
                unsigned long int body = l << (2 * (j - 1));
                unsigned long int node = head + body + tail;
                if ( node != i )
                {
                    n.emplace( node << 2 );
                }
            }
        }

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
}

/*
 * Asks neighbors for possible coverage. Returns true if a feasible answer is found.
 *
 * enc         : The binary encoding of the k-mer
 * k           : The length of the k-mer
 * d           : The maximum edit distance allowed
 * coverage    : The coverage array
 * isKMinus1mer: A boolean indicating whether the given enc is of a (k-1)-mer
 */
bool askNeighbors( const unsigned long int enc,
                   const int k,
                   const int d,
                   CoverageArray &coverage,
                   bool isKMinus1mer )
{
    unordered_set<unsigned long int> asked;   // A set to store asked neighbors
    unordered_set<unsigned long int> checked; // A set to store checked possibilities
    if ( isKMinus1mer )
    {
        for (int j = 0; j <= k - 1; ++j)
        {
            unsigned long int head = (enc >> (2 * j)) << (2 * (j + 1));
            unsigned long int tail = (enc << 1 << (63 - 2 * j)) >> (63 - 2 * j) >> 1;
            for (unsigned long int l = 0; l < 4; ++l)
            {
                unsigned long int body = l << (2 * j);
                unsigned long int node = head + body + tail;
                if ( asked.emplace(node).second )
                {
                    unsigned long int temp = coverage[node];
                    int k1;
                    if ( (temp & 3) == 0 )
                    {
                        k1 = k - 1;
                    }
                    else
                    {
                        k1 = k;
                    }
                    if ( checked.emplace(temp).second && editDist(temp >> 2, enc, k1, k-1, d) <= d )
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    for ( int j = 1; j <= k; ++j )
    {
        unsigned long int head = (enc >> (2 * j)) << (2 * j);
        unsigned long int tail = (enc << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
        for ( unsigned long int l = 0; l < 4; ++l )
        {
            unsigned long int body = l << (2 * (j - 1));
            unsigned long int node = head + body + tail;
            if ( asked.emplace(node).second )
            {
                unsigned long int temp = coverage[node];
                int k1;
                if ( (temp & 3) == 0 )
                {
                    k1 = k - 1;
                }
                else
                {
                    k1 = k;
                }
                if ( checked.emplace(temp).second && editDist(temp >> 2, enc, k1, k, d) <= d )
                {
                    coverage.setCov(enc, temp);
                    return true;
                }
            }
        }
    }
    return false;
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

void doBFS( const int k, const int d )
{
    // Initialize dist arrays for BFS
    unsigned long int num_kmers = 1;
    num_kmers = num_kmers << (2 * k);
    DistArray dist_kmer(num_kmers, d);

    unsigned long int num_kMinus1mers = 1;
    num_kMinus1mers = num_kMinus1mers << (2 * (k-1));
    DistArray dist_kMinus1mer(num_kMinus1mers, d);

    unsigned long int num_indep_cliques = 0;
    cerr << "\nList of independent cliques: " << endl;
    for ( unsigned long int i = 0; i < num_kmers; ++i )
    {
        if ( dist_kmer[i] != (d + 1)/2 - 1 )
        {
            continue;
        }

        // Find the largest clique that contains i
        unordered_set<unsigned long int> clique;
        clique.emplace( (i << 2) | 1 );
        for (int j = 1; j <= k; ++j)
        {
            unordered_set<unsigned long int> temp;
            temp.emplace( (i << 2) | 1 );
            unsigned long int head = (i >> (2 * j)) << (2 * j);
            unsigned long int tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            for (unsigned long int l = 0; l < 4; ++l)
            {
                unsigned long int body = l << (2 * (j - 1));
                unsigned long int node = head + body + tail;
                if ( dist_kmer[node] == (d + 1)/2 - 1 )
                {
                    temp.emplace( (node << 2) | 1 );
                }
            }
            head = (i >> (2 * j)) << (2 * (j - 1));
            tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            if ( dist_kMinus1mer[ head + tail ] == (d + 1)/2 - 1 )
            {
                temp.emplace( (head + tail) << 2 );
            }

            if ( temp.size() > clique.size() )
            {
                clique = temp;
            }
            if ( clique.size() == 5 )
            {
                break;
            }
        }
        num_indep_cliques++;

        // Do BFS
        vector<unsigned long int> Q; // Initialize an empty queue
        unordered_map<unsigned long int, unsigned int> hist; // Keep the search history
        for ( const unsigned long int &j : clique )
        {
            if ( (j & 3) == 1 )
            {
                printKmer( j, k );
                cerr << ' ';
                dist_kmer.setDist(j >> 2, 0);
            }
            else
            {
                printKmer( j, k - 1 );
                cerr << ' ';
                dist_kMinus1mer.setDist(j >> 2, 0);
            }
            Q.push_back( j );
            hist.emplace( j, 0 );
        }
        cerr << endl;

        while ( !Q.empty() )
        {
            auto q0 = hist.find( Q[0] );
            if ( q0->second + 1 > d )
            {
                break;
            }
            unordered_set<unsigned long int> neighbors;
            getNeighbor( Q[0], k, neighbors );

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
                else
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
            }
            Q.erase( Q.begin() );
        }
    }

    cerr << "\nThe graph has " << num_indep_cliques << " independent cliques.\n\n";
    reportPerformance();
}

void doNN( const int k, const int d )
{
    unsigned long int kmerSpaceSize = 1;
    kmerSpaceSize = kmerSpaceSize << (2 * k);
    VisitedArray visit_kmer(kmerSpaceSize);
    vector<unsigned long int> MIS;
    CoverageArray coverage(kmerSpaceSize / 4);
    unsigned long int num_indep_cliques = 0;

    cerr << "\nList of independent nodes: " << endl;
    for (unsigned long int l = 0; l < 4; ++l)
    {
        MIS.push_back( (l << 2) | 1 ); // k-mers
        coverage.setCov( l, (l << 2) | 1 );
        visit_kmer.setVisit( l );
        printKmer( (l << 2) | 1, k );
        cerr << ' ';
    }
    MIS.push_back( 0 ); // (k-1)-mer
    printKmer( 0, k-1 );
    cerr << endl;
    num_indep_cliques++;

    bool isCovered = false;
    for ( unsigned long int i = 0; i < kmerSpaceSize; ++i )
    {
        if ( visit_kmer[i] )
        {
            continue;
        }

        if ( askNeighbors(i, k, d, coverage, false) )
        {
            visit_kmer.setVisit( i );
            continue;
        }

        for ( const unsigned long int &j : MIS )
        {
            int k2;
            if ( (j & 3) == 0 )
            {
                k2 = k - 1;
            }
            else
            {
                k2 = k;
            }

            if ( editDist(i, j >> 2, k, k2, d) <= d )
            {
                coverage.setCov(i, j);
                visit_kmer.setVisit( i );
                isCovered = true;
                break;
            }
        }
        if ( isCovered )
        {
            isCovered = false;
            continue;
        }

        // Find the largest clique that contains i
        unordered_set<unsigned long int> clique;
        clique.emplace( (i << 2) | 1 );
        for (int j = 1; j <= k; ++j)
        {
            unordered_set<unsigned long int> temp;
            temp.emplace( (i << 2) | 1 );
            unsigned long int head = (i >> (2 * j)) << (2 * j);
            unsigned long int tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;
            for (unsigned long int l = 0; l < 4; ++l)
            {
                unsigned long int body = l << (2 * (j - 1));
                unsigned long int node = head + body + tail;
                if ( !visit_kmer[node] )
                {
                    if ( askNeighbors(node, k, d, coverage, false) )
                    {
                        visit_kmer.setVisit( node );
                        continue;
                    }

                    for ( const unsigned long int &j : MIS )
                    {
                        int k2;
                        if ( (j & 3) == 0 )
                        {
                            k2 = k - 1;
                        }
                        else
                        {
                            k2 = k;
                        }

                        if ( editDist(node, j >> 2, k, k2, d) <= d )
                        {
                            coverage.setCov(node, j);
                            visit_kmer.setVisit( node );
                            isCovered = true;
                            break;
                        }
                    }
                    if ( isCovered )
                    {
                        isCovered = false;
                        continue;
                    }
                    temp.emplace( (node << 2) | 1 );
                }
            }

            head = (i >> (2 * j)) << (2 * (j - 1));
            tail = (i << 1 << (63 - 2 * (j - 1))) >> (63 - 2 * (j - 1)) >> 1;

            isCovered = askNeighbors(head + tail, k, d, coverage, true);
            if ( !isCovered )
            {
                for ( const unsigned long int &j : MIS )
                {
                    int k2;
                    if ( (j & 3) == 0 )
                    {
                        k2 = k - 1;
                    }
                    else
                    {
                        k2 = k;
                    }
                    if ( editDist((head + tail), j >> 2, k-1, k2, d) <= d )
                    {
                        isCovered = true;
                        break;
                    }
                }
            }
            if ( !isCovered )
            {
                temp.emplace( (head + tail) << 2 );
            }
            else
            {
                isCovered = false;
            }

            if ( temp.size() > clique.size() )
            {
                clique = temp;
            }
            if ( clique.size() == 5 )
            {
                break;
            }
        }

        for ( const unsigned long int &j : clique )
        {
            if ( (j & 3) == 1 )
            {
                visit_kmer.setVisit(j >> 2);
                printKmer( j, k );
                cerr << ' ';
                coverage.setCov( j >> 2, j );
            }
            else
            {
                printKmer( j, k - 1 );
                cerr << ' ';
            }
            MIS.push_back( j );
        }
        cerr << endl;
        num_indep_cliques++;
    }

    cerr << "\nThe graph has " << num_indep_cliques << " independent cliques.\n\n";
    reportPerformance();
}

int main(int argc, char *argv[])
{
    int k = 15;
    int d = 1;

    for (int i = 0; i < argc; ++i)
    {
        if ( strcmp(argv[i], "-k") == 0 )
        {
            k = atoi( argv[i+1] );
        }
        else if ( strcmp(argv[i], "-d") == 0 )
        {
            d = atoi( argv[i+1] );
        }
    }

    cerr << "Please enter k: ";
    cerr << k << endl;
    cerr << "Plesae enter d: ";
    cerr << d << endl;

    if ( d < 5 )
    {
        doBFS( k, d );
    }
    else
    {
        doNN( k, d );
    }

    return 0;
}
