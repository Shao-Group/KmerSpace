# On the Maximal Independent Sets of Strings with Edit Distance

This repository provides source code that implements three algorithms
to find a maximal independent set of all kmers (all strings of fixed length k)
under the edit distance.  Formally, let S be the set of all kmers over alphabet {A, C, G, T}.
Given a parameter d, d < k, a subset M of S is defined as an independent set of S
parameterized by d if the edit distance between any two kmers in M is larger than d.
An independent set M of S is defined to be maximal (i.e. an MIS) if there does not exist
another independet set M' such that M is a strict subset of M'.
Finding MIS has critical applications in hashing/indexing/clustering of all kmers.

We implemented three algorithms for finding an MIS for a given k and d.
1. The first algorithm is an greedy algorithm that iteratively examine if the current kmer
can be added to the maintained MIS. This algorithm is similar to the
greedy algorithm that finds an MIS on a general graph.
2. The second algorithm improves the first one by 
reducing redundant comparisons through recognizing the locality properties of kmers
and estimating the bounds of the edit distance. 
3. The third algorithm transforms calculating the edit distance into finding the shortest path
in a new graph, and finding an MIS is then transformed into efficient graph traversing
together with data structures to speed up.

More details can be found in our manuscript titled "On the Maximal Independent
Sets of Strings with Edit Distance" (available soon).

## Compilation

Please use the following command to compile the code.

```bash
g++ findMIS.cpp -o findMIS -std=c++11
```

## Execution

Please use the following command to run the executable and then proceed according to the instructions displayed on the console.

```bash
./findMIS
```
