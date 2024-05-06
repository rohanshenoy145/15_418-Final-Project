# 15_418-Final-Project

Non Thread Safe Union Find:
DisjointSets.cpp
DisjointSets.h

Fine Grained:
ParallelDisjointSets.cpp
ParallelDisjointSets.h

Lock Free:
LockFreeDisjointSets.cpp
LockFreeDisjointSets.cpp


For each of the following, update the #include and
makefile to the union find implementation you would like to test.

Boruvka's Implementations:
MasterWithLockShortest.cpp : Final version
Makefile.lockshortest

LockShortestCoarse.cpp : Final version with critical regions around Union Find
Makefile.coarse

boruvka_ds.cpp : fastest sequential version
Makefile.seq

boruvka_parallel : Shared grid for finding min edges
Compile with make fine, or make lockfree, but still update the include

MasterWithCAS_short.cpp : CAS for finding min edges
Makefile.cas

MasterWithSortScan.cpp: Finding min edges using parallel sort and scan
Makefile.sortscan

MasterWithPrefixSumEnd.cpp : Delay creating MST array until end using prefix sum

MasterWithPrivateArray.cpp : Select MST edges using private arrays and prefix sum

boruvka_contract.cpp : star contraction sequential version

TESTING

graphGenerator.py: create random graphs by changing the input in the file
to the desired number of nodes and edges

validate_script.py: uses NetworkX to validate result, use by updating the file pathes

Run parallel versions with :
./boruvka_parallel -f path/to/file -n NUM_THREADS


