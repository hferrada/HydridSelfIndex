# HydridSelfIndex
This code is a 64-bits and extended version of the original Hybrid Index (HI) [1], converting the HI in a compressed full-text index for any types of text (not only repetitive datasets).

The original HI is a compressed text index designed for repetitive texts to solve pattern matching queries for a maximum pattern length M. This improved version extends its functionalities to locate patterns of any length and to reproduce any text segment (the typical extract/display operation). We continue using M as a threshold on pattern length for solving queries for pattern lengths at most M in a similar way as the original index, which is our optimal query time case, giving a trade-off between locate-time and space usage.

Authors: H. Ferrada, D. Kempa and S. J. Puglisi. 
{hferrada,dkempa,puglisi}@cs.helsinki.fi}

PREREQUISITES
=============
1.- The sdsl library, which is is available from: https://github.com/simongog/sdsl-lite
2.- The RMQ library, which is available from: https://github.com/hferrada/rmq

MAKE 
======
1.- Install the libraries sdsl and rmq.
2.- Edit the Makefile writing the variables SDSL_DIR and RMQ_DIR with the path where you installed the libraries.
3.- To make the library just execute the command 'make', this will create the file: 'hsi.a'.

COMPILE and LINKING
=========================
To use the index you must compile your program linking 'hsi.a' and include the the header "HybridSelfIndex.h". In the Makefile file there are two tag in order to create/load the index from the examples cpp files build_hsi.cpp and load_hsi.cpp included here.

PARAMETERS
============
These are detailed in the code. The trade-off of our HybridSelfIndex depends on two main parameters: S_SA, which is the sampling-size of the internal FMI and the value for M. A larger M value you will obtain a larger filtered text (i.e., a bigger FMI), but also you will reduce the size of the internal structures to report secondary occurrences (given the LZ77 parser of the text). On the other hand an small M value speeds up locate queries for short patterns (m<=M) and a larger M value speeds up queries for longer patterns (m>M) and also reduce the time to report secondary occurrences.

References: The reference for this new index will be soon.
[1]. H. Ferrada, T. Gagie, T. Hirvola, and S. J. Puglisi. Hybrid indexes for repetitive datasets. Philosophical Transactions of the Royal Society A, 327, 2014. Aricle no. 2016.
