# HydridSelfIndex
This code is a 64-bits and extended version of the original Hybrid Index (HI) [1], converting the HI in a compressed full-text index for any types of text (not only repetitive datasets).

The original HI is a compressed text index designed for repetitive texts to locate occurrences of patterns for a maximum pattern length M. This code is an improved version that extends the original functionalities to locate patterns of any length and to reproduce any text segment (the typical extract/display operation). It continues receiving M, as input parameter, as a threshold on pattern length for solving queries for pattern lengths at most M in a similar way as the original index, which is our optimal query time case; when m>M we offer a new method to locate occurrences, which is asymptotically is not optimal, giving a trade-off between locate-time and space usage.

Authors: H. Ferrada, D. Kempa and S. J. Puglisi. 
{{hferrada,dkempa,puglisi}@cs.helsinki.fi}

Requisites
=============
1.- The sdsl library, which is is available from: https://github.com/simongog/sdsl-lite <br />
2.- The RMQ library, which is available from: https://github.com/hferrada/rmq <br />

MAKE 
======
1.- Install the libraries sdsl and rmq. <br />
2.- Edit the Makefile writing the variables SDSL_DIR and RMQ_DIR with the path where you have installed the libraries. <br />
3.- To make the library just execute the command 'make', this will create the file: 'hsi.a'. <br />

Compiling & Linking
=========================
To use the index you must compile your program linking 'hsi.a' and include the header "HybridSelfIndex.h". In the Makefile file there are two tag in order to create/load the index from the examples .cpp files build_hsi.cpp and load_hsi.cpp included here.

Parameters of Construction
==========================
These are also well detailed in the code. The trade-off of our HybridSelfIndex depends on two main parameters: S_SA, which is the sampling-size of the internal FMI and the value for M. A larger M value you will obtain a larger filtered text (i.e., a bigger FMI), but also you will reduce the size of the internal structures to report secondary occurrences (given the LZ77 parser of the text). On the other hand an small M value speeds up locate queries for short patterns (m<=M) and a larger M value speeds up queries for longer patterns (m>M) and also reduce the time to report secondary occurrences.<br />
In order to index a an input text, you first have to compute its LZ77 parser with the code that we have included in the folder “computelz77_hi”. The file generated by the binary “computelz77_hi/computelz77_hi” is the new input file for our constructor in the clase “HybridSelfIndex”. The “HybridSelfIndex” receives only three parameter: <br />
-. Parser File. It is the LZ77 parser computed by computelz77_hi/computelz77_hi. <br />
-. M. The optimal pattern length M value for our index.<br />
-. Prefix of the Save Path. The prefix path where to save/load the files generated by the code.<br />

References
===========
Please, if you want to include this tool as part of a work or experiments, in your references include the paper [2] .<br />
[1]. H. Ferrada, T. Gagie, T. Hirvola, and S. J. Puglisi. Hybrid indexes for repetitive datasets. Philosophical Transactions of the Royal Society A, 327, 2014. Aricle no. 2016. <br />
[2]. H. Ferrada, D. Kempa and S. J. Puglisi. Hybrid indexing Revisited. To appear in Algorithm Engineering and Experiments (ALENEX18). 2018.
