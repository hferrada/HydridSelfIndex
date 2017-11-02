The parameter for the constructor are: 
 1.- parserFile: filename of the LZ77 parser
 2.- M: the M value for the structure
 3.- dirSaveLoad: folder where save/load the structure
  
==========================================================
Another public methods:
 .- We provide also save/load functionality:
 
 // save load structure to/from the base-name 'pathFile', typically a directory with a prefix.
 void saveStructure();
 void loadStructure();
  
 .- Locate occurrences of a pattern 'pat' of length 'm':
 
 // writes in *nOcc the number of occurrences of the pattern *pat[0..m-1] allocating these in **occ.
 void locate(uchar *pat, uint m, ulong *nOcc, ulong **occ);
 
 .- We also offer extraction operation:
 
 // extracts the 'len' characters that represent the original segment T[sp.. sp+len-1] and allocates these into A[0...len-1]
 void extract(ulong sp, ulong len, uchar **A);
  
==========================================================
We also provide two small programs: build_hsi.cpp and load_hsi.cpp; as examples of how you can use this index.
Note that with 'load.cpp' also you can execute experiments to measure the size/locate_time for the index, where the patterns will read from files. The format of the file of patterns is one pattern (string) for line, and all the patterns must to be the same length for each different file: 'patt_m', where ‘m’ is the length of each pattern in the file. The number of patterns also is given as parameter in the main function (see the documentation in these files for more details).
