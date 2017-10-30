#include "HybridSelfIndex.h"

bool TRACE = false;			// true: print all details for console
bool TEST = 0;			// true: apply exhaustive test
uint N_REP = 100;

typedef struct {
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	uint M;					// maximum pattern length
	char parserFile[400];	// file with the parser
	char prefixStore[300];	// directory to save/load the data structure

	HybridSelfIndex *index;

	ulong* patt;			// data structure are for test and experiments
} ParProg;

// Parameters example
// /home/ferrada/data/documents/alabarda/hi_simple/alabardaXL/parser.lz77 6 /home/ferrada/data/documents/alabarda/hi_simple/alabardaXL/
// /home/ferrada/data/documents/alabarda/hi_simple/alabardaL/parser.lz77 4 /home/ferrada/data/documents/alabarda/hi_simple/alabardaL/
// /home/ferrada/data/repetitive/cere/hi_simple/1MB/parser.lz77 10 /home/ferrada/data/repetitive/cere/hi_simple/1MB/
// /home/ferrada/data/repetitive/kernel/hi_simple/parser.lz77 40 /home/ferrada/data/repetitive/kernel/hi_simple/HI_M40_SA32/
// /home/ferrada/data/repetitive/kernel/hi_simple/1MB/parser.lz77 10 /home/ferrada/data/repetitive/kernel/hi_simple/1MB/HI_M10_SA32/
// /home/ferrada/data/repetitive/kernel/10MB/parser.lz77 20 /home/ferrada/data/repetitive/kernel/hi_simple/10MB/HI_M20_SA32/
// /home/ferrada/data/repetitive/cere/hi_simple/parser.lz77 40 /home/ferrada/data/repetitive/cere/hi_simple/461MiB_M40_SA32/
// /home/ferrada/data/documents/english/englishX/parser.lz77 10 /home/ferrada/data/documents/english/hi_simple/englishX/

// /home/ferrada/data/repetitive/influ/hsi_simple/parserDocs.lz77 40 /home/ferrada/data/repetitive/influ/hsi_simple/M40_SA32/

void testLocate(ParProg *par);
void testExtract(ParProg *par);
void createPatterns(ParProg *par, uint m, bool avoidPattEq);

int main(int argc, char *argv[]) {
	ParProg *par = new ParProg();
	char fileName[400];

	if(argc != 4){
		cout << "ERRORR with parameters!! " << endl;
		cout << "build_hsi_ms's usage requires 3 parameters:" << endl;
		cout << "<file parser> <optimal pattern length M> <save/load prefixStore>" << endl;
		exit(1);
	}

	HybridSelfIndex::TRACE = TRACE;
	HybridSelfIndex::CREATE_FMI_TEST = TEST;
	HybridSelfIndex::SHOW_SIZE = true;

	strcpy(par->parserFile, argv[1]);
	par->M = atoi(argv[2]);
	strcpy(par->prefixStore, argv[3]);

	cout << "build_hsi_ms' parameters..." << endl;
	cout << "Parser File               : " << par->parserFile << endl;
	cout << "Optimal pattern length, M : " << par->M << endl;
	cout << "Store Folder (prefix path): " << par->prefixStore << endl;
	cout << "FMI(T') Sampling; S_SA    : " << S_SA << ", S_ISA: " << S_ISA << endl << endl;

	// Building gthe Index ...
	par->index = new HybridSelfIndex(par->parserFile, par->M, par->prefixStore);

	// saving the index ...
	par->index->saveStructure();
	par->index->~HybridSelfIndex();

	if(TEST && HybridSelfIndex::CREATE_FMI_TEST){
		par->index = new HybridSelfIndex(par->prefixStore);
		par->n = par->index->n;

		// load FMI_TEST ...
		strcpy(fileName, "");
		strcpy(fileName, par->prefixStore);
		strcat(fileName, "fmi_test.test");
		cout << "____________________________________________________" << endl;
		cout << " Loading the FMI_TEST from: " << fileName << endl;

		load_from_file(par->index->FMI_TEST, fileName);
		ulong sizeFMI = size_in_bytes(par->index->FMI_TEST);
		cout << " **  FMI_TEST size " << sizeFMI << " Bytes = " << (float)sizeFMI/(float)par->n << "|T| = " << (float)sizeFMI*8.0/(float)par->n << " bpc" << endl;
		cout << " **  FMI length " << par->index->FMI_TEST.size() << endl;
		if (par->index->FMI_TEST.size() != par->n+1){
			cout << "ERROR. FMI length != n+1 = " << par->n+1 << endl;
			exit(1);
		}

		string substr = sdsl::extract(par->index->FMI_TEST, 0, 10);
		cout << "T[0..10] : " << substr << endl;

		cout << "Running test locate().." << endl;
		testLocate(par);
		cout << "Test locate() OK !!" << endl;

		if (S_ISA <= 1024){
			cout << "Running test extract().." << endl;
			testExtract(par);
			cout << "Test extract() OK !!" << endl;
		}else
			cout << " Take into account that Testing extract() will take a lot of time given the sapling-size of the inverse SA!" << endl;

		par->index->~HybridSelfIndex();
	}

	cout << "$$$$$$$$$$$$$$$$$$$$$" << endl;
	return 0;
}

void testExtract(ParProg *par){
	uchar *A;
	ulong i, sp, t;
	uint m, Max = 4*par->index->M;
	string substr;

	for (m=5; m<Max; m*=2){
		cout << " ... Extracting patterns for m = " << m << endl;
		for (t=0; t<N_REP; t++){
			sp = rand() % (par->n-m-2);
			substr = sdsl::extract(par->index->FMI_TEST, sp, sp+m-1);
			//cout << "t = " << t << ", m = " << m << " [" << substr << "]" << endl;

			par->index->extract(sp, m, &A);
			for(i=0; i<m; i++){
				if ((int)(substr[i]) != (int)(A[i]) && (256+(int)(substr[i]) != (int)(A[i]))){
					cout << "ERROR: m=" << m << ", test=" << t << ", sp = " << sp << endl;
					cout << "substr = [" << substr << "]" << endl;
					cout << "Different symbol in position " << i << ", substr[i] = [" << substr[i] << "] != A[i] = [" << A[i] << "]" << endl;
					cout << "(int)(substr[i]) = " << (int)(substr[i]) << ", (int)(A[i]) = " << (int)(A[i]) << endl;
					exit(-2);
				}
			}
			delete [] A;
		}
	}
}

void testLocate(ParProg *par){
	uint m;
	ulong i, t, nOccHI, *occ, nOccTest, *A = nullptr;
	//ulong phra;
	bit_vector B_TEST(par->n, 0);

	par->patt = new ulong[N_REP];
	//TRACE = 1;

	for (m=par->M/3; m<=3*par->M; m+=par->M/3){
	//for (m=2*par->M+1; m<=4*par->M; m+=5){
	//for (m=10; m<=4*par->M; m+=5){
		createPatterns(par, m, true);
		uchar* pat = new uchar[m+1];
		cout << " ... Locating patterns for m = " << m << endl;
		for (t=0; t<N_REP; t++){
			//cout << "t =" << t << ", m = " << m << endl;
			//if (t==2666)
			//	TRACE = 1;//cout<<"";

			string query = sdsl::extract(par->index->FMI_TEST, par->patt[t], par->patt[t]+m-1);
			size_t occs = sdsl::count(par->index->FMI_TEST, query.begin(), query.begin()+m);
			nOccTest=0;
			if (occs){
				A = new ulong[occs];
				auto locations = locate(par->index->FMI_TEST, query.begin(), query.begin()+m);
				if (TRACE){
					cout << "# of occurrences found with FMI_TEST(T) : " << occs << " =" << endl;
					for(i=0; i<occs; i++)
						cout << locations[i] << " ";
					cout << endl;
				}
				//cout << "Total occurrences found with FMI = " << occs << endl;
				for(i=0; i<occs; i++){
					//if (par->index->isPrimaryT(locations[i], m, &phra)){
						B_TEST[locations[i]] = 1;
						A[nOccTest] = locations[i];
						nOccTest++;
						//cout << i << ", " << locations[i] << endl;
					//}
				}
				//cout << endl;
			}

			strncpy((char*) pat, query.c_str(), m);
			//cout << "patt = [";
			//for(uint j=0; j<m; j++)
			//	cout << pat[j];
			//cout << "]" << endl;

			nOccHI=0;
			occ = nullptr;
			par->index->locate(pat, m, &nOccHI, &occ);

			if(nOccTest != nOccHI){
				cout << "ERROR: m=" << m << ", test=" << t << ", occs with FMI_TEST = " << nOccTest << " != nOcc = " << nOccHI << endl;
				cout << "patt = [";
				for(uint j=0; j<m; j++)
					cout << pat[j];
				cout << "]" << endl;
				for(i=0; i<nOccHI; i++){
					if(B_TEST[occ[i]] == 0)
						cout << occ[i] << " not found with the FMI_TEST !" << endl;
					else
						B_TEST[occ[i]] = 0;
				}

				for(i=0; i<nOccTest; i++){
					if(B_TEST[A[i]] == 1)
						cout << A[i] << " not found with the Hybrid Index !" << endl;
				}
				exit(-1);
			}

			for(i=0; i<nOccHI; i++){
				if(B_TEST[occ[i]] == 0){
					cout << "ERROR: m=" << m << ", test=" << t << endl;
					cout << "patt = [";
					for(uint j=0; j<m; j++)
						cout << pat[j];
					cout << "]" << endl;
					cout << occ[i] << " not found with the FMI_TEST !" << endl;
					exit(-1);
				}
			}

			for(i=0; i<nOccHI; i++)
				B_TEST[occ[i]] = 0;

			if (occs){
				delete [] occ;
				delete [] A;
			}
		}
	}

	delete [] (par->patt);
}

void createPatterns(ParProg *par, uint m, bool avoidPattEq){
	ulong i, j, k;
	par->patt = new ulong[N_REP];
	bool foundI;
	int LOWER = 32;
	int BIGER = 255;
	string str;

	//cout << "creating random patterns of legth " << m << endl;
	if (m==1){
		for (k=0; k<N_REP; k++){
			j = rand() % (par->n-2);
			string str = sdsl::extract(par->index->FMI_TEST, j, j+1);
			while ((int)(str[0]) < LOWER || (int)(str[0]) > BIGER)
				j = rand() % (par->n-2);

			par->patt[k] = j;
		}
		return;
	}

	if (avoidPattEq){
		for(k=0; k<N_REP; k++){
			foundI = false;
			while (!foundI){
				i = rand() % (par->n-(m+1));
				str = sdsl::extract(par->index->FMI_TEST, i, i+m);

				for (j=0; j<m; j++){ ;
					if ((int)str[j] < LOWER || (int)str[j] > BIGER)
						break;
				}

				if (j==m){
					if (m <= par->M){
						for (j=1; !foundI && j<m; j++){
							if (str[j-1] != str[j])
								foundI = true;
						}
					}else{
						for (j=1; !foundI && j<par->M; j++){
							if (str[j-1] != str[j])
								foundI = true;
						}
						if (foundI){
							foundI = false;
							for (j=m-par->M; !foundI && j<m; j++){
								if (str[j-1] != str[j])
									foundI = true;
							}
						}
					}
				}
			}
			par->patt[k] = i;
		}
	}else{
		for(k=0; k<N_REP; k++){
			i = rand() % (par->n-(m+1));
			par->patt[k] = i;
		}
	}
}
