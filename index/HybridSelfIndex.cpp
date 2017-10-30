#include "HybridSelfIndex.h"

bool HybridSelfIndex::TRACE = false;
bool HybridSelfIndex::SHOW_SIZE = true;
uint HybridSelfIndex::SAMPGC = SAMP_GC;
uint HybridSelfIndex::POT_GC = log(SAMP_GC)/log(2);

bool HybridSelfIndex::CREATE_FMI_TEST = false;
uchar HybridSelfIndex::MIN_CHAR = '\1';
uchar HybridSelfIndex::REPLACE_CHAR = '\n';

HybridSelfIndex::HybridSelfIndex(char dirSaveLoad[300]){
	double tLoad = getTime_ms();
	{
		strcpy(dirStore, dirSaveLoad);
		this->sizeDS = 0;

		loadStructure();
		//testPredecessor();
	}
	tLoad = getTime_ms() - tLoad;
	cout << "====================================================" << endl;
	cout << " ### Hybrid Self-Index Loaded in = " << tLoad/(1000.0*60.0) << " minutes." << endl;
	cout << "====================================================" << endl << endl;

	cout << "====================================================" << endl;
	cout << " ### Index size " << sizeDS << " bytes = " << sizeDS/(1024.0*1024.0) << " MiB = " << (float)sizeDS*8.0/(float)n << " bpc" << endl;
	cout << "====================================================" << endl << endl;
}

HybridSelfIndex::HybridSelfIndex(char *parserFile, uint optM, char dirSaveLoad[300]) {
	double tConst = getTime_ms();
	{
		strcpy(dirStore, dirSaveLoad);
		this->M = optM;
		this->sizeDS = 0;

		// retrieving parsing from file
		cout << " Reading... " << parserFile << endl;
		readParsing(parserFile);

		this->lgN = 1 + (uint)(log(n+1)/log(2));
		this->lgM = 1 + (uint)(log(M)/log(2));

		ulong *Source = nullptr;
		buildBasicStructures(&Source);

		if (nG)
			buildPtrRMQonG(Source);
		else
			cout << "WARNING.... Length of the grid nG = 0 !!" << endl;
	}

	tConst = getTime_ms() - tConst;

	cout << "====================================================" << endl;
	cout << " ### Hybrid Self-Index built in = " << tConst/(1000.0*60.0) << " minutes." << endl;
	cout << "====================================================" << endl << endl;

	cout << "====================================================" << endl;
	cout << " ### Index size " << sizeDS << " bytes = " << sizeDS/(1024.0*1024.0) << " MiB = " << (float)sizeDS*8.0/(float)n << " bpc" << endl;
	cout << "====================================================" << endl << endl;
}

void HybridSelfIndex::buildBasicStructures(ulong **Source){
	ulong i, j, k, cont, posX, posXF, antX, antXF, lenS, posPar, posFT, nz, lenArray;
	bool found, prevSM;
	ulong *alpha = new ulong[SIGMA];

	for (i=0; i<SIGMA; i++)
		alpha[i] = n;

	j = k = posPar = 0;
	z = nFT = 0;
	uint maxM = 2*M;
	prevSM=false;
	for(i=0; i<orig_z; i++, j+=lgPOS, k+=lgLEN){
		posX = getNum64(ARR_POS, j, lgPOS);
		lenS = getNum64(ARR_LEN, k, lgLEN);

		if(lenS){
			posPar += lenS;
			if(lenS > maxM){
				nFT+=maxM+1;
				z++;
				prevSM=false;
			}else{
				if (prevSM==false){
					z++;
					prevSM=true;
				}
				nFT+=lenS;
			}
		}else{
			posPar++;
			nFT++;
			if (prevSM==false){
				z++;
				prevSM=true;
			}
		}
	}
	if(n != posPar){
		cout << "ERROR! n = " << n << " != final position of the parser = " << posPar << endl;
		exit(1);
	}
	z++;
	this->lgZ = 1 + (uint)(log(z)/log(2));
	cout << "Length of the filtered text fTxt, nFT = " << nFT << endl;

	bit_vector BL(z, 0);

	// WE REPLACE ANY SIMBOL WHOSE ASCII CODE <= MIN_CHAR BY 'REPLACE_CHAR'
	Txt = new uchar[n+1];		// original text
	j = k = posPar = nz = 0;
	for(i=0; i<orig_z; i++, j+=lgPOS, k+=lgLEN){
		posX = getNum64(ARR_POS, j, lgPOS);
		lenS = getNum64(ARR_LEN, k, lgLEN);

		if(lenS){
			if(lenS == 1){
				if (Txt[posX] <= MIN_CHAR)
					Txt[posPar] = REPLACE_CHAR;
				else
					Txt[posPar] = Txt[posX];
				posPar++;
			}else{
				for(cont=0; cont<lenS; cont++){
					if (Txt[posX+cont] <= MIN_CHAR)
						Txt[posPar] = REPLACE_CHAR;
					else
						Txt[posPar] = Txt[posX+cont];
					posPar++;
				}
			}
		}else{
			if (posX <= MIN_CHAR)
				Txt[posPar] = REPLACE_CHAR;
			else
				Txt[posPar] = (uchar)posX;
			posPar++;
			if (i+1 > orig_z){
				cout << "ERROR.. phrase = " << i+1 << " > orig_z = " << orig_z << endl;
				exit(1);
			}
		}
	}
	Txt[n] = 0;

	// look for the best Mismatch symbol
	found = false;
	for (i=0; i<n; i++){
		if (alpha[Txt[i]] > i)
			alpha[Txt[i]] = i;
	}
	if (alpha['#'] == n){
		// first preference, symbol '#'
		fSymb = '#';
		found = true;
	}else{
		// second preference, symbol '$'
		if (alpha['$'] == n){
			fSymb = '$';
			found = true;
		}else{
			// third preference, a symbol in the ASCII range [32..126]
			found = false;
			for (i=32; (i<=126) && (found==false); i++){
				if (alpha[i] == n){
					found = true;
					fSymb = (uchar)i;
				}
			}

			// Ford preference, a symbol in the ASCII range [128..254]
			for (i=128; (i<=254) && (found==false); i++){
				if (alpha[i] == n){
					found = true;
					fSymb = (uchar)i;
				}
			}

			// Last preference, a symbol in the ASCII range [2..31]
			for (i=31; i && found==false; i--){
				if (alpha[i] == n){
					found = true;
					fSymb = (uchar)i;
				}
			}
		}
	}
	if (!found){
		cout << "Symbol separator not found !!" << endl;
		exit(0);
	}

	cout << "Original orig_z = " << orig_z << ", new z = " << z << ", lgZ = " << lgZ << endl;
	if (TRACE && false){
		cout << "Separator symbol: code = " << (uint)fSymb << ", symbol = " << fSymb << endl;
		cout << endl << ".- T[] = " << endl << Txt << endl;
		cout << "Alphabet: " << endl;
		for (i=0; i<SIGMA; i++)
			cout << alpha[i] << " ";
		cout << endl;
	}
	delete [] alpha;

	// Build the filtered text fTxt[1..nFT] from the original parsing
	ulong *AUX_PhraT, *AUX_PhraFT;
	lenArray = z*lgN / (8*sizeof(ulong));
		if ((z*lgN) % (8*sizeof(ulong)))
			lenArray++;
	AUX_PhraT = new ulong[lenArray];
	AUX_PhraFT = new ulong[lenArray];
	*Source = new ulong[lenArray];

	if (nFT){
		prevSM = false;
		fTxt = new uchar[nFT+1];
		j = k = posPar = posFT = nz = 0;
		for(i=0; i<orig_z; i++, j+=lgPOS, k+=lgLEN){
			posX = getNum64(ARR_POS, j, lgPOS);
			lenS = getNum64(ARR_LEN, k, lgLEN);

			if(lenS){
				if(lenS == 1){
					if (prevSM==false){
						setNum64(*Source, nz*lgN, lgN, posX);
						setNum64(AUX_PhraT, nz*lgN, lgN, posPar);
						setNum64(AUX_PhraFT, nz*lgN, lgN, posFT);
						BL[nz]=1;
						nz++;
						prevSM = true;
					}

					if (Txt[posX] > MIN_CHAR)
						fTxt[posFT] = Txt[posX];
					posPar++;
					posFT++;
				}else{
					if(lenS<=maxM){
						if (prevSM==false){
							setNum64(*Source, nz*lgN, lgN, posX);
							setNum64(AUX_PhraT, nz*lgN, lgN, posPar);
							setNum64(AUX_PhraFT, nz*lgN, lgN, posFT);
							BL[nz]=1;
							nz++;
							prevSM = true;
						}

						for(cont=0; cont<lenS; cont++){
							if (Txt[posX+cont] <= MIN_CHAR)
								fTxt[posFT] = REPLACE_CHAR;
							else
								fTxt[posFT] = Txt[posX+cont];
							posFT++;
						}
					}else{
						setNum64(*Source, nz*lgN, lgN, posX);
						setNum64(AUX_PhraT, nz*lgN, lgN, posPar);
						setNum64(AUX_PhraFT, nz*lgN, lgN, posFT);
						nz++;
						prevSM = false;

						for(cont=0; cont<M; cont++){			// copy M symbol + '$' + M symbols...
							if (Txt[posX+cont] <= MIN_CHAR)
								fTxt[posFT] = REPLACE_CHAR;
							else
								fTxt[posFT] = Txt[posX+cont];
							posFT++;
						}

						fTxt[posFT] = fSymb;
						posFT++;

						for(cont=lenS-M; cont<lenS; cont++){
							if (Txt[posX+cont] <= MIN_CHAR)
								fTxt[posFT] = REPLACE_CHAR;
							else
								fTxt[posFT] = Txt[posX+cont];
							posFT++;
						}
					}
					posPar+=lenS;
				}
			}else{
				if (prevSM==false){
					setNum64(*Source, nz*lgN, lgN, posPar);
					setNum64(AUX_PhraT, nz*lgN, lgN, posPar);
					setNum64(AUX_PhraFT, nz*lgN, lgN, posFT);
					BL[nz]=1;
					nz++;
					prevSM = true;
				}

				if (posX <= MIN_CHAR)
					fTxt[posFT] = REPLACE_CHAR;
				else
					fTxt[posFT]= (uchar)posX;
				posPar++;
				posFT++;
			}
		}

		fTxt[nFT] = 0;
		setNum64(*Source, nz*lgN, lgN, n-1);
		setNum64(AUX_PhraFT, nz*lgN, lgN, nFT);
		setNum64(AUX_PhraT, nz*lgN, lgN, n);
	}

	BL_il = bit_vector_il<512>(BL);
	rankBL = bit_vector_il<512>::rank_1_type(&BL_il);
	sizeDS += size_in_bytes(BL_il);

	// To Create the FMI and FMI_TEST
	buildFMI();

	/*cout << "AUX_PhraT = " << endl;
	for(i=z-20; i<z; i++)
		cout << getNum64(AUX_PhraT, i*lgN, lgN) << " ";
	cout << endl;
	cout << "BL" << endl;
	for(i=z-10; i<z; i++)
		cout << BL_il[i];
	cout << endl;
	cout << endl << ".- fTxt[] = " << endl << fTxt << endl;
	exit(0);*/
	cout << "Length of the Filtered text, posFT = " << posFT << endl;

	if (TRACE){
		cout << endl << ".- fTxt[] = " << endl << fTxt << endl;

		cout << endl << ".- Txt[] = " << endl << Txt << endl;

		cout << "AUX_PhraT = " << endl;
		for(i=k=0; i<z; i++, k+=lgN)
			cout << getNum64(AUX_PhraT, k, lgN) << " ";
		cout << endl;

		cout << "AUX_PhraFT = " << endl;
		for(i=k=0; i<z; i++, k+=lgN)
			cout << getNum64(AUX_PhraFT, k, lgN) << " ";
		cout << endl;

		cout << "Source = " << endl;
		for(i=k=0; i<z; i++, k+=lgN)
			cout << getNum64(*Source, k, lgN) << " ";
		cout << endl;

		cout << "BL" << endl;
		for(i=0; i<z; i++)
			cout << BL_il[i];
		cout << endl;

		cout << "rankBL = " << endl;
		for(i=1; i<=z; i++)
			cout << rankBL.rank(i);
		cout << endl;
	}

	delete [] Txt;
	delete [] fTxt;

	//===========================================================================================================================
	this->nG = z;
	for(i=0; i<z; i++){
		if (BL[i])
			nG--;
	}
	lgNG = 1 + (uint)(log(nG)/log(2));
	nSamP = z/SAMPGC;
	if (z % SAMPGC)
		nSamP++;

	ulong maxLenT, maxLenFT;
	maxLenT = maxLenFT = 0;
	for (i=1; i<z; i++){
		posX = getNum64(AUX_PhraT, i*lgN, lgN) - getNum64(AUX_PhraT, (i-1)*lgN, lgN);
		if (posX > maxLenT)
			maxLenT = posX;

		posXF = getNum64(AUX_PhraFT, i*lgN, lgN) - getNum64(AUX_PhraFT, (i-1)*lgN, lgN);
		if (posXF > maxLenFT)
			maxLenFT = posXF;
	}
	lgPT = 1 + (uint)(log(maxLenT)/log(2));
	lgPFT = 1 + (uint)(log(maxLenFT)/log(2));

	lenArray = z*lgPT / (8*sizeof(ulong));
	if ((z*lgPT) % (8*sizeof(ulong)))
		lenArray++;
	PhraT = new ulong[lenArray];
	sizeDS += lenArray*sizeof(ulong);

	lenArray = z*lgPFT / (8*sizeof(ulong));
	if ((z*lgPFT) % (8*sizeof(ulong)))
		lenArray++;
	PhraFT = new ulong[lenArray];
	sizeDS += lenArray*sizeof(ulong);

	lenArray = nSamP*lgN / (8*sizeof(ulong));
	if ((nSamP*lgN) % (8*sizeof(ulong)))
		lenArray++;
	SGCPT = new ulong[lenArray];
	SGCPFT = new ulong[lenArray];
	sizeDS += 2*lenArray*sizeof(ulong);

	antX = getNum64(AUX_PhraT, 0, lgN);
	antXF = getNum64(AUX_PhraFT, 0, lgN);
	for (i=k=0; i<z; i++){
		posX = getNum64(AUX_PhraT, i*lgN, lgN);
		setNum64(PhraT, i*lgPT, lgPT, posX-antX);
		antX = posX;

		posXF = getNum64(AUX_PhraFT, i*lgN, lgN);
		setNum64(PhraFT, i*lgPFT, lgPFT, posXF-antXF);
		antXF = posXF;

		if (i%SAMPGC == 0){
			setNum64(SGCPT, k, lgN, antX);
			setNum64(SGCPFT, k, lgN, antXF);
			k += lgN;
		}
	}
	delete [] AUX_PhraT;
	delete [] AUX_PhraFT;

	cout << "Length of the grid, nG = " << nG << ", lgNG = " << lgNG << endl;
	if (SHOW_SIZE){
		cout << " ** size of BL_il " << size_in_bytes(BL_il) << " Bytes = " << size_in_bytes(BL_il)*8.0/(float)this->n << " bpc" << endl;

		lenArray = z*lgPT / (8*sizeof(ulong));
		if ((z*lgPT) % (8*sizeof(ulong)))
			lenArray++;
		cout << " ** size of PhraT[0..z] " << lenArray*sizeof(ulong) << " Bytes = " << lenArray*sizeof(ulong)*8.0/(float)this->n << " bpc" << endl;

		lenArray = z*lgPFT / (8*sizeof(ulong));
		if ((z*lgPFT) % (8*sizeof(ulong)))
			lenArray++;
		cout << " ** size of PhraFT[0..z] " << lenArray*sizeof(ulong) << " Bytes = " << lenArray*sizeof(ulong)*8.0/(float)this->n << " bpc" << endl;

		lenArray = nSamP*lgN / (8*sizeof(ulong));
			if ((nSamP*lgN) % (8*sizeof(ulong)))
				lenArray++;
		cout << " ** size of table SGCPT and SGCPFT " << 2*lenArray*sizeof(ulong) << " Bytes = " << 2*lenArray*sizeof(ulong)*8.0/(float)this->n << " bpc" << endl;
	}

	if (TRACE){
		cout << "lgPT = " << lgPT << ", lgPFT = " << lgPFT << endl;
		cout << "# of sampling: nSamP = " << nSamP << endl;

		cout << "PhraT[0.." << z << "] = " << endl;
		for(i=k=0; i<z; i++, k+=lgPT)
			cout << getNum64(PhraT, k, lgPT) << " ";
		cout << endl;

		cout << "PhraFT[0.." << z << "] = " << endl;
		for(i=k=0; i<z; i++, k+=lgPFT)
			cout << getNum64(PhraFT, k, lgPFT) << " ";
		cout << endl;

		cout << "SGCPT[0.." << nSamP-1 << "] = " << endl;
		for(i=k=0; i<nSamP; i++, k+=lgN)
			cout << getNum64(SGCPT, k, lgN) << " ";
		cout << endl;

		cout << "SGCPFT[0.." << nSamP-1 << "] = " << endl;
		for(i=k=0; i<nSamP; i++, k+=lgN)
			cout << getNum64(SGCPFT, k, lgN) << " ";
		cout << endl;
	}
}

// Giving newX: source of the phrase i = Q[k], we move 'newX' go up until its correct position in the partial minimum queue Q[1..k]
void HybridSelfIndex::swimInMinQX(ulong *Q, ulong newX, ulong i, ulong k, ulong *Source){
	ulong x, pX;

	pX = getNum64(Q, (k/2)*lgZ, lgZ);
	x = getNum64(Source, pX*lgN, lgN);

	while(k > 1 && x > newX){
		setNum64(Q, k*lgZ, lgZ, pX);
		k /= 2;

		pX = getNum64(Q, (k/2)*lgZ, lgZ);
		x = getNum64(Source, pX*lgN, lgN);
	}

	setNum64(Q, k*lgZ, lgZ, i);	// 'i' represents the source newX
}

// set Q[k] as top of the queue Q[1..k-1] and move down it until its correct position
void HybridSelfIndex::setTopMinQX(ulong *Q, uint k, ulong *Source){
	ulong val, pVal, x, pX, x2, pX2, m;

	pVal = getNum64(Q, k*lgZ, lgZ);
	val = getNum64(Source, pVal*lgN, lgN);

	m = 2;
	pX = getNum64(Q, m*lgZ, lgZ);
	x = getNum64(Source, pX*lgN, lgN);

	if (m+1 < k){
		pX2 = getNum64(Q, (m+1)*lgZ, lgZ);
		x2 = getNum64(Source, pX2*lgN, lgN);

		if (x2 < x){
			m++;
			x=x2;
			pX=pX2;
		}
	}

	while(m < k && val > x){
		setNum64(Q, (m/2)*lgZ, lgZ, pX);
		m <<= 1;

		if (m < k){
			pX = getNum64(Q, m*lgZ, lgZ);
			x = getNum64(Source, pX*lgN, lgN);
			if (m+1 < k){
				pX2 = getNum64(Q, (m+1)*lgZ, lgZ);
				x2 = getNum64(Source, pX2*lgN, lgN);

				if (x2 < x){
					m++;
					x=x2;
					pX=pX2;
				}
			}
		}
	}

	setNum64(Q, (m/2)*lgZ, lgZ, pVal);
}

// It Builds the grid G(X, Y, Pr) sorted by x-coordinate using a queue.
void HybridSelfIndex::buildPtrRMQonG(ulong *Source){
	ulong i, k, c, posX, lenS, lenArray, r;

	ulong *QUEUE_G;
	lenArray = (nG+1)*lgZ / (8*sizeof(ulong));
	if (((nG+1)*lgZ) % (8*sizeof(ulong)))
		lenArray++;
	QUEUE_G = new ulong[lenArray];
	setNum64(QUEUE_G, 0, lgZ, 0);

	posX=0;
	c=1;
	for(i=k=0; i<z; i++, k+=lgN){
		//	WE ONLY INCLUDE SOURCES WHAT NOT GENERATE LITERAL PHRASES (ie. BL_il[i] == false)
		if(BL_il[i]==0){
			posX = getNum64(Source, k, lgN);
			setNum64(QUEUE_G, c*lgZ, lgZ, i);
			swimInMinQX(QUEUE_G, posX, i, c, Source);
			c++;
		}
	}

	// for the X array
	lenArray = nG*lgN / (8*sizeof(ulong));
	if ((nG*lgN) % (8*sizeof(ulong)))
		lenArray++;
	X = new ulong[lenArray];
	sizeDS += lenArray*sizeof(ulong);

	// Y array for the RMQ...
	long int *Y = new long int[nG];

	// Pr array
	lenArray = nG*lgZ / (8*sizeof(ulong));
	if ((nG*lgZ) % (8*sizeof(ulong)))
		lenArray++;
	Pr = new ulong[lenArray];
	sizeDS += lenArray*sizeof(ulong);

	//=============================================================================================
	{
		// The additional Src array to solve queries for any pattern length
		lenArray = nG*lgNG / (8*sizeof(ulong));
		if ((nG*lgNG) % (8*sizeof(ulong)))
			lenArray++;
		Src = new ulong[lenArray];
		sizeDS += lenArray*sizeof(ulong);
	}
	//=============================================================================================

	for (i=0, c=nG; i<nG; i++, c--){
		k = getNum64(QUEUE_G, lgZ, lgZ); // phrase from the top of the queue
		setNum64(Pr, i*lgZ, lgZ, k);

		r = k - rankBL.rank(k+1);
		setNum64(Src, r*lgNG, lgNG, i);

		posX = getNum64(Source, k*lgN, lgN);
		setNum64(X, i*lgN, lgN, posX);

		if (k==z-1)
			lenS=1;
		else
			lenS = getNum64(PhraT, (k+1)*lgPT, lgPT);
		Y[i] = posX + lenS - 1;

		if (c > 1)
			setTopMinQX(QUEUE_G, c, Source);	// c is the current length of the queue
	}
	delete [] QUEUE_G;
	delete [] Source;

	//	LOOKUP..
	uint blk;
	ulong blk_X = posX/N_LOOKUP;
	b_X = (uint)(log(blk_X)/log(2));
	nSampX = 1+(posX>>b_X);

	cout << "Final length lookup Table nSampX = " << nSampX << ", bits per segment b_X = " << b_X << endl;

	LOOKUP_X = new uint[nSampX];
	COUNT_X = new uint[nSampX];
	sizeDS += 2*nSampX*sizeof(uint);

	for (i=0; i<nSampX; i++)
		LOOKUP_X[i]=nG+1;

	for (i=0; i<nG; i++){
		posX = getNum64(X, i*lgN, lgN);
		blk = posX>>b_X;

		if (LOOKUP_X[blk] > i)
			LOOKUP_X[blk] = i;
	}

	for (i=0; i<nSampX; i++){
		if (LOOKUP_X[i] > nG)
			LOOKUP_X[i] = LOOKUP_X[i-1];
	}

	for (i=0; i<nSampX-1; i++){
		if (LOOKUP_X[i+1] < LOOKUP_X[i]){
			cout << "ERR.. i = " << i << ", LOOKUP_X[i+1] = " << LOOKUP_X[i+1] << " < LOOKUP_X[i] = " << LOOKUP_X[i] << endl;
			exit(0);
		}
		COUNT_X[i] = LOOKUP_X[i+1] - LOOKUP_X[i];
	}
	COUNT_X[i] = nG - LOOKUP_X[i];

	for (i=0; i<nSampX-1; i++){
		if (COUNT_X[i]==0){
			for (c=i+1; c<nSampX && COUNT_X[c]==0; c++);
			if (c<nSampX && COUNT_X[c])
				COUNT_X[i]=COUNT_X[c];
		}
	}

	float avgCount = 0;
	uint maxCC = 0;
	for (i=0; i<nSampX; i++){
		avgCount += COUNT_X[i];
		if (maxCC < COUNT_X[i]) maxCC = COUNT_X[i];
	}
	avgCount = avgCount/nSampX;
	cout << "maxCC = " << maxCC << ", avgCount = " << avgCount << endl;

	testPredecessor();

	if (TRACE){
		ulong ph;

		cout << "b_X = " << b_X << endl;
		cout << "blk_X = " << blk_X << endl;

		cout << "LOOKUP_X[0.." << nSampX-1 << "](count) = " << endl;
		for(i=0; i<nSampX; i++)
			cout << LOOKUP_X[i] << "(" << COUNT_X[i] << ") ";
		cout << endl;

		cout << "Pr[0.." << nG-1 << "] = " << endl;
		for(i=k=0; i<nG; i++, k+=lgZ)
			cout << getNum64(Pr, k, lgZ) << " ";
		cout << endl;

		cout << "(X, Y) -- > PrX " << endl;
		for(i=k=0; i<nG; i++, k+=lgN){
			ph = getNum64(Pr, i*lgZ, lgZ);
			posX = getPosPhraT(ph);
			cout << i << "(" << getNum64(X, k, lgN) << "," << Y[i] << ") --> " << posX << endl;
		}
		cout << endl;
	}

	for(i=0; i<nG; i++)
		Y[i] = -1*Y[i];

	// Build a RMQ on YG[]
	rmqY = new RMQRMM64(Y, nG);
	sizeDS += rmqY->getSize();

	delete [] Y;

	if (SHOW_SIZE){
		lenArray = nG*lgN / (8*sizeof(ulong));
		if ((nG*lgN) % (8*sizeof(ulong)))
			lenArray++;
		cout << " ** size of X[]  = " << lenArray << " Bytes = " << lenArray*8.0/(float)this->n << " bpc" << endl;

		cout << " ** size of LOOKUP_X[] + COUNT_X[]  = " << 2*nSampX*sizeof(uint) << " Bytes = " << 2.0*nSampX*sizeof(uint)*8.0/(float)this->n << " bpc" << endl;

		lenArray = nG*lgZ / (8*sizeof(ulong));
		if ((nG*lgZ) % (8*sizeof(ulong)))
			lenArray++;
		cout << " ** size of Pr[1..nG] " << lenArray*sizeof(ulong) << " Bytes = " << lenArray*sizeof(ulong)*8.0/(float)this->n << " bpc" << endl;

		lenArray = nG*lgNG / (8*sizeof(ulong));
		if ((nG*lgNG) % (8*sizeof(ulong)))
			lenArray++;
		cout << " ** size of Src[1..nG] " << lenArray*sizeof(ulong) << " Bytes = " << lenArray*sizeof(ulong)*8.0/(float)this->n << " bpc" << endl;

		cout << " ** size of rmqY: " << rmqY->getSize() << " = " << rmqY->getSize()*8.0/(float)n << " bpc" << endl;
	}

	if (TRACE){
		cout << "X[0.." << nG-1 << "] = " << endl;
		for(i=k=0; i<nG; i++, k+=lgN)
			cout << getNum64(X, k, lgN) << " ";
		cout << endl;

		if (false){
			cout << "predForX[0.." << n-1 << "] = " << endl;
			for(i=0; i<n; i++){
				if (findPredecesor(i, &k))
					cout << i << ":" << k << " ";
			}
			cout << endl;
		}

		cout << "Src[0.." << nG-1 << "] = " << endl;
		for(i=k=0; i<nG; i++, k+=lgNG)
			cout << getNum64(Src, k, lgNG) << " ";
		cout << endl;
	}
}
//======================================================================================================================================

void HybridSelfIndex::buildFMI(){
	char fileName[400];

	if (nFT){
		cout << " Build the FMI of fTxt[1..nFT]..." << endl;

		string theFilTxt = (char *)fTxt;
		construct_im(FMI, theFilTxt, 1);
		sizeDS += size_in_bytes(FMI);
		cout << "     FMI length " << FMI.size() << endl;

		// ===================
		//sdsl::write_structure<sdsl::JSON_FORMAT>(FMI, std::cout);
		// ===================

		if (FMI.size() != nFT+1){
			cout << "ERROR. FMI.size() = " << FMI.size() << " != nFT+1 = " << nFT+1 << endl;
			exit(1);
		}
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "fmi.hsi");
		if (store_to_file(FMI, fileName))
			cout << " **  FMI size " << size_in_bytes(FMI) << " bytes = " << size_in_bytes(FMI)/(1024.0*1024.0) << " MiB = " << (float)size_in_bytes(FMI)*8.0/(float)n << " bpc" << endl;
		else{
			cout << " ERROR STORING THE FMI IN: " << fileName<< endl;
			exit(0);
		}
	}

	if (CREATE_FMI_TEST){
		cout << "Creating the FMI_TEST..." << endl;
		string theTxt = (char *)Txt;
		construct_im(FMI_TEST, theTxt, 1);

		if (FMI_TEST.size() != n+1){
			cout << "ERROR. FMI_TEST.size() = " << FMI_TEST.size() << " != n+1 = " << n+1 << endl;
			exit(1);
		}
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "fmi_test.test");

		if (store_to_file(FMI_TEST, fileName)){
			if (TRACE){
				cout << " **  FMI_TEST size " << size_in_bytes(FMI_TEST) << " bytes = " << (float)size_in_bytes(FMI_TEST)/(float)n << "|T| = " << (float)size_in_bytes(FMI_TEST)*8.0/(float)n << " bpc" << endl;
				cout << "     FMI_TEST length " << FMI_TEST.size() << endl;
				cout << "saving the FMI_TEST in " << fileName << " ..." << endl;
			}
		}else{
			cout << " ERROR STORING THE FULL FMI FOR TEST!! " << endl;
			exit(0);
		}
	}
}

// Read the modified LZ77 parsing, in which we have forced to that every phrase do not split in two documents.
void HybridSelfIndex::readParsing(char *parserFile){
	ulong i, j, k, pos, len, x;
	ifstream is(parserFile, ios::binary);

	if (is.is_open()==false){
		cout << " File not found !!" << endl;
		exit(0);
	}

	is.read((char*)&n, sizeof(ulong));
	is.read((char*)&orig_z, sizeof(ulong));
	is.read((char*)&lgPOS, sizeof(uint));
	is.read((char*)&lgLEN, sizeof(uint));

	ulong size_ARR = orig_z*lgPOS/(8*sizeof(ulong));		// number of 'ulong' segments for ARR_POS
	if ((orig_z*lgPOS)%(8*sizeof(ulong)))
		size_ARR++;
	ARR_POS = new ulong[size_ARR];
	is.read((char*)ARR_POS, size_ARR*sizeof(ulong));

	size_ARR = orig_z*lgLEN/(8*sizeof(ulong));
	if ((orig_z*lgLEN)%(8*sizeof(ulong)))
		size_ARR++;
	ARR_LEN = new ulong[size_ARR];
	is.read((char*)ARR_LEN, size_ARR*sizeof(ulong));
	is.close();

	cout << " Text with n = " << n << endl;
	if(TRACE){
		cout << " Original Dictionary with orig_z = " << orig_z << endl;

		cout << "Dictionary: " << endl;
		pos = getNum64(ARR_POS, 0, lgPOS);
		len = getNum64(ARR_LEN, 0, lgLEN);

		if(len)
			cout << "0:  (" << pos << "," << len<< ") --> 0" << endl;
		else
			cout << "0:  (" << pos << ", * ) --> 0" << endl;
		j=lgPOS;
		k=lgLEN;
		x=1;
		for(i=1; i<orig_z; i++, j+=lgPOS, k+=lgLEN){
			pos = getNum64(ARR_POS, j, lgPOS);
			len = getNum64(ARR_LEN, k, lgLEN);
			if(len){
			//	if (i>=orig_z-20 || x==1048379)
				cout << i << ":  (" <<pos << "," << pos+len-1 << ") --> " << x << endl;
				x += len;
			}else{
			//	if (i>=orig_z-20 || x==1048379)
				cout << i << ":  (" <<pos << ", * ) --> " << x << endl;
				x++;
			}
		}
		cout << endl;
	}
}

//======================================================================================================
void HybridSelfIndex::testPredecessor(){
	cout << "Testing Predecessor of i ..." << endl;
	ulong j, pred, lastVal, x;

	if (nG>100){
		cout << "X[0..99] = " << endl;
		for(j=0; j<100; j++){
			cout << getNum64(X, j*lgN, lgN) << " ";
			if ((j+1) % 20 == 0)
				cout << endl;
		}
		cout << endl;
	}

	lastVal = getNum64(X, (nG-1)*lgN, lgN);
	for (j=0, x=getNum64(X, 0, lgN); x<lastVal ; x++){
		for(; getNum64(X, (j+1)*lgN, lgN) <= x; j++);

		//if (x == 65536)
		//	cout << "";

		if (findPredecesor(x, &pred)==false || pred != j){
			cout << "ERR. predecessor of " << x << " = " << pred << " Not found or it is different than its real position pred = " << j << endl;
			exit(1);
		}
	}

	findPredecesor(lastVal, &pred);
	if (pred != nG-1){
		cout << "ERR. predecessor of " << lastVal << " = " << pred << " != real pos pred = " << nG-1 << endl;
		exit(1);
	}

	cout << "Test Predecessor OK !!" << endl;
}

// return in '*pred' the greater position pos, such as X[pos] <= x. If x<A[0] then we return false
bool HybridSelfIndex::findPredecesor(ulong x, ulong *pred){
	if (x < getNum64(X, 0, lgN))
		return false;

	//if (x==17907990)
	//	cout<<"";

	uint r, j, pos = x>>b_X;

	if (pos>=nSampX){
		for (r=nG-1, j=r*lgN; getNum64(X, j, lgN) > x; j-=lgN)
			r--;
		*pred = r;
		return true;
	}

	uint c = COUNT_X[pos];
	if (c > N_SCAN){
		// Binary searching in the X-segment...
		ulong u,v;
		uint m, l = LOOKUP_X[pos];

		r = l+c-1;
		m = r>>1;
		if (l) l--;
		while (l<=r){
			u = getNum64(X, m*lgN, lgN);
			if (u<=x){
				if (m<r){
					m++;
					v = getNum64(X, m*lgN, lgN);
					if (v>x){
						*pred = m-1;
						return true;
					}else{
						l=m;
						m=(l+r)>>1;
					}
				}else{
					*pred = m;
					return true;
				}
			}else{
				if (l<m){
					m--;
					v = getNum64(X, m*lgN, lgN);
					if (v<=x){
						*pred = m;
						return true;
					}else{
						r=m;
						m=(l+r)>>1;
					}
				}else
					return false;
			}
		}

		return false;
	}else{
		// Scanning a X-segment of maximum c cells...

		r = LOOKUP_X[pos]+c-1;
		for (j=r*lgN; getNum64(X, j, lgN) > x; j-=lgN)
			r--;
		*pred = r;
	}
	return true;
}

bool HybridSelfIndex::isPrimaryLeftLarge(ulong x, uint len, ulong *phra, uint *dx){
	ulong l=0, r=nSamP-1, m=nSamP>>1, x1, x2, phr;

	while (l<=r){
		x1 = getNum64(SGCPFT, m*lgN, lgN);
		if (x1 <= x){
			if (m<r){
				x2 = getNum64(SGCPFT, (m+1)*lgN, lgN);
				if (x2 > x){
					phr = m<<POT_GC;
					break;
				}else{
					l=m+1;
					m=(l+r)>>1;
				}
			}else{
				phr = m<<POT_GC;
				break;
			}
		}else {
			m--;
			if (l<=m){
				x1 = getNum64(SGCPFT, m*lgN, lgN);
				if (x1 <= x){
					phr = m<<POT_GC;
					break;
				}else{
					r=m-1;
					m=(l+r)>>1;
				}
			}else {
				phr = m<<POT_GC;
				break;
			}
		}
	}

	// retrieve gaps from x1 to x
	for(r = (phr+1)*lgPFT; x1<=x; r+=lgPFT){
		phr++;
		x1 += getNum64(PhraFT, r, lgPFT);
	}
	*phra = phr;
	*dx = x1-x;

	if (x+len < x1)
		return false;

	return true;
}

// pEnd is the last phrase boundary in T_{x..x+len-1}, and dEnd is the distance
bool HybridSelfIndex::isPrimaryLeft_b(ulong x, uint len, ulong *pIni, uint *dIni){
	ulong l=0, r=nSamP-1, m=nSamP>>1, x1, x2, phr;

	while (l<=r){
		x1 = getNum64(SGCPFT, m*lgN, lgN);
		if (x1 <= x){
			if (m<r){
				x2 = getNum64(SGCPFT, (m+1)*lgN, lgN);
				if (x2 > x){
					phr = m<<POT_GC;
					break;
				}else{
					l=m+1;
					m=(l+r)>>1;
				}
			}else{
				phr = m<<POT_GC;
				break;
			}
		}else {
			m--;
			if (l<=m){
				x1 = getNum64(SGCPFT, m*lgN, lgN);
				if (x1 <= x){
					phr = m<<POT_GC;
					break;
				}else{
					r=m-1;
					m=(l+r)>>1;
				}
			}else {
				phr = m<<POT_GC;
				break;
			}
		}
	}

	// retrieve gaps from x1 to x
	for(r = (phr+1)*lgPFT; x1<=x; phr++, r+=lgPFT)
		x1 += getNum64(PhraFT, r, lgPFT);

	if (x+len <= x1){
		*pIni = phr;
		*dIni = x1-x;
		return false;
	}

	m = getNum64(PhraFT, r, lgPFT);
	if (x+len > x1+m){
		*pIni = phr+1;
		*dIni = x1-x+m;
	}else{
		*pIni = phr;
		*dIni = x1-x;
	}

	return true;
}

bool HybridSelfIndex::isPrimaryMidle(ulong x, uint len, ulong *phra, uint *dx){
	ulong l=0, r=nSamP-1, m=nSamP>>1, x1, x2, phr;

	while (l<=r){
		x1 = getNum64(SGCPFT, m*lgN, lgN);
		if (x1 <= x){
			if (m<r){
				x2 = getNum64(SGCPFT, (m+1)*lgN, lgN);
				if (x2 > x){
					phr = m<<POT_GC;
					break;
				}else{
					l=m+1;
					m=(l+r)>>1;
				}
			}else{
				phr = m<<POT_GC;
				break;
			}
		}else {
			m--;
			if (l<=m){
				x1 = getNum64(SGCPFT, m*lgN, lgN);
				if (x1 <= x){
					phr = m<<POT_GC;
					break;
				}else{
					r=m-1;
					m=(l+r)>>1;
				}
			}else {
				phr = m<<POT_GC;
				break;
			}
		}
	}

	// retrieve gaps from x1 to x
	for(r = (phr+1)*lgPFT; x1<x; phr++, r+=lgPFT)
		x1 += getNum64(PhraFT, r, lgPFT);

	*phra = phr;
	if (x1 == x)
		*dx = 0;
	else {
		*dx = x1-x;
		if (x+len <= x1)
			return false;
	}

	return true;
}

// return true if in T'_{x..x+len-1} there is a phrase boundary (including the position x)
bool HybridSelfIndex::isPrimaryRight(ulong x, uint len, ulong *phra, uint *dx){
	ulong l=0, r=nSamP-1, m=nSamP>>1, x1, x2, phr;

	while (l<=r){
		x1 = getNum64(SGCPFT, m*lgN, lgN);
		if (x1 <= x){
			if (m<r){
				x2 = getNum64(SGCPFT, (m+1)*lgN, lgN);
				if (x2 > x){
					phr = m<<POT_GC;
					break;
				}else{
					l=m+1;
					m=(l+r)>>1;
				}
			}else{
				phr = m<<POT_GC;
				break;
			}
		}else {
			m--;
			if (l<=m){
				x1 = getNum64(SGCPFT, m*lgN, lgN);
				if (x1 <= x){
					phr = m<<POT_GC;
					break;
				}else{
					r=m-1;
					m=(l+r)>>1;
				}
			}else {
				phr = m<<POT_GC;
				break;
			}
		}
	}

	// retrieve gaps from x1 to x
	for(r = (phr+1)*lgPFT; x1<x; phr++, r+=lgPFT)
		x1 += getNum64(PhraFT, r, lgPFT);

	*phra = phr;
	if (x==x1)
		*dx = 0;
	else{
		*dx = x1-x;
		if (x+len <= x1)
			return false;
	}

	return true;
}

// return the position of the phrase boundary in the text Txt
ulong HybridSelfIndex::getPosPhraT(ulong phra){
	ulong ph = phra>>POT_GC;
	ulong x = getNum64(SGCPT, ph*lgN, lgN);
	if (phra%SAMP_GC == 0)
		return x;

	// extract from the sampled phrase 'ph' to 'phra'...
	ph=(ph<<POT_GC)+1;
	for (ulong c = ph*lgPT; ph<=phra; ph++, c+=lgPT)
		x += getNum64(PhraT, c, lgPT);

	return x;
}

// return the position of the phrase boundary in the filtered text fTxt
ulong HybridSelfIndex::getPosPhraFT(ulong phra){
	ulong ph = phra>>POT_GC;
	ulong x = getNum64(SGCPFT, ph*lgN, lgN);
	if (phra%SAMP_GC == 0)
		return x;

	// extract from the sampled phrase 'ph' to 'phra'...
	ph=(ph<<POT_GC)+1;
	for (ulong c = ph*lgPFT; ph<=phra; ph++, c+=lgPFT)
		x += getNum64(PhraFT, c, lgPFT);

	return x;
}

// =============================================================================================================
// 'pIni' is the first Id-phrase inside of fTxt_{x..x+len-1}
// 'dIni' is the distance from fTxt_x to the 'pIni'
bool HybridSelfIndex::isPrimary(ulong x, uint len, ulong *pIni, uint *dIni){
	ulong l=0, r=nSamP-1, m=nSamP>>1, x1, x2, phr;

	while (l<=r){
		x1 = getNum64(SGCPFT, m*lgN, lgN);
		if (x1 <= x){
			if (m<r){
				x2 = getNum64(SGCPFT, (m+1)*lgN, lgN);
				if (x2 > x){
					phr = m<<POT_GC;
					break;
				}else{
					l=m+1;
					m=(l+r)>>1;
				}
			}else{
				phr = m<<POT_GC;
				break;
			}
		}else {
			m--;
			if (l<=m){
				x1 = getNum64(SGCPFT, m*lgN, lgN);
				if (x1 <= x){
					phr = m<<POT_GC;
					break;
				}else{
					r=m-1;
					m=(l+r)>>1;
				}
			}else {
				phr = m<<POT_GC;
				break;
			}
		}
	}

	// retrieve gaps from x1 to x
	for(r = (phr+1)*lgPFT; x1<x; phr++, r+=lgPFT)
		x1 += getNum64(PhraFT, r, lgPFT);

	if (x1 == x){
		*dIni = 0;
		*pIni = phr;
		x1 += getNum64(PhraFT, r, lgPFT);
	}else {
		*pIni = phr;
		*dIni = x1-x;
	}

	if (x+len <= x1)
		return false;

	return true;
}

void HybridSelfIndex::locatePryOcc(uchar *pat, uint m, ulong *nOcc, ulong **L, ulong* currN){
	string query = string((char *)pat);
	size_t nLoc = sdsl::count(FMI, query.begin(), query.begin()+m);

	if (nLoc){
		ulong i, pEnd, nn, *A;
		uint dEnd;

		*L = A = new ulong[nLoc];
		*currN = nLoc;
		auto list = sdsl::locate(FMI, query.begin(), query.begin()+m);
		for(i=nn=0; i<nLoc; i++){
			if (isPrimary(list[i], m, &pEnd, &dEnd)){
				//A[nn] = PfD_PT->extract(pEnd) - dEnd;
				A[nn] = getPosPhraT(pEnd) - dEnd;
				nn++;
				//cout << i << " --> " << A[nn-1] << endl;
			}
		}
		*nOcc = nn;
	}else
		*nOcc=0;
}

void HybridSelfIndex::locateSecOcc(ulong l, ulong r, ulong posX, uint m, ulong *nOcc, ulong **occ, ulong *currN){
	ulong pos, x, y, phra, xPr;

	pos = rmqY->queryRMQ(l,r);
	x = getNum64(X, pos*lgN, lgN);
	phra = getNum64(Pr, pos*lgZ, lgZ);	// Phrase of the pointer of the source
	y = x+getNum64(PhraT, (phra+1)*lgPT, lgPT); // = x+len

	if (y >= posX+m){
		ulong nn = *nOcc;
		ulong currLen = *currN;
		ulong *A = *occ;

		xPr = getPosPhraT(phra)+posX-x;
		//cout << " sec --> " << xPr << endl;
		if (nn < currLen)
			A[nn] = xPr;
		else{
			ulong j, *AUX;

			currLen <<= 1;
			AUX = new ulong[currLen];	// We duplicate the size of the array *occ
			for (j=0; j<nn; j++)
				AUX[j] = A[j];
			AUX[j] = xPr;
			delete [] A;
			*occ = AUX;
			A = AUX;
			*currN = currLen;
		}
		(*nOcc)++;

		 // recursive call for this new occurrence...
		if (findPredecesor(xPr, &x))
			locateSecOcc(0, x, xPr, m, nOcc, occ, currN);

		// we continue searching occurrences covering all the interval [l..r]
		if (pos>l)
			locateSecOcc(l, pos-1, posX, m, nOcc, occ, currN);
		if (pos<r)
			locateSecOcc(pos+1, r, posX, m, nOcc, occ, currN);
	}
}

ulong HybridSelfIndex::searchPhraFilTxt(ulong x, uint *dx){
	ulong l=0, r=nSamP-1, m=nSamP>>1, x1, x2, phr;

	while (l<=r){
		x1 = getNum64(SGCPFT, m*lgN, lgN);
		if (x1 <= x){
			if (m<r){
				x2 = getNum64(SGCPFT, (m+1)*lgN, lgN);
				if (x2 > x){
					phr = m<<POT_GC;
					break;
				}else{
					l=m+1;
					m=(l+r)>>1;
				}
			}else{
				phr = m<<POT_GC;
				break;
			}
		}else {
			m--;
			if (l<=m){
				x1 = getNum64(SGCPFT, m*lgN, lgN);
				if (x1 <= x){
					phr = m<<POT_GC;
					break;
				}else{
					r=m-1;
					m=(l+r)>>1;
				}
			}else {
				phr = m<<POT_GC;
				break;
			}
		}
	}

	// retrieve gaps from x1 to x
	for(r = (phr+1)*lgPFT; x1<=x; phr++, r+=lgPFT)
		x1 += getNum64(PhraFT, r, lgPFT);
	*dx = x1-x;

	return phr;
}

// return the phrase-id where x is, and writes in (pos,len) its position and length
ulong HybridSelfIndex::searchPhraTxt(ulong x, ulong *pos, uint *len){
	ulong l=0, r=nSamP-1, m=nSamP>>1, x1, x2, phr;

	while (l<=r){
		x1 = getNum64(SGCPT, m*lgN, lgN);
		if (x1 <= x){
			if (m<r){
				x2 = getNum64(SGCPT, (m+1)*lgN, lgN);
				if (x2 > x){
					phr = 1+(m<<POT_GC);
					break;
				}else{
					l=m+1;
					m=(l+r)>>1;
				}
			}else{
				phr = 1+(m<<POT_GC);
				break;
			}
		}else {
			m--;
			if (l<=m){
				x1 = getNum64(SGCPT, m*lgN, lgN);
				if (x1 <= x){
					phr = 1+(m<<POT_GC);
					break;
				}else{
					r=m-1;
					m=(l+r)>>1;
				}
			}else {
				phr = 1+(m<<POT_GC);
				break;
			}
		}
	}

	// retrieve gaps from x1 to x
	for(r = phr*lgPT; x1<=x; phr++, r+=lgPT){
		l = getNum64(PhraT, r, lgPT);
		x1 += l;
	}
	*len = l;
	*pos = x1-l;

	return phr-2;
}

// =============================================================================================================
// LOCATE METHODS
// =============================================================================================================
void HybridSelfIndex::locate(uchar *pat, uint m, ulong *nOcc, ulong **occ){
	if (m==1)
		locateAChar(pat, nOcc, occ);
	else{
		if (m<=M)
			locateUptoM(pat, m, nOcc, occ);
		else{
			if (m<=(M<<1))
				locateUpto2M(pat, m, nOcc, occ);
			else
				locateLargeM_b(pat, m, nOcc, occ);
		}
	}
}

void HybridSelfIndex::locateAChar(uchar *pat, ulong *nOcc, ulong **occ){
	pat[1]='\0';
	string query = string((char *)pat);
	size_t nLoc = sdsl::count(FMI, query.begin(), query.begin()+1);

	if (nLoc){
		ulong r, currN, i, pr, nn, *A;
		uint dx;

		*occ = A = new ulong[nLoc];
		currN = nLoc;
		auto list = sdsl::locate(FMI, query.begin(), query.begin()+1);
		for(i=nn=0; i<nLoc; i++){
			pr = searchPhraFilTxt(list[i], &dx);
			if (BL_il[pr-1]){
				A[nn] = getPosPhraT(pr) - dx;
				nn++;
				//cout << i << " p*** --> " << A[nn-1] << endl;
			}
		}

		// search secondary occurrence from literal ones...
		*nOcc = nn;
		for(i=0; i<nn; i++){
			if (findPredecesor(A[i], &r)){
				locateSecOcc(0, r, A[i], 1, nOcc, occ, &currN);
				A = *occ;
			}
		}
	}else
		*nOcc=0;
}

void HybridSelfIndex::locateUptoM(uchar *pat, uint m, ulong *nOcc, ulong **occ){
	string query = string((char *)pat);
	auto list = sdsl::locate(FMI, query.begin(), query.begin()+m);
	size_t nLoc = list.size();

	if (nLoc){
		ulong r, currN, i, pr, nn, *A;
		uint dx;

		*occ = A = new ulong[nLoc];
		currN = nLoc;
		for(i=nn=0; i<nLoc; i++){
			if (isPrimary(list[i], m, &pr, &dx)){
				A[nn] = getPosPhraT(pr) - dx;
				nn++;
				//cout << i << " pri --> " << A[nn-1] << endl;
			}else{
				if (dx){
					if (BL_il[pr-1]){
						A[nn] = getPosPhraT(pr) - dx;
						nn++;
						//cout << i << " pri_2 --> " << A[nn-1] << endl;
					}
				}else{
					if (BL_il[pr]){
						A[nn] = getPosPhraT(pr);
						nn++;
						//cout << i << " pri_3* --> " << A[nn-1] << endl;
					}
				}
			}
		}

		// search secondary occurrence from primary ones...
		*nOcc = nn;
		for(i=0; i<nn; i++){
			if (findPredecesor(A[i], &r)){
				locateSecOcc(0, r, A[i], m, nOcc, occ, &currN);
				A = *occ;
			}
		}
	}else
		*nOcc=0;
}

// put the last item in the top and sort the queue
void HybridSelfIndex::setTopMinQ(ulong *Q, ulong pos){
	ulong j, l, u;

	j=1;
	u = Q[pos];
	l=2;
	while(l<pos){
		if (l+1 < pos && Q[l] > Q[l+1])
			l++;
		if(u > Q[l]){
			Q[j] = Q[l];
			j = l;
			l <<= 1;
		}else
			break;
	}
	Q[j] = u;
}

// return true if u is in A[0..len]
bool HybridSelfIndex::isOccInArr(ulong u, ulong len, ulong *A){
	ulong l, r ,m;

	l=0;r=len-1;
	while (l<=r){
		m=(r+l)>>1;
		if (A[m]==u)
			return true;

		if (A[m]<u)
			l=m+1;
		else{
			if (m)
				r=m-1;
			else
				return false;
		}
	}

	return false;
}

// M < m <= 2M
void HybridSelfIndex::locateUpto2M(uchar *pat, uint m, ulong *nOcc, ulong **occ){
	ulong *QL, *QR, i, j, l, x;
	string query = string((char *)pat);
	size_t nLocL, nLocR;

	auto listL = sdsl::locate(FMI, query.begin(), query.begin()+M);
	nLocL = listL.size();

	if (!nLocL){
		*nOcc = 0;
		return;
	}

	query = string((char *)pat+m-M);
	auto listR = sdsl::locate(FMI, query.begin(), query.begin()+M);
	nLocR = listR.size();

	if (!nLocR){
		*nOcc = 0;
		return;
	}

	QL = new ulong[nLocL+1];
	j=1;
	//auto list = sdsl::locate(FMI, query.begin(), query.begin()+M);
	QL[1]=listL[0];
	for(i=1; i<nLocL; i++){
		j=i+1;
		x = listL[i];
		l=j>>1;
		while(l && QL[l]>x){
			QL[j] = QL[l];
			j=l;
			l>>=1;
		}
		QL[j]=x;
	}
	/*cout << "QL = ";
	for(i=1; i<=nLocL; i++)
		cout << QL[i] << " ";
	cout << endl;*/

	QR = new ulong[nLocR+1];
	j=1;
	//auto list = sdsl::locate(FMI, query.begin(), query.begin()+M);
	QR[1]=listR[0];
	for(i=1; i<nLocR; i++){
		j=i+1;
		x = listR[i];
		l=j>>1;
		while(l && QR[l]>x){
			QR[j] = QR[l];
			j=l;
			l>>=1;
		}
		QR[j]=x;
	}
	/*cout << "QR = ";
	for(i=1; i<=nLocR; i++)
		cout << QR[i] << " ";
	cout << endl;*/


	ulong *A, *L, *R, nn, nl, nr, ii, pI, y, u, pos, currN;
	uint dx, len, lenS = m-M;		// distance between x and y
	bool didX, didY, procX, procY;

	L = new ulong[nLocL];
	R = new ulong[nLocR];
	currN = nLocL + nLocR;
	*occ = A = new ulong[currN];
	nl = nr = nn = 0;

	x = QL[1];	// extract from QL
	setTopMinQ(QL, nLocL);
	y = QR[1];	// extract from QR
	setTopMinQ(QR, nLocR);
	i=ii=1;
	didX = didY = procX = procY = false;

	while(true){
		if (x+lenS==y){
			// to check if this is a primary occurrence and report it
			if (isPrimaryLeft_b(x, m, &pI, &dx) || BL_il[pI-1]){
				A[nn] = L[nl] = getPosPhraT(pI) - dx;
				R[nr] = L[nl] + lenS;
				nn++;
				//cout << x << ", eq --> " << A[nn-1] << endl;
			}else{
				L[nl] = getPosPhraT(pI) - dx;
				R[nr] = L[nl] + lenS;
			}
			nl++; nr++;
			//cout << " L " << L[nl-1] << endl;
			//cout << " R " << R[nr-1] << endl;
			didX = didY = true;
		}else{
			if (x+lenS<y)
				procX = true;
			else
				procY = true;
		}

		if (procX){
			if (isPrimaryLeft_b(x, M, &pI, &dx)){
				L[nl] = x = getPosPhraT(pI)-dx;
				nl++;
				//cout << " L " << L[nl-1] << endl;

				if (nr){
					if (dx<lenS && BL_il[pI]==false){
						// go to the source if y is in a normal phrase (ie., not in T')
						l = pI - rankBL.rank(pI);
						l = getNum64(Src, l*lgNG, lgNG);		// rank index in G of the source of the phrase j
						u = getNum64(X, l*lgN, lgN) + lenS - dx;

						pI = searchPhraTxt(u, &pos, &len);
						dx = u-pos;
						while (len >= M+dx && BL_il[pI]==false){
							l = pI - rankBL.rank(pI);
							l = getNum64(Src, l*lgNG, lgNG);		// rank index in G of the source of the phrase j
							u = getNum64(X, l*lgN, lgN) + dx;

							pI = searchPhraTxt(u, &pos, &len);	 	// phrase for the k index in G.
							dx = u-pos;
						}
					}else
						u=x+lenS;
					// search u in R ........... cambiar a busqueda binaria entre 0 y nr
					if(isOccInArr(u, nr, R)){
						A[nn] = x;
						nn++;
						//cout << "x --> " << A[nn-1] << endl;
					}
				}
			}else{
				L[nl] = x = getPosPhraT(pI)-dx;
				nl++;
				//cout << " L " << L[nl-1] << endl;

				if (nr && BL_il[pI-1]){ // x no es primaria pero esta en una literal... entonces y debe estar en R[]
					u=x+lenS;
					// search u in R ........... cambiar a busqueda binaria entre 0 y nr
					if(isOccInArr(u, nr, R)){
						A[nn] = x;
						nn++;
						//cout << "x --> " << A[nn-1] << endl;
					}
				}
			}
			didX = true;
			procX = false;
		}else{
			if (procY){
				if (isPrimaryRight(y, M, &pI, &dx)){
					R[nr] = y = getPosPhraT(pI)-dx;
					nr++;
					//cout << " R " << R[nr-1] << endl;

					if (pI && nl){
						pI--;
						if (lenS+dx>=M && BL_il[pI]==false){
							l = pI - rankBL.rank(pI);
							l = getNum64(Src, l*lgNG, lgNG);		// rank index in G of the source of the phrase j
							len = getNum64(PhraT, (pI+1)*lgPT, lgPT);//PfD_PT->get_ithGap(pI+1);
							u = getNum64(X, l*lgN, lgN)+len-dx-lenS;

							pI = searchPhraTxt(u, &pos, &len);	 	// phrase for the k index in G.
							dx = u-pos;
							while (len >= M+dx && BL_il[pI]==false){
								l = pI - rankBL.rank(pI);
								l = getNum64(Src, l*lgNG, lgNG);		// rank index in G of the source of the phrase j
								u = getNum64(X, l*lgN, lgN) + dx;

								pI = searchPhraTxt(u, &pos, &len);	 	// phrase for the k index in G.
								dx = u-pos;
							}
						}else
							u=y-lenS;
						// search u in L ..
						if(isOccInArr(u, nl, L)){
							A[nn] = y -lenS;
							nn++;
							//cout << "y --> " << A[nn-1] << endl;
						}
					}
				}else{
					R[nr] = y = getPosPhraT(pI)-dx;
					nr++;
					//cout << " R " << R[nr-1] << endl;

					if (nl && BL_il[pI-1]){
						u=y-lenS;
						if(isOccInArr(u, nl, L)){
							A[nn] = y-lenS;
							nn++;
							//cout << "y --> " << A[nn-1] << endl;
						}
					}
				}
				didY = true;
				procY = false;
			}
		}

		if (didX){
			if (i==nLocL){
				if (ii<nLocR && didY){
					y = QR[1];	// extract from QR
					didY = false;
					setTopMinQ(QR, nLocR-ii);
					ii++;
				}
				break;
			}
			x = QL[1];	// extract from QL
			didX = false;
			setTopMinQ(QL, nLocL-i);
			i++;
		}
		if (didY){
			if (ii==nLocR) break;
			y = QR[1];	// extract from QR
			didY = false;
			setTopMinQ(QR, nLocR-ii);
			ii++;
		}
	}

	if (nr){
		while (!didX){		// process only x values...
			if (isPrimaryLeft_b(x, M, &pI, &dx)){
				x = getPosPhraT(pI)-dx;
				if (dx<lenS && BL_il[pI]==false){
					// go to the source if y is in a normal phrase (ie., not in T')
					l = pI - rankBL.rank(pI);
					l = getNum64(Src, l*lgNG, lgNG);		// rank index in G of the source of the phrase j
					u = getNum64(X, l*lgN, lgN)+lenS-dx;

					pI = searchPhraTxt(u, &pos, &len);	 	// phrase for the k index in G.
					dx = u-pos;
					while (len >= M+dx && BL_il[pI]==false){
						l = pI-rankBL.rank(pI);
						l = getNum64(Src, l*lgNG, lgNG);		// rank index in G of the source of the phrase j
						u = getNum64(X, l*lgN, lgN)+dx;

						pI = searchPhraTxt(u, &pos, &len);	 	// phrase for the k index in G.
						dx = u-pos;
					}
				}else
					u=x+lenS;

				// search u in R ........... cambiar a busqueda binaria entre 0 y nr
				if(isOccInArr(u, nr, R)){
					A[nn] = x;
					nn++;
					//cout << "only x --> " << A[nn-1] << endl;
				}
			}

			if (i < nLocL){
				x = QL[1];	// extract from QL
				setTopMinQ(QL, nLocL-i);
				i++;
			}else break;
		}
	}

	if (nl){
		while (!didY){	// process only y values...
			if (isPrimaryRight(y, M, &pI, &dx)){
				if (pI){
					x=getPosPhraT(pI)-dx-lenS;
					if (lenS+dx>=M && BL_il[pI-1]==false){
						y=pI;
						pI--;
						l = pI-rankBL.rank(pI);
						l = getNum64(Src, l*lgNG, lgNG);		// rank index in G of the source of the phrase j
						len = getNum64(PhraT, (pI+1)*lgPT, lgPT);//PfD_PT->get_ithGap(pI+1);
						u = getNum64(X, l*lgN, lgN)+len-dx-lenS;

						pI = searchPhraTxt(u, &pos, &len);	 	// phrase for the k index in G.
						dx = u-pos;
						while (len >= M+dx && BL_il[pI]==false){	// ** no solo buscare primarias en las listas despues
							l = pI-rankBL.rank(pI);
							l = getNum64(Src, l*lgNG, lgNG);		// rank index in G of the source of the phrase j
							u = getNum64(X, l*lgN, lgN)+dx;

							pI = searchPhraTxt(u, &pos, &len);	 	// phrase for the k index in G.
							dx = u-pos;
						}
					}else
						u=x;

					// search u in L ..
					if(isOccInArr(u, nl, L)){
						A[nn] = x;
						nn++;
						//cout << "only y --> " << A[nn-1] << endl;
					}
				}
			}

			if (ii < nLocR){
				y = QR[1];	// extract from QR
				setTopMinQ(QR, nLocR-ii);
				ii++;
			}else break;
		}
	}

	delete [] QL;
	delete [] QR;
	delete [] L;
	delete [] R;

	// search secondary occurrence from primary ones...
	*occ = A;
	*nOcc = nn;
	for(i=0; i<nn; i++){
		x = A[i];
		if (findPredecesor(A[i], &j)){
			locateSecOcc(0, j, A[i], m, nOcc, occ, &currN);
			A = *occ;
		}
	}
}

bool HybridSelfIndex::isOccInArr(ulong u, ulong len, ulong *A, ulong *pos){
	ulong l, r ,m;

	l=0;r=len-1;
	while (l<=r){
		m=(l+r)>>1;
		if (A[m]==u){
			*pos = m;
			return true;
		}

		if (A[m]<u)
			l=m+1;
		else{
			if (m)
				r=m-1;
			else
				return false;
		}
	}

	return false;
}


void HybridSelfIndex::locateLargeM_b(uchar *pat, uint m, ulong *nOcc, ulong **occ){
	PriOcc *LL;
	uint ini, k, eps, dxI, *Dx;
	ulong i, ii, j, l, x, u, pos, pI, nn, currN, *L, *Phr, *Pos, *Q;
	size_t nLocIni, nLoc, nLocFin;
	string query = string((char *)pat);

	k = m/M;			// number of segments in which we split the pattern
	if (m%M){
		k++;
		eps = M-m%M;
	}else
		eps = 0;
	//cout << "k = " << k << ", eps = " << eps << endl;

	// looking for primary occurrences in the first segment...
	//nLoc = sdsl::count(FMI, query.begin(), query.begin()+M);
	auto listIni = sdsl::locate(FMI, query.begin(), query.begin()+M);
	nLocIni = listIni.size();
	if (!nLocIni){
		*nOcc = 0;
		return;
	}

	// looking for primary occurrences in the last segment...
	ini = m-M;
	//nLoc = sdsl::count(FMI, query.begin()+ini, query.begin()+m);
	auto listFin = sdsl::locate(FMI, query.begin()+ini, query.begin()+m);
	nLocFin = listFin.size();
	if (!nLocFin){
		*nOcc = 0;
		return;
	}

	LL = new PriOcc[k];

	// looking for primary occurrences in the last segment...
	{
		LL[0].L = L = new ulong[nLocIni];
		LL[0].Phr = Phr = new ulong[nLocIni];
		LL[0].Pos = Pos = new ulong[nLocIni];
		LL[0].Dx = Dx = new uint[nLocIni];
		LL[0].B = new bool[nLocIni];
		Q = new ulong[nLocIni+1];

		Q[1]=listIni[0];
		for(i=1; i<nLocIni; i++){
			// insertInQueue(list[i])
			j=i+1;
			x = listIni[i];
			l=j>>1;
			while(l && Q[l]>x){
				Q[j] = Q[l];
				j=l;
				l>>=1;
			}
			Q[j]=x;
		}

		/*if (TRACE){ // ELIMINAR AL FINAL
			cout << "A. Total occurrences found with the FMI : " << nLoc << " =" << endl;
			for(i=0; i<nLoc; i++)
				cout << listIni[i] << " ";
			cout << endl;

			cout << "Q =" << endl;
			for(i=1; i<=nLoc; i++)
				cout << Q[i] << " ";
			cout << endl;
		}*/

		for(i=nn=0; i<nLocIni; i++){
			x = Q[1];							// extract from Q
			setTopMinQ(Q, nLocIni-i);
			//cout << x << " ";

			if (isPrimaryLeftLarge(x, M, &pI, &dxI) || BL_il[pI-1]){
				Pos[nn] = getPosPhraT(pI);
				L[nn] =  Pos[nn] - dxI;
				Phr[nn] = pI-1;
				Dx[nn] = dxI;
				nn++;
				//cout << x << " A--> " << L[nn-1] << endl;
			}
		}
		//cout << endl;
		LL[0].n = nn;
		currN = nn;
		for (i=0; i<nn; i++)
			LL[0].B[i]=0;

		/*j=0;
		for (i=1; i<nn; i++){
			if (LL[j].L[i] <= LL[j].L[i-1]){
				cout << "Err C L[i] " << LL[j].L[i] << " <= L[i-1]" << LL[j].L[i-1] << ", i=" << i << ", j=" << j << endl;
				exit(0);
			}
		}*/

		delete [] Q;
	}

	// looking for primary occurrences in intermediate segments...
	for (j=1; j<k-1; j++){
		ini = j*M;
		auto list = sdsl::locate(FMI, query.begin()+ini, query.begin()+ini+M);
		nLoc = list.size();
		if (nLoc == 0){
			*nOcc = 0;
			return;
		}
		LL[j].L = L = new ulong[nLoc];
		LL[j].Phr = Phr = new ulong[nLoc];
		LL[j].Pos = Pos = new ulong[nLoc];
		LL[j].Dx = Dx = new uint[nLoc];
		LL[j].B = new bool[nLoc];
		Q = new ulong[nLoc+1];
		ii=1;

		/*if (TRACE){ // ELIMINAR AL FINAL
			cout << "B. Total occurrences found with the FMI : " << nLoc << " =" << endl;
			for(i=0; i<nLoc; i++)
				cout << list[i] << " ";
			cout << endl;
		}*/
		Q[1]=list[0];
		for(i=1; i<nLoc; i++){
			// insertInQueue(list[i])
			ii=i+1;
			x = list[i];
			l=ii>>1;
			while(l && Q[l]>x){
				Q[ii] = Q[l];
				ii=l;
				l>>=1;
			}
			Q[ii]=x;
		}

		for(i=nn=0; i<nLoc; i++){
			x = Q[1];							// extract from Q
			setTopMinQ(Q, nLoc-i);

			if (isPrimaryMidle(x, M, &pI, &dxI) || (pI && BL_il[pI-1])){
				Pos[nn] = getPosPhraT(pI);
				L[nn] = Pos[nn] - dxI;
				Phr[nn] = pI-1;
				Dx[nn] = dxI;
				nn++;
				//if (j==1 && nn==134)
					//cout << x << " B--> " << L[nn-1] << endl;
			}

		}
		LL[j].n = nn;
		currN += nn;
		for (i=0; i<nn; i++)
			LL[j].B[i]=0;

		for (i=1; i<nn; i++){
			/*if (LL[j].L[i] <= LL[j].L[i-1]){
				cout << "Err B L[i] " << LL[j].L[i] << " <= L[i-1]" << LL[j].L[i-1] << ", i=" << i << ", j=" << j << endl;
				exit(0);
			}*/
		}
		delete [] Q;
	}

	// looking for primary occurrences in the last segment...
	j = k-1;
	ini = m-M;
	{
		LL[j].L = L = new ulong[nLocFin];
		LL[j].Phr = Phr = new ulong[nLocFin];
		LL[j].Pos = Pos = new ulong[nLocFin];
		LL[j].Dx = Dx = new uint[nLocFin];
		LL[j].B = new bool[nLocFin];
		Q = new ulong[nLocFin+1];
		ii=1;

		/*if (TRACE){
			cout << "C. Total occurrences found with the FMI : " << nLocFin << " =" << endl;
			for(i=0; i<nLocFin; i++)
				cout << list[i] << " ";
			cout << endl;
		}*/
		Q[1]=listFin[0];
		for(i=1; i<nLocFin; i++){
			// insertInQueue(list[i])...
			ii=i+1;
			x = listFin[i];
			l=ii>>1;
			while(l && Q[l]>x){
				Q[ii] = Q[l];
				ii=l;
				l>>=1;
			}
			Q[ii]=x;
		}

		for(i=nn=0; i<nLocFin; i++){
			x = Q[1];							// extract from Q
			setTopMinQ(Q, nLocFin-i);

			if (isPrimaryRight(x, m-ini, &pI, &dxI) || (pI && BL_il[pI-1])){
				Pos[nn] = getPosPhraT(pI);
				L[nn] = Pos[nn] - dxI;
				Phr[nn] = pI-1;
				Dx[nn] = dxI;
				nn++;
				//cout << x << " C--> " << L[nn-1] << endl;
			}

		}
		LL[j].n = nn;
		currN += nn;
		for (i=0; i<nLocFin; i++)
			LL[j].B[i]=0;
		delete [] Q;

		/*for (i=1; i<nn; i++){
			if (LL[j].L[i] <= LL[j].L[i-1]){
				cout << "Err C L[i] " << LL[j].L[i] << " <= L[i-1]" << LL[j].L[i-1] << ", i=" << i << ", j=" << j << endl;
				exit(0);
			}
		}*/
	}

	/*if (TRACE){
		cout << "FRASES..." << endl;
		for (i=0; i<k; i++){
			for (ii=0; ii<LL[i].n; ii++){
				cout << "LL[" << i << "].Phr[" <<ii<< "] = " << LL[i].Phr[ii] << endl;
				if (LL[i].Phr[ii]>=z){
					cout << "Err LL[i].Phr[ii] = " << LL[i].Phr[ii] << " >= z = " << z << endl;
					cout << "i = " << i << ", ii = " << ii << endl;
					exit(0);
				}
			}
		}
		cout << endl;
	}*/

	//===================================================================
	// Process the list
	// ==================================================================
	ulong *A = new ulong[currN];		// vector of final results...
	ulong phr, src;
	uint len, delta, inc;
	bool isKerR, isNot;

	for (i=0; i<k; i++){								// processing each list L_i
		nn=0;
		//cout << " L " << i << endl;
		if (i<k-1){
			// I.- Processing kernel L_i with L_{i+1}
			//cout << " I.- Processing kernel L_i with L_{i+1} " << endl;
			isKerR = true;
			for (ii=0; ii<LL[i].n; ii++){				// processing all the points of Li and searching concatenations with L_{i+1}
				if (LL[i].B[ii]==1){					// point already included
					LL[i].B[ii]=0;
					continue;
				}

				x = LL[i].L[ii];
				if (x > n-2*M)
					continue;
				j=i+1;
				u = x+M;
				if (j==k-1)
					u -= eps;

				if (isOccInArr(u, LL[j].n, LL[j].L, &l)){	// search u in L_j
					LL[i].B[ii]=1;						// we mark this occurrence in L_i in order to search further concatenations
					LL[j].B[l]=1;						// Discard the position u in L_j to not include it again L_i[ii] with L_j[l]
					nn++;
					continue;							// evaluate next x
				}

				phr = LL[i].Phr[ii];
				len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
				pos = LL[i].Pos[ii]-len;

				while (pos+len <= u){
					phr++;
					pos += len;
					len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
				}

				delta = u - pos;
				while (BL_il[phr]==false && delta+M <= len){
					phr -= rankBL(phr);
					src = getNum64(Src, phr*lgNG, lgNG);	// k index in G of the source of the phrase j
					u = getNum64(X, src*lgN, lgN) + delta;

					if (isOccInArr(u, LL[j].n, LL[j].L)){	// search u in L_j
						LL[i].B[ii]=1;					// we mark this occurrence in L_i in order to search further concatenations
						nn++;
						break;							// evaluate next x
					}

					phr = searchPhraTxt(u, &pos, &len);	 // phrase for the k index in G.. AQUI NECESITO EL INICIO DE LA FRASE POS
					delta = u - pos;
				}
			}
			if (nn==0) continue; 						// the list L_i does not has any candidate
		}else
			isKerR = false;


		if (i){
			// II.- Processing kernel L_i with L_{i-1}
			//cout << " II.- Processing kernel L_i with L_{i-1} " << endl;
			j=i-1;
			if (isKerR){
				for (ii=0; ii<LL[i].n; ii++){				// processing all the points of Li searching concatenations with L_{i-1}
					if (LL[i].B[ii]==0) continue;
					x = LL[i].L[ii];
					if (x < M){
						LL[i].B[ii]=0;
						nn--;
						continue;
					}
					u = x-M;
					if (i==k-1)
						u += eps;

					if (isOccInArr(u, LL[j].n, LL[j].L))	// search u in L_j
						continue;							// evaluate next x. here LL[i].B[ii] is already marked

					phr = LL[i].Phr[ii];
					len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
					pos = LL[i].Pos[ii]-len;
					while (pos > u){
						phr--;
						len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
						pos -= len;
					}
					delta = pos + len - u;

					isNot = true;
					while (delta >= M && BL_il[phr]==false){
						phr -= rankBL(phr);
						src = getNum64(Src, phr*lgNG, lgNG);	// k index in G of the source of the phrase j
						u = getNum64(X, src*lgN, lgN) + len - delta;

						if (isOccInArr(u, LL[j].n, LL[j].L)){	// search u in L_j
							isNot = false;
							break;							// evaluate next x. here LL[i].B[ii] is already marked
						}

						phr = searchPhraTxt(u, &pos, &len);	// phrase that contains u.  AQUI NECESITO EL FIN DE LA FRASE (POS+LEN)
						delta = pos + len - u;
					}
					if (isNot){
						LL[i].B[ii]=0;						// occurrence not found for this x. We discarding the point x as primary occurrence of the pattern
						nn--;
					}
				}
			}else {
				for (ii=0; ii<LL[i].n; ii++){				// processing all the points of Li searching concatenations with L_{i-1}
					if (LL[i].B[ii]){						// point already marked in order to not to re-evaluate
						LL[i].B[ii]=0;
						continue;
					}
					x = LL[i].L[ii];
					if (x < M){
						LL[i].B[ii]=0;
						continue;
					}
					u = x-M;
					if (i==k-1)
						u += eps;

					if (isOccInArr(u, LL[j].n, LL[j].L)){	// search u in L_j
						LL[i].B[ii]=1;						// new occurrence !! we mark this occurrence in L_i is order to check further concatenations
						nn++;
						continue;							// evaluate next x.
					}

					phr = LL[i].Phr[ii];
					len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
					pos = LL[i].Pos[ii]-len;
					while (pos > u){
						phr--;
						len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
						pos -= len;
					}
					delta = pos + len - u;

					while (delta >= M && BL_il[phr]==false){
						phr -= rankBL(phr);
						src = getNum64(Src, phr*lgNG, lgNG);	// k index in G of the source of the phrase j
						u = getNum64(X, src*lgN, lgN) + len - delta;

						if (isOccInArr(u, LL[j].n, LL[j].L)){	// search u in L_j
							LL[i].B[ii]=1;					// we mark this occurrence in L_i
							nn++;
							break;							// evaluate next x. here LL[i].B[ii] is already marked
						}

						phr = searchPhraTxt(u, &pos, &len);	// phrase that contains u.
						delta = pos + len - u;
					}
				}
			}
		}
		if (nn==0) continue; // the list L_i does not has candidates

		if (i>=2){
			// III.- Backward processing of L_i with lists from i-2 to 0
			//cout << " III.- Backward processing of L_i with lists from i-2 to 0 " << endl;
			j = i-1;
			while (nn && j){
				j--;
				inc = i-j;

				for (ii=0; nn && ii<LL[i].n; ii++){			// processing all the points of Li and searching concatenations with Lj
					if (LL[i].B[ii]==0) continue;
					x = LL[i].L[ii];
					u = inc*M;
					if (x < u){
						LL[i].B[ii]=0;						// occurrence not found for this x. We discarding the point x as primary occurrence of the pattern
						nn--;
						continue;
					}
					u = x-u;
					if (i==k-1)
						u += eps;

					if (isOccInArr(u, LL[j].n, LL[j].L))	// search u in L_j
						continue;							// point ok... evaluate next x

					phr = LL[i].Phr[ii];
					len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
					pos = LL[i].Pos[ii]-len;
					while (pos > u){
						phr--;
						len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
						pos -= len;
					}
					delta = pos + len - u;

					isNot = true;
					while (delta >= M && BL_il[phr]==false){
						phr -= rankBL(phr);
						src = getNum64(Src, phr*lgNG, lgNG);	// k index in G of the source of the phrase j
						u = getNum64(X, src*lgN, lgN) + len - delta;

						if (isOccInArr(u, LL[j].n, LL[j].L)){	// search u in L_j
							isNot = false;
							break;							// Point ok... evaluate next x. here LL[i].B[ii] is already marked
						}
						phr = searchPhraTxt(u, &pos, &len);	// phrase that contains u.  AQUI NECESITO EL FIN DE LA FRASE (POS+LEN)
						delta = pos + len - u;
					}
					if (isNot){
						LL[i].B[ii]=0;						// occurrence not found for this x. We discarding the point x as primary occurrence of the pattern
						nn--;
					}
				}
			}
		}

		// IV.- Forward processing of L_i with lists from i+2 to k-1
		//cout << " IV.- Forward processing of L_i with lists from i+2 to k-1 " << endl;
		for (j=i+2; nn && j<k; j++){
			inc = j-i;

			for (ii=0; nn && ii<LL[i].n; ii++){			// processing all the points of Li and searching concatenations with Lj
				if (LL[i].B[ii]==0) continue;
				x = LL[i].L[ii];
				if (x > n-inc*M){
					LL[i].B[ii]=0;
					continue;
				}
				u = x+inc*M;
				if (j==k-1)
					u -= eps;

				if (isOccInArr(u, LL[j].n, LL[j].L, &l)){	// search u in L_j
					LL[j].B[l]=1;						// Discard position of u in list j, to not check it again
					continue;							// evaluate next x
				}

				phr = LL[i].Phr[ii];
				len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
				pos = LL[i].Pos[ii]-len;
				while (pos+len <= u){
					phr++;
					pos += len;
					len = getNum64(PhraT, (phr+1)*lgPT, lgPT);
				}
				delta = u-pos;

				isNot = true;
				while (BL_il[phr]==false && delta+M <= len){
					phr -= rankBL(phr);
					src = getNum64(Src, phr*lgNG, lgNG);	// k index in G of the source of the phrase j
					u = getNum64(X, src*lgN, lgN)+delta;

					if (isOccInArr(u, LL[j].n, LL[j].L)){	// search u in L_j
						isNot = false;
						break;							// evaluate next x
					}
					phr = searchPhraTxt(u, &pos, &len);	// phrase for the k index in G.. AQUI NECESITO EL INICIO DE LA FRASE POS
					delta = u-pos;
				}
				if (isNot){
					LL[i].B[ii]=0;						// occurrence not found for this x. We discarding the point x as primary occurrence of the pattern
					nn--;
					continue;
				}
			}
		}
	}

	// making the final answer of primary occurrences by the union of found points in each list...
	for (i=nn=0; i<k; i++){
		if (i==k-1)
			inc = i*M-eps;
		else
			inc = i*M;
		for (ii=0; ii<LL[i].n; ii++){
			if (LL[i].B[ii]){
				x = LL[i].L[ii] - inc;
				A[nn] = x;
				nn++;
				//cout << " Pri x = " << x << endl;
			}
		}

		delete [] LL[i].L;
		delete [] LL[i].Phr;
		delete [] LL[i].Pos;
		delete [] LL[i].Dx;
		delete [] LL[i].B;
	}
	delete [] LL;

	// search secondary occurrence from primary ones...
	*occ = A;
	*nOcc = nn;
	for(i=0; i<nn; i++){
		if (findPredecesor(A[i], &j)){
			locateSecOcc(0, j, A[i], m, nOcc, occ, &currN);
			A = *occ;
		}
	}
}

// =============================================================================================================
// Extracting text...
// =============================================================================================================
/*
 * The algorithm to extract the segment T[x..y] is:
 *
 * 1.- we maps T[x..y] to T'[x'..y'] using the array PhrT[] (which is encoded with pfrodelta)
 * 2.- We extract using I_M(T') the symbols T'[x'..y']
 * 3.- for each mismatch symbol '#' in T'[x'..y'] we go recursively to the source of the phrase that contain the mismatch symbol using the array Sr[]
 * 4.- And recursively retrieve the symbol from this source until that we do not find a symbol "#"
 *
 */

// Extracts the len characters that represent the original segment T[sp.. sp+len-1] and allocates these into A[0...len-1]
void HybridSelfIndex::extract(ulong sp, ulong lenS, uchar **A){
	*A = new uchar[lenS+1];
	extract(sp, lenS, *A);
}

void HybridSelfIndex::extract(ulong x, ulong lenS, uchar *A){
	ulong y, pI, pE, posPI, posPE, i, ini, end, cStr;
	uint lIni, sIni, l, dx, dy, lenPI, lenPE, twoM=M<<1;
	string strA;

	// make mapping from Txt[sp..ep] to fTxt[ini..end]
	y = x+lenS-1;
	pI = searchPhraTxt(x, &posPI, &lenPI);
	dx = x-posPI;
	dy = lenPI-dx;
	ini = end = getPosPhraFT(pI);

	// initial segment...
	i=l=0;
	if (BL_il[pI]==false){
		if(dx < M){
			ini += dx;
			lIni = M-dx;
			if (lIni >= lenS){
				strA = sdsl::extract(FMI, ini, ini+lenS-1);				// ** extract fTxt[ini..ini+lenS-1] from the FMI
				strcpy((char *)A, strA.substr(0, lenS).c_str());
				return;
			}else{
				sIni = lenPI-twoM;
				if (lIni+sIni >= lenS){  // ES NECESARIO, ES MUY POCO PROBABLE PARECE ????????????????
					strA = sdsl::extract(FMI, ini, ini+lIni-1);			// ** extract fTxt[ini..ini+lIni-1] from the FMI
					strcpy((char *)A, strA.c_str());
					x = pI-rankBL.rank(pI);
					x = getNum64(Src, x*lgNG, lgNG);					// index in G of the source of the phrase pI
					x = getNum64(X, x*lgN, lgN)+dx+lIni;
					lenS -= lIni;
					extract(x, lenS, A+lIni);							//retrieveTxtSrc(pI, dx, sIni, A+lIni);
					return;
				}
			}
			l = M;
		}else{
			lIni = 0;
			if (dy <= M){
				ini += twoM-dy+1;
				sIni = 0;
				l = dy;
			}else{
				sIni = dy-M;
				x = pI-rankBL.rank(pI);
				x = getNum64(Src, x*lgNG, lgNG);						// index in G of the source of the phrase pI
				x = getNum64(X, x*lgN, lgN)+dx;

				if (sIni >= lenS){
					extract(x, lenS, A);
					return;
				}
				l = M;
				ini += l+1;
				extract(x, sIni, A);									//retrieveTxtSrc(pI, dx, sIni, A);
				lenS -= sIni;
				i += sIni;
			}
		}
	}else{																// PI is completely in the filtered text
		ini += dx;
		if (lenPI-dx >= lenS){
			strA = sdsl::extract(FMI, ini, ini+lenS-1);					// ** extract fTxt[ini..ini+lIni-1] from the FMI
			strcpy((char *)A, strA.c_str());
			return;
		}
		lIni = sIni = 0;
		l = lenPI-dx;
	}

	// last segment...
	if (posPI+lenPI<=y){
		pE = searchPhraTxt(y, &posPE, &lenPE);
		end = getPosPhraFT(pE);
	}else{
		pE = pI;	// ESTO NO VA A PASAR MAS!!
		posPE = posPI;
		lenPE = lenPI;
	}
	dx = y-posPE;
	dy = lenPE-dx;
	if (BL_il[pE]==false && dx>M){
		if (dy>M)
			end += M-1;
		else
			end += twoM+1-dy;
	}else 																// PE is completely in the filtered text
		end += dx;

	strA = sdsl::extract(FMI, ini, end);								// ** extract fTxt[ini..end] from the FMI

	if (lIni){
		strcpy((char *)A+i, strA.substr(0, lIni).c_str());				// copy strA[0... lIni-1] into A[0..lIni-1]
		cStr=lIni+1;//+k													// +1 to avoid the mismatch symbol
		x = pI-rankBL.rank(pI+1);
		x = getNum64(Src, x*lgNG, lgNG);								// index in G of the source of the phrase pI
		x = getNum64(X, x*lgN, lgN)+M;
		extract(x, sIni, A+lIni);										//retrieveTxtSrc(pI, dx, sIni, A+lIni);
		i =lIni+sIni;
		lenS -= i;
	}else
		cStr=0;

	// intermediate and complete phrase segments...
	pI++;
	while(pI<=pE){
		if (BL_il[pI]){
			pI++;
			lenPI = getNum64(PhraT, pI*lgPT, lgPT);
			//lenPI = PfD_PT->get_ithGap(pI);
			l += lenPI;
		}
		l += M;
		if (l >= lenS){
			strcpy((char *)A+i, strA.substr(cStr, cStr+lenS).c_str());	// copy in A[i..i+lenS-1] <- strA[cStr... cStr+lenS-1]
			return;
		}else{
			strcpy((char *)A+i, strA.substr(cStr, cStr+l).c_str());		// copy strA[cStr... l-1] into A[i..i+l-1]
			//cout << "   A = [" << A << "]" << endl;
			cStr+=l+1;  // that is +k=1
			lenS -= l;
		}
		x = pI-rankBL.rank(pI);
		x = getNum64(Src, x*lgNG, lgNG);								// index in G of the source of the phrase pI
		x = getNum64(X, x*lgN, lgN)+M;
		lenPI = getNum64(PhraT, (pI+1)*lgPT, lgPT);
		lIni = lenPI-twoM;
		if (lIni >= lenS){
			extract(x, lenS, A+l+i);									//retrieveTxtSrc(pI, M, lenPI-twoM, A+l);
			return;
		}
		extract(x, lIni, A+l+i);										//retrieveTxtSrc(pI, M, lenPI-twoM, A+l);
		lenS -= lIni;
		i += lIni+l;
		//cout << "   A = [" << A << "]" << endl;
		l = M;
		pI++;
	}

	if (l >= lenS){
		strcpy((char *)A+i, strA.substr(cStr, cStr+lenS).c_str());		// copy in A[i..i+lenS-1] <- strA[cStr... cStr+lenS-1]
		return;
	}

	strcpy((char *)A+i, strA.substr(cStr, cStr+l).c_str());				// copy in A[i..i+l-1] <- strA[cStr... cStr+l-1]
	i+=l;
	x = pE-rankBL.rank(pE+1);
	x = getNum64(Src, x*lgNG, lgNG);									// index in G of the source of the phrase pE
	x = getNum64(X, x*lgN, lgN)+M;
	if (dy > M){														// retrieve from source of PE
		extract(x, dx-M+1, A+i);										// retrieveTxtSrc(pE, M, dx-M, A+i);
		return;
	}

	cStr +=l;
	lenS -= l;
	l = lenPE-twoM;
	extract(x, l, A+i);													//retrieveTxtSrc(pE, M, lenPE-twoM, A+i);
	i += l;
	strcpy((char *)A+i, strA.substr(cStr, cStr+lenS-l).c_str());		// copy strA[cStr... lenS-1] into A[i..i+lenS-1]
}

// return true if the occurrence T'[pos...pos+m-1] is primary in T[x...x+m-1].
// public function for testing
bool HybridSelfIndex::isPrimaryT(ulong x, uint len, ulong *phra){
	ulong x1;
	uint wide;

	*phra = searchPhraTxt(x+len-1, &x1, &wide);
	if (x1 > x)
		return true;

	if (BL_il[*phra])
		return true;

	return false;
}

// ====================================================================================================
// save/load methods
// ====================================================================================================
// save the Data Structure in the directory 'pathName'. pathName has to end with /
void HybridSelfIndex::saveStructure(){
	char *fileName = new char[400];
	ulong lenArray;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "selfHI.hsi");
	ofstream os (fileName, ios::binary);

	os.write((const char*)&n, sizeof(ulong));
	os.write((const char*)&nFT, sizeof(ulong));

	os.write((const char*)&nG, sizeof(uint));
	os.write((const char*)&z, sizeof(uint));
	os.write((const char*)&M, sizeof(uint));
	os.write((const char*)&lgPT, sizeof(uint));
	os.write((const char*)&lgPFT, sizeof(uint));
	os.write((const char*)&nSamP, sizeof(uint));
	os.write((const char*)&b_X, sizeof(uint));
	os.write((const char*)&nSampX, sizeof(uint));

	// Saving arrays o fphrases...
	lenArray = z*lgPT / (8*sizeof(ulong));
	if ((z*lgPT) % (8*sizeof(ulong)))
		lenArray++;
	os.write((const char*)PhraT, lenArray*sizeof(ulong));

	lenArray = nSamP*lgN / (8*sizeof(ulong));
	if ((nSamP*lgN) % (8*sizeof(ulong)))
		lenArray++;
	os.write((const char*)SGCPT, lenArray*sizeof(ulong));

	lenArray = z*lgPFT / (8*sizeof(ulong));
	if ((z*lgPFT) % (8*sizeof(ulong)))
		lenArray++;
	os.write((const char*)PhraFT, lenArray*sizeof(ulong));

	lenArray = nSamP*lgN / (8*sizeof(ulong));
	if ((nSamP*lgN) % (8*sizeof(ulong)))
		lenArray++;
	os.write((const char*)SGCPFT, lenArray*sizeof(ulong));

	if(nG){
		// Saving arrays of the grid...
		lenArray = nG*lgN / (8*sizeof(ulong));
		if ((nG*lgN) % (8*sizeof(ulong)))
			lenArray++;
		os.write((const char*)X, lenArray*sizeof(ulong));

		os.write((const char*)LOOKUP_X, nSampX*sizeof(uint));
		os.write((const char*)COUNT_X, nSampX*sizeof(uint));

		lenArray = nG*lgZ / (8*sizeof(ulong));
		if ((nG*lgZ) % (8*sizeof(ulong)))
			lenArray++;
		os.write((const char*)Pr, lenArray*sizeof(ulong));

		lenArray = nG*lgNG / (8*sizeof(ulong));
		if ((nG*lgNG) % (8*sizeof(ulong)))
			lenArray++;
		os.write((const char*)Src, lenArray*sizeof(ulong));
	}

	os.close();

	// ======================================================================================
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "BL_il.hsi");
	store_to_file(BL_il, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "rankBL.hsi");
	store_to_file(rankBL, fileName);

	// ======================================================================================
	if(nG){
		// Saving rmqY...
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "rmqY.hsi");
		rmqY->saveDS(fileName);
	}

	if (SHOW_SIZE || TRACE){
		cout << "Save FastHybridIndex data structures in " << fileName << endl;
		cout << "Saving structures..." << endl;

		cout << " .- BL_il[] " << size_in_bytes(BL_il) << " Bytes" << endl;

		lenArray = z*lgPT / (8*sizeof(ulong));
		if ((z*lgPT) % (8*sizeof(ulong)))
			lenArray++;
		cout << " .- PhraT[] " << lenArray*sizeof(ulong) << " Bytes" << endl;

		lenArray = nSamP*lgN / (8*sizeof(ulong));
		if ((nSamP*lgN) % (8*sizeof(ulong)))
			lenArray++;
		cout << " .- SGCPT[] " << lenArray*sizeof(ulong) << " Bytes" << endl;

		lenArray = z*lgPFT / (8*sizeof(ulong));
		if ((z*lgPFT) % (8*sizeof(ulong)))
			lenArray++;
		cout << " .- PhraFT[] " << lenArray*sizeof(ulong) << " Bytes" << endl;

		lenArray = nSamP*lgN / (8*sizeof(ulong));
		if ((nSamP*lgN) % (8*sizeof(ulong)))
			lenArray++;
		cout << " .- SGCPFT[] " << lenArray*sizeof(ulong) << " Bytes" << endl;

		if(nG){
			lenArray = nG*lgN / (8*sizeof(ulong));
			if ((nG*lgN) % (8*sizeof(ulong)))
				lenArray++;
			cout << " .- X[] " << lenArray*sizeof(ulong) << " Bytes" << endl;

			cout << " .- LOOKUP_X[]+COUNT_X[] " << 2*nSampX*sizeof(uint) << " Bytes" << endl;

			lenArray = nG*lgZ / (8*sizeof(ulong));
			if ((nG*lgZ) % (8*sizeof(ulong)))
				lenArray++;
			cout << " .- Pr[] " << lenArray*sizeof(ulong) << " Bytes" << endl;

			lenArray = nG*lgNG / (8*sizeof(ulong));
			if ((nG*lgNG) % (8*sizeof(ulong)))
				lenArray++;
			cout << " .- Src[] " << lenArray*sizeof(ulong) << " Bytes" << endl;

			cout << " .- rmqY " << rmqY->getSize() << " Bytes" << endl;
		}

		cout << "   Structures saved !" << endl;
		cout << "______________________________________________________________" << endl;
	}
}

void HybridSelfIndex::loadStructure(){
	char *fileName = new char[400];
	ulong lenArray;
	sizeDS = 0;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "selfHI.hsi");
	ifstream is(fileName, ios::binary);

	// ======================================================================================
	// Loading variables...
	is.read((char*)&n, sizeof(ulong));
	is.read((char*)&nFT, sizeof(ulong));

	is.read((char*)&nG, sizeof(uint));
	is.read((char*)&z, sizeof(uint));
	is.read((char*)&M, sizeof(uint));
	is.read((char*)&lgPT, sizeof(uint));
	is.read((char*)&lgPFT, sizeof(uint));
	is.read((char*)&nSamP, sizeof(uint));
	is.read((char*)&b_X, sizeof(uint));
	is.read((char*)&nSampX, sizeof(uint));

	this->lgN = 1 + (uint)(log(n+1)/log(2));
	this->lgNG = 1 + (uint)(log(nG+1)/log(2));
	this->lgM = 1 + (uint)(log(M)/log(2));
	this->lgZ = 1 + (uint)(log(z)/log(2));

	lenArray = z*lgPT / (8*sizeof(ulong));
	if ((z*lgPT) % (8*sizeof(ulong)))
		lenArray++;
	PhraT = new ulong[lenArray];
	is.read((char*)PhraT, lenArray*sizeof(ulong));
	sizeDS += lenArray*sizeof(ulong);

	lenArray = nSamP*lgN / (8*sizeof(ulong));
	if ((nSamP*lgN) % (8*sizeof(ulong)))
		lenArray++;
	SGCPT = new ulong[lenArray];
	is.read((char*)SGCPT, lenArray*sizeof(ulong));
	sizeDS += lenArray*sizeof(ulong);

	lenArray = z*lgPFT / (8*sizeof(ulong));
	if ((z*lgPFT) % (8*sizeof(ulong)))
		lenArray++;
	PhraFT = new ulong[lenArray];
	is.read((char*)PhraFT, lenArray*sizeof(ulong));
	sizeDS += lenArray*sizeof(ulong);

	lenArray = nSamP*lgN / (8*sizeof(ulong));
	if ((nSamP*lgN) % (8*sizeof(ulong)))
		lenArray++;
	SGCPFT = new ulong[lenArray];
	is.read((char*)SGCPFT, lenArray*sizeof(ulong));
	sizeDS += lenArray*sizeof(ulong);

	if(nG){
		// Loading arrays of the grid...
		lenArray = nG*lgN / (8*sizeof(ulong));
		if ((nG*lgN) % (8*sizeof(ulong)))
			lenArray++;
		X = new ulong[lenArray];
		is.read((char*)X, lenArray*sizeof(ulong));
		sizeDS += lenArray*sizeof(ulong);

		LOOKUP_X = new uint[nSampX];
		is.read((char*)LOOKUP_X, nSampX*sizeof(uint));
		sizeDS += nSampX*sizeof(uint);

		COUNT_X = new uint[nSampX];
		is.read((char*)COUNT_X, nSampX*sizeof(uint));
		sizeDS += nSampX*sizeof(uint);

		lenArray = nG*lgZ / (8*sizeof(ulong));
		if ((nG*lgZ) % (8*sizeof(ulong)))
			lenArray++;
		Pr = new ulong[lenArray];
		is.read((char*)Pr, lenArray*sizeof(ulong));
		sizeDS += lenArray*sizeof(ulong);

		lenArray = nG*lgNG / (8*sizeof(ulong));
		if ((nG*lgNG) % (8*sizeof(ulong)))
			lenArray++;
		Src = new ulong[lenArray];
		is.read((char*)Src, lenArray*sizeof(ulong));
		sizeDS += lenArray*sizeof(ulong);
	}
	is.close();

	// ======================================================================================
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "BL_il.hsi");
	load_from_file(BL_il, fileName);
	sizeDS += size_in_bytes(BL_il);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "rankBL.hsi");
	load_from_file(rankBL, fileName);
	util::init_support(rankBL, &BL_il);

	// ======================================================================================
	if(nG){
		// Loading rmqY...
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "rmqY.hsi");
		rmqY = new RMQRMM64(fileName);
		sizeDS += rmqY->getSize();
	}

	// ======================================================================================
	if(nFT){
		// Loading FMI
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "fmi.hsi");
		load_from_file(FMI, fileName);
		sizeDS += size_in_bytes(FMI);
	}

	if (SHOW_SIZE || TRACE){
		cout << " Load data structures from " << fileName << endl;
		cout << "BreakDown --size of the structures:" << endl;

		cout << " ** load BL_il[] of size " << size_in_bytes(BL_il) << " Bytes" << endl;

		lenArray = z*lgPT / (8*sizeof(ulong));
		if ((z*lgPT) % (8*sizeof(ulong)))
			lenArray++;
		cout << " ** load PhraT[] of size " << lenArray*sizeof(ulong) << " Bytes" << endl;

		lenArray = nSamP*lgN / (8*sizeof(ulong));
		if ((nSamP*lgN) % (8*sizeof(ulong)))
			lenArray++;
		cout << " ** load SGCPT[] of size " << lenArray*sizeof(ulong) << " Bytes" << endl;

		lenArray = z*lgPFT / (8*sizeof(ulong));
		if ((z*lgPFT) % (8*sizeof(ulong)))
			lenArray++;
		cout << " ** load PhraFT[] of size " << lenArray*sizeof(ulong) << " Bytes" << endl;

		lenArray = nSamP*lgN / (8*sizeof(ulong));
		if ((nSamP*lgN) % (8*sizeof(ulong)))
			lenArray++;
		cout << " ** load SGCPFT[] of size " << lenArray*sizeof(uint) << " Bytes" << endl;

		if(nG){
			lenArray = nG*lgN / (8*sizeof(ulong));
			if ((nG*lgN) % (8*sizeof(ulong)))
				lenArray++;
			cout << " ** load X[] of size " << lenArray*sizeof(ulong) << " Bytes" << endl;

			cout << " ** load LOOKUP_X[]+COUNT_X[] of size " << 2*nSampX*sizeof(uint) << " Bytes" << endl;

			lenArray = nG*lgZ / (8*sizeof(ulong));
			if ((nG*lgZ) % (8*sizeof(ulong)))
				lenArray++;
			cout << " ** load Pr[] of size " << lenArray*sizeof(ulong) << " Bytes" << endl;

			lenArray = nG*lgNG / (8*sizeof(ulong));
			if ((nG*lgNG) % (8*sizeof(ulong)))
				lenArray++;
			cout << " ** load Src[] of size " << lenArray*sizeof(ulong) << " Bytes" << endl;

			cout << " ** load rmqY of size " << rmqY->getSize() << " Bytes" << endl;
		}

		if (nFT)
			cout << " **  load FMI (S_SA=" << S_SA << ", S_ISA=" << S_ISA<< ") size " << size_in_bytes(FMI) << " Bytes = " << (float)size_in_bytes(FMI)/(1024.0*1024.0) << "MiB = " << (float)size_in_bytes(FMI)*8.0/(float)n << " bpc" << endl;
	}

	cout << "______________________________________________________________" << endl;
	cout << "   Index Structures loaded. sizeDS = " << sizeDS/(1024.0*1024.0) << " MiB = " << (float)sizeDS*8.0/(float)n << " bpc" << endl;
	cout << "______________________________________________________________" << endl;
}

HybridSelfIndex::~HybridSelfIndex() {
	delete [] PhraT;
	delete [] PhraFT;
	delete [] SGCPT;
	delete [] SGCPFT;

	decltype(BL_il) emptyBL;
	BL_il.swap(emptyBL);

	if(nG){
		delete [] X;
		delete [] LOOKUP_X;
		delete [] COUNT_X;

		delete [] Pr;
		delete [] Src;

		rmqY->~RMQRMM64();
	}

	if (nFT){
		decltype(FMI) emptyFMI;
		FMI.swap(emptyFMI);
	}
	if (CREATE_FMI_TEST){
		decltype(FMI_TEST) emptyFMITEST;
		FMI_TEST.swap(emptyFMITEST);
	}
}

