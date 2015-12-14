//============================================================================
// Name        : buildTopkLZ.cpp
// Author      : Hector Ferrada
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "TopkLZ.h"
#include <ConfigFile.h>

bool TRACE = false;			// true: print all details for console
bool TEST = false;			// true: apply exhaustive test
uint N_REP = 10;
uint M = 10;

// Structure with all globals parameters program
typedef struct {
	string configFile;		// properties file

	uchar *seq;				// original sequence (1 byte for symbol)
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	char cutDoc;			// symbol to separate documents
	uint g;					// minimum length of document segments to store a topk list in a node. Maximum 8 bits
	TopkLZ *index;

	char inputFile[200];	// list of files
	bool lowercase;			// 1: transform to lowercase
	bool filesInList;			// 1: list of files, 0: Unique file
	char boundSymbol;		// original symbol delimiter of documents when we read all documents in 1 file.
	char dirStore[200];		// directory to save/load the data structure (files *.tk)

	// The following data structure are only for test !!
	ulong* patterns;
} ParamProgram;

void testSearchPatterns_CSA_LZ(ParamProgram *par);
bool searchPatternInFMI(ParamProgram *par, uchar* pat, uint m, ulong *occt1, ulong *occt2, ulong *occt3, ulong *fqCSA);
void createTopK_Test(ParamProgram *par, uint *docList, ulong *frqList, ulong *fqCSA, uint kst);
void createPatterns(ParamProgram *par, uint m, ulong repeat);

int main(int argc, char *argv[]) {
	ParamProgram *par = new ParamProgram();
	char fileName[300];

	if(argc != 2){
		cout << "ERRORR !! " << endl;
		cout << "buildAppTopkLZhq's usage requires the Properties File as parameter !! " << endl;
		cout << "Example for the file 'config.txt' in the same directory: ./buildAppTopkLZhq config.txt" << endl;
		exit(1);
	}
	par->configFile = string(argv[1]);
	cout << "Congif file: " << par->configFile << endl;
	ConfigFile cf(par->configFile);

	TRACE = cf.Value("GLOBALS","TRACE");
	TEST = cf.Value("GLOBALS","TEST");
	N_REP = cf.Value("GLOBALS","N_REP");
	M = cf.Value("GLOBALS","MAX_M");

	par->g = cf.Value("TOPK","g");
	strcpy(par->inputFile, ((string)(cf.Value("TOPK","inputFile"))).c_str());
	par->filesInList = cf.Value("TOPK","filesInList");
	strcpy(par->dirStore, ((string)(cf.Value("TOPK","dirStore"))).c_str());
	par->lowercase = cf.Value("TOPK","lowercase");
	par->boundSymbol = cf.Value("TOPK","boundSymbol");
	par->cutDoc = cf.Value("TOPK","cutDoc");

	cout << "buildTopkLZ parameters..." << endl;
	cout << "  g = " << (uint)par->g << endl;
	cout << "  Input file: " << par->inputFile << endl;
	cout << "  Files in list: " << par->filesInList << endl;		// 0:archivo unico para todocl por ejemplo
	cout << "  dirStore: " << par->dirStore << endl;
	cout << "  lowercase: " << par->lowercase << endl;
	cout << "  boundSymbol code: " << (int)(par->boundSymbol) << endl;
	cout << "  cutDoc code: " << (int)(par->cutDoc) << endl;

	TopkLZ::TRACE = TRACE;
	TopkLZ::TEST = TEST;
	par->index = new TopkLZ(par->g, par->inputFile, par->filesInList, par->cutDoc, par->lowercase, par->dirStore, par->boundSymbol);
	par->index = new TopkLZ(par->dirStore, true);
	par->n = par->index->n;
	cout << "____________________________________________________" << endl;
	cout << "***  Index size " << par->index->sizeDS << " bytes = " << (float)par->index->sizeDS*8.0/(float)par->n << " bpc" << endl;
	cout << "====================================================" << endl;

	if (TEST){
		// load Sequence...
		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sequence.test");
		ifstream is(fileName, ios::binary);
		par->seq = new uchar[par->n];
		is.read((char*)par->seq, par->n*sizeof(uchar));
		is.close();

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "EndDocs.test");
		ifstream isEndDocs(fileName, ios::binary);
		par->index->EndDocs = new ulong[par->index->nDocs];
		isEndDocs.read((char*)par->index->EndDocs, par->index->nDocs*sizeof(ulong));
		isEndDocs.close();

		// load separator bit sequence and rank/select support for test!
		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sep_rrr.test");
		load_from_file(par->index->sep_rrr, fileName);

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sep_rank.test");
		load_from_file(par->index->sep_rank, fileName);
		util::init_support(par->index->sep_rank, &(par->index->sep_rrr));

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sep_sel.test");
		load_from_file(par->index->sep_sel, fileName);
		util::init_support(par->index->sep_sel, &(par->index->sep_rrr));

		// load fmi...
		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "fmi.temp");
		load_from_file(par->index->fmi, fileName);
		cout << " FMI loaded, it requires " << size_in_mega_bytes(par->index->fmi) << " MiB." << endl;

		cout << "Running test searchPattern.." << endl;
		testSearchPatterns_CSA_LZ(par);
		cout << "Test searchPattern OK !!" << endl;
	}

	par->index->~TopkLZ();

	cout << "$$$$$$$$$$$$$$$$$$$$$" << endl;
	return 0;
}

void testSearchPatterns_CSA_LZ(ParamProgram *par){
	//TRACE = true;
	ulong j, k, t, pv, x, realK;
	uint m, l;
	ulong occ1, occ2, occ3, diffK;
	bool isInRev;
	uchar *pat = new uchar[M+1];
	uint *docList = new uint[par->index->nDocs+1];
	uint *docList_test = new uint[par->index->nDocs+1];
	ulong *frqList_test = new ulong[par->index->nDocs+1];
	ulong *fqCSA = new ulong[par->index->nDocs+1];
	par->patterns = new ulong[N_REP];

	uint MAX_K = 10, PLUS_K = 1, INI_K = 1;
	if (par->index->nDocs >= 100){
		MAX_K = 100;
		PLUS_K = 20;
		INI_K = 20;
	}else{
		if (MAX_K >= par->index->nDocs)
			MAX_K = par->index->nDocs/2;
	}

	//INI_K = 5;
	for (uint kst=INI_K; kst<=MAX_K; kst+=PLUS_K){
		cout << " Test for k = " << kst << "..." << endl;
		for (m=2; m<=M; m++){
			createPatterns(par, m, N_REP);
			for (j=0; j<=par->index->nDocs; j++)
				par->index->dictD[j] = par->index->keyFq[j] = par->index->keyDoc[j] = 0;

			//cout << " ... test in patterns of length m = " << m << endl;
			for (t=0; t<N_REP; t++){
				for (l=0; l<m; l++)
					pat[l] = par->seq[par->patterns[t]+l];
				pat[m] = '\0'; // [ssian ]
				/*pat[0] = 'd';
				pat[1] = '\n';
				/*pat[2] = 'i';
				pat[3] = 'a';
				pat[4] = 'n';
				pat[5] = ' ';
				pat[6] = 'a';
				pat[7] = 'b';
				pat[8] = 'a';
				pat[9] = 'p';
				pat[10] = 'a';
				pat[11] = 'l';
				pat[12] = 'a';
				pat[13] = 'b';
				pat[14] = 'r';
				pat[15] = 'a';
				pat[16] = 'r';
				pat[17] = 'l';
				pat[18] = 'a';
				pat[19] = 'a';*/

				for(k=0; k<=par->index->nDocs; k++)
					fqCSA[k] = docList_test[k] = frqList_test[k] = 0;

				searchPatternInFMI(par, pat, m, &occ1, &occ2, &occ3, fqCSA);
				if (TRACE){
					cout << "pat=[" << pat << "], m=" << m << ", test=" << t << ", k*=" << kst << endl;
					cout << "occ1 = " << occ1 << ", occ2 = " << occ2 << ", occ3 = " << occ3 << endl;
					for(k=0; k<par->index->nDocs; k++){
						if(fqCSA[k] > 0)
							cout << k << "(" << fqCSA[k] << ") ";
					}
					cout << endl;
				}

				// TEST OCCURRENCES TYPE 1...
				x = 1;
				isInRev = par->index->searchPattern_Rev3(pat, m, &x, &pv, 0);
				if (isInRev == false && occ1){
					cout << "ERROR, [" << pat << "], m=" << m << ", k* =" << kst << ", test=" << t << ", NOT FOUND in RevTrie, but IS FOUND in the FMI with " << occ1 << " occurrences !!" << endl;
					exit(1);
				}
				if (isInRev && occ1==0){
					cout << "ERROR, [" << pat << "], m=" << m << ", k* =" << kst << ", test=" << t << ", is in the node with preorder " << pv << " but NOT FOUND int the FMI!!" << endl;
					exit(1);
				}

				diffK = 0;
				for(k=0; k<par->index->nDocs; k++){
					if(fqCSA[k] > 0)
						diffK++;
				}

				if(kst > diffK)
					realK = diffK;
				else
					realK = kst;
				createTopK_Test(par, docList_test, frqList_test, fqCSA, realK);
				if (TRACE){
					cout << "Real Docs and Frequencies..." << endl;
					for(k=0; k<realK && frqList_test[k]; k++)
						cout << k+1 << ", (" << docList_test[k] << "," << frqList_test[k] << ") " << endl;
				}

				diffK = par->index->topKDocument(pat, m, docList, realK);
				if (TRACE){
					cout << "diffK= " << diffK << endl;
					cout << "List returned..." << endl;
					for(k=0; k<diffK; k++)
						cout << k+1 << ", (" << docList[k] << ") " << endl;
				}
				if(diffK < realK){
					cout << "ERROR. In pat = [" << pat << "], t = " << t << ", kst = " << kst << ", realK = " << realK << ", App report only " << diffK << " documents !!" << endl;
					exit(1);
				}

				for(k=0; k<diffK; k++){
					if (fqCSA[docList[k]] == 0){
						cout << "ERROR. In pat = [" << pat << "], t = " << t << ", kst = " << kst << ", App report the document: " << docList[k] << " in position " << k << ", but it has not frequency!!" << endl;
						exit(1);
					}
					if (fqCSA[docList[k]] != frqList_test[k]){
						cout << "ERROR. In pat = [" << pat << "], t = " << t << ", kst = " << kst << ", fqCSA[docList[k]] = " << fqCSA[docList[k]] << " !=  frqList_test[k] = " << frqList_test[k] << endl;
						exit(1);
					}
				}
			}
		}
	}
	delete [] docList;
	delete [] docList_test;
	delete [] frqList_test;
	delete [] fqCSA;
	delete [] (par->patterns);
}

void createTopK_Test(ParamProgram *par, uint *docList, ulong *frqList, ulong *fqCSA, uint kst){
	ulong i, j, x, lenQ;
	//cout << "nDocs " << par->index->nDocs << ", kst=" << kst << endl;
	uint *maxQ = new uint[par->index->nDocs+1]; // queue for the D documents

	/*cout << "fqCSA:" << endl;
	for(i=0; i<par->index->nDocs; i++){
		cout << fqCSA[i] << " ";
	}
	cout << endl;*/

	for(i=0; fqCSA[i]==0; i++);
	maxQ[1] = i;	// first document with frequency > 0
	lenQ = 2;

	for(i++; i<par->index->nDocs; i++){
		// go down in maxQ
		if (fqCSA[i]){
			maxQ[lenQ] = i;
			j = lenQ;
			lenQ++;

			// j must to climb...
			while (j > 1 && fqCSA[maxQ[j/2]] < fqCSA[maxQ[j]]) {
				x = maxQ[j/2];
				maxQ[j/2] = maxQ[j];
				maxQ[j] = x;
				j /= 2;
			}
		}
	}

	/*cout << "maxQ:" << endl;
	for(i=1; i<lenQ; i++)
		cout << maxQ[i] << " ";
	cout << endl;*/

	// make answer...
	for(i=0; i<kst && lenQ>1; i++){
		docList[i] = maxQ[1];
		frqList[i] = fqCSA[maxQ[1]];
		lenQ--;
		j=2;
		x = maxQ[lenQ];

		// maxQ[j] must to go down...
		while (j < lenQ) {
			if(j+1<lenQ){
				if(fqCSA[maxQ[j+1]]>fqCSA[maxQ[j]])
					j++;
			}
			if(fqCSA[maxQ[j]] > fqCSA[x])
				maxQ[j/2] = maxQ[j];
			else{
				break;
			}
			j*=2;
		}
		maxQ[j/2] = x;
	}
	delete [] maxQ;
	//cout << "ok" << endl;

	/*cout << "**** docList:" << endl;
	for(i=0; i<par->kDocs && frqList[i]; i++)
		cout << docList[i] << " ";
	cout << endl;
	for(i=0; i<par->kDocs && frqList[i]; i++)
		cout << frqList[i] << " ";
	cout << endl;*/
}

void createPatterns(ParamProgram *par, uint m, ulong repeat){
	ulong i, j, k;
	bool eq;

	//cout << "Creating patterns of length m=" << m << endl;
	for(k=0; k<repeat; k++){
		eq = true;
		while(eq){
			i = (rand() % (par->n-(m+1)))+1;
			for(j=i; j<i+(m+1); j++){
				if (par->seq[j] == par->cutDoc){
					i = (rand() % (par->n-(m+1)))+1;
					j = i-1;
				}
				if(j==0) break;
				else{
					if (j>i && par->seq[j-1] != par->seq[j])
						eq = false;
				}
			}
		}
		par->patterns[k] = i;
		//cout << "["<<k<<"] " << i << endl;
	}
	//cout << "Patterns created !!" << endl;
}

// it founds all occurrences and count the frequencies by document
bool searchPatternInFMI(ParamProgram *par, uchar* pat, uint m, ulong *occt1, ulong *occt2, ulong *occt3, ulong *fqCSA){
	ulong r, o1, o2, o3, doc;
	string query = string((char *)pat);
	size_t occs = sdsl::count(par->index->fmi, query.begin(), query.begin()+m);
	auto locations = locate(par->index->fmi, query.begin(), query.begin()+m);
	if (TRACE) cout << endl << "Total occurrences found with FMI : " << occs << endl;

	o1 = o2 = o3 = 0;
	//cout << "locations..." << endl;
	for(ulong i=0; i<occs; i++){
		//cout << locations[i] << endl;
		r = par->index->sep_rank.rank(locations[i]+m) - par->index->sep_rank.rank(locations[i]);
		ulong rank = par->index->sep_rank.rank(locations[i]+1);
		doc = par->index->searchDocument(rank);
		(fqCSA[doc])++;

		if(TRACE){
			cout << par->index->sep_rrr[locations[i]] << par->index->sep_rrr[locations[i]+1] << par->index->sep_rrr[locations[i]+2] <<
					par->index->sep_rrr[locations[i]+3] << par->index->sep_rrr[locations[i]+4] << par->index->sep_rrr[locations[i]+5] <<
					par->index->sep_rrr[locations[i]+6] << par->index->sep_rrr[locations[i]+7] << par->index->sep_rrr[locations[i]+8] << endl;
			cout << "locations " << locations[i] << ",doc:" << doc << endl;
		}

		if(par->index->sep_rrr[locations[i]]){
			if (r==1){
				o1++;
				//cout << locations[i] << ",(t1)doc:" << doc << endl;
			}else{
				if (r==2){
					o2++;
					//cout << locations[i] << ",(t2)doc:" << doc << endl;
				}else{
					o3++;
					//cout << locations[i] << ",(t3)doc:" << doc << endl;
				}
			}
		}else{
			if (r==0){
				o1++;
				//cout << locations[i] << ",(t1)doc:" << doc << endl;
			}else{
				if (r==1){
					o2++;
					//cout << locations[i] << ",(t2)doc:" << doc << endl;
				}else{
					o3++;
					//cout << locations[i] << ",(t3)doc:" << doc << endl;
				}
			}
		}
	}
	// cout << endl;
	*occt1 = o1;
	*occt2 = o2;
	*occt3 = o3;

	return true;
}
