/*
 * loadAppTopkLZhq.cpp
 *
 *  Created on: 15-09-2014
 *      Author: hector
 */

#include "TopkLZ.h"
#include <ConfigFile.h>

bool PRINT = false;			// true: print all details for console
bool TEST_IND = true;			// true: apply exhaustive test
bool RUN_EXP = false;
uint REPEAT = 5;				// number of repetitions for each experiment
uint MAX_M = 10;			// maximum length pattern value to compute quality

// Structure with all globals parameters program
typedef struct {
	string configFile;		// properties file

	uchar *seq;				// original sequence (1 byte for symbol)
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	char cutDoc;			// symbol to separate documents
	bool lowerCPatt;		// 1: transform the patterns to lowercase
	uint g;

	TopkLZ *index;

	bool pattFromFile;			// 0: random patterns, 1: from file (as in todoCL)
	char dirStore[300];		// directory to save/load the data structure (files *.tk)
	char dirResult[300];	// directory to save tables

	char fileTodoCL1W[300];
	char fileTodoCL2W[300];

	ulong* patterns;
	uchar **patt;
	suint* lenPat;
} ParamProgram;

bool searchPatternInFMI_load(ParamProgram *par, uchar* pat, uint m, ulong *occt1, ulong *occt2, ulong *occt3, ulong *fqCSA);
void testSearchPatterns_CSA_LZ_load(ParamProgram *par);
void runExperimentsTwo(ParamProgram *par, uint k);
void runExperiments(ParamProgram *par);
void createPatterns_load(ParamProgram *par, uint m, ulong repeat);
void loadPatterns(ParamProgram *par, char *fileQueries);
void loadPatternsAll(ParamProgram *par, char *fileQueries);
void createTopK_Test_load(ParamProgram *par, uint *docList, ulong *frqList, ulong *fqCSA, uint kst);

int main(int argc, char *argv[]) {
	ParamProgram *par = new ParamProgram();
	char fileName[300];

	if(argc != 2){
		cout << "ERRORR !! " << endl;
		cout << "loadAppTopkLZhq's usage requires the Properties File as parameter !! " << endl;
		cout << "Example for the file 'config.txt' in the same directory: ./loadAppTopkLZhq config.txt" << endl;
		exit(1);
	}
	par->configFile = string(argv[1]);
	cout << "Congif file: " << par->configFile << endl;
	ConfigFile cf(par->configFile);

	PRINT = cf.Value("GLOBALS","TRACE");
	TEST_IND = cf.Value("GLOBALS","TEST");
	REPEAT = cf.Value("GLOBALS","N_REP");
	RUN_EXP = cf.Value("GLOBALS","RUN_EXP");
	MAX_M = cf.Value("GLOBALS","MAX_M");

	par->g = cf.Value("TOPK","g");
	strcpy(par->dirStore, ((string)(cf.Value("TOPK","dirStore"))).c_str());
	strcpy(par->dirResult, ((string)(cf.Value("TOPK","dirResult"))).c_str());
	par->pattFromFile = cf.Value("TOPK","pattFromFile");	// 0:random
	par->lowerCPatt = cf.Value("TOPK","lowercase");
	if(par->pattFromFile){
		strcpy(par->fileTodoCL1W, ((string)(cf.Value("TOPK","fileTodoCL1W"))).c_str());
		strcpy(par->fileTodoCL2W, ((string)(cf.Value("TOPK","fileTodoCL2W"))).c_str());
	}

	cout << "loadTopkLZ parameters..." << endl;
	cout << "dirStore: " << par->dirStore << endl;
	cout << "dirResult: " << par->dirResult << endl;
	cout << "patterns from file: " << par->pattFromFile << endl;
	cout << "lowercase patterns: " << par->lowerCPatt << endl;
	cout << "REPEATs: " << REPEAT << endl;
	if(par->pattFromFile){
		cout << "fileTodoCL1W: " << par->fileTodoCL1W << endl;
		cout << "fileTodoCL2W: " << par->fileTodoCL2W << endl;
	}

	par->index = new TopkLZ(par->dirStore, true);
	par->n = par->index->n;
	par->cutDoc = par->index->cutDoc;

	cout << "____________________________________________________" << endl;
	cout << "***  Index size " << par->index->sizeDS << " bytes = " << (float)par->index->sizeDS*8.0/(float)par->n << " bpc" << endl;
	cout << "====================================================" << endl;

	{
		// load Sequence to create random patterns...
		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sequence.test");
		ifstream is(fileName, ios::binary);
		par->seq = new uchar[par->n];
		is.read((char*)par->seq, par->n*sizeof(uchar));
		is.close();
	}

	if (TEST_IND){
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

		cout << "Test Index..." << endl;
		testSearchPatterns_CSA_LZ_load(par);
		cout << "Test Index OK !!" << endl;
	}

	if (RUN_EXP)
		runExperiments(par); // this store the patterns in par->patterns[]

	cout << "$$$$$$$$$$$$$$$$$$$$$" << endl;
	return 0;
}

// it founds all occurrences and count the frequencies by document
bool searchPatternInFMI_load(ParamProgram *par, uchar* pat, uint m, ulong *occt1, ulong *occt2, ulong *occt3, ulong *fqCSA){
	ulong r, o1, o2, o3, doc;
	string query = string((char *)pat);
	size_t occs = sdsl::count(par->index->fmi, query.begin(), query.begin()+m);
	auto locations = locate(par->index->fmi, query.begin(), query.begin()+m);
	if (PRINT) cout << endl << "Total occurrences found with FMI : " << occs << endl;

	o1 = o2 = o3 = 0;
	// cout << "locations..." << endl;
	for(ulong i=0; i<occs; i++){
		//cout << locations[i] << endl;
		r = par->index->sep_rank.rank(locations[i]+m) - par->index->sep_rank.rank(locations[i]);
		ulong rank = par->index->sep_rank.rank(locations[i]+1);
		doc = par->index->searchDocument(rank);
		(fqCSA[doc])++;

		//cout << "locations " << locations[i] << ",doc:" << doc << endl;

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


// ordena en 'sortedApp' la lista App por frecuencia real
void sortAppReal(uint *docList, ulong *fqCSA, uint diffK, ulong *sortedApp){
	uint i, k;
	ulong frec;

	sortedApp[0] = fqCSA[docList[0]];
	// insertion sort
	for (k=1; k<diffK; k++){
		frec = fqCSA[docList[k]];
		for (i=k; i>0 && sortedApp[i-1] < frec; i--)
			sortedApp[i] = sortedApp[i-1];
		sortedApp[i] = frec;
	}
}

void createTopK_Test_load(ParamProgram *par, uint *docList, ulong *frqList, ulong *fqCSA, uint kst){
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

void testSearchPatterns_CSA_LZ_load(ParamProgram *par){
	//PRINT = true;
	ulong j, k, t, pv, x, realK;
	uint m, l, M=20;
	ulong occ1, occ2, occ3, diffK, cWrong;
	bool isInRev;
	uchar *pat = new uchar[M+1];
	uint *docList = new uint[par->index->nDocs+1];
	uint *docList_test = new uint[par->index->nDocs+1];
	ulong *frqList_test = new ulong[par->index->nDocs+1];
	ulong *fqCSA = new ulong[par->index->nDocs+1];
	par->patterns = new ulong[REPEAT];

	uint MAX_K = 10, PLUS_K = 1, INI_K = 1;
	if (par->index->nDocs >= 100){
		MAX_K = 100;
		PLUS_K = 20;
		INI_K = 20;
	}else{
		if (MAX_K >= par->index->nDocs)
			MAX_K = par->index->nDocs/2;
	}

	//INI_K = 40;
	for (uint kst=INI_K; kst<=MAX_K; kst+=PLUS_K){
		cout << " Test for k = " << kst << "..." << endl;
		for (m=4; m<=M; m++){
			createPatterns_load(par, m, REPEAT);
			for (j=0; j<=par->index->nDocs; j++)
				par->index->dictD[j] = par->index->keyFq[j] = par->index->keyDoc[j] = 0;
			//cout << " ... test in patterns of length m = " << m << endl;
			for (t=cWrong=0; t<REPEAT; t++){
				for (l=0; l<m; l++)
					pat[l] = par->seq[par->patterns[t]+l];
				pat[m] = '\0'; // [ PW[somerv]
				/*pat[0] = '\n';
				pat[1] = 'P';
				pat[2] = 'W';
				pat[3] = '[';
				pat[4] = 's';
				pat[5] = 'o';
				pat[6] = 'm';
				pat[7] = 'e';
				pat[8] = 'r';
				pat[9] = 'v';*/

				//cout << "pat=[" << pat << "], m=" << m << ", test=" << t << ", k*=" << kst << endl;

				for(k=0; k<=par->index->nDocs; k++)
					fqCSA[k] = docList_test[k] = frqList_test[k] = 0;

				searchPatternInFMI_load(par, pat, m, &occ1, &occ2, &occ3, fqCSA);
				if (PRINT){
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
				createTopK_Test_load(par, docList_test, frqList_test, fqCSA, realK);
				if (PRINT){
					cout << "Real Docs and Frequencies..." << endl;
					for(k=0; k<realK && frqList_test[k]; k++)
						cout << k+1 << ", (" << docList_test[k] << "," << frqList_test[k] << ") " << endl;
				}

				diffK = par->index->topKDocument(pat, m, docList, realK);
				if (PRINT){
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
						cout << "ERROR. In pat = [" << pat << "], t = " << t << ", kst = " << kst << ", fqCSA["<<docList[k]<<"] = " << fqCSA[docList[k]] << " !=  frqList_test["<<k<<"] = " << frqList_test[k] << endl;
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


void createPatterns_load(ParamProgram *par, uint m, ulong repeat){
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


void loadPatternsAll(ParamProgram *par, char *fileQueries){
	uint i, j, k, len, totLines;
	string line;
	ifstream input (fileQueries);
	assert(input.good()); 				// check open
	if (input.is_open()){
		totLines = len = 0;
		while (input.good()){
			totLines++;
			getline(input,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			len++;
		}
		input.close();
	} else{
		cout << "Unable to open file " << fileQueries << endl;
		exit(0);
	}

	cout << "Total lines: " << totLines << ", patterns selected: " << len << endl;
	REPEAT = len;
	par->patt = new uchar*[len];
	par->lenPat = new suint[len];
	ifstream input2(fileQueries);// open file
	if (input2.is_open()){
		i = 0;
		while (input2.good() && i < len){
			getline(input2,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			par->patt[i] = new uchar[k+1];
			par->lenPat[i] = k;
			for(j=0; j<k; j++){
				if (par->lowerCPatt)
					par->patt[i][j] = (uchar)tolower(line.at(j));
				else
					par->patt[i][j] = (uchar)(line.at(j));
			}
			par->patt[i][k] = '\0';
			i++;
			//cout << par->patt[i] << endl;
		}
		input2.close();
	}
	cout << i << " patterns loaded !!" << endl;
}

void loadPatterns(ParamProgram *par, char *fileQueries){
	uint i, j, k, len, numspaces, totLines;
	string line;
	ifstream input (fileQueries);
	assert(input.good()); 				// check open
	if (input.is_open()){
		totLines = len =0;
		while (input.good()){
			totLines++;
			getline(input,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			for(j=numspaces=0; j<k; j++){
				if (isspace(line[j]))
					numspaces++;
			}
			if(numspaces <= 1)
				len++;
		}
		input.close();
	} else{
		cout << "Unable to open file " << fileQueries << endl;
		exit(0);
	}

	cout << "Total lines: " << totLines << ", patterns selected: " << len << endl;
	REPEAT = len;
	par->patt = new uchar*[len];
	par->lenPat = new suint[len];
	ifstream input2(fileQueries);// open file
	if (input2.is_open()){
		i = 0;
		while (input2.good() && i < len){
			getline(input2,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			for(j=numspaces=0; j<k; j++){
				if (isspace(line[j]))
					numspaces++;
			}
			if(numspaces <= 1){
				par->patt[i] = new uchar[k+1];
				par->lenPat[i] = k;
				for(j=0; j<k; j++){
					if (par->lowerCPatt)
						par->patt[i][j] = (uchar)tolower(line.at(j));
					else
						par->patt[i][j] = (uchar)(line.at(j));
				}
				par->patt[i][k] = '\0';
				i++;
				//cout << par->patt[i] << endl;
			}
		}
		input2.close();
	}
	cout << i << " patterns loaded !!" << endl;
}

void runExperimentsOne(ParamProgram *par, uint m, uint k){
	double t, avgTime=0.0;
	char aFile[300];
	uint *docList = new uint[par->index->nDocs];

	cout << "____________________________________________________" << endl;
	cout << "Start Top-k for " << REPEAT << " patterns of length m = " << m << " and k = " << k << endl;

	for (uint j=0; j<REPEAT; j++){
		t = getTime_ms();
		par->index->topKDocument(par->seq+par->patterns[j], m, docList, k);
		avgTime += getTime_ms() - t;
	}

	avgTime /= (double)REPEAT;
	cout << "Average CPU time for execution, Search type 1 : " << avgTime*1000.0 << " Microseconds" << endl;
	cout << "Size : " << par->index->sizeDS*8.0/(float)par->n << endl;
	cout << "____________________________________________________" << endl;

	strcpy(aFile, "");
	strcpy(aFile, par->dirStore);
	strcat(aFile, "resume.topk");
	FILE *fp = fopen(aFile, "a+" );
	// [m] [k] [g] [size] [time]
	fprintf(fp, "%d %d %d %f %G\n", m, k, par->index->g, par->index->sizeDS*8.0/(float)par->n, avgTime*1000.0);
	fclose(fp);
}

void runExperimentsTwo(ParamProgram *par, uint k){
	double t, avgTime=0.0;
	char aFile[300];
	uint m;
	uchar *patron;
	uint *docList = new uint[par->index->nDocs];
	float avgM = 0.0;

	cout << "____________________________________________________" << endl;
	cout << "Start Top-k for k = " << k << endl;

	for (uint j=0; j<REPEAT; j++){
		patron = par->patt[j];
		m = par->lenPat[j];
		t = getTime_ms();
		par->index->topKDocument(patron, m, docList, k);
		avgTime += getTime_ms() - t;
		avgM += m;
	}

	avgM /= (float)REPEAT;
	avgTime /= (double)REPEAT;
	cout << "Average CPU time for execution, Search type 1 : " << avgTime*1000.0 << " Microseconds" << endl;
	cout << "Average pattern length : " << avgM << endl;
	cout << "Size : " << par->index->sizeDS*8.0/(float)par->n << endl;
	cout << "____________________________________________________" << endl;

	strcpy(aFile, "");
	strcpy(aFile, par->dirResult);
	strcat(aFile, "resume.topk");
	FILE *fp = fopen(aFile, "a+" );
	// [avgM] [k] [g] [size] [time]
	fprintf(fp, "%f %d %d %f %G\n", avgM, k, par->index->g, par->index->sizeDS*8.0/(float)par->n, avgTime*1000.0);
	fclose(fp);
}

void runExperiments(ParamProgram *par){
	cout << " RUN TOP-K DOCUMENT LISTING..." << endl;
	cout << "====================================================" << endl;
	uint j, m;

	// This initialization is done only one time.
	for (j=0; j<par->index->nDocs; j++)
		par->index->dictD[j] = par->index->keyFq[j] = par->index->keyDoc[j] = 0;

	if (par->pattFromFile){
		cout << " Todocl for patterns with 1 word..." << endl;
		loadPatternsAll(par, par->fileTodoCL1W);
		for (uint kst=10; kst<=100; kst+=90)
			runExperimentsTwo(par, kst);
		delete [] par->patt;
		delete [] par->lenPat;

		cout << " Todocl for patterns with 2 words..." << endl;
		loadPatternsAll(par, par->fileTodoCL2W);
		for (uint kst=10; kst<=100; kst+=90)
			runExperimentsTwo(par, kst);
		delete [] par->patt;
		delete [] par->lenPat;

	}else{
		par->patterns = new ulong[REPEAT];
		m=6;
		uint kst=10;
		createPatterns_load(par, m, REPEAT);
		runExperimentsOne(par, m, kst);

		m=10;
		createPatterns_load(par, m, REPEAT);
		runExperimentsOne(par, m, kst);

		if (par->index->nDocs > 100){
			kst=100;
			m=6;
			createPatterns_load(par, m, REPEAT);
			runExperimentsOne(par, m, kst);
			m=10;
			createPatterns_load(par, m, REPEAT);
			runExperimentsOne(par, m, kst);
		}

		delete [] (par->patterns);
	}
}

