/*
 * TopkLZ.cpp
 *
 *  Created on: 19-08-2015
 *      Author: hector
 */

#include "TopkLZ.h"
bool TopkLZ::TRACE = false;
bool TopkLZ::TEST = true;
const uint MAX_DEP = 100000;
const uint MAX_LEN_DOCS = 40; // it is 40 bpc

TopkLZ::TopkLZ(char dirSaveLoad[200], bool showSize){
	strcpy(dirStore, dirSaveLoad);
	loadDS(showSize);

	ulong size = nDocs+1;
	this->dictD = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	if (showSize) cout << " ** size of dictD[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->keyDoc = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	if (showSize) cout << " ** size of keyDoc[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->keyFq = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	if (showSize) cout << " ** size of keyFq[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->heap = new uint[size];
	size *= sizeof(uint);
	sizeDS += size;
	if (showSize) cout << " ** size of heap[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

}

TopkLZ::TopkLZ(uint gVal, char *inputFile, bool filesInList, char cutDocCode, bool bitLowercase, char dirSaveLoad[300], char boundS) {
	ulong i, j, sizeAux;
	char fileName[400];

	// size for variables
	this->sizeDS = 9*sizeof(ulong) + 4*sizeof(uint) + sizeof(char);
	this->sizeDSav = 0;

	this->g = gVal;
	this->cutDoc = cutDocCode;
	this->boundSymbol = boundS;
	strcpy(dirStore, dirSaveLoad);
	this->nDocs = this->n = this->nEFHead = 0;

	charTf = new ulong[LZ78Tries64::SIGMA];
	for(i=0; i<LZ78Tries64::SIGMA; i++)
		charTf[i] = 0;
	if (filesInList){
		readListFiles(inputFile, bitLowercase);
	}else
		readUniqueFile(inputFile, bitLowercase);

	ulong* cPosTb = new ulong[LZ78Tries64::SIGMA]; // hasTable position for each character
	makePriorities(cPosTb);

	// create the FMI, needed at time construction, to compute the real top-k list for market nodes.
	cout << "____________________________________________________" << endl;
	cout << " Make fmi (for top-k brute force)..." << endl;
	cout << " Reading... " << inputFile << endl;
	strcpy(fileName, "");
	strcpy(fileName, inputFile);
	cout << " Reading... " << inputFile << endl;
	strcat(fileName, "_copy.txt");
	construct(fmi, fileName, 1); // generate index
	cout << " **  FMI size " << size_in_bytes(fmi) << " bytes = " << (float)size_in_bytes(fmi)/(float)n << "|T|" << endl;
	cout << " **  FMI length " << fmi.size() << endl;
	if (fmi.size() != n){
		cout << "ERROR. FMI length != n = " << n << endl;
		exit(1);
	}
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fmi.temp");
	store_to_file(fmi, fileName);

	this->EndDocs = new ulong[nDocs];
	this->EndDocsPos = new ulong[nDocs];	// used only at construction time

	//cout << "  Create Original LZTrie and RevTrie..." << endl;
	tries = new LZ78Tries64(nDocs, cutDoc, cPosTb);
	this->separator_b = bit_vector(n, 0);

	tries->createTriesFromText(seq, n, EndDocs, EndDocsPos, &separator_b);
	cout << "  LZTrie with " << tries->nPhra << " nodes (phrases)" << endl;
	this->lgD = ceilingLog64(nDocs, 2);
	this->nod = tries->nPhra;
	this->lgNod = ceilingLog64(nod, 2);

	sep_rrr = rrr_vector<127>(separator_b);
	sep_rank = rrr_vector<127>::rank_1_type(&sep_rrr);
	sep_sel = rrr_vector<127>::select_1_type(&sep_rrr);

	if (TRACE){
		cout << " separator_b[0.." << sep_rrr.size()-1 << "]:" << endl;
		cout << sep_rrr[0];
		for(i=1; i<sep_rrr.size(); i++)
			cout << sep_rrr[i];
		cout << endl;
	}

	{	// save seq[] for testing and delete it !!
		strcpy(fileName, "");
		strcpy(fileName, dirSaveLoad);
		strcat(fileName, "sequence.test");
		ofstream osSeq (fileName, ios::binary);
		osSeq.write((const char*)seq, n*sizeof(uchar));
		if(TRACE) cout << " .- Seq[] " << n*sizeof(uchar) << " Bytes saved" << endl;
		osSeq.close();
		delete []seq;

		strcpy(fileName, "");
		strcpy(fileName, dirSaveLoad);
		strcat(fileName, "EndDocs.test");
		ofstream osEndDocs (fileName, ios::binary);
		osEndDocs.write((const char*)EndDocs, nDocs*sizeof(ulong));
		osEndDocs.close();

		// save sep_rrr and sep_rank for testing and delete them !!
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "sep_rrr.test");
		store_to_file(sep_rrr, fileName);

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "sep_sel.test");
		store_to_file(sep_sel, fileName);

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "sep_rank.test");
		store_to_file(sep_rank, fileName);

		// delete separator_b and sep_rrr
		{
		    bit_vector empty_b;
		    separator_b.swap(empty_b);

		    rrr_vector<127> empty_rrr;
		    sep_rrr.swap(empty_rrr);
		}
	}

	sizeAux = tries->nPhra*lgNod/W64;
	if ((tries->nPhra*lgNod)%W64)
		sizeAux++;
	Node = new ulong[sizeAux];
	for(i=0; i<sizeAux; i++)
		Node[i] = 0;
	sizeAux *= sizeof(ulong);
	sizeDS += sizeAux;
	cout << " ** Size of Node array " << sizeAux << " = " << sizeAux*8.0/(float)n << " bpc" << endl;
	i = 0;
	cout << "To generate NodeArray..."<< endl;
	ulong *ArrRange = new ulong[tries->nPhra];	// it represents the pairs (j,i) <--> If node j in RevTrie has id=k, then node i in LZTrie has id=k+1
	genNodeArray(tries->revTrie, ArrRange, &i);
	if (tries->nPhra != i){
		cout << "total RevTrie preorder's for true nodes = " << i << " != nPhra = " << tries->nPhra << endl;
		exit(1);
	}

	if(TRACE){
		cout << "ArrRange[1.." << tries->nPhra << "]:" << endl;
		cout << ArrRange[0] << " ";
		for(i=1; i<tries->nPhra; i++){
			if(i%10 == 0)
				cout << "- ";
			cout << ArrRange[i] << " ";
		}
		cout << endl;

		cout << " Node array[0.." << nod-1 << "]:" << endl;
		cout << getNum64(Node, 0, lgNod) << " ";
		for(i=1, j=lgNod; i<nod; i++, j+=lgNod){
			if(i%10 == 0)
				cout << "- ";
			cout << getNum64(Node, j, lgNod) << " ";
		}
		cout << endl;
	}

	{	// to save and delete Node...
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "Node.tk");
		ofstream osNode (fileName, ios::binary);
		ulong sizeSaved = nod*lgNod/W64;
		if ((nod*lgNod)%W64)
			sizeSaved++;
		sizeSaved *= sizeof(ulong);
		sizeDSav += sizeSaved;
		osNode.write((const char*)Node, sizeSaved);				// save Node[]
		osNode.close();
		if(TRACE) cout << " .- Node[] " << sizeSaved << " Bytes saved" << endl;
		delete [] Node;
	}

	cout << "Creating Range structure..." << endl;
	createRange(ArrRange);
	if(TRACE){
		cout << " ## Range[0.." << h-1 << "]x[0.." << nod-1 << "] ..." << endl;
		for(i=0; i<h; i++){
			cout << Range[i].B_rrr[0];
			for(j=1; j<nod; j++){
				if(j%10 == 0) cout << "-";
				cout << Range[i].B_rrr[j];
			}
			cout << endl;
		}
		cout << endl;
	}

	{	// to save Range and delete it...
		char str[100];
		cout << "  .- Range[]'s pointers " << h*sizeof(LevelRange) << " = " << (h*sizeof(LevelRange))*8.0/(float)this->n << " bpc" << endl;
		ulong sizeBitMapsRange = 0;
		//cout << "h = " << h << endl;

		for(uint i=0; i<this->h; i++){
			strcpy(fileName, "");
			strcpy(fileName, dirStore);
			strcpy(str, "");
			sprintf(str, "Range_%d.B_rrr.tk", i);
			strcat(fileName, str);
			store_to_file(Range[i].B_rrr, fileName);
			sizeBitMapsRange += size_in_bytes(Range[i].B_rrr);
			{
				decltype(Range[i].B_rrr) empty;
				Range[i].B_rrr.swap(empty);
			}

			strcpy(fileName, "");
			strcpy(fileName, dirStore);
			strcpy(str, "");
			sprintf(str, "Range_%d.B_rank0.tk", i);
			strcat(fileName, str);
			store_to_file(Range[i].B_rank0, fileName);

			strcpy(fileName, "");
			strcpy(fileName, dirStore);
			strcpy(str, "");
			sprintf(str, "Range_%d.B_rank1.tk", i);
			strcat(fileName, str);
			store_to_file(Range[i].B_rank1, fileName);
		}
		sizeDSav += sizeBitMapsRange;
		if(TRACE) cout << " .- Range[]'s BitVector(RRR): " << sizeBitMapsRange << " Bytes saved" << endl;
		delete [] Range;
	}

	tries->countFictNodRevTrie(tries->revTrie);	// this set tries->nFictU
	cout << "  RevTrie with " << tries->nPhra << " nodes, " << tries->nFictU << " nFictU nodes and " << tries->nExtraFict << " nExtraFict nodes" << endl;

	this->nF = tries->nFictU;
	this->nEF = tries->nExtraFict;
	this->nodRev = nod+nF;
	this->nRev = 2*nodRev;
	this->nLz = 2*nod;

	//tries->listLZTrie();
	//exit(0);
	if (TRACE){
		cout << endl << " LZTrie with n'= " << tries->nPhra << " nodes..." << endl;
		tries->listLZTrie();
		cout << "====================================================" << endl;
		//cout << endl << " RevTrie with n'= " << tries->nPhra << " phrase nodes, nFictU = " << tries->nFictU << " fictitious nodes, and " << tries->nExtraFict << " extra fictitious nodes (nExtraFict)... " << endl;

		cout << " EndDocs[0.." << nDocs-1 << "]:" << endl;
		for(i=0; i<nDocs; i++)
			cout << EndDocs[i] << " ";
		cout << endl;
	}

	cout << "  DFUFS RevTrie seq with " << nRev << " bits" << endl;
	cout << "  DFUFS LZTrie seq with " << nLz << " bits" << endl;

	sizeAux = nLz/W64;
	if (nLz%W64)
		sizeAux++;
	this->PLz = new ulong[sizeAux];		// the size of PLz are included in its RMM representatio

	sizeAux = nRev/W64;
	if (nRev%W64)
		sizeAux++;
	this->PRev = new ulong[sizeAux];	// the size of PRev are included in its RMM representatio

	// set initial open parenthesis in DFUDS sequences
	setBit64(PLz, 0);
	setBit64(PRev, 0);

	this->LbRev = new uchar[nodRev];
	cout << " length of LbRev  " << nodRev << endl;
	sizeAux = nodRev*sizeof(uchar);
	sizeDS += sizeAux;
	cout << " ** size of LbRev[]: " << sizeAux << " = " << sizeAux*8.0/(float)n << " bpc" << endl;
	cout << " LbRev lenght " << nodRev << endl;

	this->LbRevF = new uchar[nEF];
	//cout << "## len of LbRevF  " << nEF << endl;
	sizeAux = nEF*sizeof(uchar);
	sizeDS += sizeAux;
	cout << " ** size of LbRevF[]: " << sizeAux << " = " << sizeAux*8.0/(float)n << " bpc" << endl;
	cout << " LbRevF lenght " << nEF << endl;

	sizeAux = nDocs+1;
	this->dictD = new ulong[sizeAux];
	sizeAux *= sizeof(ulong);
	sizeDS += sizeAux;
	cout << " ** size of dictD[1..nDocs]: " << sizeAux << " = " << sizeAux*8.0/(float)this->n << " bpc" << endl;

	sizeAux = nDocs+1;
	this->keyDoc = new ulong[sizeAux];
	sizeAux *= sizeof(ulong);
	sizeDS += sizeAux;
	cout << " ** size of keyDoc[1..nDocs]: " << sizeAux << " = " << sizeAux*8.0/(float)this->n << " bpc" << endl;

	sizeAux = nDocs+1;
	this->keyFq = new ulong[sizeAux];
	sizeAux *= sizeof(ulong);
	sizeDS += sizeAux;
	cout << " ** size of keyFq[1..nDocs]: " << sizeAux << " = " << sizeAux*8.0/(float)this->n << " bpc" << endl;

	sizeAux = nDocs+1;
	this->heap = new uint[sizeAux];
	sizeAux *= sizeof(uint);
	sizeDS += sizeAux;
	cout << " ** size of heap[1..nDocs]" << sizeAux << " = " << sizeAux*8.0/(float)this->n << " bpc" << endl;

	createStructure(dirSaveLoad);

	cout << "Saving the hole structure..." << endl;
	saveDS(true);
	destroyConstruct();

	cout << " PERFECT !!" << endl;
}

void TopkLZ::createStructure(char dirSaveLoad[300]){
	char fileName[400];
	ulong i, j, pos, posRev, posF, posExtF, posFU;
	cout << "  Create LZ Top-k list Structure from RevTrie..." << endl;

	// Fictitious nodes...
	{
		//[1] Create the DFUDS sequences of parentheses PRev, the bitstring FRev to mark the fictitious nodes, and the label vectors LbRev and LbRevF...
		fRev_b = bit_vector(nRev/2, 0);
		//cout << "lenght of fRev_b (posF) = " << fRev_b.bit_size() << endl;
		fictU_b = bit_vector(nF, 0);
		//cout << "lenght of fictU_b (posFU) = " << fictU_b.bit_size() << endl;
		listF_b = bit_vector(nEF, 0);
		//cout << "lenght of listF_b (posExtF) = " << listF_b.bit_size() << endl;

		posF = posExtF = posFU = 0;
		genSeqExtraFictRevTrie(tries->revTrie, &posExtF, &posF, &posFU);
		//cout << "posExtF " << posExtF << ", posF " << posF << ", posFU " << posFU << endl;
		//cout << "posF = " << posF << " ? nRev/2 = " << nRev/2 << endl;
		//cout << "genSeqExtraFictRevTrie OK" << endl;

		fRev_rrr = rrr_vector<127>(fRev_b);
		sizeDS += size_in_bytes(fRev_rrr);
		cout << " ** size of fRev_rrr: " << size_in_bytes(fRev_rrr) << " = " << size_in_bytes(fRev_rrr)*8.0/(float)n << " bpc" << endl;
		fRev_rank = rrr_vector<127>::rank_1_type(&fRev_rrr);
		{
		    bit_vector empty_b;
		    fRev_b.swap(empty_b);
		}

		fictU_rrr = rrr_vector<127>(fictU_b);
		sizeDS += size_in_bytes(fictU_rrr);
		cout << " ** size of fictU_rrr: " << size_in_bytes(fictU_rrr) << " = " << size_in_bytes(fictU_rrr)*8.0/(float)n << " bpc" << endl;
		fictU_rank = rrr_vector<127>::rank_1_type(&fictU_rrr);
		{
			bit_vector empty_b;
			fictU_b.swap(empty_b);
		}

		listF_rrr = rrr_vector<127>(listF_b);
		sizeDS += size_in_bytes(listF_rrr);
		cout << " ** size of listF_rrr: " << size_in_bytes(listF_rrr) << " = " << size_in_bytes(listF_rrr)*8.0/(float)n << " bpc" << endl;
		lisF_sel = rrr_vector<127>::select_1_type(&listF_rrr);
		{
			bit_vector empty_b;
			listF_b.swap(empty_b);
		}

	}

	{	// LABELS
		pos = 1; posRev = 0;
		genDFUDSSeqRevTrie(tries->revTrie, &pos, &posRev);
		//cout << "pos " << pos << ", posRev (LbRev) " << posRev << endl;

		pos = 1;
		genDFUDSSeqLZTrieNotLbl(tries->lzTrie, &pos);
		//cout << "pos (Plz) " << pos << " =? nLz = " << nLz <<endl;

		if(TEST){ // test for balanced sequences
			long int sum = 0;
			ulong r=0;
			for (; r<nLz; r++){
				if(readBit64(PLz, r))
					sum++;
				else sum--;
			}
			if(sum != 0){
				cout << " PLz DFUDS is not a balanced sequence of parentheses !! " << endl;
				exit(1);
			}
			//else cout << " Test for PLz. DFUDS is a well balanced sequence of parentheses !! " << endl;
			sum = r = 0;
			for (; r<nRev; r++){
				if(readBit64(PRev, r))
					sum++;
				else sum--;
			}
			if(sum != 0){
				cout << " PRev DFUDS is not a balanced sequence of parentheses !! " << endl;
				exit(1);
			}
			//else cout << " Test for PRev. DFUDS is a well balanced sequence of parentheses !! " << endl;
		}

		//tries->listRevTrie();
		//exit(0);
		if(TRACE){
			cout << "RevTrie with " << nod << " nodes, " << tries->nFictU << " nFictU nodes and " << tries->nExtraFict << " nExtraFict nodes" << endl;
			tries->listRevTrie();

			cout << " DFUDFS PLz[0.." << nLz-1 << "]:" << endl;
			cout << readBit64(PLz, 0);
			for(i=1; i<nLz; i++){
				if(i%10 == 0)
					cout << "-";
				cout << readBit64(PLz, i);
			}
			cout << endl;

			cout << " DFUDFS PRev[0.." << nRev-1 << "]:" << endl;
			cout << readBit64(PRev, 0);
			for(i=1; i<nRev; i++){
				if(i%10 == 0)
					cout << "-";
				cout << readBit64(PRev, i);
			}
			cout << endl;
			cout << " fRev[0.." << fRev_rrr.size()-1 << "]:" << endl;
			cout << fRev_rrr[0];
			for(i=1; i<fRev_rrr.size(); i++){
				if(i%10 == 0)
					cout << "-";
				cout << fRev_rrr[i];
			}
			cout << endl;
			cout << " fictU[0.." << fictU_rrr.size()-1 << "]:" << endl;
			cout << fictU_rrr[0];
			for(i=1; i<fictU_rrr.size(); i++){
				if(i%10 == 0)
					cout << "-";
				cout << fictU_rrr[i];
			}
			cout << endl;
			cout << " LbRev[0.." << nodRev-1 << "]:" << endl;
			for(i=0; i<nodRev; i++)
				cout << LbRev[i];

			cout << endl;
			if (nEF){
				cout << " LbRevF[0.." << nEF-1 << "]:" << endl;
				for(i=0; i<nEF; i++)
					cout << LbRevF[i];
				cout << endl;
				cout << " listF[0.." << listF_rrr.size()-1 << "]:" << endl;
				for(i=0; i<listF_rrr.size(); i++)
					cout << listF_rrr[i];
				cout << endl;
			}else
				cout << " LbRevF does not have any extra labels (unary path of fictitious node with length > 1)" << endl;
		}
	}


	{	//	TREES !!
		cout << "  Create the representation with Range min-max of PLz and PRev DFUDS..." << endl;
		treeLz = new RangeMMTree64(PLz, nLz, NULL, false);
		sizeDS += treeLz->sizeRMM;
		cout << " ** size of treeLz: " <<  treeLz->sizeRMM << " = " <<  treeLz->sizeRMM*8.0/(float)n << " bpc" << endl;

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "LzTrie.tk");
		cout << "***   Saving LzTrie... " << endl;

		treeLz->saveDS(fileName, false);
		sizeDSav += treeLz->sizeRMM;
		if(TRACE) cout << " .- treeLz: " << treeLz->sizeRMM << " Bytes saved" << endl;
		cout << "Destroying treeLZ..." << endl;
		treeLz->~RangeMMTree64();

		treeRev = new RangeMMTree64(PRev, nRev, LbRev, false);
		sizeDS += treeRev->sizeRMM;
		cout << " ** size of treeRev: " <<  treeRev->sizeRMM << " = " <<  treeRev->sizeRMM*8.0/(float)n << " bpc" << endl;

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "RevTrie.tk");
		cout << "***   Saving RevTrie... " << endl;
		treeRev->saveDS(fileName, false);
		sizeDSav += treeRev->sizeRMM;
		if(TRACE) cout << " .- treeRev: " << treeRev->sizeRMM << " Bytes saved" << endl;
	}

	{	// DOCUMENT ARRAY OF PHRASES DocLZ...
		cout << "  Create DocLZ..." << endl;
		ulong size = nod*lgD;
		if(size%W64)
			size = size/W64 + 1;
		else
			size /= W64;
		this->DocLZ = new ulong[size];
		size *= sizeof(ulong);
		sizeDS += size;
		cout << " ** size of DocLZ[1..n']: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

		pos = 0;
		// making DocLZ in LZ78Tries
		genDocArray(tries->lzTrie, DocLZ, &pos);

		{	// free memory do not need it more
			delete [] EndDocs;
		}

		cout << "  Document array[0.." << nod-1 << "] and " << pos << " bits" << endl;
		if(TRACE){
			cout << "  Document array[0.." << nod-1 << "]:" << endl;
			for(i=j=0; i<nod; i++, j+=lgD){
				if(i%10 == 0) cout << "-";
				cout << getNum64(DocLZ, j, lgD);
			}
			cout << endl;
		}

		// To save and delete document array...
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "DocLZ.tk");
		ofstream osDocLZ (fileName, ios::binary);
		osDocLZ.write((const char*)DocLZ, size);				// save DocLZ[]
		osDocLZ.close();
		sizeDSav += size;
		if(TRACE) cout << " .- DocLZ[]: " << size << " Bytes saved" << endl;
		delete [] DocLZ;
	}

	{
		cout << "Destroying Original tries with pointers..." << endl;
		tries->~LZ78Tries64();
	}

	//////////////////////////////////////////////////// TOP-K STRUCTURE ////////////////////////////////////////////////////
	//cout << " Create the top-k list in each node..." << endl;
	bTopk_b = bit_vector(nRev/2+1, 0);
	cout << "lenght of bTopk_b = " << bTopk_b.bit_size() << endl;

	uint* dNodK = new uint[nDocs];		// array of K most frequent id's
	uint* fNodK = new uint[nDocs];		// array of K frequencies
	uint* fNod = new uint[nDocs];
	ulong lenD, nDTot, len, pre, posLb;
	pos = nDTot = lenD = 0;

	len = MAX_LEN_DOCS*n;					// Maximum length allowed is MAX_LEN_DOCS*n bits
	if(len%W64)
		len = len/W64 + 1;
	else
		len /= W64;
	ulong *docTopkAux = new ulong[len];		// Auxiliary array for a temporal store of all doc's Id
	cout << "largo(W) de docTopkAux = " << len << endl;

	len = MAX_LEN_DOCS*n;
	if(len%lgD)
		len = len/lgD + 1;
	else
		len /= lgD;
	bool *limDTk_bAux = new bool[len]; 		// auxiliary bitstring to mark list boundaries
	cout << "Maximum number of doc's Id allowed for limDTk_bAux = " << len << "; lgD = " << lgD << endl;
	for(i=0; i<len; i++)
		limDTk_bAux[i] = 0;

	len = nRev/2+1;
	bool *bDL_bAux = new bool[len]; 		// auxiliary bitstring to mark DL answer
	for(i=0; i<len; i++)
		bDL_bAux[i] = 0;
	uchar *patt = new uchar[MAX_DEP];
	uint m=0;
	pre = pos = 1; // pos count only the n' true nodes
	posLb = lenD = nDTot = 0;
	createDocTopkSeq(&pos, &posLb, &pre, patt, m, fNod, dNodK, fNodK, docTopkAux, &lenD, &nDTot, bDL_bAux, limDTk_bAux);

	cout << "After createDocTopkSeq() :" << endl;
	cout << "pos :" << pos << ", posLb :" << posLb << ", pre :" << pre << ", lenD :" << lenD << ", nDTot :" << nDTot << endl;
	delete [] dNodK;
	delete [] fNodK;
	delete [] fNod;

	{	// to store the documents word by word from docTopkAux to docTopk...
		len = lenD*lgD;
		if(len%W64)
			len = len/W64 + 1;
		else
			len /= W64;
		docTopk = new ulong[len];
		for(i=0; i<len; i++)
			docTopk[i] = docTopkAux[i];
		delete [] docTopkAux;
		len *= sizeof(ulong);
		sizeDS += len;
		cout << " ** size of docTopk: " << len << " = " << len*8.0/(float)n << " bpc" << endl;

		// To save and delete sequence of document stored...
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcat(fileName, "docTopk.tk");
		ofstream osdocTopk (fileName, ios::binary);
		osdocTopk.write((const char*)docTopk, len);				// save DocLZ[]
		osdocTopk.close();
		sizeDSav += len;
		if(TRACE) cout << " .- docTopk[]: " << len << " Bytes saved" << endl;
	}

	{	// set sublists boundaries
		limDTk_b = bit_vector(lenD, 0);
		for(i=0; i<lenD; i++)
			limDTk_b[i] = limDTk_bAux[i];
		delete [] limDTk_bAux;

		limDTk_rrr = rrr_vector<127>(limDTk_b);
		limDTk_rank = rrr_vector<127>::rank_1_type(&limDTk_rrr);
		limDTk_sel = rrr_vector<127>::select_1_type(&limDTk_rrr);
		sizeDS += size_in_bytes(limDTk_rrr);
		cout << " ** size of limDTk_rrr: " << size_in_bytes(limDTk_rrr) << " = " << size_in_bytes(limDTk_rrr)*8.0/(float)n << " bpc" << endl;
		cout << " limDTk has " << limDTk_rank.rank(limDTk_rrr.size()) << " 1's " << endl;
		{
			bit_vector empty_b;
			limDTk_b.swap(empty_b);
		}
	}

	{	// RRR and rank support to marked documents bitstring
		// bTopk has nTopk 1's and z 0's, where n'=nTopk+z
		bTopk_rrr = rrr_vector<127>(bTopk_b);
		sizeDS += size_in_bytes(bTopk_rrr);
		cout << " ** size of bTopk_rrr: " << size_in_bytes(bTopk_rrr) << " = " << size_in_bytes(bTopk_rrr)*8.0/(float)n << " bpc" << endl;

		bTopk_rank = rrr_vector<127>::rank_1_type(&bTopk_rrr);
		nTopk = bTopk_rank.rank(bTopk_rrr.size());
		cout << "   bTopk has " << nTopk << " nodes marked (" << ((float)nTopk/(float)nod)*100.0 << "% of total = " << bTopk_b.bit_size() << ")" << endl;

		{
		    bit_vector empty_b;
		    bTopk_b.swap(empty_b);
		}
	}

	{	// DL marks for each sublist of documents stored
		bDL_b = bit_vector(nTopk, 0);
		for(i=1; i<nTopk; i++)
			bDL_b[i] = bDL_bAux[i];
		delete [] bDL_bAux;
		sizeDS += size_in_bytes(bDL_b);
		cout << " ** size of bDL_b: " << size_in_bytes(bDL_b) << " = " << size_in_bytes(bDL_b)*8.0/(float)n << " bpc" << endl;
	}

	if(TRACE){
		len = bTopk_rrr.size();
		cout << " bTopk[0.." << len-1 << "]:" << endl;
		cout << bTopk_rrr[0];
		for(i=1; i<len; i++){
			if(i%10 == 0)
				cout << "-";
			cout << bTopk_rrr[i];
		}
		cout << endl;

		cout << " bDL[0.." << nTopk << "]:" << endl;
		cout << bDL_b[0];
		for(i=1; i<nTopk; i++){
			if(i%10 == 0)
				cout << "-";
			cout << bDL_b[i];
		}
		cout << endl;

		cout << " limDTk[0.." << limDTk_rrr.size()-1 << "]:" << endl;
		cout << limDTk_rrr[0];
		for(i=1; i<limDTk_rrr.size(); i++){
			if(i%10 == 0)
				cout << "-";
			cout << limDTk_rrr[i];
		}
		cout << endl;

		ulong pv, r, l;
		cout << "documents stored as Topk-List:" << endl;
		for(uint pre=0; pre<bTopk_rrr.size(); pre++){
			if (bTopk_rrr[pre]){
				pv = bTopk_rank.rank(pre+1);

				l = limDTk_sel.select(pv);		// start position of this block of id documents
				if(pv<nTopk)
					r = limDTk_sel.select(pv+1);
				else
					r = limDTk_rrr.size();

				for(;l<r; l++){
					cout << getNum64(docTopk, l*lgD, lgD) << " ";
				}
				cout << endl;
			}
		}
	}
}

// create Node array
void TopkLZ::genNodeArray(LZ78Tries64::RevNode* nod, ulong *ArrRange, ulong *preRev){
	LZ78Tries64::RevNode* p = nod->fChild;
	if (nod->idNode == 0){				// mejora: no preguntar por esto, hacerlo al inicio y luego llamar al metodo
		setNum64(Node, 0, lgNod, 0);
		ArrRange[0] = tries->IdDocPreLZ[1];
		(*preRev)++;
	}else{
		if(nod->fict == false){
			setNum64(Node, (*preRev)*lgNod, lgNod, nod->nodeLZ->preorder);

			// store a new point for range structure...
			if (nod->idNode == (int)(tries->nPhra-1))
				ArrRange[*preRev] = tries->IdDocPreLZ[0];
			else
				ArrRange[*preRev] = tries->IdDocPreLZ[nod->idNode+1];

			(*preRev)++;
		}
	}

	for(uint i=0; i<nod->nChildren; i++){
		genNodeArray(p, ArrRange, preRev);
		p = p->nextSib;
	}
	if(p){
		cout << "ERROR in RevTrie: nod->nChildren = " << nod->nChildren << ". But there is another node p = " << p->idNode << endl;
		exit(1);
	}
}

// 'pos' is the osition in PRev[] topology. 'pre' is the preorder value of the current node
void TopkLZ::createDocTopkSeq(ulong *pP, ulong *pL, ulong *pre, uchar *patt, uint m, uint *fNod, uint *dNodK, uint *fNodK,
								ulong *docTopkAux, ulong *posD, ulong *nodTopk, bool *bDL_bAux, bool *limDTk_bAux){
	uint nCh, i;
	ulong sig0, nOcc, iniLb=*pL, iniP=*pP;

	if (*pre >= bTopk_b.size()){
		cout << "ERR: *pre  = " << *pre << " >= bTopk_b.size() = " << bTopk_b.size() << endl;
		exit(0);
	}

	if(readBit64(PRev, iniP)){
		sig0 = treeRev->selectNext0(iniP);
		nCh = sig0 - iniP;	// Number of children for 'pre'
		*pP = sig0+1;
	}else{
		nCh = 0;
		(*pP)++;
	}

	if(*pre > 1 && nCh > LZ78Tries64::SIGMA){
		cout << "ERR: *pre  = " << *pre << " nCh = " << nCh << " > SIGMA" << endl;
		exit(0);
	}
	//cout << "pre " << *pre << ", nCh " << nCh << ", pP " << *pP << ", pL " << *pL << ", m=" << m << ", path = [" << patt << "]" << endl;
	//if(*pre == 17)
		//cout << "";

	*pL += nCh;
	for(i=0; i<nCh; i++){
		patt[m] = treeRev->labels[iniLb+i];
		if (fRev_rrr[*pre] == false && patt[m] != cutDoc){
			// reverse pattern in auxPath[]
			uchar *auxpatt = new uchar[m+2];
			for(uint k=0; k<=m; k++)
				auxpatt[k] = patt[m-k];
			auxpatt[m+1] = '\0';
			//cout << endl << "pre = " << *pre << ", auxPath = [" << auxpatt << "], DOCS: ";
			string query = string((char *)auxpatt);
			nOcc = sdsl::count(fmi, query.begin(), query.begin()+m+1);

			if(nOcc > g){
				ulong j, K, kStar;

				bTopk_b[*pre] = 1;		// to mark it
				for (K=0; K<nDocs; K++)
					fNod[K] = 0;

				// count frequencies by document
				freqPattInFMI(auxpatt, m+1, nOcc, fNod);
				for (K=j=0; K<nDocs; K++){
					if (fNod[K])
						j++;		// number of document found
				}
				kStar = log(nOcc/g) / log(2);
				kStar = pow(2, kStar);
				if(g*kStar <= nOcc)
					kStar <<= 1;
				if (kStar > nDocs)
					kStar = nDocs;
				if (j <= kStar){
					// to mark DL for this node...
					bDL_bAux[*nodTopk] = 1;
					kStar = j;
				}
				if((*posD+kStar)*lgD >= MAX_LEN_DOCS*n){
					cout << "ERROR.... docTopk requires more than " << MAX_LEN_DOCS << " bpc !!" << endl;
					exit(0);
				}
				sortTopkDocsByFq(fNod, kStar , dNodK, fNodK);

				// to mark the first position in this sublist of documents
				limDTk_bAux[*posD] = 1;

				//cout << "pre " << *pre << ", nCh " << nCh << ", pP " << *pP << ", pL " << *pL << endl << "DOCS: ";
				for (K=0, j=(*posD)*lgD; K<kStar; K++, j+=lgD){
					if (TRACE) cout << dNodK[K] << " ";
					setNum64(docTopkAux, j, lgD, dNodK[K]);
				}
				//cout << endl;

				(*posD) += kStar;
				(*nodTopk)++;
				//cout << "kStar " << kStar << ", posD " << *posD << endl;
			}
		}

		uint mCh = m;
		if (fRev_rrr[*pre]){					// is it a fictitious node ?
			ulong w = fRev_rank.rank(*pre);
			if (fictU_rrr[w]){					// is it a header of a unary path of various fictitious nodes ?
				w = fictU_rank.rank(w+1);
				ulong r,l=lisF_sel.select(w);
				if (w < nEFHead)
					r = lisF_sel.select(w+1);
				else
					r = nEF;
				// extract the extra labels from LbRevF[l..r]
				for(; l<r; l++){
					mCh++;
					patt[mCh]=LbRevF[l];
				}
			}
		}

		(*pre)++;
		createDocTopkSeq(pP, pL, pre, patt, mCh+1, fNod, dNodK, fNodK, docTopkAux, posD, nodTopk, bDL_bAux, limDTk_bAux);
	}
}

// sort the list in increasing mode by frequency
void TopkLZ::sortTopkDocsByFq(uint* fNod, uint k, uint* dNodK, uint* fNodK){
	uint i, freq, len, items;
	int j;

	for (i=0; i<nDocs; i++){
		fNodK[i] = 0;
		dNodK[i] = 0;
	}
	for (i=0; i<nDocs && fNod[i]==0; i++);
	fNodK[0] = fNod[i];	// insert first document 'i' with freq>0
	dNodK[0] = i;
	items = 1;

	for (i++; i<nDocs; i++){
		freq = fNod[i];
		len = docLength[i];
		if(freq && items < k){
			// insertion sort for this document...
			for (j=items; j>0; j--){
				if (freq > fNodK[j-1]){
					fNodK[j] = fNodK[j-1];
					dNodK[j] = dNodK[j-1];
				}else{
					if (freq == fNodK[j-1] && len > docLength[dNodK[j-1]]){
						fNodK[j] = fNodK[j-1];
						dNodK[j] = dNodK[j-1];
					}else
						break;
				}
			}
			fNodK[j] = freq;
			dNodK[j] = i;
			items++;
		}else{
			if (freq && freq >= fNodK[k-1]){
				// insertion sort for this document...
				for (j=k-2; j>=0 && freq > fNodK[j]; j--){
					fNodK[j+1] = fNodK[j];
					dNodK[j+1] = dNodK[j];
				}
				for (; j>=0 && freq == fNodK[j] && len > docLength[dNodK[j]]; j--){
					fNodK[j+1] = fNodK[j];
					dNodK[j+1] = dNodK[j];
				}
				fNodK[j+1] = freq;
				dNodK[j+1] = i;
			}
		}
	}

	if (TRACE && false){
		cout << "sorted top-k list..." << endl;
		for (i=0; i<k; i++)
			cout << "doc: " << dNodK[i] << ", freq: " << fNodK[i] << endl;
	}
}

void TopkLZ::create_sNode(LZ78Tries64::RevNode* nod, ulong *pre, ulong *aux_sNode, ulong *posNod){
	LZ78Tries64::RevNode* p = nod->fChild;
	ulong i;

	if (bTopk_b[*pre] == false && nod->fict == false){
		aux_sNode[*posNod] = nod->nodeLZ->preorder;		// set initial document in this list
		(*posNod)++;
	}
	(*pre)++;

	for(i=0; i<nod->nChildren; i++){
		create_sNode(p, pre, aux_sNode, posNod);
		p = p->nextSib;
	}
}

void TopkLZ::freqPattInFMI(uchar *patt, uint m, ulong occs, uint *fqCSA){
	uint doc;
	string query = string((char *)patt);
	auto locations = locate(fmi, query.begin(), query.begin()+m);
	//cout << "Total occurrences found with FMI : " << occs << endl;

	//cout << "locations..." << endl;
	for(ulong i=0; i<occs; i++){
		//cout << locations[i] << endl;
		doc = searchDocPosition(locations[i]);
		(fqCSA[doc])++;
	}
	//cout << endl;

	{	// free memory...
//		decltype(locations) empty;
//		locations.swap(empty);
	}

}

// create wavelet tree for range represented as an array of BitSequences...
void TopkLZ::createRange(ulong *ArrRange){
	ulong i;
	ulong *Rcopy=NULL;
	this->h = ceilingLog64(nod,2);

	if (TEST){
		Rcopy = new ulong[nod];
		for(i=0; i<nod; i++)
			Rcopy[i] = ArrRange[i];
	}

	BaseRange *netRange = new BaseRange[h];
	for (i=0; i<h; i++)
		netRange[i].B_b = bit_vector(nod, 0);

	Range = new LevelRange[h];
	sizeDS += h*sizeof(LevelRange);
	cout << " ** size of Range[]'s pointers " << h*sizeof(LevelRange) << " = " << ( h*sizeof(LevelRange))*8.0/(float)this->n << " bpc" << endl;
	genNodoWTRange(netRange, 0, ArrRange, 0, nod, 0);

	if (TEST == false)
		delete []ArrRange;

	ulong sizeBitMapsRange = 0;
	for (i=0; i<h; i++){
		Range[i].B_rrr = rrr_vector<127>(netRange[i].B_b);
		Range[i].B_rank1 = rrr_vector<127>::rank_1_type(&(Range[i].B_rrr));
		Range[i].B_rank0 = rrr_vector<127>::rank_0_type(&(Range[i].B_rrr));
		sizeBitMapsRange += size_in_bytes(Range[i].B_rrr);
		{
			bit_vector empty_b;
			netRange[i].B_b.swap(empty_b);
		}
	}

	cout << " ** size of bitVectors (RRR) of Range: " << sizeBitMapsRange << " = " << sizeBitMapsRange*8.0/(float)this->n << " bpc" << endl;
	sizeDS += sizeBitMapsRange;

	if (TEST){
		//uint doc;
		ulong x;
		cout << " Test searchLZCoord..." << endl;
		for (i=0; i<nod; i++){
			x = searchLZCoord(0, 0, nod-1, i);
			if (Rcopy[i] != x){
				cout << "ERROR!! points[" << i << "] = " << ArrRange[i] << " != searchLZCoord = " << x << endl;
				exit(0);
			}
		}
		cout << " Test Ok!" << endl;
		delete []Rcopy;
		delete []ArrRange;
	}
}

// search the 'Y' coordinate that matches with obj. => ('Y', obj) is in Range
ulong TopkLZ::searchLZCoord(uint lev, ulong lcur, ulong rcur, ulong obj){
	if (lev == h-1){	// this is a leaf !!
		if (Range[lev].B_rrr[obj])
			return rcur;
		else
			return lcur;
	}

	ulong med = lcur + ((rcur-lcur)>>1);
	if (Range[lev].B_rrr[obj] == 0){
		// we go down to the left...
		obj = Range[lev].B_rank0.rank(obj+1);
		if (lcur)
			obj -= Range[lev].B_rank0.rank(lcur);
		return searchLZCoord(lev+1, lcur, med, lcur+obj-1);
	}else{
		//we go down to the right...
		obj = Range[lev].B_rank1.rank(obj+1);
		if (lcur)
			obj -= Range[lev].B_rank1.rank(lcur);
		return searchLZCoord(lev+1, med+1, rcur, med+obj);
	}

	return 0;
}

void TopkLZ::genNodoWTRange(BaseRange *netRange, uint level, ulong *ArrRange, ulong start, ulong length, ulong lowSymbol){
	ulong i, lenR = length>>1;
	ulong medSymbol = lowSymbol + lenR;
	if (length%2)
		medSymbol++;

	for(i=start; i<(start+length); i++){
		if(ArrRange[i] < medSymbol)
			netRange[level].B_b[i] = 0;
		else
			netRange[level].B_b[i] = 1;
	}

	if (length>2){
		ulong l, r, *auxR = new ulong[lenR];
		for(i=l=start, r=0; i<(start+length); i++){
			if (netRange[level].B_b[i]){
				auxR[r] = ArrRange[i];
				r++;
			}else{
				ArrRange[l] = ArrRange[i];
				l++;
			}
		}
		for(i=0; i<lenR; i++)
			ArrRange[start+i+length-lenR] = auxR[i];

		delete [] auxR;

		genNodoWTRange(netRange, level+1, ArrRange, start, length-lenR, lowSymbol);
		if(lenR>1)
			genNodoWTRange(netRange, level+1, ArrRange, medSymbol, lenR, medSymbol);
	}
}

// determine the doc Id for the text position
uint TopkLZ::searchDocPosition(ulong pos){
	uint m, ini=0, end=nDocs-1;

	// binary search in the interval endDocs[0..nPhra-1]
	while(ini < end){
		m = ini+(end-ini)/2;
		if (pos > EndDocsPos[m])
			ini = m+1;
		else{
			if(m){
				if (pos > EndDocsPos[m-1])
					return m;
				else
					end = m-1;
			}else
				return 0;
		}
	}

	return ini;
}

// determine the doc Id for the phrase idNode
uint TopkLZ::searchDocument(ulong idNode){
	uint m, ini=0, end=nDocs-1;

	// binary search in the interval endDocs[0..nPhra-1]
	while(ini < end){
		m = ini+(end-ini)/2;
		if (idNode > EndDocs[m])
			ini = m+1;
		else{
			if(m){
				if (idNode > EndDocs[m-1])
					return m;
				else
					end = m-1;
			}else
				return 0;
		}
	}

	return ini;
}

// create the document array for LZTrie's phrases
void TopkLZ::genDocArray(LZ78Tries64::LZNode* nod, ulong *DocLZ, ulong *ini){
	LZ78Tries64::LZNode* p = nod->fChild;
	uint doc = searchDocument(nod->idNode);

	setNum64(DocLZ, *ini, lgD, doc);
	(*ini) += lgD;
	for(uint i=0; i<nod->nChildren; i++){
		genDocArray(p, DocLZ, ini);
		p = p->nextSib;
	}
}

// set extra fictitious labels...
void TopkLZ::genSeqExtraFictRevTrie(LZ78Tries64::RevNode *node, ulong *posExtF, ulong *pfRev, ulong *posFU){
	LZ78Tries64::RevNode *child;
	ulong i;

	//if(node->idNode ==19)
		//cout << "";

	if (node->fict){
		fRev_b[*pfRev] = 1;
		(*posFU)++;
		if(node->nChildren == 1){
			child = node->fChild;
			if (child->fict){
				nEFHead++;
				fictU_b[*posFU-1] = 1;
				listF_b[*posExtF] = 1;

				// search the last node in this unary path...
				uint unary = 1;
				while (child->nChildren == 1 && child->fChild->fict){
					LbRevF[*posExtF] = child->symbol;
					unary++;
					(*posExtF)++;
					child = child->fChild;
				}
				node->nChildren = child->nChildren;
				node->fChild = child->fChild;
				node->uPath = true;
				if (unary >= MAX_DEP){
					cout << "ERROR unary = " << unary << " >= MAX_DEP = " << MAX_DEP << endl;
					exit(1);
				}

				LbRevF[*posExtF] = child->symbol;
				(*posExtF)++;
			}
		}
	}
	(*pfRev)++;
	/*if (*pfRev > fRev_b.bit_size()){
		cout << "ERRRRR: *pfRev = " << *pfRev << " > " << "fRev_b.bit_size() = " << fRev_b.bit_size() << endl;
		exit(1);
	}*/

	// recursive call for all children
	if (node->nChildren){
		child = node->fChild;
		for(i=0; i<node->nChildren; i++){
			/*if(!child){
				cout << "ERRRRR: NO child and node->nChildren = " << node->nChildren << ", for node = " << node->idNode << endl;
				exit(1);
			}*/
			genSeqExtraFictRevTrie(child, posExtF, pfRev, posFU);
			child = child->nextSib;
		}
	}
}

//	create 	PRev, LbRev and LbRevF structures
void TopkLZ::genDFUDSSeqRevTrie(LZ78Tries64::RevNode *node, ulong *pos, ulong *posRev){
	LZ78Tries64::RevNode *child;
	ulong i;

	// set open parentheses
	for(i=0, child = node->fChild; i<node->nChildren; i++, (*pos)++, child=child->nextSib, (*posRev)++){
		LbRev[*posRev] = child->symbol;
		setBit64(PRev, *pos);
		//cout << "(";
	}
	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}

	// set close parenthesis
	cleanBit64(PRev, *pos);
	(*pos)++;
	//cout << ")";

	// recursive call for all children
	child = node->fChild;
	for(i=0; i<node->nChildren; i++){
		genDFUDSSeqRevTrie(child, pos, posRev);
		child = child->nextSib;
	}
	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}
}

//	create 	PLZ	:	LZTrie's DFUDS representation without Labels
void TopkLZ::genDFUDSSeqLZTrieNotLbl(LZ78Tries64::LZNode *node, ulong *pos){
	LZ78Tries64::LZNode *child;
	ulong i;

	// set open parentheses
	for(i=0, child = node->fChild; i<node->nChildren; i++, (*pos)++, child=child->nextSib)
		setBit64(PLz, *pos);

	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}

	// set close parenthesis
	cleanBit64(PLz, *pos);
	(*pos)++;

	// recursive call for all children
	child = node->fChild;
	for(i=0; i<node->nChildren; i++){
		genDFUDSSeqLZTrieNotLbl(child, pos);
		child = child->nextSib;
	}
	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}
}

void TopkLZ::readListFiles(char *inputFile, bool lowercase){
	ulong i, len, lenText;
	char fileName[300];

	std::ifstream in(inputFile);
	string line;
	std::getline(in,line);
	while(in){
	    strcpy(fileName,line.c_str());
	    //cout << "File: " << fileName << endl;
	    std::ifstream input(fileName);
		assert(input.good());
		input.seekg(0, ios_base::end);
		len = (size_t)input.tellg();
		if(len > 1){
			n += len;
			nDocs++;
		}
		input.close();
		std::getline(in,line);
	}
	in.close();
	docLength = new ulong[nDocs];
	sizeDS += nDocs*sizeof(ulong);
	cout << " ** size of docLength[1..nDocs] " << nDocs*sizeof(ulong) << " = " << (nDocs*sizeof(ulong))*8.0/(float)this->n << " bpc" << endl;

	cout << "Length of generalize text(n): " << n << ", in " << nDocs << " Documents" << endl;
	// allocate to memory for text...
	seq = new uchar[n];
	char *aux;

	lenText = 0;

	bool *present = new bool[LZ78Tries64::SIGMA];
	for (i=0; i<LZ78Tries64::SIGMA; i++)
		present[i] = false;
	uint sigma = 0;

	std::ifstream in2(inputFile);
	for(ulong texts=0; texts < nDocs;){
		std::getline(in2,line);
		//cout << "File: " << line.c_str() << endl;
		strcpy(fileName,line.c_str());
		std::ifstream input(fileName); 			// open file
		//cout << "... File: " << fileName << endl;
		assert(input.good()); 				// check open
		input.seekg(0, ios_base::end);			// move to the end
		len = (size_t)input.tellg();			// add the final symbol (smallest)
		if(len > 1){
			docLength[texts] = len;
			aux = new char[len];
			//if (TRACE)
			//cout << "To read " << fileName << " pos: " << lenText << "..." << lenText+len-1 << endl;
			input.seekg(0, ios_base::beg);		// move back to the beginning
			if (input.good()){
				input.read(aux, len);

				len--;
				aux[len] = cutDoc;
				//cout << aux << endl;
				for (i=0; i<len; i++, lenText++){

					if (((int)aux[i]) < 0)
						aux[i] = ' ';
					if (present[(uint)aux[i]] == false){
						sigma++;
						present[(uint)aux[i]] = true;
					}
					if((uchar)(aux[i]) <= cutDoc)
						seq[lenText] = ' ';
					else{
						if (lowercase)
							seq[lenText] = (uchar)(tolower(aux[i]));
						else
							seq[lenText] = (uchar)(aux[i]);
					}
					(charTf[seq[lenText]])++;
				}
				seq[lenText] = cutDoc;
				lenText++;
				assert(!input.fail());
				input.close();
			}else{
				cout << "Can't to open the file: <" << fileName << ">";
				exit(1);
			}
			delete [] aux;
			texts++;
		}
	}
	//cout << " ## sigma original = " << sigma << endl;

	seq[n-1] = '\0';
	in2.close();
	char fileCpy[300];
	strcpy(fileCpy,inputFile);
	strcat(fileCpy, "_copy.txt");
	cout << "Saving original sequence in " << fileCpy << endl;
	ofstream myfile;
	myfile.open (fileCpy);
	myfile << seq;
	myfile.close();
	seq[n-1] = cutDoc;

	if(TEST){
		uint DD = 0;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				DD++;
		}
		if(nDocs != DD){
			cout << "ERROR with nDocs in Sequence !! " << endl;
			cout << "nDocs = " << nDocs << " != " << DD << endl;
			exit(1);
		}
	}
	if(TRACE){		// listing original sequence
		cout << endl << "T[0.." << n-1 << "]:" << endl;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				cout << "$";
			else
				cout << seq[i];
		}
		cout << endl;
		cout << "pdocLength[0.." << nDocs-1 << "]: ";
		for(i=0; i<nDocs; i++)
			cout << docLength[i] << " ";
		cout << endl;
		cout << "charTf[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++)
			if (charTf[i])
				cout << i << "(" << charTf[i] << ") ";
		cout << endl;
	}
}

void TopkLZ::readUniqueFile(char *inputFile, bool lowercase){
	ulong i, j, len;
	n = nDocs = 0;
	ifstream input(inputFile);			// open file
	assert(input.good()); 				// check open
	input.seekg(0, ios_base::end);		// move to the end
	n = (size_t)input.tellg();
	seq = new uchar[n];
	char *aux = new char[n];

	input.seekg(0, ios_base::beg);		// move back to the beginning
	if (input.good()){
		input.read(aux, n-1);
		for (i=0; i<n-1; i++){
			if((uchar)(aux[i]) == boundSymbol){
				nDocs++;
				while((uchar)(aux[i+1]) == boundSymbol){
					seq[i] = '\n';
					(charTf['\n'])++;
					i++;
				}
				seq[i] = cutDoc;
			}else{
				if((uchar)(aux[i]) <= cutDoc)
					seq[i] = '\n';
				else{
					if (lowercase)
						seq[i] = (uchar)(tolower(aux[i]));
					else
						seq[i] = (uchar)(aux[i]);
				}
			}
			(charTf[seq[i]])++;
		}
		if(seq[n-2] == cutDoc){
			seq[n-2] = '\n';
			(charTf['\n'])++;
			nDocs--;
		}
		seq[n-1] = cutDoc;
		nDocs++;
		//cout << seq << endl;
		assert(!input.fail());
		input.close();
	}else{
		cout << "Can't to open the file: <" << inputFile << ">";
		exit(1);
	}
	input.close();
	delete [] aux;

	docLength = new ulong[nDocs];
	sizeDS += nDocs*sizeof(ulong);
	cout << " ** size of docLength[1..nDocs] " << nDocs*sizeof(ulong) << " = " << (nDocs*sizeof(ulong))*8.0/(float)this->n << " bpc" << endl;
	cout << "Length of generalize text(n): " << n << ", in " << nDocs << " Documents" << endl;

	for (i=j=len=0; i<n; i++){
		if(seq[i] == cutDoc){
			docLength[j] = len;
			len=0;
			j++;
		}else
			len++;
	}
	if (j != nDocs){
		cout << "Error cutDocs Symbols = " << j << " != nDocs = " << nDocs << endl;
		exit(1);
	}

	seq[n-1] = '\0';
	char fileCpy[300];
	strcpy(fileCpy,inputFile);
	strcat(fileCpy, "_copy.txt");
	ofstream myfile;
	myfile.open (fileCpy);
	myfile << seq;
	myfile.close();
	seq[n-1] = cutDoc;

	if(TEST){
		uint DD = 0;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				DD++;
		}
		if(nDocs != DD){
			cout << "ERROR with nDocs in Sequence !! " << endl;
			cout << "nDocs = " << nDocs << " != " << DD << endl;
			exit(1);
		}
	}

	if(TRACE){		// listing original sequence
		cout << "T[0.." << n-1 << "]:" << endl;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				cout << "$";
			else
				cout << seq[i];
		}
		cout << endl;
		cout << "pdocLength[0.." << nDocs-1 << "]: ";
		for(i=0; i<nDocs; i++)
			cout << docLength[i] << " ";
		cout << endl;
		cout << "charTf[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++)
			if (charTf[i])
				cout << i << "(" << charTf[i] << ") ";
		cout << endl;
	}
}

// make priorities for hash table by probability
void TopkLZ::makePriorities(ulong* cPosTb){
	uint i, j, sigma = 0;
	uchar *prior = new uchar[LZ78Tries64::SIGMA];
	prior[0] = 0; // --> charTf[0];

	bool tr = TRACE;
	TRACE = false;

	for (i=1; i<LZ78Tries64::SIGMA; i++){
		if (charTf[i] > 0)
			sigma++;
		for (j=i; j>0 && charTf[prior[j-1]] < charTf[i]; j--){
			prior[j] = prior[j-1];
		}
		prior[j] = i; // --> charTf[i];
	}

	cout << " ## sigma = " << sigma << endl;

	if(TRACE){
		cout << "prior[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++){
			cout << i << "(" << (uint)prior[i] << ") ";
		}
		cout << endl;
	}

	uint rows = LZ78Tries64::SIGMA / LZ78Tries64::LENHASH;
	if (LZ78Tries64::SIGMA % LZ78Tries64::LENHASH)
		rows++;
	uint code = 0;
	for(uint fil=0; fil<rows; fil++){
		for(i=0; i<LZ78Tries64::LENHASH && code<LZ78Tries64::SIGMA; i++){
			if(fil%2 == 0)
				cPosTb[prior[code]] = i;
			else
				cPosTb[prior[code]] = LZ78Tries64::LENHASH-1-i;
			code++;
		}
	}

	if(TRACE){
		cout << "charPosTb[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++){
			cout << i << "(" << cPosTb[i] << ") ";
		}
		cout << endl;
	}
	delete [] prior;
	TRACE = tr;
}

// return true if the pattern is in the RevTrie and in '*pv' the respective preorder value in the tree,
// in '*x' stores its DFUDS bit position. Return false if the pattern does not exist
bool TopkLZ::searchPattern_Rev(uchar* pat, uint m, ulong *x, ulong *pv){
	ulong d, w1, w2, sig0, pre=0, pos=*x;

	while (m > 0){
		m--;
		if (fRev_rrr[pre]){						// is it a fictitious node ?
			w1 = fRev_rank.rank(pre+1)-1;

			if (fictU_rrr[w1]){					// is it a fictitious node with unary path of other fictitious nodes ?
				w2 = fictU_rank.rank(w1+1);
				ulong r,l=lisF_sel.select(w2);
				if (w2 < nEFHead)
					r = lisF_sel.select(w2+1);
				else
					r = nEF;
				if(m){
					while (l<r && LbRevF[l]==pat[m]){
						m--; l++;
						if (m==0) break;
					}
					if (l<r){
						if (m==0){
							if (LbRevF[l]!=pat[0])return false;
						}else return false;
					}else{
						pre++;
						if (fRev_rrr[pre]){
							if (fictU_rrr[w1+1]){
								sig0 = treeRev->selectNext0(pos);
								if (sig0 > n) sig0 = n;
								else if (pos>=sig0) return false;

								if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
									if(w2 > 1){
										d=1;
										treeRev->fwd_search(sig0-w2+1, &d, &pos);
										pos++;
									}else
										pos = sig0+1;
								}else return false;

								pre = pos-treeRev->rank_1(pos-1);
							}else m++;
						}else m++;
					}
				}else{
					if (LbRevF[l]==pat[0]){
						*x = pos;
						*pv = pre;
						return true;
					}
					else return false;
				}
			}else{									// normal case
				sig0 = treeRev->selectNext0(pos);
				if (sig0 > nRev) sig0 = nRev;
				else if (pos>=sig0) return false;

				if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
					if(w2 > 1){
						d=1;
						treeRev->fwd_search(sig0-w2+1, &d, &pos);
						pos++;
					}else
						pos = sig0+1;
					pre = pos-treeRev->rank_1(pos-1);
				}else
					return false;
			}
		}else{										// normal case
			sig0 = treeRev->selectNext0(pos);
			if (sig0 > nRev) sig0 = nRev;
			else if (pos>=sig0) return false;

			// w2 is the relative position for the charter pat[m]
			if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
				if(w2 > 1){
					d=1;
					treeRev->fwd_search(sig0-w2+1, &d, &pos);
					pos++;
				}else
					pos = sig0+1;
				pre = pos-treeRev->rank_1(pos-1);
			}else
				return false;
		}
	}
	*x = pos;
	*pv = pre;
	return true;
}

bool TopkLZ::searchPattern_Rev2(uchar* pat, uint m, ulong *x, ulong *pv, bool realP){
	ulong d, w1, w2, sig0, pre=0, pos=*x;
	uint len;

	while (m){
		m--;
		sig0 = treeRev->selectNext0(pos);
		if (sig0 > nRev) sig0 = nRev;
		else if (pos>=sig0) return false;

		len = sig0-pos;
		// w2 is the relative position for the symbol 'pat[m]'
		if (treeRev->thereIsChild(pos, pat[m], &w2, len)){
			if(w2 > 1){
				d=1;
				treeRev->fwd_search(sig0-w2+1, &d, &pos);
				pos++;
			}else
				pos = sig0+1;
			pre = pos-treeRev->rank_1(pos-1);

			if(m){
				if (fRev_rrr[pre]){						// is it a fictitious node ?
					w1 = fRev_rank.rank(pre+1)-1;
					if (fictU_rrr[w1]){					// is it a fictitious node with unary path of other fictitious nodes ?
						m--;
						// [1].- initial symbol of the unary path
						sig0 = treeRev->selectNext0(pos);
						if (sig0 > n) sig0 = n;
						else if (pos>=sig0) return false;

						if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
							if(w2 > 1){
								d=1;
								treeRev->fwd_search(sig0-w2+1, &d, &pos);
								pos++;
							}else
								pos = sig0+1;
						}else return false;

						pre++;
						if(m){
							m--;
							// [2].- symbols of the unary path
							w2 = fictU_rank.rank(w1+1);
							ulong r,l=lisF_sel.select(w2);
							if (w2 < nEFHead)
								r = lisF_sel.select(w2+1);
							else
								r = nEF;

							if (realP){
								if ((m+1) < (r-l))
									return false;
								while (l<r && LbRevF[l]==pat[m]){
									m--; l++;
								}
								if (l<r)
									return false;
								else
									m++;	// ??????????????????????????
							}else{
								while (m && l<r && LbRevF[l]==pat[m]){
									m--; l++;
								}
								if (m==0){
									if (l<r){
										if (LbRevF[l]!=pat[0])
											return false;
									}else
										m++;
								}else{
									if (l<r)
										return false;
									m++;
								}
							}
						}else
							if (realP)
								return false; // no hay frase que empieze con lo caracteres extra ficticios !!
					}
				}
			}
		}else
			return false;
	}
	*x = pos;
	*pv = pre;
	return true;
}

bool TopkLZ::searchPattern_Rev3(uchar* pat, uint m, ulong *x, ulong *pv, bool realP){
	ulong d, w1, w2, sig0, pre=0, pos=*x;
	uint len;

	while (m){
		m--;
		sig0 = treeRev->selectNext0(pos);
		if (sig0 > nRev) sig0 = nRev;
		else if (pos>=sig0) return false;

		len = sig0-pos;
		// w2 is the relative position for the symbol 'pat[m]'
		if (treeRev->thereIsChild(pos, pat[m], &w2, len)){
			if(w2 > 1){
				d=1;
				treeRev->fwd_search(sig0-w2+1, &d, &pos);
				pos++;
			}else
				pos = sig0+1;
			pre = pos-treeRev->rank_1(pos-1);

			if(m){
				if (fRev_rrr[pre]){						// is it a fictitious node ?
					w1 = fRev_rank.rank(pre+1)-1;
					if (fictU_rrr[w1]){					// is it a fictitious node with unary path of other fictitious nodes ?
						m--;
						// read symbols of the unary path
						w2 = fictU_rank.rank(w1+1);
						ulong r,l=lisF_sel.select(w2);
						if (w2 < nEFHead)
							r = lisF_sel.select(w2+1);
						else
							r = nEF;

						if (realP){
							if ((m+1) < (r-l))
								return false;
							while (m && l<r && LbRevF[l]==pat[m]){
								m--; l++;
							}
							if (l<r)
								return false;
							m++;
							pre++;
						}else{
							while (m && l<r && LbRevF[l]==pat[m]){
								m--; l++;
							}
							if (m==0){
								if (l<r){
									if (LbRevF[l]!=pat[0])
										return false;
									//else pos=treeRev->selectNext0(pos)+1;
								}else{
									m++;
									pre++;
								}
							}else{
								if (l<r)
									return false;
								m++;
								pre++;
							}
						}
					}
				}
			}
		}else
			return false;
	}
	*x = pos;
	*pv = pre;
	return true;
}

// It stores in docList[1..k] and frqList[1..k] the list and frequencies of the top-k most frequent documents for the pattern p[1..m]
ulong TopkLZ::topKDocument(uchar* pat, uint m, uint* docList, uint k){
	uint key = 1;
	ulong l, pv, x=1;

	if (searchPattern_Rev3(pat, m, &x, &pv, 0)){
		ulong pos, r;
		uint st=0;

		//cout << "pv: " << pv << endl;
		if (bTopk_rrr[pv]){
			// Answer stored in the node pv, from lowest to highest frequency
			ulong list = bTopk_rank.rank(pv+1);
			//cout << "list num " << list << endl;
			l = limDTk_sel.select(list);		// number of this block of id documents
			if(list<nTopk)
				r = limDTk_sel.select(list+1);
			else
				r = limDTk_rrr.size();		// 	CAMBIAR EN APP-TOPK limDTk_b POR limDTk_rrr

			st=r-l;						// number of different id's
			if (st > k)
				st = k;
			if (st == k || bDL_b[list-1]){
				//cout << "Stored with " << st << " different id's" << endl;
				//cout << "bDL->getBit(ra-1) " << bDL_b[list-1] << endl;
				//cout << "Documents, from l*lgD, l = " << l << ", lgD = " << lgD << endl;
				l *= lgD;
				for(pos=0; pos<st; l+=lgD, pos++){
					docList[pos] = getNum64(docTopk, l, lgD);
					//cout << docList[pos] << " ";
				}
				//cout << endl;
				return st;
			}
		}

		// Computes frequencies by brute force
		key=1;
		ulong i, len, fictBef, doc;

		// [1].- Occ type one...
		len = treeRev->subTreeSize(x);
		fictBef = fRev_rank.rank(pv);
		//ulong zzzz = fRev_rank.rank(pv+len);
		len -= fRev_rank.rank(pv+len) - fictBef;
		pv -= fictBef;
		for(i=pv*lgNod; len > 0; len--, i+=lgNod){
			l = getNum64(Node, i, lgNod);
			x = treeLz->select_0(l)+1;
			r = l + treeLz->subTreeSize(x);
			// scan DocLZ[l..r]. A new document doc has a key = dictD[doc], and his frequencies is in fqKey[key]
			//cout << "Docs for l(" << l << ") :  ";

			for (pos=l*lgD; l<r; l++, pos+=lgD){
				doc = getNum64(DocLZ, pos, lgD);
				//cout << doc << " ";
				x = dictD[doc];	// key
				if(x)
					(keyFq[x])++;
				else{
					// x==0 then it is a new ID
					dictD[doc] = key;
					keyDoc[key] = doc;
					keyFq[key] = 1;
					//cout << doc << "[" << key << "](" << keyFq[key] << ")" << endl;
					key++;
				}
			}
		}
		//cout << endl;
	}

	// [2].- Occ type two and three..
	if (m>1)
		extractFqType23(pat, m, &key);

	// [3].- Make final answer...
	//cout << endl;
	/*cout << "key " << key << endl;
	cout << "Freq [1...]" << endl;
	for (l=1; l<key; l++)
		cout << l << ", doc "<<keyDoc[l]<< ", Fq:" << keyFq[l] << endl;*/
	//exit(0);

	for (l=1; l<key; l++){
		heap[l] = l;
		swimInMaxQ(l);
		//for (uint a=1; a<=l; a++)
		//	cout << " " << keyFq[heap[a]];
		//cout << endl;
	}
	//for (l=1; l<key; l++)
	//	cout << " f" << (uint)keyFq[heap[l]];
	//cout << endl;

	key--;
	if (key < k)
		k = key;

	for (l=0; l<k; l++){
		x = heap[1];
		setTopMaxQ(key-l);
		//cout << " F" << (uint)keyFq[x] << endl;
		docList[l] = keyDoc[x];
	}

	// CLEAN ... to unmark documents
	for (l=1; l<=key; l++){
		dictD[keyDoc[l]] = 0;
		keyDoc[l] = keyFq[l] = 0;
	}

	return k;
}

// apply sequential search in the range [l,r]...
void TopkLZ::reportDocOcc2_levRange(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, ulong lobj, ulong robj, uint* key){
	uint auxLev, doc;
	ulong i, obj, auxL, auxR, med;

	for(i=l; i<=r; i++){
		for(obj=i, auxLev=lev, auxL=lcur, auxR=rcur; auxLev<h-1; auxLev++){
			med = auxL + ((auxR-auxL)>>1);
			if (Range[auxLev].B_rrr[obj] == 0){
				// we go down to the left...
				obj = auxL+Range[auxLev].B_rank0.rank(obj+1);
				obj -= Range[auxLev].B_rank0.rank(auxL)+1;
				auxR = med;
			}else{
				//we go down to the right...
				obj = med+Range[auxLev].B_rank1.rank(obj+1);
				obj -= Range[auxLev].B_rank1.rank(auxL);
				auxL = med+1;
			}
		}

		if (auxLev == h-1){
			// this is a leaf --> report...
			if (Range[auxLev].B_rrr[obj])
				obj = auxR;
			else
				obj = auxL;
			doc = getNum64(DocLZ, obj*lgD, lgD);
			//cout << "Occ type 2 --> doc " << doc << endl;

			obj = dictD[doc];	// key
			if(obj)
				(keyFq[obj])++;
			else{
				obj = *key;
				// x==0 then it is a new ID
				dictD[doc] = obj;
				keyDoc[obj] = doc;
				keyFq[obj] = 1;
				//cout << doc << "[" << obj << "](" << keyFq[obj] << ") ";
				(*key)++;
			}
		}
	}
}

// [lcur, rcur] is the current bitstring in the level 'lev' (B[lcur, rcur]). [lobj, robj] is the target interval to search in LZTrie
void TopkLZ::searchIntervalInRange(uint lev, ulong lcur, ulong rcur, ulong lobj, ulong robj, ulong lrmq, ulong rrmq, uint *key){
	if (lev == h-1){	// this is a leaf (or two leaves) !!
		if (lrmq < rrmq){
			if (lobj <= lcur){
				uint doc = getNum64(DocLZ, lcur*lgD, lgD);
				//cout << " hoja -> doc " << doc << endl;
				ulong x = dictD[doc];	// key
				uint k = *key;
				if(x)
					(keyFq[x])++;
				else{
					// x==0 then it is a new ID
					dictD[doc] = k;
					keyDoc[k] = doc;
					keyFq[k] = 1;
					//cout << doc << "[" << key << "](" << keyFq[key] << ") ";
					(*key)++;
				}

				if ((rcur <= robj) && (lcur!=rcur)){
					uint doc = getNum64(DocLZ, rcur*lgD, lgD);
					//cout << " hoja -> doc " << doc << endl;
					ulong x = dictD[doc];	// key
					uint k = *key;
					if(x)
						(keyFq[x])++;
					else{
						// x==0 then it is a new ID
						dictD[doc] = k;
						keyDoc[k] = doc;
						keyFq[k] = 1;
						//cout << doc << "[" << key << "](" << keyFq[key] << ") ";
						(*key)++;
					}
				}

			}else{
				if (rcur <= robj){
					uint doc = getNum64(DocLZ, rcur*lgD, lgD);
					//cout << " hoja -> doc " << doc << endl;
					ulong x = dictD[doc];	// key
					uint k = *key;
					if(x)
						(keyFq[x])++;
					else{
						// x==0 then it is a new ID
						dictD[doc] = k;
						keyDoc[k] = doc;
						keyFq[k] = 1;
						//cout << doc << "[" << key << "](" << keyFq[key] << ") ";
						(*key)++;
					}
				}
			}
		}else{ // then lcur == rcur
			if (Range[lev].B_rrr[lrmq] == 0 && lobj <= lcur){
				uint doc = getNum64(DocLZ, lcur*lgD, lgD);
				//cout << " hoja -> doc " << doc << endl;
				ulong x = dictD[doc];	// key
				uint k = *key;
				if(x)
					(keyFq[x])++;
				else{
					// x==0 then it is a new ID
					dictD[doc] = k;
					keyDoc[k] = doc;
					keyFq[k] = 1;
					//cout << doc << "[" << key << "](" << keyFq[key] << ") ";
					(*key)++;
				}
			}else{
				if (Range[lev].B_rrr[lrmq] && rcur <= robj){
					uint doc = getNum64(DocLZ, rcur*lgD, lgD);
					//cout << " hoja -> doc " << doc << endl;
					ulong x = dictD[doc];	// key
					uint k = *key;
					if(x)
						(keyFq[x])++;
					else{
						// x==0 then it is a new ID
						dictD[doc] = k;
						keyDoc[k] = doc;
						keyFq[k] = 1;
						//cout << doc << "[" << key << "](" << keyFq[key] << ") ";
						(*key)++;
					}
				}
			}
		}
	}else{
		if(lobj <= lcur && rcur <= robj)			// we consider the complete node --> the node is totally inside of objective
			reportDocOcc2_levRange(lev, lcur, rcur, lrmq, rrmq, lobj, robj, key);
		else{
			ulong med = lcur + ((rcur-lcur)>>1);

			if(robj <= med){
				// we go down only to the left...
				ulong befL0=0, befR0=0, aux;
				aux = Range[lev].B_rank0.rank(lcur);
				befL0 = Range[lev].B_rank0.rank(lrmq) - aux;
				befR0 = Range[lev].B_rank0.rank(rrmq+1) - aux;
				if (lev<h && befR0>befL0)
					searchIntervalInRange(lev+1, lcur, med, lobj, robj, lcur+befL0, lcur+befR0-1, key);
			}else{
				if(lobj > med){
					// we go down only to the right...
					ulong befL1=0, befR1=0, aux;
					aux = Range[lev].B_rank1.rank(lcur);
					befL1 = Range[lev].B_rank1.rank(lrmq) - aux;
					befR1 = Range[lev].B_rank1.rank(rrmq+1) - aux;
					if (befR1>befL1)
						searchIntervalInRange(lev+1, med+1, rcur, lobj, robj, med+befL1+1, med+befR1, key);
				}else{
					// first, we go down to the left...
					ulong befL0=0, befR0=0, aux;
					aux = Range[lev].B_rank0.rank(lcur);
					befL0 = Range[lev].B_rank0.rank(lrmq) - aux;
					befR0 = Range[lev].B_rank0.rank(rrmq+1) - aux;
					if (lev<h && befR0>befL0)
						searchIntervalInRange(lev+1, lcur, med, lobj, robj, lcur+befL0, lcur+befR0-1, key);

					// second, we go down to the right...
					ulong befL1=0, befR1=0;
					aux = Range[lev].B_rank1.rank(lcur);
					befL1 = Range[lev].B_rank1.rank(lrmq) - aux;
					befR1 = Range[lev].B_rank1.rank(rrmq+1) - aux;
					if (lev<h && befR1>befL1)
						searchIntervalInRange(lev+1, med+1, rcur, lobj, robj, med+befL1+1, med+befR1, key);
				}
			}
		}
	}
}

// check if obj (LZ) matches with some point in [lRev... rRev]
bool TopkLZ::searchRevRange(uint lev, ulong lcur, ulong rcur, ulong obj, ulong lRev, ulong rRev){
	if(lRev <= lcur && rRev >= rcur)
		return true;

	if (lev == h-1){	// this is a leaf !!
		if (Range[lev].B_rrr[lRev] == 0 && obj == lcur)
			return true;
		else{
			if (Range[lev].B_rrr[lRev] && obj == rcur)
				return true;
		}
	}else{
		ulong med = lcur + ((rcur-lcur)>>1);
		if(obj <= med){	// we go down to the left...
			ulong befL0=0, befR0=0, aux;
			aux = Range[lev].B_rank0.rank(lcur);
			befL0 = Range[lev].B_rank0.rank(lRev) - aux;
			befR0 = Range[lev].B_rank0.rank(rRev+1) - aux;
			if (lev<h && befR0>befL0)
				return (searchRevRange(lev+1, lcur, med, obj, lcur+befL0, lcur+befR0-1));
		}else{
			// we go down to the right...
			ulong befL1=0, befR1=0, aux;
			aux = Range[lev].B_rank1.rank(lcur);
			befL1 = Range[lev].B_rank1.rank(lRev) - aux;
			befR1 = Range[lev].B_rank1.rank(rRev+1) - aux;
			if (befR1>befL1)
				return (searchRevRange(lev+1, med+1, rcur, obj, med+befL1+1, med+befR1));
		}
	}
	return false;
}

// this stores in mat[i] the preorder of the node that matches with a phrase p_{i..m} (find all phrase i for this m)
void TopkLZ::searchPhase_Rev(uchar* pat, uint m, ulong *mat){
	ulong d, w1, w2, sig0, pre=0, pos=1;

	while (m > 0){
		if (fRev_rrr[pre]){						// is it a fictitious node ?
			w1 = fRev_rank.rank(pre+1)-1;
			if (fictU_rrr[w1]){					// is it a fictitious node with unary path of other fictitious nodes ?
				w2 = fictU_rank.rank(w1+1);
				uint r,l=lisF_sel.select(w2);
				if (w2 < nEFHead)
					r = lisF_sel.select(w2+1);
				else
					r = nEF;
				if ((r-l) > m)
					return;
				for(;l<r && LbRevF[l]==pat[m];m--,l++);
				if (l<r || m==0)
					return;
				else{
					pre++;
					if (fRev_rrr[pre]){
						if (fictU_rrr[w1+1]){
							sig0 = treeRev->selectNext0(pos);
							if (sig0 > nRev) sig0 = nRev;
							else if (pos>=sig0) return;
							if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
								if(w2 > 1){
									d=1;
									treeRev->fwd_search(sig0-w2+1, &d, &pos);
									pos++;
								}else
									pos = sig0+1;
							}else return;
							pre = pos-treeRev->rank_1(pos-1);
							if (fRev_rrr[pre]==0)
								mat[m] = pre;
						}else m++;
					}else m++;
				}
			}else{									// one fictitious node
				sig0 = treeRev->selectNext0(pos);
				if (sig0 > nRev) sig0 = nRev;
				else if (pos>=sig0) return;
				if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
					if(w2 > 1){
						d=1;
						treeRev->fwd_search(sig0-w2+1, &d, &pos);
						pos++;
					}else
						pos = sig0+1;
				}else return;
				pre = pos-treeRev->rank_1(pos-1);
				if (fRev_rrr[pre]==0)
					mat[m] = pre;
			}
		}else{										// phrase node
			sig0 = treeRev->selectNext0(pos);
			if (sig0 > nRev) sig0 = nRev;
			else if (pos>=sig0) return;
			if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
				if(w2 > 1){
					d=1;
					treeRev->fwd_search(sig0-w2+1, &d, &pos);
					pos++;
				}else
					pos = sig0+1;
			}else return;
			pre = pos-treeRev->rank_1(pos-1);
			if (fRev_rrr[pre]==0)
				mat[m] = pre;
		}
		m--;
	}
}
void TopkLZ::searchPhase_Rev2(uchar* pat, uint m, ulong *mat){
	ulong d, w1, w2, sig0, pre=0, pos=1;

	while (m > 0){
		sig0 = treeRev->selectNext0(pos);
		if (sig0 > nRev) sig0 = nRev;
		else if (pos>=sig0) return;
		if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
			if(w2 > 1){
				d=1;
				treeRev->fwd_search(sig0-w2+1, &d, &pos);
				pos++;
			}else
				pos = sig0+1;
		}else return;
		pre = pos-treeRev->rank_1(pos-1);
		if (fRev_rrr[pre]==0)
			mat[m] = pre;

		if(m){
			if (fRev_rrr[pre]){						// is it a fictitious node ?
				w1 = fRev_rank.rank(pre+1)-1;
				if (fictU_rrr[w1]){					// is it a fictitious node with unary path of other fictitious nodes ?
					m--;
					// [1].- initial symbol of the unary path
					sig0 = treeRev->selectNext0(pos);
					if (sig0 > n) sig0 = n;
					else if (pos>=sig0) return;

					if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
						if(w2 > 1){
							d=1;
							treeRev->fwd_search(sig0-w2+1, &d, &pos);
							pos++;
						}else
							pos = sig0+1;
					}else return;

					pre++;
					//if (fRev_rrr[pre]==0)
						//mat[m] = pre;
					if(m){
						m--;
						// [2].- symbols of the unary path
						w2 = fictU_rank.rank(w1+1);
						ulong r,l=lisF_sel.select(w2);
						if (w2 < nEFHead)
							r = lisF_sel.select(w2+1);
						else
							r = nEF;

						while (m && l<r && LbRevF[l]==pat[m]){
							m--;
							l++;
						}
						if (l==r){
							m++;
							mat[m] = pre;
						}

						if(m==0)
							break;

						/*if (l<r){
							if (m==0){
								m++;
								mat[m] = pre;
								//if (LbRevF[l]!=pat[0])return;
							}else return;
						}else{
							m++;
							mat[m] = pre;
						}*/
					}else
						break;
				}
			}
		}
		m--;
	}
}

void TopkLZ::searchPhase_Rev3(uchar* pat, uint m, ulong *mat){
	ulong d, w1, w2, sig0, pre=0, pos=1;

	while (m){
		sig0 = treeRev->selectNext0(pos);
		if (sig0 > nRev) sig0 = nRev;
		else if (pos>=sig0) return;
		if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
			if(w2 > 1){
				d=1;
				treeRev->fwd_search(sig0-w2+1, &d, &pos);
				pos++;
			}else
				pos = sig0+1;
		}else return;
		pre = pos-treeRev->rank_1(pos-1);
		if (fRev_rrr[pre]==0)
			mat[m] = pre;

		if(m){
			if (fRev_rrr[pre]){						// is it a fictitious node ?
				w1 = fRev_rank.rank(pre+1)-1;
				if (fictU_rrr[w1]){					// is it a fictitious node with unary path of other fictitious nodes ?
					m--;
					// symbols of the unary path
					w2 = fictU_rank.rank(w1+1);
					ulong r,l=lisF_sel.select(w2);
					if (w2 < nEFHead)
						r = lisF_sel.select(w2+1);
					else
						r = nEF;

					while (m && l<r && LbRevF[l]==pat[m]){
						m--;
						l++;
					}
					if (l==r){
						if(m==0)
							break;
						m++;
						//mat[m] = pre;
					}else{
						if(m)
							break;
					}

					if(m==0)
						break;
				}
			}
		}
		m--;
	}
}
void TopkLZ::extractFqType23(uchar* pat, uint m, uint *key){
	ulong pv, len, fictBef, x=1, l, r;
	uint i, j, k, doc;

	// =========================================================================================================================
	// occurrences type 2...
	// =========================================================================================================================
	// for each partition, pref[k] stores the prefix P_{0..k}. [0]=preorder-Rev, [1]=subtree size of [0]
	// suff[k] stores the suffix P_{k+1..m-1}.: [0]:preorder-LZ, [1]:subtree size of [0].
	ulong pref[m-1][2], suff[m-1][2], M[m-1][m-1][3];
	ulong *mat = new ulong[m-1];
	for (i=0; i<(m-1); i++){
		mat[i] = pref[i][0] = pref[i][1] = suff[i][0] = suff[i][1] = 0;
		for (j=0; j<(m-1); j++)
			M[i][j][0] = M[i][j][1] = 0;
	}

	// [2] occ type 2...
	for (l=1; l<m; l++){
		x=1;
		if (searchPattern_Rev3(pat, l, &x, &pv, 0)){				// search prefix P_{0..l-1}...
			len = treeRev->subTreeSize(x);
			fictBef = fRev_rank.rank(pv);
			pref[l][0] = pv - fictBef;
			pref[l][1] = len - fRev_rank.rank(pv+len) + fictBef;

			x = 1;
			if (searchPattern_Rev3(pat+l, m-l, &x, &pv, 1)){		// search suffix P_{l..m-1}...
				if (fRev_rrr[pv]){
					suff[l-1][0] = 0;
				}else{
					pv -= fRev_rank.rank(pv+1);
					pv = getNum64(Node, pv*lgNod, lgNod);					//pos in Plz of the preorder pv
					x = treeLz->select_0(pv)+1;
					len = treeLz->subTreeSize(x);
					suff[l-1][0] = pv;
					suff[l-1][1] = len;	// subtree_size(x) in LZTrie
					searchIntervalInRange(0, 0, nod-1, pv, pv+len-1, pref[l][0], pref[l][0]+pref[l][1]-1, key);
				}
			}
		}else{
			if (m>2){	// for occurrences type 3...
				x = 1;
				if (searchPattern_Rev3(pat+l, m-l, &x, &pv, 1)){		// search suffix P_{l..m-1}...
					if (fRev_rrr[pv]){
						suff[l-1][0] = 0;
					}else{
						pv -= fRev_rank.rank(pv+1);
						pv = getNum64(Node, pv*lgNod, lgNod);					//pos in Plz of the preorder pv
						x = treeLz->select_0(pv)+1;
						len = treeLz->subTreeSize(x);
						suff[l-1][0] = pv;
						suff[l-1][1] = len;
					}
				}
			}
		}
	}
	/*cout << "type2, Prefixes and Suffixes..." << endl;
	for (i=0; i<(m-1); i++)
		cout << i << ":(" << pref[i][0] << "," << pref[i][1] << "," << suff[i][0] << "," << suff[i][1]<< ") " << endl;
	 */

// =========================================================================================================================
// [3] occ type 3...
// =========================================================================================================================
	// M[i][j] represents to the phrase P_{i..j}, where:
	// M[i][j][0] stores the preorderLZ for the node that matches with p_{i..j}
	// M[i][j][1] stores the M's index r of phrase(K+1) which is stored in M[j+1][r]=P_{j+1..r} that join to the right with phrase(K)=P_{i..j}
	// M[i][j][2] stores the preorderLZ of phrase(K+1) that join to the right with phrase(K) = P_{i..j}
	if (m>2){
		for (j=m-2; j>0; j--){
			searchPhase_Rev3(pat, j, mat);
			for (i=1; i<=j; i++){
				if (mat[i]){
					pv = mat[i] - fRev_rank.rank(mat[i]);//FRev->rank1(mat[i]-1);
					M[i][j][0] = getNum64(Node, pv*lgNod, lgNod);
					//cout << "M["<<i<<"]["<<j<<"][0] = " << M[i][j][0] << endl;
					// to search preorderLZ of phrase(K+1) that join with phrase(K) = P_{i..j}
					M[i][j][2] = searchLZCoord(0, 0, nod-1, pv);
					//cout << "M["<<i<<"]["<<j<<"[2] = " << M[i][j][2] << endl;
					for (l=r=j+1; r<m-1; r++){ 	// search concatenations...
						if(M[l][r][0] == M[i][j][2]){
							M[i][j][1]=r;
							break;
						}
					}
				}
			}
		}

		/*cout << "M (preLZ, next M's index (phrase concatenated), preLZ of next phrase)...";
		for (i=1; i<m-1; i++){
			cout << endl << "[" << i << "] ";
			for (j=i; j<m-1; j++)
				cout << j << ":(" << M[i][j][0] << "," << M[i][j][1] << "," << M[i][j][2]<< ") ";
		}
		cout << endl;*/


		// *** HERE THE MODEL FOR P_{0..m-1} IS SPLIT P IN P_{0..i-1} + P_{i..j} + P_{j+1..m-1},
		// where P_{0..i-1} is in pref[i], P_{i..j} is in M[i][j] and P_{j+1..m-1} is in suff[j]
		// and we also consider to split P_{i..j} in more phrases stored in M[][]
		for (i=1; i<m-1; i++){
			if(pref[i][0]){
				for (j=i; j<m-1; j++){
					pv = M[i][j][0];
					if(pv){
						// check pref[i]+M[i][j]...
						if (searchRevRange(0, 0, nod-1, pv, pref[i][0], pref[i][0]+pref[i][1]-1)){
							doc = getNum64(DocLZ, pv*lgD, lgD);
							if(doc == 487)
								cout << "";
							if (suff[j][0] > 0){
								// can be 'M phase' and 'suffix' concatenated ?
								if (M[i][j][2] >= suff[j][0] && M[i][j][2] <= (suff[j][0]+suff[j][1]-1)){
									//cout << "Occ t3 doc " << doc << endl;
									x = dictD[doc];	// key
									k = *key;
									if(x)
										(keyFq[x])++;
									else{
										// x==0 then it is a new ID
										dictD[doc] = k;
										keyDoc[k] = doc;
										keyFq[k] = 1;
										//cout << doc << " Occ type3 [" << *key << "](" << keyFq[*key] << ") **************" << endl;
										(*key)++;
									}
								}
							}

							// to search concatenations of full phrase...
							l=i;r=j;
							while(M[l][r][1]){
								k=l;
								l = r+1;
								r = M[k][r][1];
								if (suff[r][0] > 0){
									// can be 'M phase' and 'suffix' concatenated ?
									if (M[l][r][2] >= suff[r][0] && M[l][r][2] <= (suff[r][0]+suff[r][1]-1)){
										//cout << "Occ type3 doc " << doc << endl;
										x = dictD[doc];	// key
										k = *key;
										if(x)
											(keyFq[x])++;
										else{
											// x==0 then it is a new ID
											dictD[doc] = k;
											keyDoc[k] = doc;
											keyFq[k] = 1;
											//cout << doc << " type3 [" << *key << "](" << keyFq[*key] << ") **************" << endl;
											(*key)++;
										}
									}
								}
							}

						}
					}
				}
			}
		}
	}
	//cout << "Found " << nOcc << " occ" << endl;
}

//=================================================================================================
// the doc with key dictD[k] go up until its correct position in the maximum queue
void TopkLZ::swimInMaxQ(uint k){
	uint sw=heap[k];
	ulong doc=keyDoc[k];
	float num=keyFq[sw];

	while(k > 1 && keyFq[heap[k>>1]] < num){
		heap[k] = heap[k>>1];
		k >>= 1;
	}
	while(k > 1 && keyFq[heap[k>>1]] == num && docLength[keyDoc[k>>1]] < docLength[doc]){
		heap[k] = heap[k>>1];
		k >>= 1;
	}
	heap[k] = sw;
}

// set the new item at top of the queue heap[1..k]
void TopkLZ::setTopMaxQ(uint k){
	uint j, m=2, i=1;
	ulong doc=keyDoc[k];

	float num=keyFq[heap[k]];

	while(m<k){
		j=m+1;
		if (j < k){
			if (keyFq[heap[m]] < keyFq[heap[j]])
				m=j;
			else{
				if (keyFq[heap[m]] == keyFq[heap[j]] && docLength[keyDoc[m]] < docLength[keyDoc[j]])
					m=j;
			}
		}
		if(num <= keyFq[heap[m]]){
			if(num < keyFq[heap[m]] || docLength[doc] < docLength[keyDoc[m]]){
				heap[i] = heap[m];
				i = m;
				m <<= 1;
			}else
				break;
		}else
			break;
	}
	heap[i] = heap[k];
}

// save the Data Structure in file 'fileName'
void TopkLZ::saveDS(bool showSize){
	char *fileName = new char[300];
	cout << "Save data structure in folder " << dirStore << endl;

	strcpy(fileName, dirStore);
	strcat(fileName, "dataStructures.tk");
	ofstream os (fileName, ios::binary);
	cout << "   Saving. Data structure size: " << sizeDS << endl;

	os.write((const char*)&n, sizeof(ulong));
	os.write((const char*)&nDocs, sizeof(uint));
	os.write((const char*)&lgD, sizeof(ulong));
	os.write((const char*)&g, sizeof(uint));
	os.write((const char*)&cutDoc, sizeof(char));
	os.write((const char*)&nod, sizeof(ulong));
	os.write((const char*)&lgNod, sizeof(uint));
	os.write((const char*)&nF, sizeof(ulong));
	os.write((const char*)&nEF, sizeof(ulong));
	os.write((const char*)&nEFHead, sizeof(uint));
	os.write((const char*)&nodRev, sizeof(ulong));
	os.write((const char*)&nLz, sizeof(ulong));
	os.write((const char*)&nRev, sizeof(ulong));
	os.write((const char*)&nTopk, sizeof(ulong));

	sizeDSav += 9*sizeof(ulong) + 4*sizeof(uint) + sizeof(char);	// size for variables

	os.write((const char*)docLength, nDocs*sizeof(ulong));				// save docLength[]
	sizeDSav += nDocs*sizeof(ulong);
	if(showSize) cout << " .- docLength[] " << nDocs*sizeof(ulong) << " Bytes" << endl;

	os.write((const char*)LbRev, nodRev*sizeof(uchar));					// save LbRevF[]
	sizeDSav += nodRev*sizeof(uchar);
	if(showSize) cout << " .- LbRev[] " << nodRev*sizeof(uchar) << " Bytes" << endl;

	os.write((const char*)LbRevF, nEF*sizeof(uchar));					// save LbRevF[]
	sizeDSav += nEF*sizeof(uchar);
	if(showSize) cout << " .- LbRevF[] " << nEF*sizeof(uchar) << " Bytes" << endl;
	os.close();

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rrr.tk");
	store_to_file(fRev_rrr, fileName);
	sizeDSav += size_in_bytes(fRev_rrr);
	if(showSize) cout << " .- fRev_rrr " << size_in_bytes(fRev_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rank.tk");
	store_to_file(fRev_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rrr.tk");
	store_to_file(fictU_rrr, fileName);
	sizeDSav += size_in_bytes(fictU_rrr);
	if(showSize) cout << " .- fictU_rrr " << size_in_bytes(fictU_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rank.tk");
	store_to_file(fictU_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "listF_rrr.tk");
	store_to_file(listF_rrr, fileName);
	sizeDSav += size_in_bytes(listF_rrr);
	if(showSize) cout << " .- listF_rrr " << size_in_bytes(listF_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "lisF_sel.tk");
	store_to_file(lisF_sel, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bTopk_rrr.tk");
	store_to_file(bTopk_rrr, fileName);
	sizeDSav += size_in_bytes(bTopk_rrr);
	if(showSize) cout << " .- bTopk_rrr " << size_in_bytes(bTopk_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bTopk_rank.tk");
	store_to_file(bTopk_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bDL_b.tk");
	store_to_file(bDL_b, fileName);
	sizeDSav += size_in_bytes(bDL_b);
	if(showSize) cout << " .- bDL_b " << size_in_bytes(bDL_b) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_rrr.tk");
	store_to_file(limDTk_rrr, fileName);
	sizeDSav += size_in_bytes(limDTk_rrr);
	if(showSize) cout << " .- limDTk_rrr " << size_in_bytes(limDTk_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_rank.tk");
	store_to_file(limDTk_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_sel.tk");
	store_to_file(limDTk_sel, fileName);

	cout << "   Total bytes saved from data structure: " << sizeDSav << endl;
}

// load the Data Structure from the file 'fileName'
void TopkLZ::loadDS(bool showSize){
	ulong auxSize;
	cout << " Load data structure from " << dirStore << endl;
	char *fileName = new char[400];
	sizeDS = 0;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "dataStructures.tk");
	ifstream is(fileName, ios::binary);

	is.read((char*)&n, sizeof(ulong));
	is.read((char*)&nDocs, sizeof(uint));
	is.read((char*)&lgD, sizeof(ulong));
	is.read((char*)&g, sizeof(uint));
	is.read((char*)&cutDoc, sizeof(char));
	is.read((char*)&nod, sizeof(ulong));
	is.read((char*)&lgNod, sizeof(uint));
	is.read((char*)&nF, sizeof(ulong));
	is.read((char*)&nEF, sizeof(ulong));
	is.read((char*)&nEFHead, sizeof(uint));
	is.read((char*)&nodRev, sizeof(ulong));
	is.read((char*)&nLz, sizeof(ulong));
	is.read((char*)&nRev, sizeof(ulong));
	is.read((char*)&nTopk, sizeof(ulong));
	//cout << "nDocs " << nDocs << endl;

	sizeDS += 9*sizeof(ulong) + 4*sizeof(uint) + sizeof(char);	// size for variables

	docLength = new ulong[nDocs];
	is.read((char*)docLength, nDocs*sizeof(ulong));
	sizeDS += nDocs*sizeof(ulong);
	if(showSize) cout << " ** size of docLength[] " << nDocs*sizeof(ulong) << " Bytes" << endl;

	LbRev = new uchar[nodRev];
	is.read((char*)LbRev, nodRev*sizeof(uchar));
	sizeDS += nodRev*sizeof(uchar);
	if(showSize) cout << " ** size of LbRev[] " << nodRev*sizeof(uchar) << " Bytes" << endl;

	LbRevF = new uchar[nEF];
	is.read((char*)LbRevF, nEF*sizeof(uchar));
	sizeDS += nEF*sizeof(uchar);
	if(showSize) cout << " ** size of LbRevF[] " << nEF*sizeof(uchar) << " Bytes" << endl;
	is.close();

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rrr.tk");
	load_from_file(fRev_rrr, fileName);
	sizeDS += size_in_bytes(fRev_rrr);
	if(showSize) cout << " ** size of fRev_rrr " << size_in_bytes(fRev_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rank.tk");
	load_from_file(fRev_rank, fileName);
	util::init_support(fRev_rank, &fRev_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rrr.tk");
	load_from_file(fictU_rrr, fileName);
	sizeDS += size_in_bytes(fictU_rrr);
	if(showSize) cout << " ** size of fictU_rrr " << size_in_bytes(fictU_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rank.tk");
	load_from_file(fictU_rank, fileName);
	util::init_support(fictU_rank, &fictU_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "listF_rrr.tk");
	load_from_file(listF_rrr, fileName);
	sizeDS += size_in_bytes(listF_rrr);
	if(showSize) cout << " ** size of listF_rrr " << size_in_bytes(listF_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "lisF_sel.tk");
	load_from_file(lisF_sel, fileName);
	util::init_support(lisF_sel, &listF_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bTopk_rrr.tk");
	load_from_file(bTopk_rrr, fileName);
	sizeDS += size_in_bytes(bTopk_rrr);
	if(showSize) cout << " ** size of bTopk_rrr " << size_in_bytes(bTopk_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bTopk_rank.tk");
	load_from_file(bTopk_rank, fileName);
	util::init_support(bTopk_rank, &bTopk_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bDL_b.tk");
	load_from_file(bDL_b, fileName);
	sizeDS += size_in_bytes(bDL_b);
	if(showSize) cout << " ** size of bDL_b " << size_in_bytes(bDL_b) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_rrr.tk");
	load_from_file(limDTk_rrr, fileName);
	sizeDS += size_in_bytes(limDTk_rrr);
	if(showSize) cout << " ** size of limDTk_rrr " << size_in_bytes(limDTk_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_rank.tk");
	load_from_file(limDTk_rank, fileName);
	util::init_support(limDTk_rank, &limDTk_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_sel.tk");
	load_from_file(limDTk_sel, fileName);
	util::init_support(limDTk_sel, &limDTk_rrr);

	// docTopk
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "docTopk.tk");
	ifstream isdocTopk(fileName, ios::binary);
	auxSize = limDTk_rrr.size()*lgD/W64;
	if ((limDTk_rrr.size()*lgD)%W64)
		auxSize++;
	docTopk = new ulong[auxSize];
	isdocTopk.read((char*)docTopk, auxSize*sizeof(ulong));
	sizeDS += auxSize*sizeof(ulong);
	if(showSize) cout << " ** size of docTopk[] " << auxSize*sizeof(ulong) << " Bytes" << endl;
	isdocTopk.close();

	// tries
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "LzTrie.tk");
	treeLz = new RangeMMTree64(fileName, false);
	sizeDS += treeLz->sizeRMM;
	if(showSize) cout << " ** size of treeLz " << treeLz->sizeRMM << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "RevTrie.tk");
	treeRev = new RangeMMTree64(fileName, false);
	treeRev->labels = LbRev;
	sizeDS += treeRev->sizeRMM;
	if(showSize) cout << " ** size of treeRev " << treeRev->sizeRMM << " Bytes" << endl;

	// DocLZ
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "DocLZ.tk");
	ifstream isDocLZ(fileName, ios::binary);
	auxSize = nod*lgD/W64;
	if ((nod*lgD)%W64)
		auxSize++;
	DocLZ = new ulong[auxSize];
	isDocLZ.read((char*)DocLZ, auxSize*sizeof(ulong));
	sizeDS += auxSize*sizeof(ulong);
	if(showSize) cout << " ** size of DocLZ[] " << auxSize*sizeof(ulong) << " Bytes" << endl;
	isDocLZ.close();

	// NODE
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "Node.tk");
	ifstream isNode(fileName, ios::binary);
	auxSize = nod*lgNod/W64;
	if ((nod*lgNod)%W64)
		auxSize++;
	Node = new ulong[auxSize];
	isNode.read((char*)Node, auxSize*sizeof(ulong));
	sizeDS += auxSize*sizeof(ulong);
	if(showSize) cout << " ** size of Node[] " << auxSize*sizeof(ulong) << " Bytes" << endl;
	isNode.close();

	// Range
	this->h = ceilingLog64(nod,2);
	Range = new LevelRange[h];
	sizeDS += h*sizeof(LevelRange);
	if(showSize) cout << " ** Size of Range[]'s pointers " << h*sizeof(LevelRange) << " = " << (h*sizeof(LevelRange))*8.0/(float)this->n << " bpc" << endl;

	ulong sizeBitMapsRange = 0;
	char str[100];
	for(uint i=0; i<h; i++){
		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcpy(str, "");
		sprintf(str, "Range_%d.B_rrr.tk", i);
		strcat(fileName, str);
		load_from_file(Range[i].B_rrr, fileName);
		sizeBitMapsRange += size_in_bytes(Range[i].B_rrr);

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcpy(str, "");
		sprintf(str, "Range_%d.B_rank0.tk", i);
		strcat(fileName, str);
		load_from_file(Range[i].B_rank0, fileName);
		util::init_support(Range[i].B_rank0, &(Range[i].B_rrr));

		strcpy(fileName, "");
		strcpy(fileName, dirStore);
		strcpy(str, "");
		sprintf(str, "Range_%d.B_rank1.tk", i);
		strcat(fileName, str);
		load_from_file(Range[i].B_rank1, fileName);
		util::init_support(Range[i].B_rank1, &(Range[i].B_rrr));
	}
	if(showSize)
		cout << " ** Range[]'s BitVector(RRR) " << sizeBitMapsRange << " Bytes" << endl;
	sizeDS += sizeBitMapsRange;

	cout << " Data Structure loaded !!" << endl;
}

void TopkLZ::destroyConstruct() {
	sizeDS = n = nDocs = lgD = g = cutDoc = nod = lgNod = 0;
	nF = nEF = nEFHead = nodRev = nLz = nRev = nTopk = 0;

	PLz = PRev = NULL;

	delete []docTopk;
	delete []docLength;
	delete []EndDocsPos;
	delete []charTf;
	delete []treeLz->P;
	delete []treeRev->P;
	delete []LbRev;
	delete []LbRevF;
	delete []dictD;
	delete []keyDoc;
	delete []keyFq;
	delete []heap;

	cout << "destroyeing treeRev..." << endl;
	treeRev->~RangeMMTree64();

	// to delete sdsl objects...
	{
	    bit_vector empty_b;
	    bDL_b.swap(empty_b);

	    rrr_vector<127> empty_rrr;
	    fRev_rrr.swap(empty_rrr);

	    rrr_vector<127> empty_rrr2;
	    fictU_rrr.swap(empty_rrr2);

	    rrr_vector<127> empty_rrr3;
	    listF_rrr.swap(empty_rrr3);

	    rrr_vector<127> empty_rrr4;
	    bTopk_rrr.swap(empty_rrr4);

	    rrr_vector<127> empty_rrr5;
	    limDTk_rrr.swap(empty_rrr5);
	}

}

TopkLZ::~TopkLZ() {

	// to delete Range
	for(uint i=0; i<this->h; i++){
		decltype(Range[i].B_rrr) empty;
		Range[i].B_rrr.swap(empty);

		decltype(Range[i].B_rank0) empty2;
		Range[i].B_rank0.swap(empty2);

		decltype(Range[i].B_rank1) empty3;
		Range[i].B_rank1.swap(empty3);
	}

	delete [] Range;

	sizeDS = n = nDocs = lgD = g = cutDoc = nod = lgNod = 0;
	nF = nEF = nEFHead = nodRev = nLz = nRev = nTopk = 0;

	delete []docLength;
	delete []treeLz->P;
	delete []treeRev->P;
	delete []LbRev;
	delete []LbRevF;
	delete []DocLZ;
	delete []Node;
	delete []docTopk;
	delete []dictD;
	delete []keyDoc;
	delete []keyFq;
	delete []heap;

	cout << "destroyeing treeLZ..." << endl;
	treeLz->~RangeMMTree64();
	cout << "destroyeing treeRev..." << endl;
	treeRev->~RangeMMTree64();

	// to delete sdsl objects...
	{
	    bit_vector empty_b;
	    bDL_b.swap(empty_b);

	    rrr_vector<127> empty_rrr;
	    fRev_rrr.swap(empty_rrr);

	    rrr_vector<127> empty_rrr2;
	    fictU_rrr.swap(empty_rrr2);

	    rrr_vector<127> empty_rrr3;
	    listF_rrr.swap(empty_rrr3);

	    rrr_vector<127> empty_rrr4;
	    bTopk_rrr.swap(empty_rrr4);

	    rrr_vector<127> empty_rrr5;
	    limDTk_rrr.swap(empty_rrr5);
	}
}

