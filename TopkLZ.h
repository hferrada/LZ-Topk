/*
 * TopkLZ.h
 *
 *  Created on: 19-08-2015
 *      Author: hector
 */

#ifndef TOPKLZ_H_
#define TOPKLZ_H_

#include "LZ78Tries64.h"
#include "RangeMMTree64.h"

using namespace sdsl;
using namespace std;
using namespace drf64;

class TopkLZ {

	// structure for wt of Range, each level has a biststring B_rrr (of legth n'-1) and a rmq for its respective pointer of documents
	typedef struct {
		rrr_vector<127> B_rrr;
		rrr_vector<127>::rank_1_type B_rank1;
		rrr_vector<127>::rank_0_type B_rank0;
	} LevelRange;

	typedef struct {
		bit_vector B_b;
	} BaseRange;

private:
	ulong lgD;				// Binary logarithm (ceiling) of nDocs
	char boundSymbol;		// original symbol delimiter of documents when we read all documents in 1 file.
	ulong* docLength;		// store the length for each document
	ulong* charTf;			// store the tf for each symbol. This is no include in the final structure

	ulong nod;				// number of phrase nodes nod. Here nod = n'+1, n' of LZIndex
	uint lgNod;				// ceiling binary log of nod
	ulong nF;				// number of fictitious nodes in RevTrie, nF only include one fictitious node for each unary path of fictitious nodes
	ulong nEF;				// this is the difference between the total fictitious nodes in the original RevTrie an the RevTrie represented where we join the fictitious nodes in unary path
	uint nEFHead;			// number of fictitious nodes are head of a unary path
	ulong nodRev;			// nod + nF
	ulong *PLz;				// LZTrie's DFUDS sequence representation in 2*nod bits
	ulong nLz;				// length of PLz DFUDS sequences, these have exactly nLz = 2*nod nodes
	ulong nRev;				// length of PRev DFUDS sequences, these have exactly nRev = 2(nod+nF) nodes
	ulong *PRev;			// RevTrie's DFUDS sequence representation in 2(nod+nF) bits

	// Labels...
	uchar *LbRev;			// RevTrie's labels, for the n' phrase nodes, in DFUDS representation
	bit_vector fRev_b;		// this marks the fictitious nodes in LbRev
	rrr_vector<127> fRev_rrr;
	rrr_vector<127>::rank_1_type fRev_rank;
	uchar *LbRevF;			// RevTrie's labels for fictitious nodes with unary paths
	bit_vector fictU_b;		// this marks the fictitious nodes with unary path in FRev
	rrr_vector<127> fictU_rrr;
	rrr_vector<127>::rank_1_type fictU_rank;
	bit_vector listF_b;		// Delimiter for each unary sequence in LbRevF
	rrr_vector<127> listF_rrr;
	rrr_vector<127>::select_1_type lisF_sel;

	RangeMMTree64 *treeLz;	// range min-max tree representations for PLz
	RangeMMTree64 *treeRev;	// range min-max tree representations for PRev

	ulong *DocLZ;			// Document array to phrases with length nod-1
	ulong *Node;			// This array give the preorder value in LZTrie, that is: Range[j] = i <-> the node with preorder j in RevTrie corresponds to the node with preorder i in LZTrie

	// MARKS NODES...
	ulong nTopk;			// number of nodes with topk list
	bit_vector bTopk_b;		// marks all node that stores a topk List, length: [1..n']. This has nTopk 1's and z 0's
	rrr_vector<127> bTopk_rrr;
	rrr_vector<127>::rank_1_type bTopk_rank;

	bit_vector bDL_b;		// for the bTopk nodes marked, it marks the nodes that store all its different documents

	// TOP-K LIST...
	ulong *docTopk;			// this are the nTopk lists of id documents, with log(D) bits per item, maximum k id per list
	bit_vector limDTk_b;	// this marks the boundaries for each list in docTopk
	rrr_vector<127> limDTk_rrr;
	rrr_vector<127>::rank_1_type limDTk_rank;
	rrr_vector<127>::select_1_type limDTk_sel;

	// sort the list in increasing mode by frequency
	void sortTopkDocsByFq(uint* fNod, uint k, uint* dNodK, uint* fNodK);

	void create_sNode(LZ78Tries64::RevNode* nod, ulong *pre, ulong *aux_sNode, ulong *posNod);

	// search the current pattern (given the cross of the RevTrie) in order to construct the real top-k answer
	bool searchPatternInFMI(uchar* pat, uint m, uint *fqCSA);

	// Compute the document frequencies for the 'query' and stored it in 'fqCSA'
	void freqPattInFMI(uchar *patt, uint m, ulong occs, uint *fqCSA);

	// determine the doc Id for the text position
	uint searchDocPosition(ulong pos);

	// marks the nodes that store a topk list in bTopk[]. Set the bits in bBlockLZ[], and filled the array fatherLZ[1..n']
	void marksTopkNodes(LZ78Tries64::RevNode *node, ulong *pre, uint childFather, uchar *patt, uint m);

	// create the document array for LZTrie's phrases
	void genDocArray(LZ78Tries64::LZNode* nod, ulong *DocLZ, ulong *ini);

	//	create 	LbRevF	:	RevTrie's labels for fictitious nodes with unary paths
	// 			ListF	:	Delimiter for each unary sequence in LbRevF
	// 			fRev	:	this marks the fictitious nodes in LbRev (RavTrie's labels)
	// 			FictU	:	this marks the fictitious nodes with unary path in FRev
	void genSeqExtraFictRevTrie(LZ78Tries64::RevNode *node, ulong *posExtF, ulong *pfRev, ulong *posFU);

	//	create 	PRev, LbRev and LbRevF structures
	void genDFUDSSeqRevTrie(LZ78Tries64::RevNode *node, ulong *pos, ulong *posRev);

	//	create 	PLZ	:	LZTrie's DFUDS representation without Labels
	void genDFUDSSeqLZTrieNotLbl(LZ78Tries64::LZNode *node, ulong *pos);

	void genNodeArray(LZ78Tries64::RevNode* nod, ulong *ArrRange, ulong *preRev);
	void createStructure(char dirSaveLoad[300]);

	// *** CREATE ID-DOCUMENTS
	void createDocTopkSeq(ulong *pP, ulong *pL, ulong *pre, uchar *patt, uint m,
							uint *fNod, uint *dNodK, uint *fNodK,
							ulong *docTopkAux, ulong *posD, ulong *nodTopk,
							bool *bDL_bAux, bool *limDTk_bAux);

	// create wavelet tree for range represented as an array of BitSequences...
	void createRange(ulong *ArrRange);
	void genNodoWTRange(BaseRange *netRange, uint level, ulong *ArrRange, ulong start, ulong length, ulong lowSymbol);
	ulong searchLZCoord(uint lev, ulong lcur, ulong rcur, ulong obj);
	void readListFiles(char *inputFile, bool bitLowercase);
	void readListFilesTrec(char *inputFile, ulong max_size);			// ELIMINAR
	void readUniqueFile(char *inputFile, bool bitLowercase);
	void makePriorities(ulong* cPosTb);

	// the doc with key dictD[k] go up until its correct position in the maximum queue
	void swimInMaxQ(uint k);

	// set the new item at top of the queue heap[1..k]
	void setTopMaxQ(uint k);

public:
	static bool TRACE;		// true: print all details for console
	static bool TEST;		// true: print all details for console

	ulong sizeDS;			// size in bytes for all data structure
	ulong sizeDSav;			// size in bytes for data structure saved in the system file;

	uchar *seq;				// original sequence (1 byte for symbol)
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	uint nDocs;				// Number of distinct documents D
	char cutDoc;			// symbol to separate documents
	uint g;					// minimum length of document segments to store a topk list in a node
	LZ78Tries64 *tries;		// structure to store the original LZTrie and RevTrie.

	ulong* EndDocs;			// this stores the final phrase number for each document. This is no include in the final structure
	ulong *EndDocsPos;		// this stores the final text position for each document. This is no include in the final structure

	// BRUTE FORCE FOR COUNT AND SORT FREQUENCIES...
	// ** IMPORTAT: dictD only can store keys of length of uint (16 bits)
	ulong *dictD;			// dictionary of keys, length: [1..D] with log(k*) bits each id-Doc
	ulong *keyDoc;			// document of each key, length: [1..k*] with lgD bits
	ulong *keyFq;			// frequencies for each key, length: [1..k*] with 8 bits
	uint *heap;				// maximum priority queue for keyFq, length: [1..k*]

	uint h;					// h = lg(n'-1) = lg(nod-2). Binary logarithm. Range's height
	LevelRange *Range;	    // Wavelet tree represented as an array of (h-1) Bitmaps of length (n'-1), with support rank/select. Range[1..h-1][i..n'-1]
							// it represents the pairs (j,i) <--> If node j in RevTrie has id=k, then node i in LZTrie has id=k+1

	char dirStore[200];		// directory to save/load the data structure (files *.tk)


	// this bit_vector is only for test search of patterns and not belong to the index!!
	bit_vector separator_b;		// This marks the beginning of each LZ-phrase
	rrr_vector<127> sep_rrr;
	rrr_vector<127>::rank_1_type sep_rank;
	rrr_vector<127>::select_1_type sep_sel;
	sdsl::csa_wt<wt_huff<rrr_vector<127> >, 4, 4> fmi;

	TopkLZ(uint gVal, char *inputFile, bool filesInList, char cutDocCode, bool bitLowercase, char dirSaveLoad[300], char boundS);
	TopkLZ(char dirSaveLoad[200], bool showSize);
	virtual ~TopkLZ();
	void destroyConstruct(); 	// destroying object did not salve in the construction

	// It stores in docList[1..k] and frqList[1..k] the list and frequencies of the top-k most frequent documents for the pattern p[1..m]
	ulong topKDocument(uchar* pat, uint m, uint* docList, uint k);
	void reportDocOcc2_levRange(uint lev, ulong lcur, ulong rcur, ulong l, ulong r, ulong lobj, ulong robj, uint* key);
	void searchIntervalInRange(uint lev, ulong lcur, ulong rcur, ulong lobj, ulong robj, ulong lrmq, ulong rrmq, uint *key);
	bool searchRevRange(uint lev, ulong lcur, ulong rcur, ulong obj, ulong lRev, ulong rRev);
	void searchPhase_Rev(uchar* pat, uint m, ulong *mat);
	void searchPhase_Rev3(uchar* pat, uint m, ulong *mat);
	void searchPhase_Rev2(uchar* pat, uint m, ulong *mat);
	void extractFqType23(uchar* pat, uint m, uint *key);

	// save the Data Structure in folder 'dirStore'
	void saveDS(bool showSize);

	// load the Data Structure from the folder 'fileName'
	void loadDS(bool showSize);

	// determine the doc Id for the phrase idNode
	uint searchDocument(ulong idNode);

	// return true if the pattern is in the RevTrie and in '*pv' the respective preorder value in the tree,
	// in '*x' stores its DFUDS bit position. Return false if the pattern does not exist
	bool searchPattern_Rev(uchar* pat, uint m, ulong *x, ulong *pv);
	bool searchPattern_Rev2(uchar* pat, uint m, ulong *x, ulong *pv, bool realP);
	bool searchPattern_Rev3(uchar* pat, uint m, ulong *x, ulong *pv, bool realP);
};

#endif /* TOPKLZ_H_ */
