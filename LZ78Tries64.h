/*
 * LZ78Tries64.h
 *
 *  Created on: 01-07-2014
 *      Author: hector
 */

#ifndef LZ78TRIES64_H_
#define LZ78TRIES64_H_

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/int_vector.hpp>
#include <assert.h>
#include <Basic_drf64.h>

using namespace sdsl;
using namespace std;
using namespace drf64;

class LZ78Tries64 {
public:

	class LZNode {
	public:
		uchar symbol;		// edge symbol from its father to this node
		long int idNode;	// node's ID
		ulong preorder;		// preorder value
		uint nChildren;		// number of children
		LZNode** hTab;
		LZNode* fChild;
		LZNode* nextSib;
	};

	class RevNode{
	public:
		uchar symbol;		// edge symbol from its father to this node
		long int idNode;	// node's ID
		ulong preorder;		// preorder value
		uint nChildren;		// number of children
		bool fict;
		bool uPath;
		LZNode* nodeLZ;
		RevNode** hTab;
		RevNode* fChild;
		RevNode* nextSib;
	};

	LZNode* lzTrie;
	RevNode* revTrie;

	LZ78Tries64(uint nDocs, uchar endSymbolDoc, ulong* symbolHTable);
	virtual ~LZ78Tries64();
	void destroyLZNode(LZNode* nod);
	void destroyRevNode(RevNode* nod);

	void destroyNodeTrie(LZNode* nod);												// destroy a NodeTrie
	void createTriesFromText(uchar* text, ulong n, ulong *endDocs, ulong *endDocsPos, bit_vector *separator_b);
	void insertPhrase(uchar* text, ulong* textPos, LZNode* newNode);				// insert a new phrase in LZtrie
	void insertRevPhrase(uchar* text, ulong posIni, ulong posFin, LZNode* nodeLZ);	// insert the reverse phrase in RevTrie

	void sortLZNodes(LZNode* nod);
	void sortRevNodes(RevNode* nod);

	// this assigns the preorder value for each LZTrie's node
	void setPreorderValues(LZNode* nod, ulong *pre);

	// create Node array
	void genNodeArray(RevNode* nod, ulong *preRev, ulong *Node);

	// determine the number of fictitious nodes in RevTrie, this only include one fictitious node for each unary path of fictitious nodes
	void countFictNodRevTrie(RevNode *node);


	//	create 	LbRevF	:	RevTrie's labels for fictitious nodes with unary paths
	// 			ListF	:	Delimiter for each unary sequence in LbRevF
	// 			fRev	:	this marks the fictitious nodes in LbRev (RavTrie's labels)
	// 			FictU	:	this marks the fictitious nodes with unary path in FRev
	void genSeqExtraFictRevTrie(RevNode *node, uchar *LbRevF, uint *posExtF, bit_vector *ListF, bit_vector *fRev, uint *pfRev, bit_vector *FictU, uint *posFU);

	//	create 	PRev, LbRev and LbRevF structures
	void genDFUDSSeqRevTrie(RevNode *node, ulong* PRev, ulong *pos, uchar *LbRev, ulong *pLbRev);

	// return the preorder value of the most right leaf of this subtree in a LZTrie
	ulong getLastLZLeaf(LZNode *node);

	// it returns the length of the subtree rooted at node
	void getLengthSubTree(RevNode *nod, ulong *len);

	void listNodeLZ(LZNode* nod, ulong* preorder, uint deph);		// list a NodeTrie
	void listNodeRev(RevNode* nod, ulong* preorder, uint deph);		// list a NodeTrie
	void listLZTrie();												// list the LZTrie
	void listRevTrie();												// list the RevTrie

	ulong *cPosTb;
	uchar cutDoc;		// last symbol for each doc
	uint nDocs;			// number of documents
	ulong nPhra;		// number of LZTrie's phrases (n' in the paper)
	uint nFictU;		// number of LZTrie's fictitious nodes in RevTrie, this include only one fictitious node for each unary path of fictitious nodes
	uint nExtraFict;	// this is the difference between the total fictitious nodes in the original RevTrie an the RevTrie represented where we join the fictitious nodes in unary path
	long int bad;		// partial number of fictitious nodes
	uchar lgPhr;		// ceiling binary log of nPhra
	uchar lgD;			// ceiling binary log of nDocs
	ulong *IdDocPreLZ;	// array used to construct range structure. For each node in LZTrie we stores its respective preorder value.

	static bool TRACE;	// true: print all details for console
	static uint SIGMA;
	static uint LENHASH;

};

#endif /* LZ78TRIES64_H_ */
