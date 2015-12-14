/*
 * LZ78Tries64.cpp
 *
 *  Created on: 01-07-2014
 *      Author: hector
 */

#include "LZ78Tries64.h"
bool LZ78Tries64::TRACE = false;
uint LZ78Tries64::SIGMA = 256;
uint LZ78Tries64::LENHASH = 16;

LZ78Tries64::LZ78Tries64(uint nDocs, uchar endSymbolDoc, ulong* symbolHTable) {
	lzTrie = new LZNode();
	revTrie = new RevNode();

	lzTrie->idNode = 0;
	lzTrie->nChildren = 0;
	lzTrie->preorder = 0;
	lzTrie->nextSib = lzTrie->fChild = NULL;

	revTrie->idNode = 0;
	revTrie->nChildren = 0;
	revTrie->preorder = 0;
	revTrie->fict = 0;
	revTrie->nodeLZ = lzTrie;
	revTrie->nextSib = revTrie->fChild = NULL;

	this->cPosTb = symbolHTable;
	this->cutDoc = endSymbolDoc;
	this->nPhra = 0;
	this->nFictU = 0;
	this->nExtraFict = 0;
	this->nDocs = nDocs;
	this->bad = -1;
	this->lgD = ceilingLog64(nDocs, 2);
}

void LZ78Tries64::createTriesFromText(uchar* text, ulong n, ulong *endDocs, ulong *endDocsPos, bit_vector *separator_b){
	ulong posIni, textPos;
	uint nDoc;
	LZNode* nodLZ;

	nDoc = textPos = 0;
	while (textPos < n){
		nPhra++;
		(*separator_b)[textPos] = 1;

		//if(textPos >= 1600) //"...ada$7"
			//cout << "";

		if (text[textPos] == cutDoc){
			//nPhra--; //we do not omite this id because we can create a wrong occurrence type two!! ?????
			endDocs[nDoc] = nPhra;
			endDocsPos[nDoc] = textPos;
			//(*separator_b)[textPos] = 0;
			if (nDoc % 10000 == 0) cout << "nDoc = " << nDoc << ", textPos = " << textPos << ", nPhra = " << nPhra << endl;
			nDoc++;
			posIni = textPos;
			nodLZ = new LZNode();
			nodLZ->idNode = nPhra;
			nodLZ->nChildren = 0;
			nodLZ->nextSib = nodLZ->fChild = NULL;

			// create an insert the new child nodLZ in lzTrie
			{
				LZNode* father = lzTrie;
				nodLZ->symbol = cutDoc;
				nodLZ->symbol = cutDoc;
				(father->nChildren)++;
				if(father->fChild == NULL)
					father->fChild = nodLZ;
				else{
					nodLZ->nextSib = father->fChild;
					father->fChild = nodLZ;
				}
				if (father->hTab == NULL){
					father->hTab = new LZNode*[LENHASH];
					for(uint i=0; i<LENHASH; i++)
						father->hTab[i] = NULL;
				}
			}
			{
				// we insert a special phrase for each terminal symbol
				RevNode *newNode = new RevNode();
				newNode->idNode = nPhra;
				newNode->symbol = cutDoc;
				newNode->nChildren = 0;
				newNode->fict = newNode->uPath = false;
				newNode->nodeLZ = nodLZ;
				newNode->nextSib = revTrie->fChild;
				revTrie->fChild = newNode;
				(revTrie->nChildren)++;
			}
		}else{
			posIni = textPos;
			nodLZ = new LZNode();
			nodLZ->idNode = nPhra;
			nodLZ->nChildren = 0;
			nodLZ->nextSib = nodLZ->fChild = NULL;

			// create an insert the new child
			insertPhrase(text, &textPos, nodLZ);
			if (text[textPos] == cutDoc){
				endDocs[nDoc] = nPhra;
				endDocsPos[nDoc] = textPos;
				if (nDoc % 10000 == 0)
					cout << "nDoc = " << nDoc << ", textPos = " << textPos << ", nPhra = " << nPhra << endl;
				nDoc++;

				// we insert a special phrase for each terminal symbol
				RevNode *newNode = new RevNode();
				newNode->idNode = nPhra;
				newNode->symbol = cutDoc;
				newNode->nChildren = 0;
				newNode->fict = newNode->uPath = false;
				newNode->nodeLZ = nodLZ;
				newNode->nextSib = revTrie->fChild;
				revTrie->fChild = newNode;
				(revTrie->nChildren)++;

				if(textPos >= n)
					break;
			}else
				insertRevPhrase(text, posIni, textPos, nodLZ);
		}
		textPos++;
	}
	cout << "Final nDoc = " << nDoc << ", textPos = " << textPos << ", nPhra = " << nPhra << endl;
	if(nDoc != nDocs){
		cout << "ERROR. Different number of documents, nDocs = " << nDocs << " != nDoc = "  << nDoc << endl;
		exit(1);
	}
	nPhra++;
	cout << "LZTrie OK ! (n'=" << nPhra << ")" << endl;
	lgPhr = ceilingLog64(nPhra, 2);

	// sort nodes in each trie...
	sortLZNodes(this->lzTrie);
	cout << "LZNodes sorted !"<< endl;

	sortRevNodes(this->revTrie);
	cout << "RevNodes sorted !"<< endl;

	// set preorder values for each node in LZTrie...
	ulong pre = 0;
	cout << "To set preorder values for each node in LZTrie..."<< endl;
	IdDocPreLZ = new ulong[nPhra];
	setPreorderValues(this->lzTrie, &pre);
	if (nPhra != pre){
		cout << "total LZTrie preorder's = " << pre << " != nPhra = " << nPhra << endl;
		exit(1);
	}
}

void LZ78Tries64::insertPhrase(uchar* text, ulong* textPos, LZNode* newNode){
	uchar letter = text[*textPos];
	LZNode* father = lzTrie;
	bool doesntInsert = true;

	while(doesntInsert){
		if(father->nChildren == 0){						// 1.- father without any children
			father->hTab = new LZNode*[LENHASH];
			for(uint i=0; i<LENHASH; i++)
				father->hTab[i] = NULL;
			father->nChildren = 1;
			newNode->symbol = letter;
			father->hTab[cPosTb[letter]] = newNode;
			doesntInsert = false;
		}else{
			suint pos = cPosTb[letter];

			if(father->hTab[pos] == NULL){			// 2.- father without any child for this symbol 'letter'
				newNode->symbol = letter;
				father->hTab[pos] = newNode;
				(father->nChildren)++;
				doesntInsert = false;
			}else{
				LZNode* child = father->hTab[pos];

				if(child->symbol == letter){			// 3.- There is a child for this symbol 'letter'
					father = child;
					(*textPos)++;
					letter = text[*textPos];
					if (letter == cutDoc)
						break;
				}else{									// 4.- There are other children in this slot but with different symbol
					while(child->nextSib && child->nextSib->symbol != letter)
						child = child->nextSib;

					if(child->nextSib && child->nextSib->symbol == letter){
						father = child->nextSib;
						(*textPos)++;
						letter = text[*textPos];
						if (letter == cutDoc)
							break;
					}else{
						(father->nChildren)++;
						newNode->symbol = letter;
						child->nextSib = newNode;
						doesntInsert = false;
					}
				}
			}
		}
	}

	if (letter == cutDoc){	// special case, end of document --> we will insert at the beginning
		newNode->symbol = cutDoc;
		(father->nChildren)++;
		if(father->fChild == NULL)
			father->fChild = newNode;
		else{
			newNode->nextSib = father->fChild;
			father->fChild = newNode;
		}

		if (father->hTab == NULL){
			father->hTab = new LZNode*[LENHASH];
			for(uint i=0; i<LENHASH; i++)
				father->hTab[i] = NULL;
		}
	}
}

void LZ78Tries64::insertRevPhrase(uchar* text, ulong posIni, ulong posFin, LZNode* nodeLZ){
	ulong cont, len = posFin - posIni + 1;
	uchar letter;
	RevNode *auxNode;
	RevNode *father = revTrie;
	RevNode *newNode = new RevNode();
	newNode->idNode = nPhra;
	newNode->nChildren = 0;
	newNode->fict = newNode->uPath = false;
	newNode->nodeLZ = nodeLZ;
	newNode->nextSib = newNode->fChild = NULL;

	for(cont=0; cont < len; cont++){
		letter = text[posFin-cont];

		if(father->nChildren == 0 || father->hTab[cPosTb[letter]] == NULL){
			if(father->nChildren == 0){
				father->hTab = new RevNode*[LENHASH];
				for(uint i=0; i<LENHASH; i++)
					father->hTab[i] = NULL;
			}
			for(;cont<len-1; cont++){
				(father->nChildren)++;
				auxNode = new RevNode();
				auxNode->idNode = bad--;
				auxNode->uPath = false;
				auxNode->nChildren = 0;
				auxNode->fict = true;
				auxNode->symbol = letter;
				auxNode->nextSib = auxNode->fChild = NULL;
				father->hTab[cPosTb[letter]] = auxNode;
				father = auxNode;
				father->hTab = new RevNode*[LENHASH];
				for(uint i=0; i<LENHASH; i++)
					father->hTab[i] = NULL;
				letter = text[posFin-cont-1];
			}
			(father->nChildren)++;
			newNode->symbol = letter;
			father->hTab[cPosTb[letter]] = newNode;
		}else{
			RevNode* child = father->hTab[cPosTb[letter]];

			if(child->symbol == letter)			// 3.- There is a child for this symbol 'letter'
				if (child->fict && cont == len-1){
					delete newNode;				// child is transformed in true node
					child->fict = false;
					child->idNode = nPhra;
					child->nodeLZ = nodeLZ;
				}else{
					father = child;
					if (cont == len-1){
						newNode->symbol = letter;
						father->nextSib = newNode;
					}
				}
			else{									// 4.- There are other children in this slot but with different symbol
				while(child->nextSib && child->nextSib->symbol != letter)
					child = child->nextSib;

				if(child->nextSib && child->nextSib->symbol == letter){
					if(child->nextSib->fict && cont == len-1){
						delete newNode;				// child->nextSib is transformed in true node
						child->nextSib->fict = false;
						child->nextSib->idNode = nPhra;
						child->nextSib->nodeLZ = nodeLZ;
					}else
						father = child->nextSib;
				}else{
					(father->nChildren)++;
					if(cont<len-1){
						auxNode = new RevNode();
						auxNode->idNode = bad--;
						auxNode->uPath = false;
						auxNode->nChildren = 0;
						auxNode->fict = true;
						auxNode->symbol = letter;
						auxNode->nextSib = auxNode->fChild = NULL;
						child->nextSib = auxNode;
						father = auxNode;
					}else{
						newNode->symbol = letter;
						child->nextSib = newNode;
					}
				}
			}
		}
	}
}

void LZ78Tries64::sortLZNodes(LZNode* nod){
	LZNode** lNodesLZ = new LZNode*[SIGMA];
	LZNode* aChild;
	LZNode* aux;
	uint i, cH0;
	for (i=0; i<SIGMA; i++)
		lNodesLZ[i] = NULL;

	cH0=0;
	if (nod->hTab){ // 1. cuenta hijos distintos a $
		for (i=0; i<LENHASH; i++){
			if(nod->hTab[i]){
				aChild = nod->hTab[i];
				while(aChild){
					if(aChild->symbol != cutDoc){
						cH0++;
						lNodesLZ[aChild->symbol] = aChild;
					}
					aChild = aChild->nextSib;
				}
			}
		}
	}

	// set the first child
	if (nod->fChild){ // If nod->fChild is NOT NULL then nod->fChild->symbol == cutDoc, and also all its sibling
		aChild = nod->fChild;
		cH0++;
		if(aChild->symbol != cutDoc){
			cout << "EEEEEEERRRRRRRRRRR 0, fChild != cutDoc" << endl;
			exit(1);
		}
		while(aChild->nextSib){
			if(aChild->symbol != cutDoc){
				cout << "EEEEEEERRRRRRRRRRR 0, fChild != cutDoc" << endl;
				exit(1);
			}
			cH0++;
			aChild = aChild->nextSib;
		}
		i=0;
	}else{
		for (i=0; i<SIGMA && lNodesLZ[i] == NULL; i++);
		nod->fChild = lNodesLZ[i];
		aChild = nod->fChild;
		i++;
	}

	if (cH0 != nod->nChildren){
		cout << "ERR. cH0 = " << cH0 << " != nod->nChildren = " << nod->nChildren << ", Id: " << nod->idNode << endl;
		exit(0);
	}

	aChild->nextSib = NULL;
	for (; i<SIGMA; i++){
		if (lNodesLZ[i]){
			aux = lNodesLZ[i];
			aux->nextSib = NULL;
			aChild->nextSib = aux;
			aChild = aux;
		}
	}
	aChild->nextSib = NULL;
	if (nod->hTab)
		delete [] nod->hTab;

	if (nod->nChildren){
		aChild = nod->fChild;
		while(aChild && aChild->symbol == cutDoc)
			aChild = aChild->nextSib;

		while(aChild){
			if(aChild->nChildren)
				sortLZNodes(aChild);
			aChild = aChild->nextSib;
		}
	}
	delete [] lNodesLZ;
}

void LZ78Tries64::sortRevNodes(RevNode* nod){
	//cout << "Rev sort " << nod->idNode << " (" << nod->nChildren << "):" << endl;
	RevNode** lNodesRev = new RevNode*[SIGMA];
	RevNode* child;
	uint i;
	for (i=0; i<SIGMA; i++)
		lNodesRev[i] = NULL;

	uint cH0 = 0;
	if (nod->hTab){
		for (i=0; i<LENHASH; i++){
			if(nod->hTab[i]){
				child = nod->hTab[i];
				while(child){
					cH0++;
					lNodesRev[child->symbol] = child;
					child = child->nextSib;
				}
			}
		}
	}
	if (nod->fChild){ // If nod->fChild is NOT NULL then nod->fChild->symbol == cutDoc, and also all its sibling
		child = nod->fChild;
		if(child->symbol != cutDoc){
			cout << "EEEEEEERRRRRRRRRRR 1" << endl;
			exit(1);
		}
		cH0++;
		while(child->nextSib){
			cH0++;
			if(child->symbol != cutDoc){
				cout << "EEEEEEERRRRRRRRRRR 1" << endl;
				exit(1);
			}
			child = child->nextSib;
		}
		i=0;
	}else{
		for (i=0; i<SIGMA && lNodesRev[i] == NULL; i++);
		nod->fChild = lNodesRev[i];
		child=nod->fChild;
		i++;
	}

	if (cH0 != nod->nChildren){
		cout << "ERROR in sortRevNodes cH0 = " << cH0 << " != nod->nChildren = " << nod->nChildren << ", Id: " << nod->idNode << endl;
		exit(0);
	}

	for (; i<SIGMA; i++){
		if (lNodesRev[i]){
			child->nextSib = lNodesRev[i];
			child = lNodesRev[i];
		}
	}
	child->nextSib = NULL;
	if (nod->hTab)
		delete [] nod->hTab;

	if (nod->nChildren){
		child = nod->fChild;
		while(child->symbol == cutDoc)
			child = child->nextSib;

		while(child){
			if (child->nChildren)
				sortRevNodes(child);
			child = child->nextSib;
		}
	}
	delete [] lNodesRev;
}

// this assigns the preorder value for each LZTrie's node
void LZ78Tries64::setPreorderValues(LZNode* nod, ulong *pre){
	if (nod){
		//cout << "nod " << nod->idNode << " , pre: " << *pre << endl;
		LZNode* p = nod->fChild;
		nod->preorder = *pre;
		IdDocPreLZ[nod->idNode] = *pre;
		(*pre)++;

		for(uint i=0; i<nod->nChildren; i++){
			setPreorderValues(p, pre);
			p = p->nextSib;
		}
		if(p){
			cout << "ERROR in LzTrie: nod->nChildren = " << nod->nChildren << ". But there is another node p = " << p->idNode << endl;
			exit(1);
		}
	}
}

//	count the fictitious nodes and extra fictitious nodes
void LZ78Tries64::countFictNodRevTrie(RevNode *node){
	RevNode *child;
	uint i;

	if (node->fict){
		// if I and my unique child are fictitious nodes --> count as extra fictitious node
		if (node->nChildren == 1 ){
			if (node->fChild->fict)
				nExtraFict++;
			else
				nFictU++;
		}else
			nFictU++;
	}

	if (node->nChildren){
		child = node->fChild;
		for(i=0; i<node->nChildren; i++){
			if(child == NULL){
				cout << "EEEEEEEEEEEEEEERRRRRRRRRRRRR child == NULL" << endl;
				cout << "i = " << i << ", node->nChildren = " << node->nChildren << ", id = " << node->idNode << endl;
				exit(1);
			}
			countFictNodRevTrie(child);
			child = child->nextSib;
		}
	}
}

//	create 	PRev, LbRev and LbRevF structures
void LZ78Tries64::genDFUDSSeqRevTrie(RevNode *node, ulong* PRev, ulong *pos, uchar *LbRev, ulong *pLbRev){
	RevNode *child;
	uint i;

	// set open parentheses
	for(i=0, child = node->fChild; i<node->nChildren; i++, (*pos)++, child=child->nextSib, (*pLbRev)++){
		LbRev[*pLbRev] = child->symbol;
		setBit64(PRev, *pos);
		//cout << "(";
	}

	// set close parenthesis
	cleanBit64(PRev, *pos);
	(*pos)++;
	//cout << ")";

	// recursive call for all children
	child = node->fChild;
	for(i=0; i<node->nChildren; i++){
		genDFUDSSeqRevTrie(child, PRev, pos, LbRev, pLbRev);
		child = child->nextSib;
	}
}

// return the preorder value of the most right leaf of this subtree in a LZTrie
ulong LZ78Tries64::getLastLZLeaf(LZNode *node){
	while(node->nChildren){
		node = node->fChild;
		while(node->nextSib)
			node = node->nextSib;
	}

	return node->preorder;
}

// it returns the length of the subtree rooted at node
void LZ78Tries64::getLengthSubTree(RevNode *nod, ulong *len){
	(*len) += nod->nChildren;

	RevNode* p = nod->fChild;
	for(uint i=0; i<nod->nChildren; i++){
		getLengthSubTree(p, len);
		p = p->nextSib;
	}
}

//////////////////////////////////////////////////////////////////////////////
void LZ78Tries64::listNodeLZ(LZNode* nod, ulong* preorder, uint deph){
	LZNode* p = nod->fChild;
	uint i;
	uchar ss = nod->symbol;
	if(ss == cutDoc) ss = '$';
	if (*preorder){
		for (i=0; i<deph;i++) cout << "-";
		cout << "pre " << *preorder << ", '" << ss << "' Id: " << nod->idNode << ", children " << nod->nChildren << endl;
	}else
		cout << "pre 0, '#' Id: " << nod->idNode << ", children " << nod->nChildren << endl;
	(*preorder)++;
	for(i=0; i<nod->nChildren; i++){
		listNodeLZ(p, preorder, deph+1);
		p = p->nextSib;
	}
}

void LZ78Tries64::listNodeRev(RevNode* nod, ulong* preorder, uint deph){
	RevNode* p = nod->fChild;
	uint i;
	uchar ss = nod->symbol;
	if(ss == cutDoc) ss = '$';
	if (*preorder){
		for (i=0; i<deph;i++) cout << "-";
		cout << "pre " << *preorder << ", '" << ss << "' Id: " << nod->idNode << ", children " << nod->nChildren << ", fict " << nod->fict << endl;
	}else
		cout << "pre 0, '#' Id: " << nod->idNode << ", children " << nod->nChildren << ", fict " << nod->fict << endl;
	(*preorder)++;
	for(i=0; i<nod->nChildren; i++){
		listNodeRev(p, preorder, deph+1);
		p = p->nextSib;
	}
}

// list LZTrie in preorder way
void LZ78Tries64::listLZTrie(){
	ulong preorder = 0;
	listNodeLZ(this->lzTrie, &preorder, 0);
}

// list RevTrie in preorder way
void LZ78Tries64::listRevTrie(){
	ulong preorder = 0;
	listNodeRev(this->revTrie, &preorder, 0);
}

void LZ78Tries64::destroyLZNode(LZNode* nod){
	LZNode* p = nod->fChild;
	LZNode* q;
	uint i;

	for(i=0; i<nod->nChildren; i++){
		q = p->nextSib;
		destroyLZNode(p);
		p = q;
	}

	p = nod->fChild;
	for(i=0; i<nod->nChildren; i++){
		q = p->nextSib;
		delete p;
		p = q;
	}
}

void LZ78Tries64::destroyRevNode(RevNode* nod){
	RevNode* p = nod->fChild;
	RevNode* q;
	uint i;

	for(i=0; i<nod->nChildren; i++){
		q = p->nextSib;
		destroyRevNode(p);
		p = q;
	}

	p = nod->fChild;
	for(i=0; i<nod->nChildren; i++){
		q = p->nextSib;
		delete p;
		p = q;
	}
}

LZ78Tries64::~LZ78Tries64() {
	destroyLZNode(this->lzTrie);
	destroyRevNode(this->revTrie);

	delete this->lzTrie;
	delete this->revTrie;

	delete [] IdDocPreLZ;
}

void LZ78Tries64::destroyNodeTrie(LZNode* nod){
	LZNode* p = nod->fChild;
	LZNode* q = NULL;

	// delete all children of p
	for(uint i=0; i<nod->nChildren; i++){
		q = p->nextSib;
		destroyNodeTrie(p);
		p = q;
	}

	// now destroy nod
	delete nod;
}


