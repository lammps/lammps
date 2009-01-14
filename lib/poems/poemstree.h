/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: poemstree.h                                             *
 *      AUTHORS: See Author List                                           * 
 *      GRANTS: See Grants List                                            *
 *      COPYRIGHT: (C) 2005 by Authors as listed in Author's List          *
 *      LICENSE: Please see License Agreement                              *
 *      DOWNLOAD: Free at www.rpi.edu/~anderk5                             *
 *      ADMINISTRATOR: Prof. Kurt Anderson                                 *
 *                     Computational Dynamics Lab                          *
 *                     Rensselaer Polytechnic Institute                    *
 *                     110 8th St. Troy NY 12180                           * 
 *      CONTACT:        anderk5@rpi.edu                                    *
 *_________________________________________________________________________*/

#ifndef TREE_H
#define TREE_H

#include "poemstreenode.h"
#include "poemsnodelib.h"


// constants to indicate the balance factor of a node
const int leftheavy = -1;
const int balanced = 0;
const int rightheavy = 1;



class Tree{
protected:
	// pointer to tree root and node most recently accessed
	TreeNode *root;
	TreeNode *current;

	// number of elements in the tree
	int size;
	// used by the copy constructor and assignment operator
	TreeNode *CopyTree(TreeNode *t);

	// used by insert and delete method to re-establish
	// the avl conditions after a node is added or deleted
	// from a subtree
	void SingleRotateLeft (TreeNode* &p);
	void SingleRotateRight (TreeNode* &p);
	void DoubleRotateLeft (TreeNode* &p);
	void DoubleRotateRight (TreeNode* &p);
	void UpdateLeftTree (TreeNode* &p, int &reviseBalanceFactor);
	void UpdateRightTree (TreeNode* &p, int &reviseBalanceFactor);


	// used by destructor, assignment operator and ClearList
	void DeleteTree(TreeNode *t);
	void ClearTree(TreeNode * &t);

	// locate a node with data item and its parent in tree
	// used by Find and Delete
	TreeNode *FindNode(const int& item, TreeNode* & parent) const;

public:
	// constructor, destructor
	Tree(void);
	~Tree(void)
	{
		ClearTree(root);
	}; 

	// assignment operator
	Tree& operator= (const Tree& rhs);

	// standard list handling methods
	void * Find(int& item);
	void * GetAuxData(int item) { return (void *)(FindNode(item, root)->GetAuxData());}
	void Insert(const int& item, const int& data, void * AuxData = NULL);
	void Delete(const int& item);
	void AVLInsert(TreeNode* &tree, TreeNode* newNode, int &reviseBalanceFactor);	
	void ClearList(void);
	
	// tree specific methods
	void Update(const int& item);
	TreeNode *GetRoot(void) const;
};


// constructor
Tree::Tree(void)
{
	root = 0; 
	current = 0; 
	size = 0;
}



// return root pointer
TreeNode *Tree::GetRoot(void) const
{
	return root;
}


// assignment operator
Tree& Tree::operator = (const Tree& rhs)
{
	// can't copy a tree to itself
	if (this == &rhs)
		return *this;

	// clear current tree. copy new tree into current object
	ClearList();
	root = CopyTree(rhs.root);

	// assign current to root and set the tree size
	current = root;
	size = rhs.size;

	// return reference to current object
	return *this;
}

// search for data item in the tree. if found, return its node
// address and a pointer to its parent; otherwise, return NULL
TreeNode *Tree::FindNode(const int& item,
								   TreeNode* & parent) const
{
	// cycle t through the tree starting with root
	TreeNode *t = root;

	// the parent of the root is NULL
	parent = NULL;

	// terminate on empty subtree
	while(t != NULL)
	{
		// stop on a match
		if (item == t->data)
			break;
		else
		{
			// update the parent pointer and move right of left
			parent = t;
			if (item < t->data)
				t = t->left;
			else
				t = t->right;
		}
	}

	// return pointer to node; NULL if not found
	return t;
}

// search for item. if found, assign the node data to item
void * Tree::Find(int& item)
{
	// we use FindNode, which requires a parent parameter
	TreeNode *parent;

	// search tree for item. assign matching node to current
	current = FindNode (item, parent);

	// if item found, assign data to item and return True
	if (current != NULL)
	{
		item = current->data;
        return current->GetAuxData();
	}
	else
		// item not found in the tree. return False
		return NULL;
}


void Tree::Insert(const int& item, const int& data, void * AuxData)
{
	// declare AVL tree node pointer; using base class method
	// GetRoot. cast to larger node and assign root pointer
	TreeNode *treeRoot, *newNode;
	 treeRoot = GetRoot();

	// flag used by AVLInsert to rebalance nodes
	int reviseBalanceFactor = 0;

	// get a new AVL tree node with empty pointer fields
	newNode = GetTreeNode(item,NULL,NULL);
	newNode->data = data;
	newNode->SetAuxData(AuxData);
	// call recursive routine to actually insert the element
	AVLInsert(treeRoot, newNode, reviseBalanceFactor);

	// assign new values to data members in the base class
	root = treeRoot;
	current = newNode;
	size++;

}

void Tree::AVLInsert(TreeNode *&tree, TreeNode *newNode, int &reviseBalanceFactor)
{
	// flag indicates change node's balanceFactor will occur
	int rebalanceCurrNode;

	// scan reaches an empty tree; time to insert the new node
	if (tree == NULL)
	{
		// update the parent to point at newNode
		tree = newNode;

		// assign balanceFactor = 0 to new node
		tree->balanceFactor = balanced;
		// broadcast message; balanceFactor value is modified
		reviseBalanceFactor = 1;
	}
	// recursively move left if new data < current data
	else if (newNode->data < tree->data)
	{
		AVLInsert(tree->left,newNode,rebalanceCurrNode);
		// check if balanceFactor must be updated.
		if (rebalanceCurrNode)
		{
			// went left from node that is left heavy. will
			// violate AVL condition; use rotation (case 3)
			if (tree->balanceFactor == leftheavy)
				UpdateLeftTree(tree,reviseBalanceFactor);

			// went left from balanced node. will create
			// node left on the left. AVL condition OK (case 1)
			else if (tree->balanceFactor == balanced)
			{
				tree->balanceFactor = leftheavy;
				reviseBalanceFactor = 1;
			}
			// went left from node that is right heavy. will
			// balance the node. AVL condition OK (case 2)
			else
			{
				tree->balanceFactor = balanced;
				reviseBalanceFactor = 0;
			}
		}
        else 
            // no balancing occurs; do not ask previous nodes
			reviseBalanceFactor = 0;
	}
    // otherwise recursively move right
    else
	{
		AVLInsert(tree->right, newNode, rebalanceCurrNode);
		// check if balanceFactor must be updated.
		if (rebalanceCurrNode)
		{
			// went right from node that is left heavy. wil;
			// balance the node. AVL condition OK (case 2)
			if (tree->balanceFactor == leftheavy)
			{
				// scanning right subtree. node heavy on left.
				// the node will become balanced
				tree->balanceFactor = balanced;
				reviseBalanceFactor = 0;
			}
			// went right from balanced node. will create
			// node heavy on the right. AVL condition OK (case 1)
			else if (tree->balanceFactor == balanced)
			{
				// node is balanced; will become heavy on right
				tree->balanceFactor = rightheavy;
				reviseBalanceFactor = 1;
			}
			// went right from node that is right heavy. will
			// violate AVL condition; use rotation (case 3)
			else
				UpdateRightTree(tree, reviseBalanceFactor);
		}
		else
			reviseBalanceFactor = 0;
	}
}


void Tree::UpdateLeftTree (TreeNode* &p, int &reviseBalanceFactor)
{
	TreeNode *lc;

	lc = p->Left();			// left subtree is also heavy
	if (lc->balanceFactor == leftheavy)
	{
		SingleRotateRight(p);
		reviseBalanceFactor = 0;
	}
	// is right subtree heavy?
	else if (lc->balanceFactor == rightheavy)
	{
		// make a double rotation
		DoubleRotateRight(p);
		// root is now balance
		reviseBalanceFactor = 0;
	}
}

void Tree::UpdateRightTree (TreeNode* &p, int &reviseBalanceFactor)
{
	TreeNode *lc;

	lc = p->Right();			// right subtree is also heavy
	if (lc->balanceFactor == rightheavy)
	{
		SingleRotateLeft(p);
		reviseBalanceFactor = 0;
	}
	// is left subtree heavy?
	else if (lc->balanceFactor == leftheavy)
	{
		// make a double rotation
		DoubleRotateLeft(p);
		// root is now balance
		reviseBalanceFactor = 0;
	}
}

void Tree::SingleRotateRight (TreeNode* &p)
{
	// the left subtree of p is heavy
	TreeNode *lc;

	// assign the left subtree to lc
	lc = p->Left();

	// update the balance factor for parent and left child
	p->balanceFactor = balanced;
	lc->balanceFactor = balanced;

	// any right subtree st of lc must continue as right
	// subtree of lc. do by making it a left subtree of p
	p->left = lc->Right();

	// rotate p (larger node) into right subtree of lc
	// make lc the pivot node
	lc->right = p;
	p = lc;
}

void Tree::SingleRotateLeft (TreeNode* &p)
{
	// the right subtree of p is heavy
	TreeNode *lc;

	// assign the left subtree to lc
	lc = p->Right();

	// update the balance factor for parent and left child
	p->balanceFactor = balanced;
	lc->balanceFactor = balanced;

	// any right subtree st of lc must continue as right
	// subtree of lc. do by making it a left subtree of p
	p->right = lc->Left();

	// rotate p (larger node) into right subtree of lc
	// make lc the pivot node
	lc->left = p;
	p = lc;
}

// double rotation right about node p
void Tree::DoubleRotateRight (TreeNode* &p)
{
	// two subtrees that are rotated
	TreeNode *lc, *np;

	// in the tree, node(lc) <= node(np) < node(p)
	lc = p->Left();			// lc is left child of parent
	np = lc->Right();		// np is right child of lc
	
	// update balance factors for p, lc, and np
	if (np->balanceFactor == rightheavy)
	{
		p->balanceFactor = balanced;
		lc->balanceFactor = rightheavy;
	}
	else if (np->balanceFactor == balanced)
	{
		p->balanceFactor = balanced;
		lc->balanceFactor = balanced;
	}
	else
	{
		p->balanceFactor = rightheavy;
		lc->balanceFactor = balanced;
	}	
	np->balanceFactor = balanced;

	// before np replaces the parent p, take care of subtrees
	// detach old children and attach new children
	lc->right = np->Left();
	np->left = lc;
	p->left = np->Right();
	np->right = p;
	p = np;
}

void Tree::DoubleRotateLeft (TreeNode* &p)
{
	// two subtrees that are rotated
	TreeNode *lc, *np;

	// in the tree, node(lc) <= node(np) < node(p)
	lc = p->Right();			// lc is right child of parent
	np = lc->Left();		// np is left child of lc
	
	// update balance factors for p, lc, and np
	if (np->balanceFactor == leftheavy)
	{
		p->balanceFactor = balanced;
		lc->balanceFactor = leftheavy;
	}
	else if (np->balanceFactor == balanced)
	{
		p->balanceFactor = balanced;
		lc->balanceFactor = balanced;
	}
	else
	{
		p->balanceFactor = leftheavy;
		lc->balanceFactor = balanced;
	}	
	np->balanceFactor = balanced;

	// before np replaces the parent p, take care of subtrees
	// detach old children and attach new children
	lc->left = np->Right();
	np->right = lc;
	p->right = np->Left();
	np->left = p;
	p = np;
}

// if item is in the tree, delete it
void Tree::Delete(const int& item)
{
	// DNodePtr = pointer to node D that is deleted
	// PNodePtr = pointer to parent P of node D
	// RNodePtr = pointer to node R that replaces D
	TreeNode *DNodePtr, *PNodePtr, *RNodePtr;

	// search for a node containing data value item. obtain its
	// node adress and that of its parent
	if ((DNodePtr = FindNode (item, PNodePtr)) == NULL)
		return;

	// If D has NULL pointer, the 
	// replacement node is the one on the other branch
	if (DNodePtr->right == NULL)
		RNodePtr = DNodePtr->left;
	else if (DNodePtr->left == NULL)
		RNodePtr = DNodePtr->right;
	// Both pointers of DNodePtr are non-NULL
	else
	{
		// Find and unlink replacement node for D
		// Starting on the left branch of node D,
		// find node whose data value is the largest of all
		// nodes whose values are less than the value in D
		// Unlink the node from the tree

		// PofRNodePtr = pointer to parent of replacement node
		TreeNode *PofRNodePtr = DNodePtr;

		// frist possible replacement is left child D
		RNodePtr = DNodePtr->left;

		// descend down right subtree of the left child of D
		// keeping a record of current node and its parent.
		// when we stop, we have found the replacement
		while (RNodePtr->right != NULL)
		{
			PofRNodePtr = RNodePtr;
			RNodePtr = RNodePtr;
		}

		if (PofRNodePtr == DNodePtr)
			// left child of deleted node is the replacement
			// assign right subtree of D to R
			RNodePtr->right =  DNodePtr->right;
		else
		{
			// we moved at least one node down a right brance
			// delete replacement node from tree by assigning
			// its left branc to its parent
			PofRNodePtr->right = RNodePtr->left;

			// put replacement node in place of DNodePtr.
			RNodePtr->left = DNodePtr->left;
			RNodePtr->right = DNodePtr->right;
		}
	}

	// complete the link to the parent node
	// deleting the root node. assign new root
	if (PNodePtr == NULL)
		root = RNodePtr;
	// attach R to the correct branch of P
	else if (DNodePtr->data < PNodePtr->data)
		PNodePtr->left = RNodePtr;
	else
		PNodePtr->right = RNodePtr;

	// delete the node from memory and decrement list size
	FreeTreeNode(DNodePtr);  // this says FirstTreeNode in the book, should be a typo
	size--;
}





// if current node is defined and its data value matches item,
// assign node value to item; otherwise, insert item in tree
void Tree::Update(const int& item)
{
	if (current !=NULL && current->data == item)
		current->data = item;
	else
		Insert(item, item);
}

// create duplicate of tree t; return the new root
TreeNode *Tree::CopyTree(TreeNode *t)
{
	// variable newnode points at each new node that is
	// created by a call to GetTreeNode and later attached to
	// the new tree. newlptr and newrptr point to the child of
	// newnode and are passed as parameters to GetTreeNode
	TreeNode *newlptr, *newrptr, *newnode;

	// stop the recursive scan when we arrive at an empty tree
	if (t == NULL)
		return NULL;

	// CopyTree builds a new tree by scanning the nodes of t.
	// At each node in t, CopyTree checks for a left child. if
	// present it makes a copy of left child or returns NULL.
	// the algorithm similarly checks for a right child.
	// CopyTree builds a copy of node using GetTreeNode and
	// appends copy of left and right children to node.

	if (t->Left() !=NULL)
		newlptr = CopyTree(t->Left());
	else
		newlptr = NULL;

	if (t->Right() !=NULL)
		newrptr = CopyTree(t->Right());
	else
		newrptr = NULL;


	// Build new tree from the bottom up by building the two
	// children and then building the parent
	newnode = GetTreeNode(t->data, newlptr, newrptr);

	// return a pointer to the newly created node
	return newnode;
}


// us the postorder scanning algorithm to traverse the nodes in
// the tree and delete each node as the vist operation
void Tree::DeleteTree(TreeNode *t)
{
	if (t != NULL)
	{
		DeleteTree(t->Left());
		DeleteTree(t->Right());
		if (t->GetAuxData() != NULL)
                    delete (TreeNode *) t->GetAuxData();
		FreeTreeNode(t);
	}
}

// call the function DeleteTree to deallocate the nodes. then
// set the root pointer back to NULL
void Tree::ClearTree(TreeNode * &t)
{
	DeleteTree(t);
	t = NULL;		// root now NULL
}

// delete all nodes in list
void Tree::ClearList(void)
{
	delete root;
	delete current;
	size = 0;
}

#endif
