/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: poemsnodelib.h                                          *
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

#ifndef NODELIB_H
#define NODELIB_H

#include <iostream>

using namespace std;


TreeNode *GetTreeNode(int item,TreeNode *lptr = nullptr,TreeNode *rptr =nullptr);

void FreeTreeNode(TreeNode *p);

void Postorder(TreeNode *t, void visit(TreeNode* &t));

void Preorder (TreeNode *t, void visit(TreeNode* &t));

void CountLeaf (TreeNode *t, int& count);

int Depth (TreeNode *t);

void IndentBlanks(int num);

void PrintTree (TreeNode *t, int level);


// ---------------- Global functions-----------------//

// postorder recursive scan of the nodes in a tree
void Postorder (TreeNode *t, void visit(TreeNode* &t))
{
        // the recursive scan terminates on a empty subtree
        if (t != nullptr)
        {
                Postorder(t->Left(), visit);    // descend left
                Postorder(t->Right(), visit);   // descend right
                visit(t);                                       // visit the node
        }
}


// preorder recursive scan of the nodes in a tree
void Preorder (TreeNode *t, void visit(TreeNode* &t))
{
        // the recursive scan terminates on a empty subtree
        if (t != nullptr)
        {
                visit(t);                               // visit the node
                Preorder(t->Left(), visit);             // descend left
                Preorder(t->Right(), visit);    // descend right
        }
}


//create TreeNode object with pointer fields lptr and rptr
// The pointers have default value nullptr
TreeNode *GetTreeNode(int item,TreeNode *lptr,TreeNode *rptr)
{
        TreeNode *p;

        // call new to allocate the new node
        // pass parameters lptr and rptr to the function
        p = new TreeNode(item, lptr, rptr);

        // if insufficient memory, terminatewith an error message
        if (p == nullptr)
        {
                cerr << "Memory allocation failure!\n";
                exit(1);
        }

        // return the pointer to the system generated memory
        return p;
}


// deallocate dynamic memory associated with the node

void FreeTreeNode(TreeNode *p)
{
        delete p;
}


// the function uses the postorder scan. a visit
// tests whether the node is a leaf node
void CountLeaf (TreeNode *t, int& count)
{
        //use postorder descent
        if(t !=nullptr)
        {
                CountLeaf(t->Left(), count); // descend left
                CountLeaf(t->Right(), count); // descend right

                // check if node t is a leaf node (no descendants)
                // if so, increment the variable count
                if (t->Left() == nullptr && t->Right() == nullptr)
                        count++;
        }
}


// the function uses the postorder scan. it computes the
// depth of the left and right subtrees of a node and
// returns the depth as 1 + max(depthLeft,depthRight).
// the depth of an empty tree is -1
int Depth (TreeNode *t)
{
        int depthLeft, depthRight, depthval;

        if (t == nullptr)
                depthval = -1;
        else
        {
                depthLeft = Depth(t->Left());
                depthRight = Depth(t->Right());
                depthval = 1+(depthLeft > depthRight?depthLeft:depthRight);
        }
        return depthval;
}

void IndentBlanks(int num)
{
//      const int indentblock = 6;

        for(int i = 0; i < num; i++)
                cout << " ";
}

void PrintTree (TreeNode *t, int level)
{
        //print tree with root t, as long as t!=nullptr
        if (t != nullptr)
        {
                int indentUnit = 5;
                // print right branch of tree t
                PrintTree(t->Right(),level + 1);
                // indent to current level; output node data
                IndentBlanks(indentUnit*level);
                cout << t->GetData() << endl;
                // print left branch of tree t
                PrintTree(t->Left(),level + 1);
        }
}
#endif

