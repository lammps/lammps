/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: poemstreenode.h                                         *
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

#ifndef TREENODE_H
#define TREENODE_H

//#define NULL 0


//Tree depends on TreeNode
class Tree;

// declares a tree node object for a binary tree
class TreeNode{

private:
// points to the left and right children of the node
        TreeNode *left;
        TreeNode *right;

        int balanceFactor;
        int data;
        void * aux_data;
public:
        // make Tree a friend because it needs access to left and right pointer fields of a node
        friend class Tree;
        TreeNode * Left();
        TreeNode * Right();
        int GetData();
        void * GetAuxData() {return aux_data;};
        void SetAuxData(void * AuxData) {aux_data = AuxData;};
        int GetBalanceFactor();
        TreeNode(const int &item, TreeNode *lptr, TreeNode *rptr, int balfac = 0);
        //friend class DCASolver;
};

#endif

