/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: poemstreenode.cpp                                       *
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

#include "poemstreenode.h"

// constructor; initialize the data and pointer fields
// The pointer value NULL assigns a empty subtree
TreeNode::TreeNode (const int & item, TreeNode *lptr,TreeNode *rptr,
					int balfac):data(item), left(lptr), right(rptr), balanceFactor(balfac)
{
}



// return left
TreeNode* TreeNode::Left(void)
{
	return left;
}

// return right
TreeNode* TreeNode::Right(void)
{
	return right;
}

int TreeNode::GetBalanceFactor(void)
{
	return balanceFactor;
}

int TreeNode::GetData(void)
{
	return data;
}
