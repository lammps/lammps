/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: PoemsChain.h                                            *
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

#ifndef POEMSCHAIN_H_
#define POEMSCHAIN_H_

#include "poemslist.h"

struct ChildRingData {
	List<int> * childRing;
	int entranceNodeId;
};

struct POEMSChain{
	~POEMSChain(){
		for(int i = 0; i < childChains.GetNumElements(); i++)
		{
			delete childChains(i);
		}
        listOfNodes.DeleteValues();
	}
	//void printTreeStructure(int tabs);
	//void getTreeAsList(List<int> * temp);
	List<int> listOfNodes;
	List<POEMSChain> childChains;
	POEMSChain * parentChain;
	List<ChildRingData> childRings;
	
	
	void printTreeStructure(int tabs){
		for(int i = 0; i < tabs; i++)
		{
			cout << "\t";
		}
		cout << "Chain: ";
		for(int i = 0; i < listOfNodes.GetNumElements(); i++)
		{
			cout << *(listOfNodes(i)) << " ";
		}
		cout << endl;
		for(int i = 0; i < childChains.GetNumElements(); i++)
		{
			childChains(i)->printTreeStructure(tabs + 1);
		}
	}
	void getTreeAsList(List<int> * temp)
	{
		for(int i = 0; i < listOfNodes.GetNumElements(); i++)
		{
			int * integer = new int;
			*integer = *(listOfNodes(i));
			temp->Append(integer);
		}
		for(int i = 0; i < childChains.GetNumElements(); i++)
		{
			childChains(i)->getTreeAsList(temp);
		}
	}
};
#endif
