/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: SystemProcessor.h                                       *
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

#ifndef _SYS_PROCESSOR_H_
#define _SYS_PROCESSOR_H_
#include "poemslist.h"
#include "poemstree.h"
#include "POEMSChain.h"


struct POEMSNode {
	List<POEMSNode> links;
	List<bool> taken;
	int idNumber;
	bool visited;
	
	~POEMSNode(){
		for(int i = 0; i < taken.GetNumElements(); i++)
		{
			delete taken(i);
		}
	};	
};


class SystemProcessor{
private:
	Tree nodes;
	static void POEMSNodeDelete_cb(void *node) {
		delete (POEMSNode *) node;
	}
	List<POEMSChain> headsOfSystems;
	List<List<int> > ringsInSystem;
	POEMSNode * findSingleLink(TreeNode * aNode);
	POEMSChain * AddNewChain(POEMSNode * currentNode);
	bool setLinkVisited(POEMSNode * firstNode, POEMSNode * secondNode);
public:
	SystemProcessor(void);
	
	~SystemProcessor(void)	{
		headsOfSystems.DeleteValues();
		for(int i = 0; i < ringsInSystem.GetNumElements(); i++)
		{
			for(int k = 0; k < ringsInSystem(i)->GetNumElements(); i++)
			{
				delete (*ringsInSystem(i))(k);
			}
		}
	};
	void processArray(int** links, int numLinks);
	List<POEMSChain> * getSystemData();
	int getNumberOfHeadChains();
};

SystemProcessor::SystemProcessor(void){
  // register callback for deleting auxiliary data from tree nodes.
  nodes.SetDeleteAuxData(&POEMSNodeDelete_cb);
}

void SystemProcessor::processArray(int** links, int numLinks)
{
	bool * false_var;					//holds the value false; needed because a constant cannot be put into a list; the list requires a
										//reference.
	for(int i = 0; i < numLinks; i++)	//go through all the links in the input array
	{
		if(!nodes.Find(links[i][0]))	//if the first node in the pair is not found in the storage tree
		{
			POEMSNode * newNode = new POEMSNode;	//make a new node
//			forDeletion.Append(newNode);
			newNode->idNumber = links[i][0];		//set its ID to the value
			newNode->visited = false;				//set it to be unvisited
			nodes.Insert(links[i][0], links[i][0], (void *) newNode);	//and add it to the tree storage structure
		}
		if(!nodes.Find(links[i][1]))				//repeat process for the other half of each link
		{
			POEMSNode * newNode = new POEMSNode;
//			forDeletion.Append(newNode);
			newNode->idNumber = links[i][1];
			newNode->visited = false;
			nodes.Insert(links[i][1], links[i][1], (void *) newNode);
		}
		POEMSNode * firstNode = (POEMSNode *)nodes.Find(links[i][0]);	//now that we are sure both nodes exist,
		POEMSNode * secondNode = (POEMSNode *)nodes.Find(links[i][1]);	//we can get both of them out of the tree
		firstNode->links.Append(secondNode);							//and add the link from the first to the second...
		false_var = new bool;
		*false_var = false;												//make a new false boolean to note that the link between these two
		firstNode->taken.Append(false_var);								//has not already been taken, and append it to the taken list
		secondNode->links.Append(firstNode);							//repeat process for link from second node to first
		false_var = new bool;
		*false_var = false;
		secondNode->taken.Append(false_var);
	}	
	
	TreeNode * temp = nodes.GetRoot();							//get the root node of the node storage tree
	POEMSNode * currentNode;
	do
	{
		currentNode = findSingleLink(temp);						//find the start of the next available chain
		if(currentNode != NULL)
		{
			headsOfSystems.Append(AddNewChain(currentNode));							//and add it to the headsOfSystems list of chains
		}
	} 
	while(currentNode != NULL);									//repeat this until all chains have been added		
}

POEMSChain * SystemProcessor::AddNewChain(POEMSNode * currentNode){
	if(currentNode == NULL)	//Termination condition; if the currentNode is null, then return null
	{
		return NULL;
	}
	int * tmp;
	POEMSNode * nextNode = NULL;	//nextNode stores the proposed next node to add to the chain.  this will be checked to make sure no backtracking is occuring before being assigned as the current node.
	POEMSChain * newChain = new POEMSChain;	//make a new POEMSChain object.  This will be the object returned

	if(currentNode->links.GetNumElements() == 0)	//if we have no links from this node, then the whole chain is only one node.  Add this node to the chain and return it; mark node as visited for future reference
	{
		currentNode->visited = true;
		tmp = new int;
		*tmp = currentNode->idNumber;
		newChain->listOfNodes.Append(tmp);
		return newChain;
	}
	while(currentNode->links.GetNumElements() <= 2)	//we go until we get to a node that branches, or both branches have already been taken both branches can already be taken if a loop with no spurs is found in the input data
	{
		currentNode->visited = true;
		tmp = new int;
		*tmp = currentNode->idNumber;
		newChain->listOfNodes.Append(tmp);	//append the current node to the chain & mark as visited
		//cout << "Appending node " << currentNode->idNumber << " to chain" << endl;
		nextNode = currentNode->links.GetHeadElement()->value;	//the next node is the first or second value stored in the links array
																//of the current node.  We get the first value...
		if(!setLinkVisited(currentNode, nextNode))					//...and see if it points back to where we came from. If it does...
		{														//either way, we set this link as visited
			if(currentNode->links.GetNumElements() == 1)		//if it does, then if that is the only link to this node, we're done with the chain, so append the chain to the list and return the newly created chain
			{
//				headsOfSystems.Append(newChain);
				return newChain;
			}
			nextNode = currentNode->links.GetHeadElement()->next->value;//follow the other link if there is one, so we go down the chain
			if(!setLinkVisited(currentNode, nextNode))				//mark link as followed, so we know not to backtrack			
			{
	//			headsOfSystems.Append(newChain);
				return newChain;								//This condition, where no branches have occurred but both links have already
																//been taken can only occur in a loop with no spurs; add this loop to the
																//system (currently added as a chain for consistency), and return.
			}
		}
		currentNode = nextNode;									//set the current node to be the next node in the chain
	}
	currentNode->visited = true;
	tmp = new int;
	*tmp = currentNode->idNumber;
	newChain->listOfNodes.Append(tmp);		//append the last node before branch (node shared jointly with branch chains)
																//re-mark as visited, just to make sure
	ListElement<POEMSNode> * tempNode = currentNode->links.GetHeadElement();	//go through all of the links, one at a time that branch
	POEMSChain * tempChain = NULL;								//temporary variable to hold data 
	while(tempNode != NULL)										//when we have followed all links, stop
	{
		if(setLinkVisited(tempNode->value, currentNode))		//dont backtrack, or create closed loops
		{
			tempChain = AddNewChain(tempNode->value);			//Add a new chain created out of the next node down that link
			tempChain->parentChain = newChain;					//set the parent to be this chain
			newChain->childChains.Append(tempChain);			//append the chain to this chain's list of child chains
		}
		tempNode = tempNode->next;								//go to process the next chain
	}
	//headsOfSystems.Append(newChain);							//append this chain to the system list		
	return newChain;	
}

POEMSNode * SystemProcessor::findSingleLink(TreeNode * aNode)
//This function takes the root of a search tree containing POEMSNodes and returns a POEMSNode corresponding to the start of a chain in the
//system.  It finds a node that has not been visited before, and only has one link; this node will be used as the head of the chain.
{
	if(aNode == NULL)	
	{
		return NULL;
	}
	POEMSNode * returnVal =  (POEMSNode *)aNode->GetAuxData();	//get the poemsnode data out of the treenode
	POEMSNode * detectLoneLoops = NULL;							//is used to handle a loop that has no protruding chains
	if(returnVal->visited == false)
	{
		detectLoneLoops = returnVal;							//if we find any node that has not been visited yet, save it
	}
	if(returnVal->links.GetNumElements() == 1 && returnVal->visited == false)	//see if it has one element and hasnt been visited already
	{
		return returnVal;														//return the node is it meets this criteria
	}
	returnVal = findSingleLink(aNode->Left());									//otherwise, check the left subtree
	if(returnVal == NULL)														//and if we find nothing...
	{
		returnVal = findSingleLink(aNode->Right());								//check the right subtree
	}
	if(returnVal == NULL)														//if we could not find any chains
	{
		returnVal = detectLoneLoops;											//see if we found any nodes at all that havent been processed
	}
	return returnVal;															//return what we find (will be NULL if no new chains are 
																				//found)
}

bool SystemProcessor::setLinkVisited(POEMSNode * firstNode, POEMSNode * secondNode)
//setLinkVisited sets the links between these two nodes as visited. If they are already visited, it returns false.  Otherwise, it sets
//them as visited and returns true.  This function is used to see whether a certain path has been taken already in the graph structure.
//If it has been, then we need to know so we dont follow it again; this prevents infinite recursion when there is a loop, and prevents
//backtracking up a chain that has already been made.  The list of booleans denoting if a link has been visited is called 'taken' and is
//part of the POEMSNode struct.  The list is the same size as the list of pointers to other nodes, and stores the boolean visited/unvisited
//value for that particular link.  Because each link is represented twice, (once at each node in the link), both of the boolean values need
//to be set in the event that the link has to be set as visited.
{
	//cout << "Checking link between nodes " << firstNode->idNumber << " and " << secondNode->idNumber << "... ";
	ListElement<POEMSNode> * tmp = firstNode->links.GetHeadElement();	//get the head element of the list of pointers for node 1
	ListElement<bool> * tmp2 = firstNode->taken.GetHeadElement();		//get the head element of the list of bool isVisited flags for node 1
	while(tmp->value != NULL || tmp2->value != NULL)					//go through untill we reach the end of the lists
	{
		if(tmp->value == secondNode)							//if we find the link to the other node
		{
			if(*(tmp2->value) == true)							//if the link has already been visited
			{
				 //cout << "visited already" << endl;
				return false;									//return false to indicate that the link has been visited before this attempt
			}
			else												//otherwise, visit it
			{
				*tmp2->value = true;
			}
			break;
		}
		tmp = tmp->next;										//go check next link
		tmp2 = tmp2->next;
	}

	tmp = secondNode->links.GetHeadElement();			//now, if the link was unvisited, we need to go set the other node's list such that
														//it also knows this link is being visited
	tmp2 = secondNode->taken.GetHeadElement();
	while(tmp->value != NULL || tmp2->value != NULL)	//go through the list
	{
		if(tmp->value == firstNode)						//if we find the link
		{
			if(*(tmp2->value) == true)					//and it has already been visited, then signal an error; this shouldnt ever happen
			{
				cout << "Error in parsing structure! Should never reach this condition! \n" << 
						"Record of visited links out of synch between two adjacent nodes.\n";
				return false;
			}
			else
			{
				*tmp2->value = true;					//set the appropriate value to true to indicate this link has been visited
			}
			break;
		}
		tmp = tmp->next;
		tmp2 = tmp2->next;
	}
	//cout << "not visited" << endl;
	return true;										//return true to indicate that this is the first time the link has been visited
}

List<POEMSChain> * SystemProcessor::getSystemData(void)	//Gets the list of POEMSChains that comprise the system.  Might eventually only
														//return chains linked to the reference plane, but currently returns every chain
														//in the system.
{
	return &headsOfSystems;
}

int SystemProcessor::getNumberOfHeadChains(void) //This function isnt implemented yet, and might be taken out entirely; this was a holdover
												//from when I intended to return an array of chain pointers, rather than a list of chains
												//It will probably be deleted once I finish figuring out exactly what needs to be returned
{
	return 0;
}
#endif
