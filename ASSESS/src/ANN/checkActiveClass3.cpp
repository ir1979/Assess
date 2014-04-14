#include <checkActiveClass3.h>
#include <iostream>

checkActiveClass3::checkActiveClass3(void)
{
	minParents = 0;
	Inactive = 0;
	int curNode = -1;
	lastNodeIdx = -1;
}


checkActiveClass3::checkActiveClass3(DisjointSets *minParents, bool *InActive, int curNode, double lastDist, int lastNodeIdx)
{
	this->minParents = minParents;
	this->Inactive = InActive;
	this->curNode = curNode; // can be removed
	this->lastDist = lastDist;
	this->lastNodeIdx = lastNodeIdx;
}


bool checkActiveClass3::isActiveFast(int testNodeIdx) 
{
	//return (!(minParents->FindSet(curNode)==minParents->FindSet(testNodeIdx) || Inactive[testNodeIdx])); 
	return minParents->FindSet(curNode)!=minParents->FindSet(testNodeIdx); 
}

/*bool checkActiveClass3::isActiveFast(int testNodeIdx) 
{
	int redirectedNodeIdx = redirectIdx[testNodeIdx];
	if (!lastNodeSeen && redirectedNodeIdx == lastNodeIdx) 
	{
		//std::cout << " " << redirectedNodeIdx << " is skipped!\n";
		lastNodeSeen = true;
		return false;
	}
	if (Inactive[redirectedNodeIdx] || minParents->FindSet(curNode)==minParents->FindSet(redirectedNodeIdx)) {
		
		return false;
	}
	return true;
}
*/


/*
bool checkActiveClass3::isActiveLazy(int testNodeIdx, double testDist) 
{
	if (testDist>lastDist || (testDist==lastDist && lastNodeSeen))
		return true;
	
	//std::cout << " " << redirectIdx[testNodeIdx] << " is skipped based on dist: " << testDist;
	return false;
}
*/

bool checkActiveClass3::isActiveLazy(double testDist) 
{
	return (testDist>=lastDist);
}


/*
bool checkActiveClass3::isActiveLazy(double testDist) 
{
	if (testDist>lastDist || (testDist==lastDist && lastNodeSeen))
		return true;
	
	//std::cout << " " << redirectIdx[testNodeIdx] << " is skipped based on dist: " << testDist;
	return false;
}
*/


checkActiveClass3::~checkActiveClass3(void)
{
	/*
	minParents = 0;
	Inactive = 0;
	redirectIdx = 0;
	curNode = -1;
	lastDist = -1;
	lastNodeIdx = -1;
	*/
}
