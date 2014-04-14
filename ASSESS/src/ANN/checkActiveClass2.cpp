#include <checkActiveClass2.h>
#include <iostream>

checkActiveClass2::checkActiveClass2(void)
{
	minParents = 0;
	Inactive = 0;
	redirectIdx = 0;
	int curNode = -1;
	lastNodeIdx = -1;
}


checkActiveClass2::checkActiveClass2(DisjointSets *minParents, bool *InActive, int *redirectIdx, int curNode, double lastDist, int lastNodeIdx)
{
	this->minParents = minParents;
	this->Inactive = InActive;
	this->redirectIdx = redirectIdx;
	this->curNode = curNode; // can be removed
	this->lastDist = lastDist;
	this->lastNodeIdx = lastNodeIdx;
}


/*bool checkActiveClass2::isActiveFast(int testNodeIdx) 
{
	//return (!(minParents->FindSet(curNode)==minParents->FindSet(redirectIdx[testNodeIdx]) || Inactive[redirectIdx[testNodeIdx]])); 
	return (minParents->FindSet(curNode)!=minParents->FindSet(redirectIdx[testNodeIdx]));// || Inactive[redirectIdx[testNodeIdx]])); 
}*/

bool checkActiveClass2::isActiveFast(int testNodeIdx) 
{
	int redirectedNodeIdx = redirectIdx[testNodeIdx];
	if (!lastNodeSeen && redirectedNodeIdx == lastNodeIdx) 
	{
		//std::cout << " " << redirectedNodeIdx << " is skipped!\n";
		lastNodeSeen = true;
		return false;
	}
	if (Inactive[redirectedNodeIdx] || minParents->FindSet(curNode)==minParents->FindSet(redirectedNodeIdx)) {
		/*
		if (curNode==3) {
		if (minParents->FindSet(curNode)==minParents->FindSet(redirectedNodeIdx))
		std::cout << " " << redirectedNodeIdx << " is skipped because of minParent\n";
		else if (Inactive[redirectedNodeIdx])
		std::cout << " " << redirectedNodeIdx << " is skipped because of inActive\n";
		}
		*/
		return false;
	}
	return true;
}



bool checkActiveClass2::isActiveLazy(int testNodeIdx, double testDist) 
{
	if (testDist>lastDist || (testDist==lastDist && lastNodeSeen))
		return true;

	//std::cout << " " << redirectIdx[testNodeIdx] << " is skipped based on dist: " << testDist;
	return false;
}


bool checkActiveClass2::isActiveLazy(double testDist) 
{
	//	return (testDist>=lastDist);


	if (testDist>lastDist || (testDist==lastDist && lastNodeSeen))
		return true;

	//std::cout << " " << redirectIdx[testNodeIdx] << " is skipped based on dist: " << testDist;
	return false;

}

/*
bool checkActiveClass2::isActiveLazy(double testDist) 
{
if (testDist>lastDist || (testDist==lastDist && lastNodeSeen))
return true;

//std::cout << " " << redirectIdx[testNodeIdx] << " is skipped based on dist: " << testDist;
return false;
}
*/


checkActiveClass2::~checkActiveClass2(void)
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
