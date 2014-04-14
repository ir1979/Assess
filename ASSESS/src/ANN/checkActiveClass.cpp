#include <ANN/checkActiveClass.h>

checkActiveClass::checkActiveClass(void)
{
	minParents = 0;
	Inactive = 0;
	redirectIdx = 0;
	int curNode = -1;
	double minDist = -1;

}


checkActiveClass::checkActiveClass(int *minParents, bool *InActive, int *redirectIdx, int curNode, double minDist)
{
	this->minParents = minParents;
	this->Inactive = InActive;
	this->redirectIdx = redirectIdx;
	this->curNode = curNode;
	this->minDist = minDist;
}


bool checkActiveClass::isActive(int testNodeIdx, double testDist) 
{
	int redirectedNodeIdx = redirectIdx[testNodeIdx];
	if (minParents[curNode]==minParents[redirectedNodeIdx] || Inactive[redirectedNodeIdx] || redirectedNodeIdx == curNode || (minDist>=0 && testDist<=minDist))
		return false;
	return true;
}


checkActiveClass::~checkActiveClass(void)
{
	minParents = 0;
	Inactive = 0;
	redirectIdx = 0;
	curNode = -1;
	minDist = -1;
}
