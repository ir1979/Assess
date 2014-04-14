#ifndef __checkActiveClass3__
#define __checkActiveClass3__

#include <DisjointSets.h>

class checkActiveClass3
{
public:
	DisjointSets *minParents;
	bool *Inactive;
	int *redirectIdx;
	int curNode; // cen be removed
	double lastDist;
	int lastNodeIdx;
	bool lastNodeSeen;
	checkActiveClass3(void);
	checkActiveClass3(DisjointSets *minParents, bool *InActive, int curNode = 0, double lastDist = 0, int lastNodeIdx = 0);
	~checkActiveClass3(void);
	
	bool isActiveFast(int testNode);
	bool isActiveLazy(double testDist);
	//bool isActiveLazy(int testNode, double testDist);
};

#endif