#ifndef __checkActiveClass2__
#define __checkActiveClass2__

#include <DisjointSets.h>

class checkActiveClass2
{
public:
	DisjointSets *minParents;
	bool *Inactive;
	int *redirectIdx;
	int curNode; // cen be removed
	double lastDist;
	int lastNodeIdx;
	bool lastNodeSeen;
	checkActiveClass2(void);
	checkActiveClass2(DisjointSets *minParents, bool *InActive, int *redirectIdx, int curNode = 0, double lastDist = 0, int lastNodeIdx = 0);
	~checkActiveClass2(void);
	
	bool isActiveFast(int testNode);
	bool isActiveLazy(double testDist);
	bool isActiveLazy(int testNodeIdx, double testDist);
	//bool isActiveLazy(int testNode, double testDist);
};

#endif