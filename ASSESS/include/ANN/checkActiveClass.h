#pragma once

class checkActiveClass
{
public:
	int *minParents;
	bool *Inactive;
	int *redirectIdx;
	int curNode;
	double minDist;
	checkActiveClass(void);
	checkActiveClass(int *minParents, bool *InActive, int *redirectIdx, int curNode = -1, double minDist=-1);
	~checkActiveClass(void);
	bool isActive(int testNode, double testDist);
};

