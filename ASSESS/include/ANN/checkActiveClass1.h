#ifndef __checkActiveClass1__
#define __checkActiveClass1__

class checkActiveClass1
{
public:
	bool *Inactive;
	checkActiveClass1(void);
	checkActiveClass1(bool *InActive);
	~checkActiveClass1(void);
	
	bool isActiveFast(int testNode);
};

#endif