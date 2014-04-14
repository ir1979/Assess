#include <ANN/checkActiveClass1.h>

checkActiveClass1::checkActiveClass1(void)
{
	Inactive = 0;
}


checkActiveClass1::checkActiveClass1(bool *InActive)
{
	this->Inactive = InActive;
}


bool checkActiveClass1::isActiveFast(int testNodeIdx) 
{
	if (Inactive[testNodeIdx]) {
		return false;
	}
	return true;
}


checkActiveClass1::~checkActiveClass1(void)
{
	this->Inactive = 0;
}
