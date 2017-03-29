#pragma once
#include "BWT.h"
class BWT2 :
	public BWT
{
public:
	BWT2();
	~BWT2();
	void makebwts2();
	void search2(string sub);
	int Occ2(int r, string c);
	int LFC2(int r, string c);
	int getC2(string c);
	string toNext2(string c);
	void run();
private:
	vector<string> BWTS2;
	vector<string> Index2;
};

