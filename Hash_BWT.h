#pragma once
#include "BWT.h"
class Hash_BWT :
	public BWT
{
public:
	Hash_BWT();
	~Hash_BWT();
	void makeHash(int k);
	//void Hash_search(string c, float e);
	void Hash_search(string sub, int k, float e);
	vector<vector<int>> Check(vector<vector<int>> &sub_res, int length);
	string IntToS(int i, int k);
	bool HasPre(vector<int> s, int num);
	int Hamming(string c1, string c2);
	void run();
	bool indic(string s);
private:
	map <string, vector<int>> HashIndex;
	vector<string> dic;
	vector<string> nt;
};

