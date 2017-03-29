#pragma once

class BWT
{
public:
	BWT();
	~BWT();
	string Read_Reference(string filename);
	vector<string> Read_Subs(string filename);
	void preprocess();
	vector<int> search(string sub);
	void unexactsearch(string sub, float e);
	int Occ(int r, char c);
	int LFC(int r, char c);
	int getC(char c);
	int editDis(string c1, string c2);
	void run();
	string T;
	vector<int> SA;
	unordered_map<char, int> C;
	unordered_map<char, char> toNext;
	vector<string> Matrix;
private:
	string BWTS;
	//vector<string> Index2;
};