#pragma once

class BWT
{
public:
	BWT();
	~BWT();
	void Read_Reference(string filename);
	vector<string> Read_Subs(string filename);
	void preprocess();
	void makebwts2();
	void search(string sub);
	void search2(string sub);
	int Occ(int r, char c);
	int Occ2(int r, string c);
	int LFC(int r, char c);
	int LFC2(int r, string c);
	int getC2(string c);
	string toNext2(string c);
	void run();
private:
	string T;
	string BWTS;
	vector<string> BWTS2;
	vector<string> Matrix;
	vector<int> SA;
	unordered_map<char, int> C;
	unordered_map<char, char> toNext;
	vector<string> Index2;
};