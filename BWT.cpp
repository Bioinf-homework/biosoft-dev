#include "stdafx.h"
#include "BWT.h"



BWT::BWT()
{
	toNext['A'] = 'C';
	toNext['C'] = 'G';
	toNext['G'] = 'T';
	string str[] = { "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "UU" };
	for (int i = 0; i < 17; i++)
	{
		Index2.push_back(str[i]);
	}
}

string BWT::toNext2(string c)
{
	int res = 0;
	for (int i = 0; i < 16; i++)
	{
		if (Index2[i] <= c)
			res++;
	}
	//cout << Index2[res];
	return Index2[res];
}

void BWT::Read_Reference(string filename)
{
	ifstream rfp;
	rfp.open(filename);
	string data;
	string res;
	//string str((istreambuf_iterator<char>(rfp)),istreambuf_iterator<char>());
	while (getline(rfp, data, '\n'))
	{
		//cout << data;
		res += data;
	}
	T = res;
	return;
}

vector<string> BWT::Read_Subs(string filename)
{
	ifstream rfp;
	rfp.open(filename);
	string data;
	vector<string> res;
	while (getline(rfp, data, '\n'))
	{
		res.push_back(data);
	}
	cout << "read subs" << endl;
	for (int i = 0; i < res.size(); i++)
	{
		cout << res[i] << endl;
	}
	return res;
}
void BWT::makebwts2()
{
	int n = T.length();
	// 构造BWT2(S)
	for (int i = 0; i < n; i++)
	{
		BWTS2.push_back(Matrix[i].substr(n - 2, 2));
		// 调试输出
		cout << Matrix[i] << "\t" << SA[i] << "\t" << Matrix[i].substr(n - 2, 2) << endl;;

	}
	//cout << BWTS << endl;
	return;

}
void BWT::preprocess()
{
	const int n = T.length();
	T.insert(n, "$");

	// 轮转矩阵
	for (int i = 0; i <= n; i++)
	{
		Matrix.push_back(T.substr(i) + T.substr(0, i));
		SA.push_back(i);
		string tmp;
		tmp = T.substr(i) + T.substr(0, i);
		//cout << tmp << endl;;
	}
	// 排序
	for (int i = 0; i <= n; i++)
	{
		for (int j = i + 1; j <= n; j++)
		{
			if (Matrix[i] > Matrix[j])
			{
				swap(Matrix[i], Matrix[j]);
				swap(SA[i], SA[j]);
			}
		}
	}
	cout << "\n\n\n";

	// 构造BWT(S)
	for (int i = 0; i <= n; i++)
	{
		// 调试输出
		//cout << Matrix[i] << "\t" << SA[i] << "\t" << Matrix[i][n] << endl;;
		BWTS += Matrix[i][n];
	}
	//cout << BWTS << endl;

	// 构造C
	for (auto c : T)
	{
		//c 为 char
		//cout << typeid(c).name();
		switch (c)
		{
		case 'A':
			C['C'] += 1;
		case 'C':
			C['G'] += 1;
		case 'G':
			C['T'] += 1;
		}
	}
	return;
}

int BWT::Occ2(int r, string c)
{
	int res = 0;
	for (int i = 0; i < r; i++)
	{
		if (BWTS2[i] == c)
			res += 1;
	}
	return res;
}

int BWT::Occ(int r, char c)
{
	int res = 0;
	for (int i = 0; i < r; i++)
	{
		if (BWTS[i] == c)
			res += 1;
	}
	return res;
}

int BWT::getC2(string c)
{
	if (c == "UU")
		return T.length();
	int n = Matrix.size(), res = 0;
	for (int i = 0; i < n; i++)
	{
		if (Matrix[i].substr(0, 2) < c)
		{
			res++;
		}
	}
	return res;
}
int BWT::LFC(int r, char c)
{
	return C[c] + Occ(r, c) + 1;
}
int BWT::LFC2(int r, string c)
{
	return getC2(c) + Occ2(r, c);
}
void BWT::search2(string sub)
{
	int n = sub.length();
	string c = sub.substr(n - 2, 2);
	int sp, ep;
	sp = getC2(c);
	//cout << toNext2(c) << endl;
	ep = getC2(toNext2(c));
	//cout << sp << "\t" << ep << endl;
	int i = n - 2;
	//cout << i << endl;
	while (sp < ep && i>1)
	{
		//c = sub[i];
		c = sub.substr(i - 2, 2);
		//cout << c
		sp = LFC2(sp, c);
		ep = LFC2(ep, c);
		i -= 2;
		//cout << i << "\t" << c << "\t" << sp << "\t" << ep << endl;
	}
	if (i == 1)
	{
		// todo:
		sp = LFC(sp, sub[0]);
		ep = LFC(ep, sub[0]);
		//cout << sp << "\t" << ep << endl;
	}

	if (sp == ep)
	{
		cout << sub << "\tNooo...\n";
	}
	else{
		cout << sub << "\t";
		for (int i = sp; i < ep; i++)
		{
			cout << SA[i] << "\t";
		}
		cout << endl;
	}
}
void BWT::search(string sub)
{
	int n = sub.length();
	char c = sub[n - 1];
	int sp, ep;
	sp = C[c] + 1;
	if (c != 'T')
	{
		ep = C[toNext[c]] + 1;
	}
	else{
		ep = T.length();
	}
	//cout << sp << "\t" << ep << endl;
	int i = n - 2;
	while (sp < ep && i > -1)
	{
		c = sub[i];
		sp = LFC(sp, c);
		ep = LFC(ep, c);
		i--;
	}
	//cout << sp << "\t" << ep << endl;

	if (sp == ep)
	{
		cout << sub << "\tNooo...";
	}
	else{
		cout << sub << "\t";
		for (int i = sp; i < ep; i++)
		{
			cout << SA[i] << "\t";
		}
		cout << endl;
	}
}


void BWT::run()
{
	//T = "GGCTTCCTAC";
	//C['G'] = 0;
	//C['A'] = 0;
	//C['T'] = 0;
	//C['C'] = 0;
	//preprocess();

	//cout << C['A'] << "\t" << C['C'] << "\t" << C['G'] << "\t" << C['T'] << "\t" << endl;;

	//search("CTTTT");

	Read_Reference("test.fa");
	vector<string> res = Read_Subs("sub.fa");
	preprocess();
	makebwts2();
	//cout << Occ2(8, "GT") << endl;
	//cout << Occ2(10, "GT") << endl;
	//cout << Occ2(19, "GG") << endl;
	//cout << Occ2(20, "GG") << endl;
	//search2("CTACTTCAGGGTCA");
	for (auto s : res){
		search(s);
		search2(s);
	}

}
BWT::~BWT()
{
}
