#include "stdafx.h"
#include "BWT2.h"


BWT2::BWT2()
{
	string str[] = { "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "UU" };
	for (int i = 0; i < 17; i++)
	{
		Index2.push_back(str[i]);
	}
}
string BWT2::toNext2(string c)
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

void BWT2::makebwts2()
{
	int n = T.length();
	// 构造BWT2(S)
	for (int i = 0; i < n; i++)
	{
		BWTS2.push_back(Matrix[i].substr(n - 2, 2));
		// 调试输出
		//cout << Matrix[i] << "\t" << SA[i] << "\t" << Matrix[i].substr(n - 2, 2) << endl;;

	}
	//cout << BWTS << endl;
	return;
}

int BWT2::Occ2(int r, string c)
{
	int res = 0;
	for (int i = 0; i < r; i++)
	{
		if (BWTS2[i] == c)
			res += 1;
	}
	return res;
}

int BWT2::getC2(string c)
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

int BWT2::LFC2(int r, string c)
{
	return getC2(c) + Occ2(r, c);
}

void BWT2::search2(string sub)
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

void BWT2::run()
{
	Read_Reference("test1.fa");
	vector<string> res = Read_Subs("sub1.fa");
	preprocess();
	makebwts2();
	for (auto s : res){
		clock_t start, mid, finish;
		start = clock();
		//search(s);
		//mid = clock();
		search2(s);
		finish = clock();
		//float d1 = (double)(mid - start) / CLOCKS_PER_SEC;
		//float d2 = (double)(finish - mid) / CLOCKS_PER_SEC;
		float d2 = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f\n", d2);
	}
	return;
}

BWT2::~BWT2()
{
}
