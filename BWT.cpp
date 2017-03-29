#include "stdafx.h"
#include "BWT.h"



BWT::BWT()
{
	toNext['A'] = 'C';
	toNext['C'] = 'G';
	toNext['G'] = 'T';
	//string str[] = { "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "UU" };
	//for (int i = 0; i < 17; i++)
	//{
	//	Index2.push_back(str[i]);
	//}
}


string BWT::Read_Reference(string filename)
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
	return T;
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
int BWT::getC(char c)
{
	int n = Matrix.size(), res = 0;
	for (int i = 0; i < n; i++)
	{
		if (Matrix[i][0] < c)
		{
			res++;
		}
	}
	return res;
}

int BWT::LFC(int r, char c)
{
	return C[c] + Occ(r, c) + 1;
	//return getC(c) + Occ(r, c);
}

void BWT::unexactsearch(string sub, float e)
{

}

vector<int> BWT::search(string sub)
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
	vector<int> res(2,-1);
	if (sp == ep)
	{
		cout << sub << "\tNooo...\n";
		return res;
	}
	else{
		
		cout << sub << "\t";
		for (int i = sp; i < ep; i++)
		{
			cout << SA[i] << "\t";
		}
		cout << endl;
		res[0] = sp;
		res[1] = ep;
		return res;
	}
}


void BWT::run()
{
	Read_Reference("test1.fa");
	vector<string> res = Read_Subs("sub1.fa");
	preprocess();
	makebwts2();
	for (auto s : res){
		clock_t start, mid, finish;
		start = clock();
		search(s);
		mid = clock();
		search2(s);
		finish = clock();
		float d1 = (double)(mid - start) / CLOCKS_PER_SEC;
		float d2 = (double)(finish - mid) / CLOCKS_PER_SEC;
		printf("%f\t%f\n", d1, d2);
	}
	return;
}
int BWT::editDis(string c1, string c2)
{
	int n, m;
	n = c1.length();
	m = c2.length();

	vector< vector<int> > M(m + 1, vector<int>(n + 1, 0));
	//边界设置
	for (int i = 0; i < M.size(); i++)
	{
		M[i][0] = i;
	}
	for (int j = 0; j < M[0].size(); j++)
	{
		M[0][j] = j;
	}
	for (int i = 1; i < M.size(); i++)
	{
		for (int j = 1; j < M[0].size(); j++)
		{
			// 状态转移..
			if (c1[j - 1] == c2[i - 1])
			{
				M[i][j] = M[i - 1][j - 1];
			}
			else
			{
				M[i][j] = __min(M[i - 1][j - 1], __min(M[i - 1][j], M[i][j - 1])) + 1;
			}
		}
	}
	return M[m][n];
}

BWT::~BWT()
{
}
