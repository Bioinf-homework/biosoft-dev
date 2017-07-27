#include "stdafx.h"
#include "BWT.h"


//初始化
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

//Summary:  读取参考基因组

//Parameters:

//       filename: 读取的文件名

//Return : 整体参考基因组串
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
		transform(data.begin(), data.end(), data.begin(), ::toupper);
		res += data;
	}
	T = res;
	return T;
}

//Summary:  读取待查序列

//Parameters:

//       filename: 读取的文件名

//Return : 待查序列数组
vector<string> BWT::Read_Subs(string filename)
{
	ifstream rfp;
	rfp.open(filename);
	string data;
	vector<string> res;
	while (getline(rfp, data, '\n'))
	{
		transform(data.begin(), data.end(),data.begin(), ::toupper);
		for (int i = 0; i < data.length(); i++)
		{
			if (data[i] != 'A' && data[i] != 'C' && data[i] != 'T' && data[i] != 'G')
			{
				cout << data[i] << endl;
			}
		}
		res.push_back(data);
	}
	cout << "read subs" << endl;
	for (int i = 0; i < res.size(); i++)
	{
		cout << res[i] << endl;
	}
	return res;
}
//Summary: 预处理，轮转矩阵，排序，构造BWT串，建C索引

//Return : 无
void BWT::preprocess2()
{
	clock_t start, finish;
	start = clock();
	const int n = T.length();
	T.insert(n, "$");
	//vector<int> SA;
	for (int i = 0; i < T.length(); i++)
	{
		SA.push_back(i);
	}
	
	finish = clock();
	float d1 = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("构造轮转矩阵时间: %f\n", d1);
	start = clock();

	QuickSort(0, T.length()-1);

	cout << "\n\n\n";

	finish = clock();
	d1 = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("排序时间: %f\n", d1);
	start = clock();
	// 构造BWT(S)
	string FirstS;
	for (int i = 0; i < SA.size(); i++)
	{
		if (SA[i]==0)
		{
			BWTS += "$";
		}
		else
		{
			BWTS += T[SA[i] - 1];
		}
		FirstS += T[SA[i]];
	}
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


	//构造checkpoint
	int now = 0;
	unordered_map<char, int> cp;
	for (int i = 0; i < BWTS.length(); i++)
	{
		if (i / 200 >= 1 && i % 200 == 0)
		{
			//cout << cp['A'] << cp['C'] << cp['G'] << cp['T'] << endl;
			Checkpoint.push_back(cp);
		}
		switch (BWTS[i])
		{
		case 'A':
			cp['A'] += 1;
			break;
		case 'C':
			cp['C'] += 1;
			break;
		case 'G':
			cp['G'] += 1;
			break;
		case 'T':
			cp['T'] += 1;
			break;
		}
	}
	d1 = (double)(finish - start) / CLOCKS_PER_SEC;

	finish = clock();

	printf("构建索引时间: %f\n", d1);

	cout << "hello" << endl;
}

//Summary: 预处理，轮转矩阵，排序，构造BWT串，建C索引

//Return : 无
void BWT::preprocess()
{}
//Summary:  计算从第1行到第r行，c字符出现的次数

//Parameters:

//       r: 截止行数

//		 c: 计算的字符

//Return : 出现次数
int BWT::Occ(int r, char c)
{
	int res = 0;
	int k = r / 200, m = r % 200;
	if (k>=1)
	{
		res += Checkpoint[k - 1][c];
	}
	for (int i = k*200; i < r; i++)
	{
		if (BWTS[i] == c)
			res += 1;
	}
	return res;
}


//Summary:  计算字符c的起始行数

//Parameters:

//		 c: 计算的字符

//Return : 行数
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

//Summary:  计算要跳转的行数

//Parameters:

//		 r: 当前行数

//		 c: 当前字符

//Return : 行数
int BWT::LFC(int r, char c)
{
	return C[c] + Occ(r, c) + 1;
	//return getC(c) + Occ(r, c);
}

//Summary:  非精确匹配条件下，在参考串中搜索待查序列所在位置

//Parameters:

//		 sub:待查序列

//		 e: 允许的错误率

//Return : 匹配位置数组
void BWT::unexactsearch(string sub, float e)
{

}

//Summary:  精切匹配条件下，在参考串中搜索待查序列所在位置

//Parameters:

//		 sub:待查序列

//Return : 匹配位置数组
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
	clock_t start, finish;
	while (sp < ep && i > -1)
	{
		start = clock();
		c = sub[i];
		sp = LFC(sp, c);
		ep = LFC(ep, c);
		i--;
		finish = clock();
		//printf("%f\n", (double)(finish - start) / CLOCKS_PER_SEC);
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

//Summary:  主函数，用于输出结果
void BWT::run()
{
	Read_Reference("test.fa");
	vector<string> res = Read_Subs("sub.fa");

	preprocess2();

	for (auto s : res){
		clock_t start, mid, finish;
		start = clock();
		search(s);
		finish = clock();
		float d1 = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f\n", d1);
	}
	return;
}

//Summary:  计算两个序列之间的编辑距离

//Parameters:

//		 c1:序列1

//		 c2:序列2

//Return : 编辑距离长度
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
