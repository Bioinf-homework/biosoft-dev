#include "stdafx.h"
#include "Hash_BWT.h"


Hash_BWT::Hash_BWT()
{
	string str[] = { "A", "C", "G", "T" };
	for (auto s:str)
	{
		nt.push_back(s);
	}
}


void Hash_BWT::makeHash(int k)
{
	vector<int> sep;
	for (int i = 0; i < pow(4, k); i++)
	{
		string s = IntToS(i, k);
		sep = search(s);
		HashIndex[s].push_back(sep[0]);
		HashIndex[s].push_back(sep[1]);
		//cout << sep[0] << "\t" << sep[1] << endl;
	}

	//cout << HashIndex["AAA"][0] << "\t" << HashIndex["AAA"][1];
	return;
}
string Hash_BWT::IntToS(int s, int k)
{
	string res;
	for (int i = 0; i < k; i++)
	{
		int x = s & 3;
		res += nt[x];
		s >>= 2;
	}
	//cout << res << endl;
	return res;
}

void Hash_BWT::Hash_search(string sub, int k, float e)
{

	int n = sub.length();
	//cout << n - 1 << endl;

	// 抽屉原理，将sub分成sn段
	float sn = n*e;
	sn = floor(sn) + 1;
	
	int each_ln = floor(n / sn)+1;
	if (sn == 1) {
		each_ln--;
	}
	if (each_ln < k)
	{
		cout << "分段后长度小于k，无效，退出。";
		return;
	}

	cout << "待查序列：" << sub << endl;
	cout << "总长度：" << n << endl;
	cout << "段：" << sn << endl;
	cout << "每段长度：" << each_ln << endl;
	
	int i = 0;
	vector<string> sub_arr;
	vector<vector<int>> sub_res;
	while (i<sn)
	{
		//分段问题。38bp。0.3错误率。运行错误11.6bp。分成12段。每段4个。报错。
		if (i*each_ln>n)
		{
			sn = i;
			cout << "实际分段：" << sn << endl;
			break;
		}
		string s = sub.substr(i*each_ln, each_ln);
		sub_arr.push_back(s);

		int sp, ep, j;
		//cout << s << s.length() << endl;
		// Si最后一个片段长度不一定
		if (s.length()!=k){

			if (s.length() < k)
			{
				j = s.length() - 1;
				//cout << "sp" << endl;
				char c = s[j];
				sp = C[c] + 1;
				if (c != 'T')
				{
					ep = C[toNext[c]] + 1;
				}
				else{
					ep = T.length();
				}
				j--;
			}
			else{
				// length > k
				string lastK = s.substr(s.length() - k, k);

				sp = HashIndex[lastK][0];
				ep = HashIndex[lastK][1];
				if (sp == -1)	{
					i++;
					continue;
				}
				//倒数k+1个字符起继续 BWT搜索
				j = s.length() - k - 1;
			}

			//cout << s << "\t" << s.substr(s.length() - k, k) << endl;
			char c = s[j];

			while (sp < ep && j > -1)
			{
				c = s[j];
				//cout << c << endl;
				sp = LFC(sp, c);
				ep = LFC(ep, c);
				j--;
			}
		}
		else{
			// length == k的情况
			string lastK = s.substr(s.length() - k, k);

			sp = HashIndex[lastK][0];
			ep = HashIndex[lastK][1];
		}

		vector<int> res;
		if (sp == ep)
		{
			cout << s << "\tNooo...\n";
			res.push_back(-1);
		}
		else{

			cout << s << "\t";
			for (int i = sp; i < ep; i++)
			{
				cout << SA[i] << "\t";
				res.push_back(SA[i]);
			}
			cout << endl;
		}

		//cout << "这一段" << sp << ep << endl;
		
		sub_res.push_back(res);
		i++;
	}

	//sub_res[0].push_back(55);
	//去重
	sub_res = Check(sub_res, each_ln);

	// 调试输出
	cout << "去重后结果" << endl;
	for (int i = 0; i < sub_res.size(); i++)
	{
		for (int j = 0; j < sub_res[i].size(); j++)
		{
			cout << sub_res[i][j] << "\t";
		}
		cout << endl;
	}

	//取相应长度字符串，计算hamming
	vector<int> Final_res;
	for (int i = 0; i < sub_res.size(); i++)
	{
		for (int j = 0; j < sub_res[i].size(); j++)
		{
			int or_index = sub_res[i][j] - i*each_ln;
			if (or_index <0)
			{
				continue;
			}
			string or = T.substr(or_index, n);

			int d = Hamming(or, sub);

			if (d < n*e)
			{
				Final_res.push_back(or_index);
			}
			cout << "原始串：" << T.substr(or_index, n) << endl;
			cout << "距离：" << d << endl;
			//cout << or_index << endl;
		}
		cout << endl;
	}
	// todo:判重？
	cout << "result:" << endl;
	for (int i = 0; i < Final_res.size(); i++)
	{
		cout << Final_res[i] << "\t";
	}
	cout << endl;
	return;
}

bool Hash_BWT::HasPre(vector<int> s, int num)
{
	for (auto i : s)
	{
		if (i == num)
			return true;
	}
	return false;
}

vector<vector<int>> Hash_BWT::Check(vector<vector<int>> &sub_res, int length)
{
	// 
	if (sub_res[0][0]==-1)
	{
		sub_res[0].erase(sub_res[0].begin());
	}
	for (int i = sub_res.size() - 1; i > 0; i--)
	{
		string f;
		for (int j = 0; j < sub_res[i].size(); j++)
		{
			if (sub_res[i][j]==-1)
			{
				f += "X";
				//continue;
				break;
			}
			if (HasPre(sub_res[i - 1], sub_res[i][j] - length))
			{
				f += "X";
			}
			else{
				f += "O";
			}

		}
		for (int j = sub_res[i].size()-1; j>=0; j--)
		{
			if (f[j] == 'O')
				continue;
			else
				sub_res[i].erase(sub_res[i].begin() + j);
		}
	}
	return sub_res;
}

int Hash_BWT::Hamming(string c1, string c2)
{
	int res = 0;
	for (int i = 0; i < c1.length(); i++)
	{
		if (c1[i]!=c2[i])
		{
			res++;
		}
	}
	return res;
}

void Hash_BWT::run()
{
	Read_Reference("test1.fa");
	vector<string> res = Read_Subs("sub1.fa");
	preprocess();
	int k = 4;
	makeHash(k);
	//Hash_search("AACCCCGTCTCTACTGAAAAATACAAAAAAAAATTAGC", k, 0.1);
	//Hash_search("AACCCCGTCTCTACTGAAAAATACAAAAAAAAATTAGCCG", 3, 0.1);
	for (auto s : res){
		clock_t start, finish;
		start = clock();
		Hash_search(s, k, 0.1);
		finish = clock();
		float d2 = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f\n",  d2);
	}
	cout << endl;
}

Hash_BWT::~Hash_BWT()
{
}
