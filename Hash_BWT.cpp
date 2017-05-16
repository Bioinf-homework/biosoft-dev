#include "stdafx.h"
#include "Hash_BWT.h"

//初始化
Hash_BWT::Hash_BWT()
{
	string str[] = { "A", "C", "G", "T" };
	for (auto s : str)
	{
		nt.push_back(s);
	}
}

bool Hash_BWT::indic(string s)
{
	for (int i = 0; i < dic.size(); i++)
	{
		if (s == dic[i])
		{
			return true;
		}
	}
	return false;
}

//Summary:  构建k-mer索引

//Parameters:

//       k: 索引序列的长度
void Hash_BWT::makeHash(int k)
{
	vector<int> sep;

	int n = 1;
	HashIndex[IntToS(0, k)].push_back(1);
	
	for (int i = 0; i < pow(4, k); i++)
	{
		string s = IntToS(i, k);
		//cout << s;
		dic.push_back(s);
		for (; n < T.length(); n++)
		{
			if (Matrix[n].substr(0, k)==s)
			{
				continue;
			}
			else
			{
				if (Matrix[n].substr(0, k).find("$") != -1)
				{
					n++;
				}
				HashIndex[s].push_back(n);
				HashIndex[IntToS(i + 1, k)].push_back(n);
				break;
			}
			
		}
		
		//sep = search(s);
		//HashIndex[s].push_back(sep[0]);
		//HashIndex[s].push_back(sep[1]);
		//cout << sep[0] << "\t" << sep[1] << endl;
	}
	HashIndex[IntToS(pow(4, k)-1, k)].push_back(T.length());
	//cout << HashIndex["AAA"][0] << "\t" << HashIndex["AAA"][1];
	return;
}

//Summary:  用二进制数映射到序列

//Parameters:

//       s: 当前数

//		 k: 序列长度

//Return : 对应的序列
string Hash_BWT::IntToS(int s, int k)
{
	string res;
	for (int i = 0; i < k; i++)
	{
		// 位操作：取最后两位
		int x = s & 3;
		res += nt[x];
		s >>= 2;
	}
	//cout << res << endl;
	string r(res.rbegin(), res.rend());
	return r;
}

//Summary:  抽屉原理划分待查序列，然后验证所有位置，放回符合的位置

//Parameters:

//       s: 待查训练

//		 k: hash索引序列长度

//		 e: 允许的错误率

//Return : 符合匹配的位置
void Hash_BWT::Hash_search(string sub, int k, float e)
{

	int n = sub.length();
	//cout << n - 1 << endl;

	// 抽屉原理，将sub分成sn段
	float sn = n*e;
	sn = floor(sn) + 1;

	int each_ln = floor(n / sn);

	int yu = n - each_ln*sn;

	//if (each_ln*sn != n)
	//{
	//	sn = ceil(n / each_ln);
	//}

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
		//分段问题。38bp。0.3错误率。允许错误11.6bp。分成12段。每段4个。报错。
		if (i*each_ln>n)
		{
			sn = i;
			cout << "实际分段：" << sn << endl;
			break;
		}
		string s;
		if (i*each_ln + each_ln + yu == n)
		{
			s = sub.substr(i*each_ln);
		}
		else
		{
			s = sub.substr(i*each_ln, each_ln);
		}

		sub_arr.push_back(s);

		vector<int> res;

		int sp, ep, j;
		//cout << s << s.length() << endl;
		// Si最后一个片段长度不一定
		if (s.length() != k){
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


				if (indic(lastK) == true)
				{
					sp = HashIndex[lastK][0];
					ep = HashIndex[lastK][1];
					if (sp == -1)	{
						i++;
						//res.push_back(-1);
						continue;
					}
					//倒数k+1个字符起继续 BWT搜索
					j = s.length() - k - 1;
				}
				else
				{
					j = -1;
					sp = 0;
					ep = 0;
				}

			}

			//cout << s << "\t" << s.substr(s.length() - k, k) << endl;
			if (j >= 0)
			{
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
		}
		else{
			// length == k的情况
			string lastK = s.substr(s.length() - k, k);
			if (indic(lastK) == true)
			{
				sp = HashIndex[lastK][0];
				ep = HashIndex[lastK][1];
			}
			else
			{
				j = -1;sp = 0;ep = 0;
			}
		}

		
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
	Check2(sub_res, each_ln);
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

			int d2 = Hamming(or, sub);
			int d = editDis(or, sub);

			if (d <= n*e)
			{
				Final_res.push_back(or_index);
			}
			cout << "原始串：" << or_index <<"\t"<< T.substr(or_index, n) << endl;
			cout << "edit距离：" << d << "\thamming 距离" << d2 << endl;
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

//Summary:  判断数组中是否已有某个数

//Parameters:

//		 s: 数组

//		 num: 待查数

//Return : bool
bool Hash_BWT::HasPre(vector<int> s, int num)
{
	for (auto i : s)
	{
		if (i == num)
			return true;
	}
	return false;
}
vector<vector<int>> Hash_BWT::Check2(vector<vector<int>> &sub_res, int length)
{
	for (int i = 0; i < sub_res.size()-1; i++)
	{
		for (int ii = 0; ii < sub_res[i].size(); ii++)
		{
			for (int j = i + 1; j < sub_res.size(); j++)
			{
				for (int k = 0; k < sub_res[j].size(); k++)
				{
					if (sub_res[j][k] == sub_res[i][ii] + length*(j - i))
					{
						sub_res[j].erase(sub_res[j].begin() + k);
					}
				}
			}
		}
	}
	return sub_res;
}


//Summary:  事先判断分段的待查序列对应的原始序列位置是否有重复，减少计算量。

//Parameters:

//		 sub_res: 各分段序列对应位置数组

//		 length: 分段长度

//Return : 去重后的对应位置数组
vector<vector<int>> Hash_BWT::Check(vector<vector<int>> &sub_res, int length)
{
	// 
	if (sub_res[0][0] == -1)
	{
		sub_res[0].erase(sub_res[0].begin());
	}
	for (int i = sub_res.size() - 1; i > 0; i--)
	{
		string f;
		for (int j = 0; j < sub_res[i].size(); j++)
		{
			if (sub_res[i][j] == -1)
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
		for (int j = sub_res[i].size() - 1; j >= 0; j--)
		{
			if (f[j] == 'O')
				continue;
			else
				sub_res[i].erase(sub_res[i].begin() + j);
		}
	}
	return sub_res;
}

//Summary:  计算两个序列之间的haming距离

//Parameters:

//		 c1:序列1

//		 c2:序列2

//Return : hamming距离长度
int Hash_BWT::Hamming(string c1, string c2)
{
	int res = 0;
	for (int i = 0; i < c1.length(); i++)
	{
		if (c1[i] != c2[i])
		{
			res++;
		}
	}
	return res;
}

//Summary:  主函数，用于输出结果
void Hash_BWT::run()
{
	clock_t start, finish;
	Read_Reference("test.fa");
	vector<string> res = Read_Subs("sub.fa");
	preprocess();
	int k = 6;
	start = clock();
	makeHash(k);
	finish = clock();
	float d2 = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("hash索引：%f\n", d2);
	//Hash_search("AACCCCGTCTCTACTGAAAAATACAAAAAAAAATTAGC", k, 0.1);
	//Hash_search("AACCCCGTCTCTACTGAAAAATACAAAAAAAAATTAGCCG", 3, 0.1);
	for (auto s : res){

 		start = clock();
		Hash_search(s, k, 0.1);
		finish = clock();
		d2 = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f\n", d2);
	}
	cout << endl;
}

Hash_BWT::~Hash_BWT()
{
}
