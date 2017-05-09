#include "stdafx.h"
#include "BWT.h"


//��ʼ��
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

//Summary:  ��ȡ�ο�������

//Parameters:

//       filename: ��ȡ���ļ���

//Return : ����ο������鴮
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

//Summary:  ��ȡ��������

//Parameters:

//       filename: ��ȡ���ļ���

//Return : ������������
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

//Summary: Ԥ������ת�������򣬹���BWT������C����

//Return : ��
void BWT::preprocess()
{
	clock_t start, finish;
	start = clock();
	const int n = T.length();
	T.insert(n, "$");

	// ��ת����
	for (int i = 0; i <= n; i++)
	{
		Matrix.push_back(T.substr(i) + T.substr(0, i));
		SA.push_back(i);
		string tmp;
		tmp = T.substr(i) + T.substr(0, i);
		//cout << tmp << endl;;
	}

	finish = clock();
	float d1 = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("������ת����ʱ��: %f\n", d1);
	start = clock();

	


	//sort(Matrix.begin(), Matrix.end());

	// ����
	QuickSort(Matrix, 0, Matrix.size() - 1);
	//for (int i = 0; i <= n; i++)
	//{
	//	for (int j = i + 1; j <= n; j++)
	//	{
	//		if (Matrix[i] > Matrix[j])
	//		{
	//			swap(Matrix[i], Matrix[j]);
	//			swap(SA[i], SA[j]);
	//		}
	//	}
	//}

	cout << "\n\n\n";

	finish = clock();
	d1 = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("����ʱ��: %f\n", d1);
	start = clock();
	// ����BWT(S)
	for (int i = 0; i <= n; i++)
	{
		// �������
		//cout << Matrix[i] << "\t" << SA[i] << "\t" << Matrix[i][n] << endl;;
		BWTS += Matrix[i][n];
	}
	//cout << BWTS << endl;

	// ����C
	for (auto c : T)
	{
		//c Ϊ char
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

	finish = clock();
	d1 = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("��������ʱ��: %f\n", d1);

	return;
}
//Summary:  ����ӵ�1�е���r�У�c�ַ����ֵĴ���

//Parameters:

//       r: ��ֹ����

//		 c: ������ַ�

//Return : ���ִ���
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


//Summary:  �����ַ�c����ʼ����

//Parameters:

//		 c: ������ַ�

//Return : ����
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

//Summary:  ����Ҫ��ת������

//Parameters:

//		 r: ��ǰ����

//		 c: ��ǰ�ַ�

//Return : ����
int BWT::LFC(int r, char c)
{
	//return C[c] + Occ(r, c) + 1;
	return getC(c) + Occ(r, c);
}

//Summary:  �Ǿ�ȷƥ�������£��ڲο���������������������λ��

//Parameters:

//		 sub:��������

//		 e: ����Ĵ�����

//Return : ƥ��λ������
void BWT::unexactsearch(string sub, float e)
{

}

//Summary:  ����ƥ�������£��ڲο���������������������λ��

//Parameters:

//		 sub:��������

//Return : ƥ��λ������
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

//Summary:  ������������������
void BWT::run()
{
	Read_Reference("test.fa");
	vector<string> res = Read_Subs("sub.fa");

	preprocess();

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

//Summary:  ������������֮��ı༭����

//Parameters:

//		 c1:����1

//		 c2:����2

//Return : �༭���볤��
int BWT::editDis(string c1, string c2)
{
	int n, m;
	n = c1.length();
	m = c2.length();

	vector< vector<int> > M(m + 1, vector<int>(n + 1, 0));
	//�߽�����
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
			// ״̬ת��..
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
