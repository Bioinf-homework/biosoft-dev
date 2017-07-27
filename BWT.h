#pragma once

class BWT
{
public:
	BWT();
	~BWT();
	string Read_Reference(string filename);
	vector<string> Read_Subs(string filename);
	void preprocess();
	void preprocess2();
	vector<int> search(string sub);
	void unexactsearch(string sub, float e);
	int Occ(int r, char c);
	int LFC(int r, char c);
	int getC(char c);
	int editDis(string c1, string c2);
	void run();
	string T;
	vector<int> SA;
	vector<unordered_map<char, int>> Checkpoint;
	unordered_map<char, int> C;
	unordered_map<char, char> toNext;
	vector<string> Matrix;

	string SAstring(int i)
	{
		return T.substr(SA[i], SA.size() - SA[i]);
	}
	//ÿһ��ѡһ����׼�����ҳ�������ĺͱ���С�ĸ���������ȷ������λ�ã�
	int partition(int low, int high) {
		string base = SAstring(high);
		int tail = low - 1;
		for (int i = low; i < high; i++)   // ������׼���������Ԫ��
		{
			if (SAstring(i) <= base)            // ��С�ڵ��ڻ�׼��Ԫ�طŵ�ǰһ����������
			{
				tail++;
				if (tail == i) continue;
				else
				{
					swap(SA[tail], SA[i]);
				}
			}
		}
		swap(SA[tail + 1], SA[high]);
		return tail + 1;
	}

	void QuickSort(int low, int high) {
		if (high < low) return;
		int Index = partition(low, high); //����һ��Ϊ��
		QuickSort(low, Index - 1);        //�ݹ�Ե��ӱ�ݹ�����
		QuickSort(Index + 1, high);
	}


private:
	string BWTS;
	//vector<string> Index2;
};