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
	vector<unordered_map<char, int>> Checkpoint;
	unordered_map<char, int> C;
	unordered_map<char, char> toNext;
	vector<string> Matrix;


	//每一趟选一个基准数，找出比它大的和比它小的个数，就能确定它的位置！
	int partition(vector<string> &arr, int low, int high) {
		//    if (low >= high) return
		string base = arr[high];
		int tail = low - 1;
		for (int i = low; i < high; i++)   // 遍历基准以外的其他元素
		{
			if (arr[i] <= base)            // 把小于等于基准的元素放到前一个子数组中
			{
				tail++;
				if (tail == i) continue;
				else 
				{
					swap(arr[tail], arr[i]);
					swap(SA[tail], SA[i]);
				}
			}
		}
		swap(arr[tail + 1], arr[high]);
		swap(SA[tail+1], SA[high]);
		return tail + 1;
	}

	void QuickSort(vector<string> &arr, int low, int high) {
		if (high < low) return;
		int Index = partition(arr, low, high); //将表一分为二
		QuickSort(arr, low, Index - 1);        //递归对低子表递归排序
		QuickSort(arr, Index + 1, high);
	}

private:
	string BWTS;
	//vector<string> Index2;
};