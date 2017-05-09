// biosoft-dev2.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"



#include "BWT.h"
#include "BWT2.h"
#include "Hash_BWT.h"


int _tmain(int argc, _TCHAR* argv[])
{	
	//BWT demo;
	//demo.run();

	//BWT2 demo;
	//demo.run();
	//string strA("abasdjkaldkjalsd");
	//transform(strA.begin(), strA.end(), strA.begin(), ::toupper);
	//cout << strA << endl;

	Hash_BWT demo2;
	demo2.run();

	system("pause");
	return 0;
}

