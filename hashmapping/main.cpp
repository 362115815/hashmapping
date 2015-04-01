#include<iostream>
#include<fstream>
#include "hash.h"
#include<cmath>
#include<ctime>
#define N 16385
#define M  1000
using namespace std;
int main()
{
	CHashTable ht[N];
	char *  filepath = "D://kmers.txt";
	//createHashTable(filepath, ht,N);
	long start, finish;
	double totaltime;
	start = clock();


	sort(filepath);

	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n运行时间为:" << totaltime << "秒" << endl;
	system("pause");
	return 0;
}
