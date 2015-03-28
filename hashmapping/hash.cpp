#pragma once
#include "hash.h"
#define N 16385
#define M  500
#include <fstream>
#include <iostream>
#include<cmath>
#include<cstring>
using namespace std;

int createHashTable(const char * filename, CHashTable ht[], long length)
{
	ifstream fin;
	fin.open(filename);

	//去除首行
	char 	line[20];
	fin.getline(line, 80);
	std::cout << line << endl;

	//创建哈希表，先读取字符，用队列处理
	long curpos = 0;//bp位置
	char seqQue[22];//kmer队列
	seqQue[21] = '\0';
	for (int i = 0; i < 20; i++)
	{
		seqQue[i] = fin.get();
	}

	/*std::cout << seqQue << endl;*/

	int head = 20;//队列头部
	int tail = 0;//队列尾部
	char kmer[22];
	kmer[21] = '\0';
	char ch;
	long flag = 21;
	
	ofstream fout("D:\\result.txt");



	while (fin >> ch)
	{
		flag++;
		curpos++;//当前位置
		if (ch == 'N')
		{
			flag = 0;
		}

		//更新队列

		head = head % 21;
		seqQue[head++] = ch;


		if (flag >= 21)
		{
			//更新kmer			
			for (int i = 0, j = tail; i < 21; j++, i++)
			{
				j %= 21;
				//kmer[i] = seqQue[j];
				fout << seqQue[j];
			}
			

			//输入kmer
				

				fout << ' ' << curpos << endl;

			//int flag = kmerInsert(kmer, curpos, ht);
		}
		tail = (tail + 1) % 21;
	}
	fout.close();
	

	return 0;
}

int  loadHashTable(const char * seqFile, const char * posFile, CHashTable ht[])
{
	return 0;
}

int hashMapping(const char * queryFile, const CHashTable ht[])
{
	return 0;
}

//插入kmer到hash表中，成功返回 1，失败返回0
int kmerInsert(char kmer[], long curpos, CHashTable ht[])
{


	//先处理前面7位
	char seq[8];
	seq[7] = '\0';
	for (int i = 0; i < 7; i++)
	{
		if (kmer[i] == 'N')
		{
			return -1;//如果是N ，则返回
		}
		seq[i] = kmer[i];
	}

	long pos = calculatePos(seq, 7);//计算前7位存储位置

	if (pos == -1)
	{
		return 0;
	}

	if (!strcmp(ht[pos].seq, ""))
	{
		strcpy_s(ht[pos].seq, seq);//将前7位存储在相应位置
	}




	//处理中间7位
	for (int i = 7; i < 14; i++)//获取中间7位序列
	{
		if (kmer[i] == 'N')
		{
			return -1;//如果是N ，则返回
		}
		seq[i - 7] = kmer[i];
	}

	CSubHT * subHT = NULL;
	if (ht[pos].psubHT == NULL)//如果第一张的表的子表为空，则申请一张子表
	{
		subHT = new CSubHT[M];
		ht[pos].psubHT = subHT;
	}
	else
	{
		subHT = ht[pos].psubHT;
	}


	//获取中间7位存储地址
	subHT = getpPos(subHT, seq, 7);
	if (!strcmp(subHT->seq, ""))//若没有存储序列，则将序列存进去
	{
		strcpy_s(subHT->seq, seq);
	}



	//处理后面7位
	
	for (int i = 14; i < 21; i++)
	{
	if (kmer[i] == 'N')
	{
	return -1;//如果是N ，则返回
	}
	seq[i - 14] = kmer[i];
	}


	
	CSubHT * subHT1 = NULL;

	if (subHT->pSubHT == NULL)
	{
		subHT->pSubHT = new CSubHT[M];
	}

	subHT1 = subHT->pSubHT;


	//获取后面7位存储地址
	subHT1 = getpPos(subHT1, seq, 7);
	if (!strcmp(subHT1->seq, ""))//若没有存储序列，则将序列存进去
	{
		strcpy_s(subHT1->seq, seq);
	}
	
	//存储bp位置
	CPosition * pPos = new CPosition;
	pPos->pos = curpos;
	pPos->pnextPos = subHT1->pPos;
	subHT1->pPos = pPos;

	return 0;
}

//hash函数，计算存储位置，用四进制,错误返回-1
long calculatePos(char seq[], int length)
{
	long pos = 0;
	for (int i = 0; i < length; i++)
	{
		switch (seq[i])
		{
		case 'A':
		{
					pos += 0;
					break;
		}
		case 'G':
		{
					pos += (long)pow(4, i);
					break;
		}
		case 'C':
		{
					pos += (long)pow(4, i) * 2;
					break;
		}
		case 'T':
		{
					pos += (long)pow(4, i) * 3;
					break;
		}
		default:
		{
				   return -1;
		}

		}
	}

	return pos;
}



unsigned int BKDRHash(char *str)
{
	unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
	unsigned int hash = 0;

	while (*str)
	{
		hash = hash * seed + (*str++);
	}

	return (hash & 0x7FFFFFFF);
}
CSubHT * getpPos(CSubHT * subHT, char seq[], int length)
{
	unsigned pos = BKDRHash(seq)%(M-1);//计算存储位置
	int flag = 0;
	
	//如果该位置序列等于"",或等于seq，则返回该位置
	while (1)
	{
		if (!strcmp((subHT + pos)->seq, "") || !strcmp((subHT + pos)->seq, seq))
		{
			return subHT + pos;
		}
		//否则说明发生冲突，继续计算
		else
		{
			flag++;
			if (flag == 30)//说明已经冲突30次
			{
				flag = 0;
				if ((subHT + M - 1)->pSubHT != NULL)//若subHT有下一张表，则进入下一张表寻找，否则新建一张表
				{
					subHT = (subHT + M - 1)->pSubHT;
				}
				else
				{

					(subHT + M - 1)->pSubHT = new CSubHT[M];
					subHT = (subHT + M - 1)->pSubHT;
				}
				pos = BKDRHash(seq) % (M - 1);

			}
			else
			{
				pos = (pos + (flag + 1) * 33) % (M - 1);
			}
			
		}
	}

	


	return	NULL;
}