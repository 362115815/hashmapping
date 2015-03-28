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

	//ȥ������
	char 	line[20];
	fin.getline(line, 80);
	std::cout << line << endl;

	//������ϣ�����ȶ�ȡ�ַ����ö��д���
	long curpos = 0;//bpλ��
	char seqQue[22];//kmer����
	seqQue[21] = '\0';
	for (int i = 0; i < 20; i++)
	{
		seqQue[i] = fin.get();
	}

	/*std::cout << seqQue << endl;*/

	int head = 20;//����ͷ��
	int tail = 0;//����β��
	char kmer[22];
	kmer[21] = '\0';
	char ch;
	long flag = 21;
	
	ofstream fout("D:\\result.txt");



	while (fin >> ch)
	{
		flag++;
		curpos++;//��ǰλ��
		if (ch == 'N')
		{
			flag = 0;
		}

		//���¶���

		head = head % 21;
		seqQue[head++] = ch;


		if (flag >= 21)
		{
			//����kmer			
			for (int i = 0, j = tail; i < 21; j++, i++)
			{
				j %= 21;
				//kmer[i] = seqQue[j];
				fout << seqQue[j];
			}
			

			//����kmer
				

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

//����kmer��hash���У��ɹ����� 1��ʧ�ܷ���0
int kmerInsert(char kmer[], long curpos, CHashTable ht[])
{


	//�ȴ���ǰ��7λ
	char seq[8];
	seq[7] = '\0';
	for (int i = 0; i < 7; i++)
	{
		if (kmer[i] == 'N')
		{
			return -1;//�����N ���򷵻�
		}
		seq[i] = kmer[i];
	}

	long pos = calculatePos(seq, 7);//����ǰ7λ�洢λ��

	if (pos == -1)
	{
		return 0;
	}

	if (!strcmp(ht[pos].seq, ""))
	{
		strcpy_s(ht[pos].seq, seq);//��ǰ7λ�洢����Ӧλ��
	}




	//�����м�7λ
	for (int i = 7; i < 14; i++)//��ȡ�м�7λ����
	{
		if (kmer[i] == 'N')
		{
			return -1;//�����N ���򷵻�
		}
		seq[i - 7] = kmer[i];
	}

	CSubHT * subHT = NULL;
	if (ht[pos].psubHT == NULL)//�����һ�ŵı����ӱ�Ϊ�գ�������һ���ӱ�
	{
		subHT = new CSubHT[M];
		ht[pos].psubHT = subHT;
	}
	else
	{
		subHT = ht[pos].psubHT;
	}


	//��ȡ�м�7λ�洢��ַ
	subHT = getpPos(subHT, seq, 7);
	if (!strcmp(subHT->seq, ""))//��û�д洢���У������д��ȥ
	{
		strcpy_s(subHT->seq, seq);
	}



	//��������7λ
	
	for (int i = 14; i < 21; i++)
	{
	if (kmer[i] == 'N')
	{
	return -1;//�����N ���򷵻�
	}
	seq[i - 14] = kmer[i];
	}


	
	CSubHT * subHT1 = NULL;

	if (subHT->pSubHT == NULL)
	{
		subHT->pSubHT = new CSubHT[M];
	}

	subHT1 = subHT->pSubHT;


	//��ȡ����7λ�洢��ַ
	subHT1 = getpPos(subHT1, seq, 7);
	if (!strcmp(subHT1->seq, ""))//��û�д洢���У������д��ȥ
	{
		strcpy_s(subHT1->seq, seq);
	}
	
	//�洢bpλ��
	CPosition * pPos = new CPosition;
	pPos->pos = curpos;
	pPos->pnextPos = subHT1->pPos;
	subHT1->pPos = pPos;

	return 0;
}

//hash����������洢λ�ã����Ľ���,���󷵻�-1
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
	unsigned pos = BKDRHash(seq)%(M-1);//����洢λ��
	int flag = 0;
	
	//�����λ�����е���"",�����seq���򷵻ظ�λ��
	while (1)
	{
		if (!strcmp((subHT + pos)->seq, "") || !strcmp((subHT + pos)->seq, seq))
		{
			return subHT + pos;
		}
		//����˵��������ͻ����������
		else
		{
			flag++;
			if (flag == 30)//˵���Ѿ���ͻ30��
			{
				flag = 0;
				if ((subHT + M - 1)->pSubHT != NULL)//��subHT����һ�ű����������һ�ű�Ѱ�ң������½�һ�ű�
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