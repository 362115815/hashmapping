#pragma once
#include "hash.h"
#define N 16384
#define M  500
#include <fstream>
#include <iostream>
#include<cmath>
#include<cstring>
#include<vector>
using namespace std;

int createHashTable(const char * filename, CHashTable ht[], long length)
{
	ifstream fin;
	fin.open(filename);

	//ȥ������
	char 	line[20];
	fin.getline(line, 80);
	std::cout << line << endl;

	//������ϣ���ȶ�ȡ�ַ����ö��д���
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
	if (ht[pos].psubHT == NULL)//�����һ�ŵı���ӱ�Ϊ�գ�������һ���ӱ�
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



	//�������7λ
	
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
				if ((subHT + M - 1)->pSubHT != NULL)//��subHT����һ�ű��������һ�ű�Ѱ�ң������½�һ�ű�
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



//���ļ��е�kmer��������
fstream * sort(char *filepath)
{
	
	/* radix sort BEGIN*/
	vector<string> kmers;
	vector<string>kmers1;
	ifstream fin(filepath);
	

	for (int i = 0; i <8000000; i++)
	{


		char ch[60];
		
		//fout<<fin.tellg()<<endl;

		fin.getline(ch, 60);

		kmers.push_back(ch);
	}


	//ifstream fin1("D:\\kmers.txt");
	//ofstream fout("D:\\1.txt");
	//for (int i = 0; i < 20000000; i++)
	//{
	//	char ch[100];
	//	fin.getline(ch, 100);
	//	fout << ch << endl;

	//}

	/*for (int i = 0; i <700; i++)
	{
		char ch[60];
		fin.getline(ch, 60);

		kmers1.push_back(ch);
	}*/

   radixSort(kmers);

	 /*ifstream fin1("D:\\radixsort.txt");
	 string ss[2];
	 int head = 0;
	 char ch[100];
	 fin1.getline(ch, 100);
	 ss[0] = ch;
	 for (int i = 0; i < 400000-1; i++)
	 {
		 int tail = head;
		 head = (head + 1) % 2;
		 fin1.getline(ch, 100);
		 ss[head] = ch;
		 if (ss[tail].compare(0, 21, ss[head], 0, 21)>0)
		 {
			 cout << ss[tail] << endl << ss[head] << endl << endl << endl;
		 }
	 }*/
	//��֤radixsort��ȷ��






	



























	/*radix sort END*/

	/* quick sort  BEGIN*/
	/*


	//long pos[N];

	vector<string> kmers;
	vector<string>kmers1;
	vector<string>kmers2;

    ifstream fin(filepath);
	for (int i = 0; i <700; i++)
	{
		char ch[60];
		fin.getline(ch, 60);
		
		kmers.push_back(ch);
	}
	quicksort(kmers, 0, kmers.size());

	//ofstream fout("D:\\result.txt");

	//for (int i = 0; i< kmers.size(); i++)
	//{

	//	fout << kmers[i] << endl;

	//}
	//fout.close();


	for (int i = 0; i <700; i++)
	{
		char ch[60];
		fin.getline(ch, 60);

		kmers1.push_back(ch);
	}

	quicksort(kmers1, 0, kmers1.size());
	

	for (int i = 0; i <700; i++)
	{
		char ch[60];
		fin.getline(ch, 60);

		kmers2.push_back(ch);
	}

	quicksort(kmers2, 0, kmers2.size());

	merge(kmers, kmers1,kmers2);
 	cout << kmers.size();



	*/

	/* quick sort END*/



	return NULL;
}


//��������

vector<string>  quicksort(vector<string>  &kmers, int begin, int end)
{
	if (begin < end)
	{
		
		int n = rand() % (end-begin) + begin;
		string temp=kmers[n];

		int i = begin;
		int j = end - 1;
		while (i<j)
		{
			//�ȴ��Ҳ࿪ʼ��
			while (i < j&& temp.compare(0,21,kmers[j],0,21)<=0)
			{
				j--;
			}
			//���i<j,����
			if (i < j)
			{
				kmers[n].assign(kmers[j]);
				n = j;
			}

			while (i<j&&temp.compare(0,21,kmers[i],0,21)>0)
			{
				i++;
			}
			if (i < j)
			{
				kmers[n].assign(kmers[i]);
				n = i;
			}

		}
		kmers[i].assign(temp);

		quicksort(kmers, begin, i - 1);
		quicksort(kmers, i + 1, end);
	}



	return kmers;
}

void merge(vector<string> &kmers1, vector<string> &kmers2, vector<string> &kmers3)
{
	ofstream fout("D:\\result1.txt");
	int count1 = kmers1.size();
	int count2 = kmers2.size();
	int count3 = kmers3.size();
	int i = 0, j = 0,k=0;
	while (i < count1&&j < count2&&k<count3)
	{
		if (kmers1[i].compare(0, 21, kmers2[j], 0, 21) <= 0 && kmers1[i].compare(0, 21, kmers3[k], 0, 21) <= 0)
		{
			fout << kmers1[i++] << endl;
		}
		else if (kmers2[j].compare(0,21,kmers3[k],0,21)<=0)
		{
			fout << kmers2[j++] << endl;
		}
		else
		{
			fout << kmers3[k++] << endl;
		}
	}
	while (i < count1&&j<count2)
	{
		if (kmers1[i].compare(0, 21, kmers2[j], 0, 21) <= 0)
		{
			fout << kmers1[i++] << endl;
		}
		else 
		{
			fout << kmers2[j++] << endl;
		}
		
	}
	while (j < count2&&k<count3)
	{

			 if (kmers2[j].compare(0, 21, kmers3[k], 0, 21) <= 0)
			{
				fout << kmers2[j++] << endl;
			}
			else
			{
				fout << kmers3[k++] << endl;
			}
		
	}

	while (i<count1&&k < count3)
	{
		if (kmers1[i].compare(0, 21, kmers3[k], 0, 21) <= 0)
		{
			fout << kmers1[i++] << endl;
		}
		else
		{
			fout << kmers3[k++] << endl;
		}
	}

	while (i < count1)
	{
		fout << kmers1[i++] << endl;
	}

	while (j < count2)
	{
		fout << kmers2[j++] << endl;
	}
	while (k < count3)
	{
		fout << kmers3[k++] << endl;
	}
	fout.close();
}

//��������
vector<string> radixSort(vector<string>& kmers)
{
	int length = kmers.size();
	vector<string>kmers1(length);
	long  radixcount[N];
	long *radix = new long[length];

	//����λ����
	for (int i = 0; i <N; i++)
	{
		radixcount[i] = 0;
	}

	for (int i = 0; i < length; i++)
	{
	
		radix[i] = getradix(kmers[i].substr(14, 7));
		radixcount[radix[i]]++;
	}

	for (int i = 1; i < N; i++)
	{
		radixcount[i] += radixcount[i - 1];
	}

	
	
	for (int i = length - 1; i >= 0; i--)
	{
		kmers1[radixcount[radix[i]]-1] = kmers[i];
		radixcount[radix[i]]--;

	}

	//�м�7λ����
	for (int i = 0; i <N; i++)
	{
		radixcount[i] = 0;
	}
	for (int i = 0; i < length; i++)
	{
		radix[i] = getradix(kmers1[i].substr(7, 7));
		radixcount[radix[i]]++;
	}

	for (int i = 1; i < N; i++)
	{
		radixcount[i] += radixcount[i - 1];
	}

	for (int i = length - 1; i >= 0; i--)
	{
		kmers[radixcount[radix[i]] - 1] = kmers1[i];
		radixcount[radix[i]]--;
	}

	//ǰ7λ����
	for (int i = 0; i <N; i++)
	{
		radixcount[i] = 0;
	}
	for (int i = 0; i < length; i++)
	{
		radix[i] = getradix(kmers[i].substr(0, 7));
		radixcount[radix[i]]++;
	}

	for (int i = 1; i < N; i++)
	{
		radixcount[i] += radixcount[i - 1];
	}

	for (int i = length - 1; i >= 0; i--)
	{
		kmers1[radixcount[radix[i]] - 1] = kmers[i];
		radixcount[radix[i]]--;
	}





	ofstream fout("D:\\radixsort.txt");
	for (int i = 0; i < length; i++)
	{
		fout << kmers1[i] << endl;
	}
		
	fout.close();





	delete[]radix;


	return kmers1;
}


//�������ж�Ӧ�Ļ���
long getradix(string s )
{
	long radix=0;
	int length = s.length();
	for (int i = 0; i < length; i++)
	{
		string temp= s.substr(i, 1);
		if (temp == "A")
		{
			radix += 0;
		}
		else if (temp == "C")
		{
			radix += (long)pow(4, length-1-i );
		}
		else if (temp == "G")
		{
			radix += (long)pow(4, length - 1 - i) * 2;
		}
		else if (temp == "T")
		{
			radix += (long)pow(4, length - 1 - i) * 3;
		}

	}
	return radix;
}