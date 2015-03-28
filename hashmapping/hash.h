#pragma once
#include<string>
#include<vector>
using namespace std;
class CPosition//kmerλ����
{
public:
	int pos;// kmerλ��
	CPosition * pnextPos;//��һ��λ��
	CPosition()
	{
		pos = 0;
		pnextPos = NULL;
	}
};
class CSubHT//��ϣ���ӱ�
{
public:
	char seq[8];//�洢����
	CPosition * pPos;//ָ������λ�ü���
	CSubHT * pSubHT;//��һ����ϣ�ӱ�
	CSubHT()
	{
		seq[0] = '\0';
		pPos = NULL;
		pSubHT = NULL;
	}

};
class CHashTable//��ϣ��ĸ��
{

public:
	char seq[8];//�洢����
	CSubHT * psubHT;//ָ���ϣ�ӱ�
	CHashTable()
	{
		seq[0] = '\0';
		psubHT = NULL;
	}

};

int createHashTable(const char * filename, CHashTable ht[],long length);


int  loadHashTable(const char * seqFile, const char * posFile, CHashTable ht[]);

int hashMapping(const char * queryFile, const CHashTable ht[]);
int kmerInsert(char kemr[], long curpos, CHashTable ht[]);
long calculatePos(char seq[], int length);
unsigned int BKDRHash(char *str);
CSubHT * getpPos(CSubHT * subHT, char seq[], int length);
fstream * sort(char *filepath);
vector<string>  quicksort(vector<string>  &kmers,  int begin,  int end);
void merge(vector<string> &kmers1, vector<string> &kmers2, vector<string> &kmers3);
