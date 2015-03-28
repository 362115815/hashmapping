#pragma once
#include<string>
#include<vector>
using namespace std;
class CPosition//kmer位置类
{
public:
	int pos;// kmer位置
	CPosition * pnextPos;//下一个位置
	CPosition()
	{
		pos = 0;
		pnextPos = NULL;
	}
};
class CSubHT//哈希表子表
{
public:
	char seq[8];//存储序列
	CPosition * pPos;//指向序列位置集合
	CSubHT * pSubHT;//下一个哈希子表
	CSubHT()
	{
		seq[0] = '\0';
		pPos = NULL;
		pSubHT = NULL;
	}

};
class CHashTable//哈希表母表
{

public:
	char seq[8];//存储序列
	CSubHT * psubHT;//指向哈希子表
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
