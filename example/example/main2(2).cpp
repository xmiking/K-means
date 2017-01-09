#include<iostream>
#include<math.h>
#include<fstream>
#include<vector>   
#include <io.h>
#include <direct.h>
//#include <stdlib.h>
#include"stdio.h"
#include "ctime"
using namespace std;

//
struct PictureFeature
{
	vector<float> m_value;
	string name;
};

// 清空每个聚类中心上的点的集合
int ClearClusterVector(vector<vector<int> > &cluster_vectors)//cluster_vectors:外层向量，存放每堆点的，向量外的向量
{
	for (size_t i = 0; i < cluster_vectors.size(); i++)
	{
		cluster_vectors[i].clear();
	}
	return 0;
}

//欧式距离：每个距离的平方之和，平方之和再开根号
double CaculateEuclideanDistance(vector<float> &v1, vector<float> &v2)
{
	if (v2.size() != v1.size())
	{
		cout << "error" << "\n";
		return -1.0;
	}

	vector<float> v3;
	v3.resize(v2.size());

	for (size_t i = 0; i < v2.size(); i++)
	{
		v3[i] = v2[i] - v1[i];
	}

	double sq_distance = 0.0;

	for (size_t i = 0; i < v3.size(); i++)
	{
		double tmp_value = v3[i] * v3[i];
		sq_distance += tmp_value;  //*****************应该为之和：内积为对应坐标的乘积之和
	}

	double distance = sqrt(sq_distance);//欧式距离值

	return distance;
}

// 计算所有点到簇中心点的欧式距离，并归类
int CaculateMatchValueSum(const std::vector<int> &current_cluster,                       //current_cluster：所有的点向量的索引
	vector<PictureFeature> &centers,                                                //centers：存放簇中心点的向量
	vector<vector<int>> &cluster_vectors,                                           //cluster_vectors：存放每个簇向量的向量
	double &match_value_sum,
	vector<PictureFeature>* m_ppicture_features)                                    //指向feature的指针

{

	size_t points_count = current_cluster.size();      //所有点的索引向量的大小，也即点的个数
	size_t centers_count = centers.size();             //簇的个数，即K           
	double tmp_distance = 0;                           //用于判断准则用的

	for (size_t i = 0; i < points_count; i++)          //遍历点 找到每个点的最近center
	{
		double smallest_distance = 0;
		int	   near_cluster_id = 0;
		// 计算两个向量的欧式距离
		//std::cout<<centers[0].m_value<<std::endl;
		PictureFeature tmp_feature1;
		tmp_feature1.m_value = (*m_ppicture_features)[current_cluster[i]].m_value;  //表示结构体类型的，指向current_cluster[i]的指针，
		//tmp_feature1.m_value = current_cluster[i].m_value;                        //current_cluster[i]表示点向量的索引
		smallest_distance = CaculateEuclideanDistance(centers[0].m_value,
			tmp_feature1.m_value);
		for (size_t j = 1; j < centers_count; j++)
		{
			//double current_distance = 0;
			PictureFeature tmp_feature2;
			tmp_feature2.m_value = (*m_ppicture_features)[current_cluster[i]].m_value;  //vector容器中存放结构体类型的变量
			//current_cluster[i]整体为索引

			double current_distance = CaculateEuclideanDistance(centers[j].m_value,
				tmp_feature2.m_value);
			if (smallest_distance > current_distance)
			{
				smallest_distance = current_distance;
				near_cluster_id = j;
			}
		}
		cluster_vectors[near_cluster_id].push_back(current_cluster[i]);
		tmp_distance += smallest_distance;   //计算所有点离它最近簇的中心点的距离之和

	}

	match_value_sum = tmp_distance / points_count;//计算所有点离它最近簇的中心点的距离之和的均值
	return 0;
}

// 归完类之后，每个簇内中心点的更新
int UpdateClusterCenter(vector<PictureFeature> &centers,
	vector<vector<int> > &cluster_vectors,                        //每个簇内中心点的更新，也即簇内均值
	float sum,
	vector<PictureFeature>* m_ppicture_features)
{

	int center_counts = cluster_vectors.size();    //存放所有簇的容器大小，K

	for (int i = 0; i < (int)center_counts; i++)  //遍历每个簇，簇的个数也即是K
	{
		for (int j = 0; j < (int)centers[i].m_value.size(); j++)
		{
			centers[i].m_value[j] = 0;                      //center置为0.0
		}

		for (int k = 0; k < (int)cluster_vectors[i].size(); k++)//第i个簇的大小 有多少个点聚到这个center
		{
			//centers[i].m_value += (*m_ppicture_features)[cluster_vectors[i][k]].m_value;
			for (size_t s = 0; s < centers[i].m_value.size(); s++)//第i个簇内的所有点的索引所在的向量的大小
			{
				centers[i].m_value[s] += (*m_ppicture_features)[cluster_vectors[i][k]].m_value[s];
			} //第i个簇内的第K个索引向量的第s个值
		}
		for (size_t r = 0; r < centers[i].m_value.size(); r++)
		{
			centers[i].m_value[r] /= cluster_vectors[i].size();       //对新的center里的每个向量，除以这个簇的点的个数
		}

	}
	return 0;
}

void main()

{

	int sum = 1;
	unsigned int FileCount = 0;                               //统计fts文件数量
	int m_MaxIterateTimes = 15;                          //最多反复次数
	vector<PictureFeature> centers;
	vector<vector<int> > cluster_vectors;
	double old_match_value_sum = 0;
	vector<PictureFeature> m_picture_features;

	//**************************读入文件
	_finddata_t file;
	long longf;
	string tempName;  //_findfirst返回的是long型; long __cdecl _findfirst(const char *, struct _finddata_t *)    
	if ((longf = _findfirst("d:\\files\\test\\*.*", &file)) == -1)
	{
		cout << "Cannot find file!" << endl;
		return;
	}
	do
	{
		tempName = file.name;
		if (tempName[0] == '.')
			continue;

		ifstream ifile;
		ifile.open("d:\\files\\test\\" + tempName, ios_base::binary | ios_base::in);        //以二进制方式打开文件用于读
		if (ifile.fail())
		{
			cout << "Open error!" << endl;
		}
		else
		{
			float temp;
			PictureFeature tmpFeature;
			while (!ifile.eof())                  //m_picture_features是特征,从文件中读取
			{
				ifile.read((char*)&temp, sizeof(float));
				tmpFeature.m_value.clear();
				tmpFeature.m_value.push_back(temp);
			}
			tmpFeature.name = tempName;
			m_picture_features.push_back(tmpFeature);
			ifile.close();
			FileCount++;
		}

	} while (_findnext(longf, &file) == 0);
	_findclose(longf);

	//***************************current_cluster为索引
	vector<int> current_cluster;
	current_cluster.clear();
	for (size_t i = 0; i < FileCount; i++)
	{
		current_cluster.push_back(i);
	}
	vector<PictureFeature> *m_ppicture_features = &m_picture_features;    //指针指向特征头

	//***************************随机产生中心点
	const int N = 2;
	for (int i = 0; i < N; i++)                                  //！存在问题：如果两个随机值一样
	{
		if (FileCount != 0)
		{
			//srand(time(NULL));
			int temprand = rand() % FileCount;
			centers.push_back(m_picture_features[temprand]);
		}
		else break;
	}

	//**************************
	cluster_vectors.resize(N);

	//**************************反复计算
	for (int i = 0; i < m_MaxIterateTimes; i++)
	{
		double current_match_value = 0;
		// 清空每个聚类中心
		ClearClusterVector(cluster_vectors);
		// 计算所有点到中心点的欧式距离

		CaculateMatchValueSum(current_cluster, centers, cluster_vectors, current_match_value, m_ppicture_features);

		if (abs(old_match_value_sum - current_match_value)>  1E-10)
		{
			old_match_value_sum = current_match_value;
		}
		else
		{
			break;
		}
		// 更新聚类中心以及集合
		UpdateClusterCenter(centers, cluster_vectors, sum, m_ppicture_features);
	}
	//******************************保存最好结果




	//******************************更改同一簇的文件名
	for (int i = 0; i < cluster_vectors.size(); i++)
	{
		string filePath = "D:\\files\\";
		char it = 0x31 + i;                       //0x31为“1”
		_mkdir((filePath + it).c_str());

		for (int k = 0; k < cluster_vectors[i].size(); k++)
		{
			rename((filePath + "test\\" + (*m_ppicture_features)[cluster_vectors[i][k]].name).c_str(),
				(filePath + it + "\\" + (*m_ppicture_features)[cluster_vectors[i][k]].name).c_str());
		}
	}

}



