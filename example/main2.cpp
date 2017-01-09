#include<iostream>
#include<math.h>
#include<fstream>
#include<vector>   
#include <io.h>
#include <direct.h>
#include"stdio.h"
#include "ctime"

#include"feature.pb.h"
#include<string>
#include<cstdlib>

using namespace std;

//
struct PictureFeature  
{
	vector<float> m_value;
	string name;
};

// ���ÿ�����������ϵĵ�ļ���
int ClearClusterVector(vector<vector<int> > &cluster_vectors)//cluster_vectors:������������ÿ�ѵ�ģ������������
{
	for (size_t i = 0; i < cluster_vectors.size(); i++)
	{
		cluster_vectors[i].clear();
	}
	return 0;
}

//ŷʽ���룺ÿ�������ƽ��֮�ͣ�ƽ��֮���ٿ�����
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
		sq_distance += tmp_value;  //*****************Ӧ��Ϊ֮�ͣ��ڻ�Ϊ��Ӧ����ĳ˻�֮��
	}

	double distance = sqrt(sq_distance);//ŷʽ����ֵ

	return distance;
}

// �������е㵽�����ĵ��ŷʽ���룬������
int CaculateMatchValueSum(const std::vector<int> &current_cluster,                       //current_cluster�����еĵ�����������
	vector<PictureFeature> &centers,                                                //centers����Ŵ����ĵ������
	vector<vector<int>> &cluster_vectors,                                           //cluster_vectors�����ÿ��������������
	double &match_value_sum, 
	vector<PictureFeature>* m_ppicture_features)                                    //ָ��feature��ָ��
	
{

	size_t points_count = current_cluster.size();      //���е�����������Ĵ�С��Ҳ����ĸ���
	size_t centers_count = centers.size();             //�صĸ�������K           
	double tmp_distance = 0;                           //�����ж�׼���õ�

	for (size_t i = 0; i < points_count; i++)          //������ �ҵ�ÿ��������center
	{
		double smallest_distance = 0;
		int	   near_cluster_id = 0;
		// ��������������ŷʽ����
		//std::cout<<centers[0].m_value<<std::endl;
		PictureFeature tmp_feature1;
		tmp_feature1.m_value = (*m_ppicture_features)[current_cluster[i]].m_value;  //��ʾ�ṹ�����͵ģ�ָ��current_cluster[i]��ָ�룬
		//tmp_feature1.m_value = current_cluster[i].m_value;                        //current_cluster[i]��ʾ������������
		smallest_distance = CaculateEuclideanDistance(centers[0].m_value,
			tmp_feature1.m_value);
		for (size_t j = 1; j < centers_count; j++)
		{
			//double current_distance = 0;
			PictureFeature tmp_feature2;
			tmp_feature2.m_value = (*m_ppicture_features)[current_cluster[i]].m_value;  //vector�����д�Žṹ�����͵ı���
			//current_cluster[i]����Ϊ����

			double current_distance = CaculateEuclideanDistance(centers[j].m_value, 
				tmp_feature2.m_value);
			if (smallest_distance > current_distance)
			{
				smallest_distance = current_distance;
				near_cluster_id = j;
			}
		}
		cluster_vectors[near_cluster_id].push_back(current_cluster[i]);
		tmp_distance += smallest_distance;   //�������е���������ص����ĵ�ľ���֮��

	}

	match_value_sum = tmp_distance / points_count;//�������е���������ص����ĵ�ľ���֮�͵ľ�ֵ
	return 0;
}

// ������֮��ÿ���������ĵ�ĸ���
int UpdateClusterCenter(vector<PictureFeature> &centers, 
	vector<vector<int> > &cluster_vectors,                        //ÿ���������ĵ�ĸ��£�Ҳ�����ھ�ֵ
	float sum,
	vector<PictureFeature>* m_ppicture_features)
{

	int center_counts = cluster_vectors.size();    //������дص�������С��K

	for (int i = 0; i < (int)center_counts; i++)  //����ÿ���أ��صĸ���Ҳ����K
	{
		for (int j = 0; j < (int)centers[i].m_value.size(); j++)
		{
			centers[i].m_value[j] = 0;                      //center��Ϊ0.0
		}

		for (int k = 0; k < (int)cluster_vectors[i].size(); k++)//��i���صĴ�С �ж��ٸ���۵����center
		{
			//centers[i].m_value += (*m_ppicture_features)[cluster_vectors[i][k]].m_value;
			for (size_t s = 0; s < centers[i].m_value.size(); s++)//��i�����ڵ����е���������ڵ������Ĵ�С
			{
				centers[i].m_value[s] += (*m_ppicture_features)[cluster_vectors[i][k]].m_value[s];
			} //��i�����ڵĵ�K�����������ĵ�s��ֵ
		}
		for (size_t r = 0; r < centers[i].m_value.size(); r++)
		{
			centers[i].m_value[r] /= cluster_vectors[i].size();       //���µ�center���ÿ����������������صĵ�ĸ���
		}

	}
	return 0;
}

void main()
{

	int sum = 1;
	unsigned int FileCount = 0;                               //ͳ��fts�ļ�����
	int m_MaxIterateTimes = 10;                            //��෴������
	vector<PictureFeature> centers;
	vector<vector<int> > cluster_vectors;
	double old_match_value_sum = 0;
	vector<PictureFeature> m_picture_features;

	//**************************�����ļ�
	_finddata_t file;
	long longf;
	string tempName;  //_findfirst���ص���long��; long __cdecl _findfirst(const char *, struct _finddata_t *)    
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
		string readfromfilestring;
		int64_t flenth = -1;
		fstream input_file("d:\\files\\test\\" + tempName, ios::in | ios::binary);
		float temp;
		input_file.read((char*)&flenth, sizeof(int64_t));//�����ļ���С����ֵ��flenth����int64_t�ķ�ʽ����
		readfromfilestring.resize(flenth);//Ϊreadfromfilestring�����flenthһ����С�Ŀռ�ֵ
		input_file.read(const_cast<char*> (readfromfilestring.data()), flenth);//��ȡ�ļ���readfromfilestring���ļ�����Ϊflenth
		Facefeature facefeature1;
		PictureFeature tmpFeature;
		facefeature1.ParseFromString(readfromfilestring);//�����ݴ�readfromfilestring�ж����ŵ�facefeature1(�õ�protobuff��ĺ�������������ת��)
		input_file.close();

		//cout << temp << endl;
		tmpFeature.m_value.clear();
		for (int j = 0; j < facefeature1.dwfacefeatsize(); j++)
		{
			tmpFeature.m_value.push_back(facefeature1.pffacefeat(j));
			//cout << "j= " << j << ": " << facefeature1.pffacefeat(j) << endl;
		}
		tmpFeature.name = tempName;
		m_picture_features.push_back(tmpFeature);
		FileCount++;

	} while (_findnext(longf, &file) == 0);
	_findclose(longf);

	//***************************current_clusterΪ����
	vector<int> current_cluster;
	current_cluster.clear();
	for (size_t i = 0; i < FileCount; i++)
	{
		current_cluster.push_back(i);
	}
	vector<PictureFeature> *m_ppicture_features = &m_picture_features;    //ָ��ָ������ͷ


	//**************************
	int N = 2;
	double min_match_value=0xFFFF;
	vector<vector<int> > best_cluster_vectors;
	cluster_vectors.resize(N);
	//**************************��������
	for (int i = 0; i < m_MaxIterateTimes; i++)
	{
		double current_match_value = 0;
		// ���ÿ����������

		for (int j = 0; j < N; j++)                                  //����������ĵ�
		{
			if (FileCount != 0)
			{
				//srand(time(NULL));
				int temprand = rand() % FileCount;
				centers.push_back(m_picture_features[temprand]);
			}
			else break;
		}

		ClearClusterVector(cluster_vectors);
		// �������е㵽���ĵ��ŷʽ����
		CaculateMatchValueSum(current_cluster, centers, cluster_vectors, current_match_value, m_ppicture_features);

		if (abs(old_match_value_sum - current_match_value) > 1E-10)      //����
		{
			old_match_value_sum = current_match_value;
		}
		else
		{
			break;
		}
		// ���¾��������Լ�����
		UpdateClusterCenter(centers, cluster_vectors, sum, m_ppicture_features);

		if (old_match_value_sum < min_match_value)
		{
			min_match_value = old_match_value_sum;
			best_cluster_vectors = cluster_vectors;
		}
		centers.erase(centers.begin(),centers.end()); //ɾ��֮ǰ��centers		
	}
	//******************************����ͬһ�ص��ļ���
	for (int i = 0; i < best_cluster_vectors.size(); i++)
	{
		string filePath = "D:\\files\\";
		char it = 0x31 + i;                       //0x31Ϊ��1��
		_mkdir((filePath + it).c_str());

		for (int k = 0; k < best_cluster_vectors[i].size(); k++)
		{
			rename((filePath + "test\\" + (*m_ppicture_features)[best_cluster_vectors[i][k]].name).c_str(),
				(filePath + it + "\\" + (*m_ppicture_features)[best_cluster_vectors[i][k]].name).c_str());
		}
	}

}



