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
int CaculateMatchValueSum(const vector<int> &current_cluster,                       //current_cluster�����еĵ�����������
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
		PictureFeature tmp_feature1;
		tmp_feature1.m_value = (*m_ppicture_features)[current_cluster[i]].m_value;  //��ʾ�ṹ�����͵ģ�ָ��current_cluster[i]��ָ�룬
		smallest_distance = CaculateEuclideanDistance(centers[0].m_value, tmp_feature1.m_value);
		for (size_t j = 0; j < centers_count; j++)
		{
			PictureFeature tmp_feature2;
			tmp_feature2.m_value = (*m_ppicture_features)[current_cluster[i]].m_value;  //vector�����д�Žṹ�����͵ı���
			//current_cluster[i]����Ϊ����
			double current_distance = CaculateEuclideanDistance(centers[j].m_value, tmp_feature2.m_value);
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

//���ļ�bi
void readbi(string Path, string tempName, unsigned int &FileCount, _finddata_t &file, long &longf, vector<PictureFeature> &m_picture_features)
{
		string readfromfilestring;
		int64_t flenth = -1;
		fstream input_file(Path + tempName, ios::in | ios::binary);
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

	
}
//���ļ�fts
void readfts(string Path, string tempName, unsigned int &FileCount, _finddata_t &file, long &longf, vector<PictureFeature> &m_picture_features)
{
	ifstream ifile;
	ifile.open(Path + tempName, ios_base::binary | ios_base::in);        //�Զ����Ʒ�ʽ���ļ����ڶ�
	if (ifile.fail())
	{
		cout << "Open error!" << endl;
	}
	else
	{

		PictureFeature tmpFeature;
		tmpFeature.m_value.clear();

		for (size_t i = 0; i < 512; i++)
		{
			float temp;
			ifile.read((char*)&temp, sizeof(float));
			tmpFeature.m_value.push_back(temp);            //m_picture_features������,���ļ��ж�ȡ
		}                 

		tmpFeature.name = tempName;
		m_picture_features.push_back(tmpFeature);
		ifile.close();
		FileCount++;
	}

}

void uncertainN(string Path, string FeatName, int &index)
{
	int NN = 512;
	//**************************���޸ı���
	//�޸�·��ʱ������1.readbi()  2.readfts()  3.������
	int m_MaxIterateTimes = 10;             //������ʱ����෴������
	float div = 1.60;                       //��ֵ����������أ�������0��
	
	//**************************

	unsigned int FileCount = 0;                             //ͳ���ļ�����
	vector<vector<PictureFeature>> centers;
	vector<vector<int>> cluster_vectors;
	double old_match_value_sum = 0;
	vector<PictureFeature> m_picture_features;

	//**************************�����ļ�
	_finddata_t file;
	long longf;    
	if ((longf = _findfirst((Path+"*.*").c_str(), &file)) == -1)
	{
		cout << "Cannot find file!" << endl;
		return;
	}
		do
		{
			char* fileName = file.name;		
			if (fileName[0] == '.')
				continue;
			for (int i = 1; i < 50; i++)
			{
				if (fileName[i] == '.')
				{
					if (fileName[i + 1] == 'f')
					{
						readfts(Path,fileName, FileCount, file, longf, m_picture_features);
						break;
					}
					else if (fileName[i + 1] == 'b')
					{
						readbi(Path,fileName, FileCount, file, longf, m_picture_features);
						break;
					}
				}
			}
		} while (_findnext(longf, &file) == 0);
	_findclose(longf);
	//***********************************

	int N = FileCount * 3 / 4;                       //�ײ�centers����Ŀ��ҪС���ļ�����
	int MaxCenters = FileCount/2;                    //centers�ϲ��Ĵ�������Ҫ���ڿ��ܵ�����centers��������Ϊ���ļ���/2
	
	//***************************current_clusterΪ����
	vector<int> current_cluster;
	current_cluster.clear();
	for (size_t i = 0; i < FileCount; i++)
	{
		current_cluster.push_back(i);
	}
	vector<PictureFeature> *m_ppicture_features = &m_picture_features;    //ָ��ָ������ͷ

	//*******************************����������ĵ�
	vector<PictureFeature > centers0;                    //centers���ݴ�
	vector<int> temprand(N,-1);                          //���ĵ��ݴ棬ȫ����ʼ��Ϊ-1    
	for (int j = 0; j < N; )                                  
	{
		if (FileCount != 0)
		{
			int randn = rand() % FileCount;
			for (int i = 0; (i < N); i++)                 //��ֹ���ĵ��ظ�
			{
				if (randn == temprand[i])
					break;
				else
				{
					if (i == N - 1)
						temprand[j] = randn;				
				}
			}
			if (temprand[j] != -1)
			{
				centers0.push_back(m_picture_features[temprand[j]]);
				j++;
			}
		}
		else
		{
			cout << "���ļ�" << endl;
			break;
		}
	}
	centers.push_back(centers0);
	cluster_vectors.resize(N);

	//**************************��������
	for (int i = 0; i < m_MaxIterateTimes; i++)
	{
		double current_match_value = 0;
		// ���ÿ����������
		ClearClusterVector(cluster_vectors);
		// �������е㵽���ĵ��ŷʽ����
		CaculateMatchValueSum(current_cluster, centers[0], cluster_vectors, current_match_value, m_ppicture_features);

		if (abs(old_match_value_sum - current_match_value) > 1E-10)      //����
		{
			old_match_value_sum = current_match_value;
		}
		else
		{
			break;
		}
		// ���¾��������Լ�����
		UpdateClusterCenter(centers[0], cluster_vectors, m_ppicture_features);
	}

	//********************************
	double SumDistance = 0;
	//vector<double> Di;
	int centers_final;                      //����centers�ĸ���           
	vector<PictureFeature > centersPick;

	//��������centers����������
	int centersCount = centers[0].size();
	for (int i = 0; i < centersCount; i++)
	{
		for (int j = 0; j < centersCount; j++)
		{
			double sumtemp;
			sumtemp = CaculateEuclideanDistance(centers[0][i].m_value, centers[0][j].m_value);               
			SumDistance += sumtemp;
		}
	}
	double Evendistance = (SumDistance / (centersCount*(centersCount-1)))/div;    //��ֵ

	for (int k = 0; k < MaxCenters; k++)            //centers�ظ����ٴ�
	{
		centersCount = centers[k].size();
		centers0.erase(centers0.begin(), centers0.end());
		centersPick.erase(centersPick.begin(), centersPick.end());
		centers0.push_back(centers[k][0]);                        //��ʼ��centers0
		int j = 0;					
		for (int m = 1; m < centersCount; m++)                   //����centers                  
		{
			if (CaculateEuclideanDistance(centers0[j].m_value, centers[k][m].m_value) < Evendistance)    //С�ھ�ֵ�����centers���ϲ�
			{
				PictureFeature tmpFeature;
				for (int r = 0; r < NN; r++)                            //����ɹ�Ϊһ�������centers�ľ�ֵΪ�µ�center
				{
					float tempValue = (centers0[j].m_value[r] + centers[k][m].m_value[r]) / 2;      
					tmpFeature.m_value.push_back(tempValue);
				}
				centers0.push_back(tmpFeature);                         //�µ�center�浽centers0�У����µ�center��centers��һ���Ƚ�
				j++;
			}
			else                                               //���������ͬһ�࣬���ݴ浽centersPick��
			{
				centersPick.push_back(centers[k][m]);
			}
		}
		centersPick.push_back(centers0[j]);                //centers0���ļ�����Ҳ�浽Pick
		centers.push_back(centersPick);                   //��Ϊ��һ��centers����centers[1]
	}//�ظ��ϲ�MaxCenters�Σ�����󲻱�
	centers_final = centers[MaxCenters].size();             //���centers�ĸ���

	//******************************�Ժϲ���ɵĵ㣬���е����¹���
	vector<vector<int>> cluster_vectors_final;
	double match_value_sum_final=0;
	cluster_vectors_final.resize(centers_final);
	CaculateMatchValueSum(current_cluster, centers[MaxCenters], cluster_vectors_final, match_value_sum_final, m_ppicture_features);

	//*******************************����ͬһ�ص��ļ���
	for (int i = 0; i < cluster_vectors_final.size(); i++)
	{
		string filePath = "D:\\files\\";
		char ge = 0x30 + i % 10;                //0x30Ϊ��0��
		char shi = 0x30 + i / 10;
		_mkdir((filePath + shi+ ge).c_str());
		for (int k = 0; k < cluster_vectors_final[i].size(); k++)
		{
			//string oldName = filePath + "test\\" + (*m_ppicture_features)[cluster_vectors_final[i][k]].name;
			//string newName = filePath + shi+ ge + "\\" + (*m_ppicture_features)[cluster_vectors_final[i][k]].name;
			rename((filePath + "test\\" + (*m_ppicture_features)[cluster_vectors_final[i][k]].name).c_str(),
				(filePath + shi+ ge + "\\" + (*m_ppicture_features)[cluster_vectors_final[i][k]].name).c_str());
		}
	}
	//********************************�ҵ�ָ��feature����
	for (int i = 0; i < cluster_vectors_final.size(); i++)
	{
		for (int k = 0; k < cluster_vectors_final[i].size(); k++)
		{
			if (FeatName == (*m_ppicture_features)[cluster_vectors_final[i][k]].name)
			{
				index = i;
			}
		}
	}
	if (index == -1)
		cout << "δ�ҵ�ָ���ļ���" << endl;
}

void certainN(string Path, string FeatName, int N, int &index)
{
	unsigned int FileCount = 0;                               //ͳ���ļ�����
	int m_MaxIterateTimes = 50;                               //��෴������
	vector<PictureFeature> centers;
	vector<vector<int> > cluster_vectors;
	double old_match_value_sum = 0;
	vector<PictureFeature> m_picture_features;

	//**************************�����ļ�
	_finddata_t file;
	long longf;
	if ((longf = _findfirst((Path + "*.*").c_str(), &file)) == -1)
	{
		cout << "Cannot find file!" << endl;
		return;
	}
	do
	{
		char* fileName = file.name;
		if (fileName[0] == '.')
			continue;
		for (int i = 1; i < 50; i++)
		{
			if (fileName[i] == '.')
			{
				if (fileName[i + 1] == 'f')
				{
					readfts(Path, fileName, FileCount, file, longf, m_picture_features);
					break;
				}
				else if (fileName[i + 1] == 'b')
				{
					readbi(Path, fileName, FileCount, file, longf, m_picture_features);
					break;
				}
			}
		}
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

	//***************************����������ĵ�

	vector<int> temprand(N, -1);                          //���ĵ��ݴ棬ȫ����ʼ��Ϊ-1    
	for (int j = 0; j < N;)
	{
		if (FileCount != 0)
		{
			int randn = rand() % FileCount;
			for (int i = 0; (i < N); i++)                 //��ֹ���ĵ��ظ�
			{
				if (randn == temprand[i])
					break;
				else
				{
					if (i == N - 1)
						temprand[j] = randn;
				}
			}
			if (temprand[j] != -1)
			{
				centers.push_back(m_picture_features[temprand[j]]);
				j++;
			}
		}
		else
		{
			cout << "���ļ�" << endl;
			break;
		}
	}
	cluster_vectors.resize(N);

	//**************************��������
	for (int i = 0; i < m_MaxIterateTimes; i++)
	{
		double current_match_value = 0;
		// ���ÿ����������
		ClearClusterVector(cluster_vectors);
		// �������е㵽���ĵ��ŷʽ����
		CaculateMatchValueSum(current_cluster, centers, cluster_vectors, current_match_value, m_ppicture_features);
		if (abs(old_match_value_sum - current_match_value)>  1E-10)
		{
			old_match_value_sum = current_match_value;
		}
		else
		{
			break;
		}
		// ���¾��������Լ�����
		UpdateClusterCenter(centers, cluster_vectors, m_ppicture_features);
	}
	//******************************����ͬһ�ص��ļ���
	for (int i = 0; i < cluster_vectors.size(); i++)
	{
		string filePath = "D:\\files\\";
		char ge = 0x30 + i % 10;                //0x30Ϊ��0��
		char shi = 0x30 + i / 10;
		_mkdir((filePath + shi + ge).c_str());
		for (int k = 0; k < cluster_vectors[i].size(); k++)
		{
			rename((filePath + "test\\" + (*m_ppicture_features)[cluster_vectors[i][k]].name).c_str(),
				(filePath + shi + ge + "\\" + (*m_ppicture_features)[cluster_vectors[i][k]].name).c_str());
		}
	}
	//********************************�ҵ�ָ��feature����
	for (int i = 0; i < cluster_vectors.size(); i++)
	{
		for (int k = 0; k < cluster_vectors[i].size(); k++)
		{
			if (FeatName == (*m_ppicture_features)[cluster_vectors[i][k]].name)
			{
				index = i;
			}
		}
	}
	if (index == -1)
		cout << "δ�ҵ�ָ���ļ���" << endl;
}

void main()
{
	string Path = "D:\\files\\test\\"; 
	string FeatName;                //Ҫ���ҵ������ļ���������׺��
	int index=-1;							                         //��ı��
	int N;                                                           //ȷ���Ĵظ���
	int tag=-1;
	cout << "���룺1-ȷ��������  2-��ȷ��������" << endl;
	cin >> tag;
	cin.get();
	if (tag == 1)
	{
		cout << "���������";
		cin >> N;
		cin.get();
		cout << "����Ҫ���ҵ��ļ�����";
		getline(cin,FeatName);
		certainN(Path, FeatName, N, index);
	}
	else if (tag == 2)
	{
		cout << "����Ҫ���ҵ��ļ�����";
		getline(cin, FeatName);
		uncertainN(Path, FeatName, index);
	}
	else
	{
		cout << "�������" << endl;
	}
	cout << FeatName << " �������ͣ�" << index << endl;
		
}

	



