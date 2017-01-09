#include<iostream>
#include<vector>
#include<math.h>
//#include <queue>
using namespace std;

struct PictureFeature
{
	std::vector<float> m_value;
};
// ���ÿ�����������ϵĵ�ļ���
int ClearClusterVector(std::vector<std::vector<int> > &cluster_vectors)//cluster_vectors:������������ÿ�ѵ�ģ������������
{
	for (std::size_t i = 0; i < cluster_vectors.size(); i++)
	{
		cluster_vectors[i].clear();
	}
	return 0;
}

double CaculateEuclideanDistance(std::vector<float> &v1, std::vector<float> &v2)//ŷʽ���룺ÿ�������ƽ��֮�ͣ�ƽ��֮���ٿ�����
{

	if (v2.size() != v1.size())
	{
		std::cout << "error" << "\n";
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
		sq_distance += tmp_value;  //1*****************Ӧ��Ϊ֮�ͣ��ڻ�Ϊ��Ӧ����ĳ˻�֮��
	}

	double distance = sqrt(sq_distance);//ŷʽ����ֵ

	return distance;
}
//rent_cluster�����еĵ�����������,������Ϊ����

// �������е㵽�����ĵ��ŷʽ���룬������
int CaculateMatchValueSum(const std::vector<int> &current_cluster,                       //current_cluster�����еĵ�����������
std::vector<PictureFeature> &centers,                                                //centers����Ŵ����ĵ������
std::vector<std::vector<int> > &cluster_vectors,                                     //cluster_vectors�����ÿ��������������
double &match_value_sum, std::vector<PictureFeature>* m_ppicture_features)          //vector�����д�Žṹ�����͵ı���

{

	std::size_t points_count = current_cluster.size();      //���е�����������Ĵ�С��Ҳ����ĸ���
	std::size_t centers_count = centers.size();            //�صĸ�������K           
	double tmp_distance = 0;                              //�����ж�׼���õ�

	for (std::size_t i = 0; i < points_count; i++)      //������
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
		for (std::size_t j = 1; j < centers_count; j++)
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

double inner_prod(std::vector<double> v1, std::vector<double> v2)//�������ڻ���.������������Ӧ����˻�֮��
{
	double value_sum = 0.0;
	for (size_t i = 0; i < v1.size(); i++)
	{
		double value = v1[i] * v2[i];
		value_sum += value;
	}
	return value_sum;
}

// ������֮��ÿ���������ĵ�ĸ���
int UpdateClusterCenter(std::vector<PictureFeature> &centers, std::vector<std::vector<int> > &cluster_vectors, //ÿ���������ĵ�ĸ��£�Ҳ�����ھ�ֵ
float sum, std::vector<PictureFeature>* m_ppicture_features)
{
	int center_counts = cluster_vectors.size();    //������дص�������С

	for (int i = 0; i < (int)center_counts; i++)  //����ÿ���أ��صĸ���Ҳ����K
	{


		for (int j = 0; j < (int)centers[i].m_value.size(); j++)//�صĸ���Ҳ����K,�����ĵ�ĸ���
		{
			centers[i].m_value[j] = 0;
		}

		for (int k = 0; k < (int)cluster_vectors[i].size(); k++)//��i���صĴ�С
		{
			//centers[i].m_value += (*m_ppicture_features)[cluster_vectors[i][k]].m_value;
			for (size_t s = 0; s < centers[i].m_value.size(); s++)//��i�����ڵ����е���������ڵ������Ĵ�С
			{
				centers[i].m_value[s] += (*m_ppicture_features)[cluster_vectors[i][k]].m_value[s];
			} //��i�����ڵĵ�K�����������ĵ�s��ֵ
		}
		for (size_t r = 0; r < centers[i].m_value.size(); r++)
		{
			centers[i].m_value[r] /= cluster_vectors[i].size();
		}


	}
	return 0;
}


int main()
{
	int sum = 1;
	int m_MaxIterateTimes = 10;

	std::vector<PictureFeature> centers;
	std::vector<std::vector<int> > cluster_vectors;
	double old_match_value_sum = 0;
	std::vector<PictureFeature> m_picture_features;

	for (size_t i = 1; i < 11; i++)                  //m_picture_features������
	{
		PictureFeature tmpFeature;
		tmpFeature.m_value.clear();
		tmpFeature.m_value.push_back(i);
		m_picture_features.push_back(tmpFeature);
	}
	for (size_t i = 1; i < 3; i++)                    //centers�����ѡ�������ĵ�
	{
		PictureFeature tmpFeature;
		tmpFeature.m_value.clear();
		tmpFeature.m_value.push_back(i);
		centers.push_back(tmpFeature);
	}
	cluster_vectors.resize(2);
	std::vector<int> current_cluster;
	current_cluster.clear();
	for (size_t i = 0; i < 10; i++)
	{
		current_cluster.push_back(i);                  //current_cluster������������       
	}

	std::vector<PictureFeature> *m_ppicture_features = &m_picture_features;
	for (int i = 0; i < m_MaxIterateTimes; i++)
	{
		double current_match_value = 0;
		// ���ÿ���������ĵ��ϵ�ļ���
		ClearClusterVector(cluster_vectors);

		// �������е㵽���ĵ��ŷʽ����
		CaculateMatchValueSum(current_cluster, centers, cluster_vectors, current_match_value, m_ppicture_features);


		if (std::abs(old_match_value_sum - current_match_value)>  1E-10)
		{
			old_match_value_sum = current_match_value;
		}
		else
		{
			break;
		}

		// ���¾��������Լ�����
		UpdateClusterCenter(centers, cluster_vectors, sum, m_ppicture_features);


	}
	system("pause");
	return 0;
}
