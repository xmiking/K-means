#include<iostream>
#include<vector>
#include<math.h>
//#include <queue>
using namespace std;

struct PictureFeature
{
	std::vector<float> m_value;
};
// 清空每个聚类中心上的点的集合
int ClearClusterVector(std::vector<std::vector<int> > &cluster_vectors)//cluster_vectors:外层向量，存放每堆点的，向量外的向量
{
	for (std::size_t i = 0; i < cluster_vectors.size(); i++)
	{
		cluster_vectors[i].clear();
	}
	return 0;
}

double CaculateEuclideanDistance(std::vector<float> &v1, std::vector<float> &v2)//欧式距离：每个距离的平方之和，平方之和再开根号
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
		sq_distance += tmp_value;  //1*****************应该为之和：内积为对应坐标的乘积之和
	}

	double distance = sqrt(sq_distance);//欧式距离值

	return distance;
}
//rent_cluster：所有的点向量的索引,且索引为向量

// 计算所有点到簇中心点的欧式距离，并归类
int CaculateMatchValueSum(const std::vector<int> &current_cluster,                       //current_cluster：所有的点向量的索引
std::vector<PictureFeature> &centers,                                                //centers：存放簇中心点的向量
std::vector<std::vector<int> > &cluster_vectors,                                     //cluster_vectors：存放每个簇向量的向量
double &match_value_sum, std::vector<PictureFeature>* m_ppicture_features)          //vector容器中存放结构体类型的变量

{

	std::size_t points_count = current_cluster.size();      //所有点的索引向量的大小，也即点的个数
	std::size_t centers_count = centers.size();            //簇的个数，即K           
	double tmp_distance = 0;                              //用于判断准则用的

	for (std::size_t i = 0; i < points_count; i++)      //遍历点
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
		for (std::size_t j = 1; j < centers_count; j++)
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

double inner_prod(std::vector<double> v1, std::vector<double> v2)//向量的内积（.）：两向量对应坐标乘积之和
{
	double value_sum = 0.0;
	for (size_t i = 0; i < v1.size(); i++)
	{
		double value = v1[i] * v2[i];
		value_sum += value;
	}
	return value_sum;
}

// 归完类之后，每个簇内中心点的更新
int UpdateClusterCenter(std::vector<PictureFeature> &centers, std::vector<std::vector<int> > &cluster_vectors, //每个簇内中心点的更新，也即簇内均值
float sum, std::vector<PictureFeature>* m_ppicture_features)
{
	int center_counts = cluster_vectors.size();    //存放所有簇的容器大小

	for (int i = 0; i < (int)center_counts; i++)  //遍历每个簇，簇的个数也即是K
	{


		for (int j = 0; j < (int)centers[i].m_value.size(); j++)//簇的个数也即是K,簇中心点的个数
		{
			centers[i].m_value[j] = 0;
		}

		for (int k = 0; k < (int)cluster_vectors[i].size(); k++)//第i个簇的大小
		{
			//centers[i].m_value += (*m_ppicture_features)[cluster_vectors[i][k]].m_value;
			for (size_t s = 0; s < centers[i].m_value.size(); s++)//第i个簇内的所有点的索引所在的向量的大小
			{
				centers[i].m_value[s] += (*m_ppicture_features)[cluster_vectors[i][k]].m_value[s];
			} //第i个簇内的第K个索引向量的第s个值
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

	for (size_t i = 1; i < 11; i++)                  //m_picture_features是特征
	{
		PictureFeature tmpFeature;
		tmpFeature.m_value.clear();
		tmpFeature.m_value.push_back(i);
		m_picture_features.push_back(tmpFeature);
	}
	for (size_t i = 1; i < 3; i++)                    //centers是随机选出的中心点
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
		current_cluster.push_back(i);                  //current_cluster是特征点的序号       
	}

	std::vector<PictureFeature> *m_ppicture_features = &m_picture_features;
	for (int i = 0; i < m_MaxIterateTimes; i++)
	{
		double current_match_value = 0;
		// 清空每个聚类中心的上点的集合
		ClearClusterVector(cluster_vectors);

		// 计算所有点到中心点的欧式距离
		CaculateMatchValueSum(current_cluster, centers, cluster_vectors, current_match_value, m_ppicture_features);


		if (std::abs(old_match_value_sum - current_match_value)>  1E-10)
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
	system("pause");
	return 0;
}
