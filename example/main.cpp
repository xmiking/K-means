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
int CaculateMatchValueSum(const vector<int> &current_cluster,                       //current_cluster：所有的点向量的索引
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
		PictureFeature tmp_feature1;
		tmp_feature1.m_value = (*m_ppicture_features)[current_cluster[i]].m_value;  //表示结构体类型的，指向current_cluster[i]的指针，
		smallest_distance = CaculateEuclideanDistance(centers[0].m_value, tmp_feature1.m_value);
		for (size_t j = 0; j < centers_count; j++)
		{
			PictureFeature tmp_feature2;
			tmp_feature2.m_value = (*m_ppicture_features)[current_cluster[i]].m_value;  //vector容器中存放结构体类型的变量
			//current_cluster[i]整体为索引
			double current_distance = CaculateEuclideanDistance(centers[j].m_value, tmp_feature2.m_value);
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

//读文件bi
void readbi(string Path, string tempName, unsigned int &FileCount, _finddata_t &file, long &longf, vector<PictureFeature> &m_picture_features)
{
		string readfromfilestring;
		int64_t flenth = -1;
		fstream input_file(Path + tempName, ios::in | ios::binary);
		float temp;
		input_file.read((char*)&flenth, sizeof(int64_t));//读出文件大小并赋值给flenth（以int64_t的方式读）
		readfromfilestring.resize(flenth);//为readfromfilestring申请和flenth一样大小的空间值
		input_file.read(const_cast<char*> (readfromfilestring.data()), flenth);//读取文件到readfromfilestring，文件长度为flenth
		Facefeature facefeature1;
		PictureFeature tmpFeature;
		facefeature1.ParseFromString(readfromfilestring);//把数据从readfromfilestring中读出放到facefeature1(用到protobuff库的函数，数据类型转换)
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
//读文件fts
void readfts(string Path, string tempName, unsigned int &FileCount, _finddata_t &file, long &longf, vector<PictureFeature> &m_picture_features)
{
	ifstream ifile;
	ifile.open(Path + tempName, ios_base::binary | ios_base::in);        //以二进制方式打开文件用于读
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
			tmpFeature.m_value.push_back(temp);            //m_picture_features是特征,从文件中读取
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
	//**************************可修改变量
	//修改路径时三处：1.readbi()  2.readfts()  3.主函数
	int m_MaxIterateTimes = 10;             //收敛的时候最多反复次数
	float div = 1.60;                       //阈值调整（负相关）（大于0）
	
	//**************************

	unsigned int FileCount = 0;                             //统计文件数量
	vector<vector<PictureFeature>> centers;
	vector<vector<int>> cluster_vectors;
	double old_match_value_sum = 0;
	vector<PictureFeature> m_picture_features;

	//**************************读入文件
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

	int N = FileCount * 3 / 4;                       //底层centers的数目（要小于文件数）
	int MaxCenters = FileCount/2;                    //centers合并的次数，需要大于可能的最终centers数，初设为：文件数/2
	
	//***************************current_cluster为索引
	vector<int> current_cluster;
	current_cluster.clear();
	for (size_t i = 0; i < FileCount; i++)
	{
		current_cluster.push_back(i);
	}
	vector<PictureFeature> *m_ppicture_features = &m_picture_features;    //指针指向特征头

	//*******************************产生随机中心点
	vector<PictureFeature > centers0;                    //centers的暂存
	vector<int> temprand(N,-1);                          //中心点暂存，全部初始化为-1    
	for (int j = 0; j < N; )                                  
	{
		if (FileCount != 0)
		{
			int randn = rand() % FileCount;
			for (int i = 0; (i < N); i++)                 //防止中心点重复
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
			cout << "无文件" << endl;
			break;
		}
	}
	centers.push_back(centers0);
	cluster_vectors.resize(N);

	//**************************收敛反复
	for (int i = 0; i < m_MaxIterateTimes; i++)
	{
		double current_match_value = 0;
		// 清空每个聚类中心
		ClearClusterVector(cluster_vectors);
		// 计算所有点到中心点的欧式距离
		CaculateMatchValueSum(current_cluster, centers[0], cluster_vectors, current_match_value, m_ppicture_features);

		if (abs(old_match_value_sum - current_match_value) > 1E-10)      //收敛
		{
			old_match_value_sum = current_match_value;
		}
		else
		{
			break;
		}
		// 更新聚类中心以及集合
		UpdateClusterCenter(centers[0], cluster_vectors, m_ppicture_features);
	}

	//********************************
	double SumDistance = 0;
	//vector<double> Di;
	int centers_final;                      //最终centers的个数           
	vector<PictureFeature > centersPick;

	//计算所有centers的两两距离
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
	double Evendistance = (SumDistance / (centersCount*(centersCount-1)))/div;    //均值

	for (int k = 0; k < MaxCenters; k++)            //centers重复多少次
	{
		centersCount = centers[k].size();
		centers0.erase(centers0.begin(), centers0.end());
		centersPick.erase(centersPick.begin(), centersPick.end());
		centers0.push_back(centers[k][0]);                        //初始化centers0
		int j = 0;					
		for (int m = 1; m < centersCount; m++)                   //遍历centers                  
		{
			if (CaculateEuclideanDistance(centers0[j].m_value, centers[k][m].m_value) < Evendistance)    //小于均值距离的centers，合并
			{
				PictureFeature tmpFeature;
				for (int r = 0; r < NN; r++)                            //求出可归为一类的两个centers的均值为新的center
				{
					float tempValue = (centers0[j].m_value[r] + centers[k][m].m_value[r]) / 2;      
					tmpFeature.m_value.push_back(tempValue);
				}
				centers0.push_back(tmpFeature);                         //新的center存到centers0中，用新的center和centers下一个比较
				j++;
			}
			else                                               //如果不属于同一类，就暂存到centersPick中
			{
				centersPick.push_back(centers[k][m]);
			}
		}
		centersPick.push_back(centers0[j]);                //centers0最后的计算结果也存到Pick
		centers.push_back(centersPick);                   //作为下一层centers，即centers[1]
	}//重复合并MaxCenters次，到最后不变
	centers_final = centers[MaxCenters].size();             //最后centers的个数

	//******************************对合并完成的点，所有点重新归类
	vector<vector<int>> cluster_vectors_final;
	double match_value_sum_final=0;
	cluster_vectors_final.resize(centers_final);
	CaculateMatchValueSum(current_cluster, centers[MaxCenters], cluster_vectors_final, match_value_sum_final, m_ppicture_features);

	//*******************************更改同一簇的文件名
	for (int i = 0; i < cluster_vectors_final.size(); i++)
	{
		string filePath = "D:\\files\\";
		char ge = 0x30 + i % 10;                //0x30为“0”
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
	//********************************找到指定feature的类
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
		cout << "未找到指定文件！" << endl;
}

void certainN(string Path, string FeatName, int N, int &index)
{
	unsigned int FileCount = 0;                               //统计文件数量
	int m_MaxIterateTimes = 50;                               //最多反复次数
	vector<PictureFeature> centers;
	vector<vector<int> > cluster_vectors;
	double old_match_value_sum = 0;
	vector<PictureFeature> m_picture_features;

	//**************************读入文件
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

	//***************************current_cluster为索引
	vector<int> current_cluster;
	current_cluster.clear();
	for (size_t i = 0; i < FileCount; i++)
	{
		current_cluster.push_back(i);
	}
	vector<PictureFeature> *m_ppicture_features = &m_picture_features;    //指针指向特征头

	//***************************随机产生中心点

	vector<int> temprand(N, -1);                          //中心点暂存，全部初始化为-1    
	for (int j = 0; j < N;)
	{
		if (FileCount != 0)
		{
			int randn = rand() % FileCount;
			for (int i = 0; (i < N); i++)                 //防止中心点重复
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
			cout << "无文件" << endl;
			break;
		}
	}
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
		UpdateClusterCenter(centers, cluster_vectors, m_ppicture_features);
	}
	//******************************更改同一簇的文件名
	for (int i = 0; i < cluster_vectors.size(); i++)
	{
		string filePath = "D:\\files\\";
		char ge = 0x30 + i % 10;                //0x30为“0”
		char shi = 0x30 + i / 10;
		_mkdir((filePath + shi + ge).c_str());
		for (int k = 0; k < cluster_vectors[i].size(); k++)
		{
			rename((filePath + "test\\" + (*m_ppicture_features)[cluster_vectors[i][k]].name).c_str(),
				(filePath + shi + ge + "\\" + (*m_ppicture_features)[cluster_vectors[i][k]].name).c_str());
		}
	}
	//********************************找到指定feature的类
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
		cout << "未找到指定文件！" << endl;
}

void main()
{
	string Path = "D:\\files\\test\\"; 
	string FeatName;                //要查找的特征文件名（带后缀）
	int index=-1;							                         //类的编号
	int N;                                                           //确定的簇个数
	int tag=-1;
	cout << "输入：1-确定类别个数  2-不确定类别个数" << endl;
	cin >> tag;
	cin.get();
	if (tag == 1)
	{
		cout << "输入个数：";
		cin >> N;
		cin.get();
		cout << "输入要查找的文件名：";
		getline(cin,FeatName);
		certainN(Path, FeatName, N, index);
	}
	else if (tag == 2)
	{
		cout << "输入要查找的文件名：";
		getline(cin, FeatName);
		uncertainN(Path, FeatName, index);
	}
	else
	{
		cout << "输入错误！" << endl;
	}
	cout << FeatName << " 属于类型：" << index << endl;
		
}

	



