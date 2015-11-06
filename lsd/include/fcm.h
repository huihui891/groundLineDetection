#define MAX_DATA_POINTS 2000
#define MAX_CLUSTER 10
#define MAX_DATA_DIMENSION 5
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

using namespace std;

class CFCM
{
public:
	CFCM(int _num_points, int _num_dim, int _num_clusters
		,double _epsilon, double _fuzziness):
	    num_data_points(_num_points)
		,num_dimensions(_num_dim)
		,num_clusters(_num_clusters)
		,epsilon(_epsilon),fuzziness(_fuzziness)
		{
			//内存分配
			//1
			degree_of_memb = new double*[num_data_points];
			for (int i = 0; i < num_data_points; i++)
			{
				degree_of_memb[i] = new double[num_clusters];
				memset(degree_of_memb[i],0.0,num_clusters*sizeof(double));
			}

			//2
			low_high = new double*[num_dimensions];
			for (int i = 0; i < num_dimensions; i++)
			{
				low_high[i] = new double[2];
				memset(low_high[i],0.0,2*sizeof(double));
			}

			//3
			data_point = new double*[num_data_points];
			for (int i = 0; i < num_data_points; i++)
			{
				data_point[i] = new double[num_dimensions];
				memset(data_point[i],0.0,num_dimensions*sizeof(double));
			}

			//4
			cluster_centre = new double*[num_clusters];
			for (int i = 0; i < num_clusters; i++)
			{
				cluster_centre[i] = new double[num_dimensions];
				memset(cluster_centre[i],0.0,sizeof(double)*num_dimensions);
			}
		}

	~CFCM()
	{
		//内存释放
		//1 
		for (int i = 0; i < num_data_points; i++)
			delete[] degree_of_memb[i];
		delete[] degree_of_memb;

		//2
		for (int i = 0; i < num_dimensions; i++)
			delete[] low_high[i];
		delete[] low_high;

		//3
		for (int i = 0; i < num_data_points; i++)
			delete[] data_point[i];
		delete[] data_point;

		//4
		for (int i = 0; i < num_clusters; i++)
			delete[] cluster_centre[i];
		delete[] cluster_centre;
	}
	
	int fcm(vector<vector<double>>& database) {
		double max_diff;
		init(database);
		do {
			calculate_centre_vectors();
			max_diff = update_degree_of_membership();
		} while (max_diff > epsilon);
		return 0;
	}

	int fcm(double* database) {
		double max_diff;
		init(database);
		do {
			calculate_centre_vectors();
			max_diff = update_degree_of_membership();
		} while (max_diff > epsilon);
		return 0;
	}

private:
	int init(vector<vector<double>>& dataBase) {
		for (int i = 0; i < num_data_points; i++) {
			for (int j = 0; j < num_dimensions; j++) {
				//填充数据
				data_point[i][j] = dataBase[i][j];
				if (data_point[i][j] < low_high[j][0])
					low_high[j][0] = data_point[i][j];
				if (data_point[i][j] > low_high[j][1])
					low_high[j][1] = data_point[i][j];
			}
		}

		//数据进行归一化
		for (int i = 0; i < num_data_points; i++) {
			for (int j = 0; j < num_dimensions; j++) {
				data_point[i][j] = (data_point[i][j] - low_high[j][0]) / 
					(low_high[j][1] - low_high[j][0]);
			}
		}

		double s;
		int rval,r;
		for (int i = 0; i < num_data_points; i++) {
			s = 0.0;
			r = 100;
			for (int j = 1; j < num_clusters; j++) {
				rval = rand() % (r + 1);
				r -= rval;
				degree_of_memb[i][j] = rval / 100.0;
				s += degree_of_memb[i][j];
			}
			degree_of_memb[i][0] = 1.0 - s;
		}
		return 0;
		// 	failure:
		// 			exit(1);
	}

	int init(double* dataBase) {
		// 			if (num_clusters > MAX_CLUSTER) {
		// 				printf("Number of clusters should be < %d\n", MAX_CLUSTER);
		// 				goto failure;
		// 			}
		// 			if (num_data_points > MAX_DATA_POINTS) {
		// 				printf("Number of data points should be < %d\n", MAX_DATA_POINTS);
		// 				goto failure;
		// 			}
		// 			if (num_dimensions > MAX_DATA_DIMENSION) {
		// 				printf("Number of dimensions should be >= 1.0 and < %d\n",
		// 					MAX_DATA_DIMENSION);
		// 				goto failure;
		// 			}
		// 
		// 			if (fuzziness <= 1.0) {
		// 				printf("Fuzzyness coefficient should be > 1.0\n");
		// 				goto failure;
		// 			}
		// 
		// 			if (epsilon <= 0.0 || epsilon > 1.0) {
		// 				printf("Termination criterion should be > 0.0 and <= 1.0\n");
		// 				goto failure;
		// 			}

		for (int i = 0; i < num_data_points; i++) {
			for (int j = 0; j < num_dimensions; j++) {
				//填充数据
				data_point[i][j] = dataBase[i*num_dimensions+j];
				if (data_point[i][j] < low_high[j][0])
					low_high[j][0] = data_point[i][j];
				if (data_point[i][j] > low_high[j][1])
					low_high[j][1] = data_point[i][j];
			}
		}

		double s;
		int rval,r;
		for (int i = 0; i < num_data_points; i++) {
			s = 0.0;
			r = 100;
			for (int j = 1; j < num_clusters; j++) {
				rval = rand() % (r + 1);
				r -= rval;
				degree_of_memb[i][j] = rval / 100.0;
				s += degree_of_memb[i][j];
			}
			degree_of_memb[i][0] = 1.0 - s;
		}
		return 0;
		// 	failure:
		// 			exit(1);
	}

	int calculate_centre_vectors() {
			double numerator, denominator;
			//分配内存 并初始化
			double **t = new double *[num_data_points];
			for (int i = 0; i < num_data_points; i++)
				t[i] = new double[num_clusters];
			for (int i = 0; i < num_data_points; i++)
			{
				for (int j = 0; j < num_clusters; j++)
				t[i][j] = 0.0;				
			}
			
			for (int i = 0; i < num_data_points; i++) {
				for (int j = 0; j < num_clusters; j++) {
					//uij的m次方
					t[i][j] = pow(degree_of_memb[i][j], fuzziness);
				}
			}
			for (int j = 0; j < num_clusters; j++) {
				for (int k = 0; k < num_dimensions; k++) {
					numerator = 0.0;
					denominator = 0.0;
					for (int i = 0; i < num_data_points; i++) {
						numerator += t[i][j] * data_point[i][k];
						denominator += t[i][j];
					}
					//每一类的每一维的中心
					cluster_centre[j][k] = numerator / denominator;
				}
			}


			//释放内存
			for (int i = 0; i < num_data_points; i++)
				delete[] t[i];
			delete[] t;
			
			return 0;
	}

	double get_norm(int i, int j) {
			int k;
			double sum = 0.0;
			for (k = 0; k < num_dimensions; k++) {
				sum += pow(data_point[i][k] - cluster_centre[j][k], 2);
			}
			return sqrt(sum);
	}

	double get_new_value(int i, int j) {
			int k;
			double t, p, sum;
			sum = 0.0;
			p = 2 / (fuzziness - 1);
			for (k = 0; k < num_clusters; k++) {
				t = get_norm(i, j) / get_norm(i, k);
				t = pow(t, p);
				sum += t;
			}
			return 1.0 / sum;
	}

	double update_degree_of_membership() {
			int i, j;
			double new_uij;
			double max_diff = 0.0, diff;
			for (j = 0; j < num_clusters; j++) {
				for (i = 0; i < num_data_points; i++) {
					new_uij = get_new_value(i, j);
					diff = new_uij - degree_of_memb[i][j];
					if (diff > max_diff)
						max_diff = diff;
					degree_of_memb[i][j] = new_uij;
				}
			}
			return max_diff;
	}

public:
	//每个点隶属某个类的程度
	double **degree_of_memb;

private:
	int num_data_points;
	int num_clusters;
	int num_dimensions;
	//每一维的最小值和最大值
	double **low_high;

	//迭代精度 停止条件
	double epsilon;
	//模糊参数 程度
	double fuzziness;
	//数据点
	double **data_point;
	//数据点中心 每一类中每一维的中心
	double **cluster_centre;
};

