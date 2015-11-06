#include <stdio.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>

#include "../include/lsd_opencv.hpp"
#include "../include/k-means.h"
#include "../include/fcm.h"

using namespace std;
using namespace cv;

const float minA = 30;
const float maxA = 150;

const float MIN_LEN = 30;

const float epsilon = 0.05;
const float fuzziness = 3;

// int main()
// {
// 	int a = -1;
// 	int b = -2;
// 
// 	std::cout << fastAtan2(a,b) << endl;
// 
// 	return 0;
// 
// }

void linefit(std::vector<Point2f>& l_point, int n_point, Point2f& line) //��Ԫ������
{
	float av_x,av_y; //��������
	float L_xx,L_yy,L_xy;
	//������ʼ��
	av_x=0; //X��ƽ��ֵ
	av_y=0; //Y��ƽ��ֵ
	L_xx=0; //Lxx
	L_yy=0; //Lyy
	L_xy=0; //Lxy
	for(int i=0;i<n_point;i++) //����X��Y��ƽ��ֵ
	{
		Point2f p = l_point[i];
		av_x+=p.x/n_point;
		av_y+=p.y/n_point;
	}

	for(int i=0;i<n_point;i++) //����Lxx��Lyy��Lxy
	{
		L_xx+=(l_point[i].x-av_x)*(l_point[i].y-av_x);
		L_yy+=(l_point[i].x-av_y)*(l_point[i].y-av_y);
		L_xy+=(l_point[i].x-av_x)*(l_point[i].x-av_y);
	}

	//This line can be fitted by y=ax+b
	line.x = L_xy/L_xx;
	line.y = av_y-L_xy*av_x/L_xx;
}

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		std::cout << "lsd_filter [in_image]" << std::endl
			<< "\tin - input image" << std::endl;
		return false;
	}

	std::string in = argv[1];

	Mat image = imread(in, CV_LOAD_IMAGE_GRAYSCALE);

	//
	// LSD call
	//
	std::vector<Vec4i> lines, filtered_lines, retained_lines, long_lines;
	std::vector<double> width, prec, nfa;
	Ptr<LineSegmentDetector> ls = createLineSegmentDetectorPtr(LSD_REFINE_STD);

	double start = double(getTickCount());
	ls->detect(image, lines);
	ls->filterSize(lines, lines, MIN_LEN, LSD_NO_SIZE_LIMIT);    // Remove all lines smaller than MIN_LEN pixels
	//ls->retainAngle(lines, retained_lines, minA, maxA);       // 
	//double duration_ms = (double(getTickCount()) - start) * 1000 / getTickFrequency();

	//cout << "It took " << duration_ms << " ms." << endl;
	//
	// Show difference
	//
	Mat drawnLines(image);
	ls->drawSegments(drawnLines, lines);
	imshow("Drawing segments", drawnLines);

// 	Mat vertical(image);
// 	ls->drawSegments(vertical, retained_lines);
// 	imshow("Retained lines", vertical);


	//����׼��
	const int dim = 3;   //Dimension of feature
	const int cluster_num = 2; //Cluster number
	//const int size = retained_lines.size(); //Number of samples
	const int size = lines.size(); //Number of samples

	double * data = new double[size*dim];
	for (int i = 0; i < size; i++)
	{
		//const Vec4i& v = retained_lines[i];
		const Vec4i& v = lines[i];
		double x = (v[0] + v[2])/2.0;
		double y = (v[1] + v[3])/2.0;

		//б�ʼ��� ����ǶȽ���[89.5,90)�İ�89.5���� (90-90.5]�İ�90.5����
		//Ҳ���� k>114.5ʱ k = 114.5; k<-114.5ʱ k = -114.5;
		double k = (v[3]-v[1]) * 1.0 / (v[2]-v[0]+0.01); //��ֹ��ĸΪ��
		if (k > 114.5)
			k = 114.5;
		else if (k < -114.5)
			k = -114.5;

		data[i*dim+0] = x;
		data[i*dim+1] = y;
		data[i*dim+2] = k;
	}

	//KMeans����
	int* labels = new int[size];
	KMeans* kmeans = new KMeans(dim,cluster_num);
	kmeans->SetInitMode(KMeans::InitRandom);
	kmeans->Cluster(data,size,labels);

// 	//FCM����
// 	CFCM f(size,dim,cluster_num,epsilon,fuzziness);
// 	f.fcm(data);
// 	for (int i = 0; i < size; i++)
// 	{
// 		if (f.degree_of_memb[i][0] > 0.5)
// 			labels[i] = 1;
// 		else
// 			labels[i] = 0;
// 	}

	//��ʾ�����Ч��
	vector<Vec4i> lh;
	lh.reserve(20);
	vector<Vec4i> rh;
	rh.reserve(20);
	for (int i = 0; i < size; i++)
	{
		if (labels[i] == 0)
			lh.push_back(lines[i]);
		else
			rh.push_back(lines[i]);
	}

	Mat lhImg(image);
	ls->drawSegments(lhImg, lh);
	imshow("lh", lhImg);

	Mat rhImg(image);
	ls->drawSegments(rhImg, rh);
	imshow("rh", rhImg);
	

// 	//����󲿷�������������� ��ô����ɸѡ
// 	//ɸѡ�󲿷ֵ� ������� ȷ����0����1
// 	int total = 0;
// 	int flag = 0;
// 	for (int i = 0; i < size; i++)
// 		total += labels[i];
// 	if (total > size / 2)
// 		flag = 1;
// 
// 
// 	//ת�ɵ�
//  	Point2f center(0,0); //���ĵ�
//  	float k = 0.0; //ƽ��б��
// 
// 	//ͳ��Ŀ���ߵ����ĵ�����
// 	for (int i = 0; i < size; i++)
// 	{
// 		if (labels[i] == flag)
// 		{
// 			center.x += data[0];
// 			center.y += data[1];
// 			k        += data[2];
// 		}
// 	}
// 	delete[] data;
// 
// 	center.x /= (flag == 1) ? total : (size-total);
// 	center.y /= (flag == 1) ? total : (size-total);
// 	k        /= (flag == 1) ? total : (size-total);
// 	
// 
// 	int y1 = 20;
// 	int y2 = image.rows - 20;
//  	int x1 = (y1-center.y) / k + center.x;
//  	int x2 = (y2-center.y) / k + center.x;
// // 
// // 	double durat = (double(getTickCount())-start)*1000/getTickFrequency();
// // 	std::cout << "time is " << durat << "ms" <<endl;
// // 
// // 
// // 
// 	//��������
// 	Mat texture(image);
// 	std::vector<Vec4i> t;
// 	t.reserve(1);
// 	Vec4i v;
// 	v[0] = x1;
// 	v[1] = y1;
// 	v[2] = x2;
// 	v[3] = y2;
// 	t.push_back(v);
// 	circle(texture,center,5,Scalar(0),-1);
// 	ls->drawSegments(texture,t);
// 	imshow("texture lines",texture);



	waitKey();
	return 0;
}
