// IICorrection.cpp : Defines the entry point for the console application.

#pragma warning(disable : 4996)
#include "stdafx.h"
#include <vector>
#include <math.h>
#include <algorithm>
#include <fstream>
#include "..\..\ImageFileReader\MRImageFileLoader.h"
#include "DicomImageLoad.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include "InhomoCorrect.h"

using namespace std;

 const double PI = 3.141592653589793;
 //double epsilon  = 1;

 #define random(x) ((rand()%x)/10);

 void normalize01(vector<double>& Img,int a)
 {
	 double min = 0.0, max = 1.0;
	 for (auto iter = Img.begin(); iter != Img.end(); iter++)
	 {
		 min = min<(*iter)?min:(*iter);
		 max = max>(*iter)?max:(*iter);
	 }

	 for (auto iter = Img.begin();iter != Img.end(); iter++)
	 {
		 (*iter) = a*(*iter - min)/(max - min);
	 }
 }

int _tmain(int argc, _TCHAR* argv[])

	{
	//int m1 = 207, n1 = 169, m2 = 17, n2 = 17;

	//算法的参数值
	clock_t starttime, endtime;
	starttime = clock();

	//读取dcm格式的数据
	DicomImageFile::CMRImageFileLoader file_loader;
	DicomInfoFor3DPtr info;
	file_loader.GetDicomInfoFromDicomFile(CString("D:\\数据\\MRIMAGES-淮海医院2015.12.16\\201509\\15\\12\\09\\012\\5\\5.dcm"), info);
	//D:\\数据\\MRIMAGES-淮海医院2015.12.16\\201509\\15\\12\\09\\012\\5\\5.dcm  (256*256)
	//D:\\数据\\MRIMAGES-淮海医院2015.12.16\\201509\\15\\12\\08\\033\\14\\11.dcm (512*512)
	//D:\\数据\\007_fl3d_tra_QSM_monopolar_20160527\\00063.dcm  (176*256)
	if (info->pixel_data == nullptr)
	{
		return 0;
	}

	vector<double> img_data_double(info->row * info->columns, 0);	
	{
		vector<short> img_data(info->row * info->columns, 0);
		std::memcpy(img_data.data(), info->pixel_data, info->row * info->columns * sizeof(short));
		for (int i = 0; i < info->row * info->columns; i++)
		{
			img_data_double[i] = static_cast<double>(img_data[i]);
		}
	}
	normalize01(img_data_double,255);

	vector<double> output;
	CInhomoCorrect corrector;
	corrector.Correct(output, img_data_double, info->row, info->columns);
	endtime = clock();
	double correction_time = (endtime - starttime) * 1.0 / CLOCKS_PER_SEC;
	cout <<"correction time: " << correction_time << "s" << endl;
	
	//保存成.txt文件。
	ofstream ofs("d:\\matlabcode\\InhomoCorrect_img.txt");
	int line = 1;
	for (auto iter = output.begin(); iter!= output.end(); iter++)
	{
		ofs<< *iter<<",";
		if (line % info->columns == 0)
		{
			ofs <<";"<<endl;
		}
		line++;
	}
	ofs.close(); 

    

	return 0;
}
