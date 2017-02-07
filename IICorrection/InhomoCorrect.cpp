#include "stdafx.h"
#include "InhomoCorrect.h"
#include <math.h>
#include <algorithm>
#include <numeric>
using namespace std;

const double PI = 3.141592653589793;

CInhomoCorrect::CInhomoCorrect(void)
{
}


CInhomoCorrect::~CInhomoCorrect(void)
{
}
namespace 
{
	double multiply(double a, double b)
	{
		return a*b;
	}

	double divide(double a, double b)
	{
		return a/b;
	}

	double plus(double a, double b)
	{
		return a+b;
	} 
}
void CInhomoCorrect::conv2(const vector<double>& image,
		   const vector<double>& kernel,
		   vector<double>& output,
		   int image_height,
		   int image_width,
		   int kernel_height, 
		   int kernel_width)
{
	//do convolution
	// image 256*256 kernel 3*3
	const double* begin = kernel.data() + kernel_height / 2 * kernel_width;
	vector<double> kernel_one(begin, begin + kernel_width);
	double sum = std::accumulate(kernel_one.begin(), kernel_one.end(), 0.0);
	for_each(kernel_one.begin(), kernel_one.end(), [sum](double& v){ v /= sum;});
	conv22(image, kernel, output, image_height, image_width, kernel_width);
	return;
	output.resize(image_height * image_width);
	double* full_output = output.data();
	int kernel_radius_height = kernel_height / 2;
	int kernel_radius_width = kernel_width / 2;
	//if((i - k) >= 0 && (i - k)<n2&&(j - l) >= 0 && (j - l) < m2)
	for(int i = kernel_radius_height; i < image_height - kernel_radius_height; i++)
	{
		for(int j = kernel_radius_width; j < image_width - kernel_radius_width; j++)
		{
			double temp = 0.0;
			for (int k = 0; k < kernel_height; k++)
			{
				int kernel_idx = k * kernel_width;
				const double* kernel_cursor = kernel.data() + kernel_idx;
				int image_idx = (i + kernel_radius_height - k) * image_width + j + kernel_radius_width; 
				const double* image_cursor = image.data() + image_idx;
				for(int l = 0; l < kernel_width; l++)
				{
					//k * n2 + l
					//(i + kernel_radius_height - k) * n1  + j - l + kernel_radius_width
					temp += (*kernel_cursor++) * (*image_cursor--);
				}
			}
			full_output[i * image_width + j] = temp;			
		}
	}

	for (int i = 0; i < kernel_radius_height; ++i)
	{
		for (int j = 0; j < image_width; ++j)
		{
			//int center = 1;
			//left = max(-kernel_width, -j);
			//j + left >= 0;

			//// top, bottom, [-kenerl_width, kernel_width]
			////[-kenerl_height, kernel_height]
			////int left, right,
			//double* image_center;
			//double* kernel_center;
			//temp = 0;
			//for (int k = top; k <= bottom; ++k)
			//{
			//	for (int l = left; l <= right; ++l)
			//	{
			//		temp += image_center[top * width + l] * kernel_center[-top * kernel_width - l]
			//	}
			//}
			//full_output[i * image_width + j] = temp;	
			full_output[i * image_width + j] = 1;
		}
	}

	for (int i = kernel_radius_height; i < image_height - kernel_radius_height; ++i)
	{
		for (int j = 0; j < kernel_radius_width; ++j)
		{
			full_output[i * image_width + j] = 1;			
		}
		for (int j = image_width - kernel_radius_width; j < image_width; ++j)
		{
			full_output[i * image_width + j] = 1;			
		}
	}

	for (int i = image_height - kernel_radius_height; i < image_height; ++i)
	{
		for (int j = 0; j < image_width; ++j)
		{
			full_output[i * image_width + j] = 1;			
		}
	}
}

void CInhomoCorrect::conv1(const vector<double>& image,
		   const vector<double>& kernel,
		   vector<double>& output,
		   int image_height,
		   int image_width,
		   int kernel_width)
{
	output.resize(image_height * image_width);
	double* full_output = output.data();
	int kernel_radius_width = kernel_width / 2;
	//if((i - k) >= 0 && (i - k)<n2&&(j - l) >= 0 && (j - l) < m2)
	for(int i = 0; i < image_height; i++)
	{
		for (int j = 0; j < kernel_radius_width; ++j)
		{
			full_output[i * image_width + j] = 1;			
		}
		for(int j = kernel_radius_width; j < image_width - kernel_radius_width; j++)
		{
			double temp = 0.0;
			const double* kernel_cursor = kernel.data();
			int image_idx = i * image_width + j + kernel_radius_width; 
			const double* image_cursor = image.data() + image_idx;
			for(int l = 0; l < kernel_width; l++)
			{
				temp += (*kernel_cursor++) * (*image_cursor--);
			}
			full_output[i * image_width + j] = temp;			
		}
		for (int j = image_width - kernel_radius_width; j < image_width; ++j)
		{
			full_output[i * image_width + j] = 1;			
		}
	}
}

void CInhomoCorrect::conv22(const vector<double>& image,
			const vector<double>& kernel,
			vector<double>& output,
			int image_height,
			int image_width,
			int kernel_width)
{
	conv1(image, kernel, output, image_height, image_width, kernel_width);
	vector<double> transpose(image_width * image_height);
	Transpose(transpose, output, image_height, image_width);
	conv1(transpose, kernel, output, image_height, image_width, kernel_width);
	Transpose(transpose, output, image_height, image_width);
	output.swap(transpose);
}

void CInhomoCorrect::CalcKernelConvOneImage(vector<double>& kernel_conv_one_image, const vector<double> &img, const vector<double> &kernel_conv1)
{
	kernel_conv_one_image.resize(img.size());
	transform(img.begin(), img.end(), img.begin(), kernel_conv_one_image.begin(), multiply);
	transform(kernel_conv_one_image.begin(), kernel_conv_one_image.end(), kernel_conv1.begin(), kernel_conv_one_image.begin(), multiply);
}

void CInhomoCorrect::Clic(double (&cluster_center)[3], //C 不再输出
		  vector<double>& outline1, vector<double>& outline2, vector<double>& bias,
		  const vector<double>& img, const vector<double>& kernel, const vector<double>& kernel_conv_one_image,
		  int image_height, int image_width, int kernel_height, int kernel_width,
		  double nu, double timestep, double mu, double epsilon, int Iter)
{
	int N_class = 3;
	//int m1,m2,n1,n2;
	vector<double> bias_conv_kernel;
	conv2(bias, kernel, bias_conv_kernel, image_height, image_width, kernel_height, kernel_width);
	vector<double> bias_square(bias.size());
	transform(bias.begin(), bias.end(), bias.begin(), bias_square.begin(), multiply);

	vector<double> bias_square_conv_kernel;
	conv2(bias_square, kernel, bias_square_conv_kernel, image_height, image_width, kernel_height, kernel_width);

	vector<double> heaviside1, heaviside2;
	Heaviside(heaviside1, outline1, epsilon);
	Heaviside(heaviside2, outline2, epsilon);

	ASSERT(sizeof(cluster_center) / sizeof(double) == 3);
	vector<double> M[3];
	for (int i = 0; i < 3; i++)
	{
		M[i].resize(img.size());
	}
	auto M0_cursor = M[0].data();
	auto M1_cursor = M[1].data();
	auto M2_cursor = M[2].data();
	for(auto iter1 = heaviside1.begin(), iter2 = heaviside2.begin(); iter1 != heaviside1.end(); ++iter1, ++iter2)
	{
		double temp = *iter1 * *iter2;
		*M0_cursor++ = temp;
		*M1_cursor++ = *iter1 - temp;
		*M2_cursor++ = 1 - *iter1; 
	}

	UpdateC(*cluster_center, img, bias_conv_kernel, bias_square_conv_kernel, M);
	UpdateLSF(outline1, outline2, img, cluster_center, N_class, kernel_conv_one_image, bias_conv_kernel, bias_square_conv_kernel, 
		mu, nu, timestep, epsilon, Iter, image_height, image_width);
	UpdateBias(bias, img, cluster_center, M, kernel, image_height, image_width,  kernel_height,  kernel_width);
}

void CInhomoCorrect::GradX(double* output, const double* input, int height, int width)
{
	for(int i = 1; i < height-1;i++)
	{
		for(int j = 0; j < width; j++)
		{
			double* center = const_cast<double*>(&input[i * width + j]);
			output[i * width + j] = (center[width] - center[-width]) / 2;
		}
	}	 
	for (int j = 0; j < width ; j++)
	{
		//tx[0][j] = array[1][j]-array[0][j];
		double* center = const_cast<double*>(&input[j]);
		output[j] = center[width] + center[0];
		center = const_cast<double*>(&input[(height - 1) * width + j]);
		output[(height - 1) * width + j] = center[0] - center[-width];
	}
}

void CInhomoCorrect::GradY(double* output, const double* input, int height, int width)
{
	for(int i = 0; i < height; i++)
	{
		for(int j = 1; j < width-1; j++)
		{
			double* center = const_cast<double*>(&input[i * width + j]);
			output[i * width + j] = (center[1] - center[-1]) / 2;
		}
	}
	for (int i = 0; i < height; i++)
	{		
		double* center = const_cast<double*>(&input[i * width]);
		output[i * width] = center[1] - center[0];
		output[i * width + width - 1] = center[width -1] - center[width - 2];
	}
}

void CInhomoCorrect::CurvatureCentral(vector<double>& output, const vector<double>& image, int height, int width)
{
	output.resize(image.size());//
	double* tx = new double[height * width];
	double* txx = new double[height * width];
	double* ty = new double[height * width];
	double* tyy = new double[height * width];

	GradX(tx, image.data(), height, width);
	GradY(ty, image.data(), height, width);

	for (int i = 0; i < height * width; i++)
	{
		double norm = sqrt(tx[i] * tx[i] + (ty[i] * ty[i])+1e-20);
		tx[i] = tx[i] / norm;
		ty[i] = ty[i] / norm;
	}

	GradX(txx, tx, height, width);
	GradY(tyy, ty, height, width);

	for(int i = 0; i < height * width; i++)
	{
		output[i] = txx[i] + tyy[i];//??
	}

	delete []tx;
	delete []txx;
	delete []ty;
	delete []tyy;
}

void CInhomoCorrect::Heaviside(vector<double>& h, const vector<double>& outline, double epsilon)
{
	double pi_inv = 1 / PI;
	double epsilon_inv = 1 / epsilon;
	h.resize(outline.size());
	auto iter_h = h.begin();
	for (auto iter = outline.begin(); iter != outline.end(); iter++, iter_h++)
	{
		*iter_h = (0.5 + pi_inv * atan(*iter * epsilon_inv));
	}
}

void CInhomoCorrect::UpdateLSF( vector<double>& outline1,vector<double>& outline2,
			   const vector<double>& img, double (&C)[3], 
			   int N_class, const vector<double>& KONE_Img, 
			   const vector<double>& KB1, const vector<double>& KB2,
			   double mu, double nu, double timestep, double epsilon, int Iter,
			   int m1, int n1)
{
	NeumannBoundCond(outline1, m1, n1);

	vector<double> Curv1;
	CurvatureCentral(Curv1, outline1, m1, n1);
	vector<double> H1; 
	Heaviside(H1, outline1, epsilon);
	//H1 = Heaviside(u1,epsilon );
	vector<double> Delta1;
	Dirac(Delta1, outline1, epsilon);
	//vector<double> Delta1 = Dirac(u1,epsilon);

	NeumannBoundCond(outline2, m1, n1);
	//u2 = NeumannBoundCond(u2);
	vector<double> Curv2;
	CurvatureCentral(Curv2, outline2, m1, n1);
	// = curvature_central(u2);
	vector<double> H2;
	Heaviside(H2, outline2, epsilon);
	//H2 = Heaviside(u2,epsilon );
	vector<double> Delta2;
	Dirac(Delta2, outline2, epsilon);
	//= Dirac(u2,epsilon);

	vector<double> DEL1 = del2(outline1,m1,n1);//?
	auto DEL2 = del2(outline2,m1,n1);
	for(int i = 0; i != img.size(); i++)
	{
		double e1 = KONE_Img[i] - 2 * img[i] * C[0] * KB1[i] + C[0] * C[0] * KB2[i];  
		double e2 = KONE_Img[i] - 2 * img[i] * C[1] * KB1[i] + C[1] * C[1] * KB2[i]; 
		double e3 = KONE_Img[i] - 2 * img[i] * C[2] * KB1[i] + C[2] * C[2] * KB2[i];  

		double A1 = - Delta1[i]*(e1 *H2[i] + e2 * (1-H2[i]) - e3);
		double P1 = mu * (4 * DEL1[i]-Curv1[i]);
		double L1 = nu*Delta1[i]*Curv1[i];
		outline1[i] = (outline1)[i] + timestep*(L1 + P1 + A1); 

		double A2 = - Delta2[i]*H1[i]*(e1 - e2);
		double P2 = mu*(4*DEL2[i]-Curv2[i]);
		double L2 = nu*Delta2[i]*Curv2[i];
		outline2[i] = (outline2)[i] + timestep*(L2 + P2 + A2);   
	}
}

void CInhomoCorrect::Dirac(vector<double>& Delta1, const vector<double>& u1, double epsilon)
 {
	 Delta1.resize(u1.size());//
	 double eps_devide_pi = epsilon / PI;
	 double epsilon_square = epsilon * epsilon;
	 for (unsigned int i = 0; i < u1.size(); i++)
	 {
		 Delta1[i] = eps_devide_pi / (epsilon_square + u1[i] * u1[i]);
	 }	 
 }

 vector<double> CInhomoCorrect::del2(const vector<double>& image, int height, int width)
 {
	 vector<double> output(width*height);

	 double* temp = output.data();

	 for(int i = 1; i<height-1;i++)
	 {
		 for (int j = 1;j < width - 1;j++)
		 {
			 double* center = const_cast<double*>(&image[i * width + j]);
			 temp[i * width + j] = (center[width] + center[-width] + center[-1] + center[1]) / 4 - center[0];
		 }
	 }
	 for(int i = 1; i < height - 1; i++)
	 {
		 double* center = const_cast<double*>(&image[i * width]);
		 temp[i * width + 0] = (center[-width] + center[width] + center[2]*4 - center[1]*5 - center[3])/4;
		 
		 center = const_cast<double*>(&image[i * width + width - 1]);
		 temp[i * width + width - 1] = (center[-width] + center[width] + center[-2]*4
			 - center[-1]*5-center[-3])/4;
		/* temp[i][0]=(array[i+1][0]+array[i-1][0]+array[i][0]+array[i][2])/4-(array[i][0]+array[i][1])/2;
		 temp[i][n-1]=(array[i+1][n-1]+array[i-1][n-1]+array[i][n-1]+array[i][n-3])/4-(array[i][n-1]+array[i][n-2])/2;*/
	 }
	 for(int j = 1; j < width - 1; j++)
	 {
		 double* center = const_cast<double*>(&image[j]);
		 temp[j] = (center[-1] + center[1] + center[2*width] - center[width] - center[3*width]) / 4;
		 center = const_cast<double*>(&image[(height-1) * width+ j]);
		 temp[(height-1) * width+ j] = (center[-1] + center[1] + center[-2 * width] * 4 - center[-width] * 5 - center[-3 * width]) / 4;
	 } 
	 //temp[0][0] = (array[0][0]*4+array[2][0]*4+array[0][2]*4-array[0][1]*5-array[1][0]*5-array[0][3]-array[3][0])/4;
	 double* center = const_cast<double*>(&image[0]);
	 temp[0] = (center[0] * 4 + center[2 * width] * 4 + center[2] * 4 - center[1] * 5 - center[width] * 5 - center[3] - center[3 * width]) / 4;
	 //temp[0][width-1] = (array[0][width-1]*4+array[0][width-3]*4+array[2][width-1]*4-array[0][width-2]*5-array[1][width-1]*5-array[0][width-4]-array[3][width-1])/4;
	 center = const_cast<double*>(&image[width-1]);
	 temp[width-1] = (center[0] * 4 + center[-2] * 4 + center[2 * width] * 4 - center[-1] * 5 - center[width] * 5 - center[-3] - center[3 * width]) / 4;
	 center = const_cast<double*>(&image[(height-1) * width]); 
	 temp[(height-1) * width] = (center[0] * 4 + center[3] * 4 + center[-2 * width] * 4 - center[1] * 5 - center[-width] * 5 - center[3] - center[-3 * width]) / 4;
	 
	 //temp[height-1][width-1] = (array[height-1][width-1]*4+array[height-1][width-3]*4+array[height-3][width-1]*4-array[height-1][width-2]*5-array[height-2][width-1]*5-array[height-1][width-4]-array[width-4][height-1])/4;
	 center = const_cast<double*>(&image[(height-1) * width + width - 1]);
	 temp[(height-1) * width + width - 1] = (center[0] * 4 + center[-2] * 4 + center[-2 * width] *4 - center[-1] * 5 - center[-width] * 5 - center[-3] - center[-3 * width]) / 4;
	 
	 return output;
 }

 void CInhomoCorrect::UpdateC(double& C,
	 const vector<double>& Img,
	 const vector<double>& KB1,
	 const vector<double>& KB2, 
	 const vector<double> (&M)[3])//?
 {
	 vector<double> N2(Img.size());

	 transform(KB1.begin(), KB1.end(), Img.begin(), N2.begin(), multiply);
	 for (int i = 0; i < 3; i++)
	 {
		 double sN2 = 0.0;
		 double sD2 = 0.0;
		 auto iter_kb2 = KB2.begin();
		 auto iter_m = M[i].begin();
		 for (auto iter_n2 = N2.begin(); iter_n2 != N2.end(); ++iter_n2, ++iter_kb2, ++iter_m)
		 {
			 sN2 += (*iter_n2) * (*iter_m);
			 sD2 += (*iter_kb2) * (*iter_m);
		 }
		 (&C)[i] = sN2 / (sD2 + (sD2 == 0));
	 }
 }

 void CInhomoCorrect::NeumannBoundCond(vector<double>& image, int height,int width)
 {
	 // col2 copy to col0; col[width-3] copy to col[width-1]
	 // without corner
	 for (int i = 1; i < height - 1; i++)
	 {
		 image[i * width] = image[i * width + 2];
		 image[i * width + width - 1] = image[i * width + width-3];
	 }

	 // row2 copy to row0; row[height-3] copy to row[height-1] 
	 // without corner
	 for (int j = 1; j < width-1; j++)
	 {
		 image[j] = image[2 * width + j];
		 image[(height-1) * width + j] = image[(height-3) * width + j];
	 }

	 // corner
	 image[0] = image[2 * width + 2];
	 image[width-1] = image[2 * width + width - 3];
	 image[(height-1) * width] = image[(height-3) * width + 2];
	 image[(height-1) * width + width-1] = image[(height-3) * width + width-3];
 }

 void CInhomoCorrect::Transpose(vector<double> &transpose, vector<double>& input, int image_height, int image_width)
 {
	 for (int i = 0; i < image_height; ++i)
	 {
		 for (int j = 0; j < image_width; ++j)
		 {
			 transpose[j * image_height + i] = input[i * image_width + j];
		 }
	 }
 }

 void CInhomoCorrect::UpdateBias(vector<double>& bias,
	 const vector<double>& image,
	 double (&C)[3],
	 vector<double> (&M)[3], 
	 const vector<double>& Ksigma,
	 int m1, 
	 int n1, 
	 int m2, 
	 int n2)
 {
	 vector<double> PC1(image.size(), 0);
	 vector<double> PC2(image.size(), 0);
	 for (unsigned int i = 0; i < 3; i ++)
	 {
		 double* m_cursor = M[i].data();
		 const double* img_cursor = image.data();
		 for(auto iter = PC1.begin(), iter_pc2 = PC2.begin(); iter != PC1.end(); ++iter, ++iter_pc2)
		 {
			 (*iter) += C[i]*(*m_cursor) * (*img_cursor);
			 (*iter_pc2) += C[i] * C[i] * (*m_cursor);
			 ++m_cursor;
			 ++img_cursor;
		 }
	 }

	 vector<double> KNm1(image.size()), KDn1(image.size());
	 conv2(PC1, Ksigma, KNm1, m1, n1, m2, n2);
	 conv2(PC2, Ksigma, KDn1, m1, n1, m2, n2);
	 transform(KNm1.begin(), KNm1.end(), KDn1.begin(), bias.begin(), divide);
 }

 void CInhomoCorrect::Correct(vector<double>& output, const vector<double>& input, int image_height, int image_width)
 {
	 //算法的参数值
	 double epsilon = 1;
	 int sigma = 4;
	 int Iter_outer = 10;
	 int Iter_inner = 10;//没用到，可考虑去掉
	 double timestep = 0.1;
	 double mu = 0.1/timestep;
	 int A = 255;
	 double nu = pow(A,2)*0.001;
	 double C[3];

	 //设置偏场与水平集函数等的初值
	 vector<double> b(image_height * image_width, 1);
	 vector<double> u1, u2;
	 u1.resize(input.size());
	 u2.resize(input.size());
	 if ((image_width - image_width / 2)>0)
	 {
		 u1[image_height*(image_width - image_width / 2)] = 1;
	 }
	 for (int i =0; i < image_height * image_width; i++) 
	 {
		 u2[i] = rand()%10/10;
	 }

	 int m2 = 17, n2 = 17;
	 vector<double> Ksigma = GaussianFilter();
	 vector<double> ONE(image_height * image_width, 1);

	 vector<double> KONE;
	 conv2(ONE, Ksigma, KONE, image_height, image_width, m2, n2); 

	 vector<double> kernel_conv_one_image;
	 CalcKernelConvOneImage(kernel_conv_one_image, input, KONE);
	 for(int i = 0; i < Iter_outer; i++)
	 {
		 Clic(C, u1, u2, b, input, Ksigma, kernel_conv_one_image, image_height, image_width, m2, n2,
			 nu, timestep, mu, epsilon, Iter_inner);
	 }

	 //灰度校正
	 output.resize(input.size());
	 for (int i = 0; i < image_height * image_width; i++)
	 {
		 output[i] = input[i]/b[i];
	 }
 }

 vector<double> CInhomoCorrect::GaussianFilter()
 {
	 double temp[17*17] =
	 {0.000194886873895782,0.000311428337740562,0.000467509306907072,0.000659293796170428,
	 0.000873422372691078,0.00108699170538682,0.00127082205360098,0.00139572497744543,
	 0.00144003004416116,0.00139572497744543,0.00127082205360098,0.00108699170538682,
	 0.000873422372691078,0.000659293796170428,0.000467509306907072,0.000311428337740562,
	 0.000194886873895782,0.000311428337740562,0.000497661066695108,
	 0.000747077745247076,0.00105354848646102,0.00139572497744543,0.00173700779934222,
	 0.00203076785935128,0.00223036216334025,0.00230116145836174,0.00223036216334025,
	 0.00203076785935128,0.00173700779934222,0.00139572497744543,0.00105354848646102,
	 0.000747077745247076,0.000497661066695108,0.000311428337740562,0.000467509306907072,
	 0.000747077745247076,0.00112149652603906,0.00158156359909902,0.00209523135104675,
	 0.00260755754679960,0.00304854362741192,0.00334817016556663,0.00345445249550854,
	 0.00334817016556663,0.00304854362741192,0.00260755754679960,0.00209523135104675,
	 0.00158156359909902,0.00112149652603906,0.000747077745247076,0.000467509306907072,
	 0.000659293796170428,0.00105354848646102,0.00158156359909902,0.00223036216334025,
	 0.00295474980043870,0.00367724554006382,0.00429913559198315,0.00472167673684358,
	 0.00487155884558024,0.00472167673684358,0.00429913559198315,0.00367724554006382,
	 0.00295474980043870,0.00223036216334025,0.00158156359909902,0.00105354848646102,
	 0.000659293796170428,0.000873422372691078,0.00139572497744543,0.00209523135104675,
	 0.00295474980043870,0.00391440750147834,0.00487155884558024,0.00569542930796808,
	 0.00625520537661502,0.00645376690987549,0.00625520537661502,0.00569542930796808,
	 0.00487155884558024,0.00391440750147834,0.00295474980043870,0.00209523135104675,
	 0.00139572497744543,0.000873422372691078,0.00108699170538682,0.00173700779934222,
	 0.00260755754679960,0.00367724554006382,0.00487155884558024,0.00606275293949040,
	 0.00708807629612661,0.00778472886940398,0.00803184269017550,0.00778472886940398,
	 0.00708807629612661,0.00606275293949040,0.00487155884558024,0.00367724554006382,
	 0.00260755754679960,0.00173700779934222,0.00108699170538682,0.00127082205360098,
	 0.00203076785935128,0.00304854362741192,0.00429913559198315,0.00569542930796808,
	 0.00708807629612661,0.00828680074565843,0.00910127011964847,0.00939017544581588,
	 0.00910127011964847,0.00828680074565843,0.00708807629612661,0.00569542930796808,
	 0.00429913559198315,0.00304854362741192,0.00203076785935128,0.00127082205360098,
	 0.00139572497744543,0.00223036216334025,0.00334817016556663,0.00472167673684358,
	 0.00625520537661502,0.00778472886940398,0.00910127011964847,0.00999578973033754,
	 0.0103130901570230,0.00999578973033754,0.00910127011964847,0.00778472886940398,
	 0.00625520537661502,0.00472167673684358,0.00334817016556663,0.00223036216334025,
	 0.00139572497744543,0.00144003004416116,0.00230116145836174,0.00345445249550854,
	 0.00487155884558024,0.00645376690987549,0.00803184269017550,0.00939017544581588,
	 0.0103130901570230,0.0106404627804524,0.0103130901570230,0.00939017544581588,
	 0.00803184269017550,0.00645376690987549,0.00487155884558024,0.00345445249550854,
	 0.00230116145836174,0.00144003004416116,0.00139572497744543,0.00223036216334025,
	 0.00334817016556663,0.00472167673684358,0.00625520537661502,0.00778472886940398,
	 0.00910127011964847,0.00999578973033754,0.0103130901570230,0.00999578973033754,
	 0.00910127011964847,0.00778472886940398,0.00625520537661502,0.00472167673684358,
	 0.00334817016556663,0.00223036216334025,0.00139572497744543,0.00127082205360098,
	 0.00203076785935128,0.00304854362741192,0.00429913559198315,0.00569542930796808,
	 0.00708807629612661,0.00828680074565843,0.00910127011964847,0.00939017544581588,
	 0.00910127011964847,0.00828680074565843,0.00708807629612661,0.00569542930796808,
	 0.00429913559198315,0.00304854362741192,0.00203076785935128,0.00127082205360098,
	 0.00108699170538682,0.00173700779934222,0.00260755754679960,0.00367724554006382,
	 0.00487155884558024,0.00606275293949040,0.00708807629612661,0.00778472886940398,
	 0.00803184269017550,0.00778472886940398,0.00708807629612661,0.00606275293949040,
	 0.00487155884558024,0.00367724554006382,0.00260755754679960,0.00173700779934222,
	 0.00108699170538682,0.000873422372691078,0.00139572497744543,0.00209523135104675,
	 0.00295474980043870,0.00391440750147834,0.00487155884558024,0.00569542930796808,
	 0.00625520537661502,0.00645376690987549,0.00625520537661502,0.00569542930796808,
	 0.00487155884558024,0.00391440750147834,0.00295474980043870,0.00209523135104675,
	 0.00139572497744543,0.000873422372691078,0.000659293796170428,0.00105354848646102,
	 0.00158156359909902,0.00223036216334025,0.00295474980043870,0.00367724554006382,
	 0.00429913559198315,0.00472167673684358,0.00487155884558024,0.00472167673684358,
	 0.00429913559198315,0.00367724554006382,0.00295474980043870,0.00223036216334025,
	 0.00158156359909902,0.00105354848646102,0.000659293796170428,0.000467509306907072,
	 0.000747077745247076,0.00112149652603906,0.00158156359909902,0.00209523135104675,
	 0.00260755754679960,0.00304854362741192,0.00334817016556663,0.00345445249550854,
	 0.00334817016556663,0.00304854362741192,0.00260755754679960,0.00209523135104675,
	 0.00158156359909902,0.00112149652603906,0.000747077745247076,0.000467509306907072,
	 0.000311428337740562,0.000497661066695108,0.000747077745247076,0.00105354848646102,
	 0.00139572497744543,0.00173700779934222,0.00203076785935128,0.00223036216334025,
	 0.00230116145836174,0.00223036216334025,0.00203076785935128,0.00173700779934222,
	 0.00139572497744543,0.00105354848646102,0.000747077745247076,0.000497661066695108,
	 0.000311428337740562,0.000194886873895782,0.000311428337740562,0.000467509306907072,
	 0.000659293796170428,0.000873422372691078,0.00108699170538682,0.00127082205360098,
	 0.00139572497744543,0.00144003004416116,0.00139572497744543,0.00127082205360098,
	 0.00108699170538682,0.000873422372691078,0.000659293796170428,0.000467509306907072,
	 0.000311428337740562,0.000194886873895782
	 };

	 return vector<double>(temp, temp + 17 * 17);
 }