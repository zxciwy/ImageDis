#pragma once
#include <vector>

class CInhomoCorrect
{
public:
	CInhomoCorrect(void);
	~CInhomoCorrect(void);

	void Correct(std::vector<double>& output, const std::vector<double>& input, int image_height, int image_width);
	void conv2(const std::vector<double>& image, const std::vector<double>& kernel, std::vector<double>& output, int image_height, int image_width, int kernel_height, int kernel_width);
	void CalcKernelConvOneImage(std::vector<double>& kernel_conv_one_image, const std::vector<double> &img, const std::vector<double> &kernel_conv1);
	void Clic(double (&cluster_center)[3], /*C ≤ª‘Ÿ ‰≥ˆ */ std::vector<double>& outline1, std::vector<double>& outline2, std::vector<double>& bias, const std::vector<double>& img, const std::vector<double>& kernel, const std::vector<double>& kernel_conv_one_image, int image_height, int image_width, int kernel_height, int kernel_width, double nu, double timestep, double mu, double epsilon, int Iter);
private:	
	void conv1(const std::vector<double>& image, const std::vector<double>& kernel, std::vector<double>& output, int image_height, int image_width, int kernel_width);
	void conv22(const std::vector<double>& image, const std::vector<double>& kernel, std::vector<double>& output, int image_height, int image_width, int kernel_width);
	void GradY(double* output, const double* input, int height, int width);
	void GradX(double* output, const double* input, int height, int width);
	void CurvatureCentral(std::vector<double>& output, const std::vector<double>& image, int height, int width);
	void Heaviside(std::vector<double>& h, const std::vector<double>& outline, double epsilon);
	void UpdateLSF( std::vector<double>& outline1,std::vector<double>& outline2, const std::vector<double>& img, double (&C)[3], int N_class, const std::vector<double>& KONE_Img, const std::vector<double>& KB1, const std::vector<double>& KB2, double mu, double nu, double timestep, double epsilon, int Iter, int m1, int n1);
	void NeumannBoundCond(std::vector<double>& image, int height,int width);
	void Dirac(std::vector<double>& Delta1, const std::vector<double>& u1, double epsilon);
	void UpdateC(double& C, const std::vector<double>& Img, const std::vector<double>& KB1, const std::vector<double>& KB2, const std::vector<double> (&M)[3]);//?;
	void Transpose(std::vector<double> &transpose, std::vector<double>& input, int image_height, int image_width);
	void UpdateBias(std::vector<double>& bias, const std::vector<double>& image, double (&C)[3], std::vector<double> (&M)[3], const std::vector<double>& Ksigma, int m1, int n1, int m2, int n2);
	std::vector<double> del2(const std::vector<double>& image, int height, int width);
	std::vector<double> GaussianFilter();
};

