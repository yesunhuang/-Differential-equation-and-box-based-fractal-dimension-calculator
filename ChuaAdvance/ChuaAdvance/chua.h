#pragma once
#include<stdio.h>
#include<math.h>
class ChuaF//标准蔡氏电路类
{
public:
	int N;
	double alpha,
		m0, m1,
		belta;
	void SetParameter(double Alpha, double Belta, double M1, double M0)
	{
		N = 3;
		alpha = Alpha; belta = Belta;
		m0 = M0; m1 = M1;
	}
	//double mL, double nC2, double nC1, double kR, double Ga, double Gb
	void SimuSetparameter(double* parameter)
	{
		N = 3;
		alpha = parameter[1] / parameter[2]; 
		belta = parameter[1] * parameter[3]*parameter[3] / parameter[0];
		m0 = parameter[3] * parameter[4]; m1 = parameter[3] * parameter[5];
	}
	double f(double* y, double t,int i)
	{
		switch (i)
		{
		case 0:return (alpha * (y[1] - y[0] - (m1*y[0] + 0.5*(m0 - m1)*(fabs(y[0] + 1) - fabs(y[0] - 1))))); break;
		case 1:return (y[0] - y[1] + y[2]); break;
		case 2:return  ((-1)*belta*y[1]); break;
		default:return 0;
			break;
		}
	}
	ChuaF() { N = 3; };
	~ChuaF() { };
private:

};


