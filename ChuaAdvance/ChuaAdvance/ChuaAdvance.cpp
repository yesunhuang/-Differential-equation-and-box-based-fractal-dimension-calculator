// ChuaAdavance.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<math.h>
#include<iostream>
#include<fstream>
using namespace std;
//double mL, double nC2, double nC1, double kR, double Ga, double Gb
/*double mL0=15,
nC10=10,
nC20=100,
nkR0=1.5,
Ga0=-0.8,
Gb0=-0.4;
*/
double parameter0[6] = { 15,100,10,1.5,-0.8,-0.4 };
double Yin[3] = { 0.1,0.1,0.1 }, tspan[2] = { 0,300 };
int min(int a, int b)
{
	if (a > b)
		return b;
	else
		return a;
}
void VariationalFractal(double rangeUp, double rangeDown, double accuracy, int code)
{
	double Fractal, parameter[6]; int i;
	ofstream data_file;
	switch (code)
	{
	case 0:data_file.open("Ldata.dat", ios::app); break;
	case 1:data_file.open("C2data.dat", ios::app); break;
	case 2:data_file.open("C1data.dat", ios::app); break;
	case 3:data_file.open("Rdata.dat", ios::app); break;
	case 4:data_file.open("GAdata.dat", ios::app); break;
	default:data_file.open("GBdata.dat", ios::app); break;
	}
	char* KProducer(int, chaos<ChuaF>*, double, char*);
	ChuaF Chua;
	chaos<ChuaF> Ckey(&Chua);
	for (i = 0; i < 6; i++) parameter[i] = parameter0[i];
	for (i = -(int)(floor(rangeUp / accuracy)); i < (int)(ceil(rangeDown / accuracy)); i++)
	{
		parameter[code] = parameter0[code] + i * accuracy;
		Chua.SimuSetparameter(parameter);
		Ckey.SetY0(Yin);
		Fractal = Ckey.ARK4FractalSolver(tspan, 0.01, 0.001, 1000);
		data_file << parameter[code] << ' ' << Fractal << endl;
	}
	data_file.close();
}
void VariationalFractal(double* rangeUp, double* rangeDown, double* accuracy, int* code)
{
	double Fractal, parameter[6]; int i, j;
	ofstream data_file;
	char* KProducer(int, chaos<ChuaF>*, double, char*);
	ChuaF Chua;
	chaos<ChuaF> Ckey(&Chua);
	for (i = 0; i < 6; i++) parameter[i] = parameter0[i];
	data_file.open("2Ddata.dat", ios::out|ios::trunc);
	data_file << (int)(floor(rangeUp[0] / accuracy[0])) + (int)(ceil(rangeDown[0] / accuracy[0]));
	data_file << ' ' << (int)(floor(rangeUp[1] / accuracy[1]))+(int)(ceil(rangeDown[1] / accuracy[1]));
	data_file << " 0" << endl;
	for (i = -(int)(floor(rangeUp[0] / accuracy[0])); i < (int)(ceil(rangeDown[0] / accuracy[0])); i++)
	{
		parameter[code[0]] = parameter0[code[0]] + i * accuracy[0];
		for (j = -(int)(floor(rangeUp[1] / accuracy[1])); j < (int)(ceil(rangeDown[1] / accuracy[1])); j++)
		{
			parameter[code[1]] = parameter0[code[1]] + j * accuracy[1];
			Chua.SimuSetparameter(parameter);
			Ckey.SetY0(Yin);
			Fractal = Ckey.ARK4FractalSolver(tspan, 0.01, 0.001, 1000);
			data_file << parameter[code[0]] << ' ' << parameter[code[1]] << ' ' << Fractal << endl;
		}
	}
	data_file.close();
}
int main()
{
	double up[2] = { 1.5,9}; double down[2] = { 2, 25 };
	double accuracy[2] = { 0.01,0.1 }; int code[2] = { 3,2 };
	VariationalFractal(up,down,accuracy,code);
}
