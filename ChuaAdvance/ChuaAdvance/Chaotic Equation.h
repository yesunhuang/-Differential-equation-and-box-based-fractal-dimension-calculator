#pragma once
#include<math.h>
#include<iostream>
#include<fstream>
#include"chua.h"
#define SETY for(i=0;i<N;i++)
using namespace std;
template<class Type>
class chaos
{
private:
	double Fractal(int Ne, int Neb, int beta)
	{
		double ChaosFractal;
		//cout << Ne << ' ' << Neb << ' ' << endl;
		ChaosFractal = (log(Ne)-log(Neb)) / log(beta);
		return ChaosFractal;
	}
	void clearM(int** r, int range)
	{
		for (int i = 0; i < range; i++)
			delete[] r[i];
		delete[] r;
	}
	bool Out(int x, int y, double*y0, int range, double Accuracy)
	{
		if ((int(ceil(fabs(y0[x] / Accuracy))) - 1 >= range) || (int(ceil(fabs(y0[y] / Accuracy))) - 1 >= range))
			return true;
		else
			return false;
	}
	bool check(int x, int y, double* y0, int** map, double Accuracy)
	{
		int x1 = int(ceil(fabs(y0[x] / Accuracy))) - 1;
		int y1 = int(ceil(fabs(y0[y] / Accuracy))) - 1;
		int Root;
		if (y0[x] > 0)
		{
			if (y0[y] > 0)
				Root = 2;
			else
				Root = 3;
		}
		else
		{
			if (y0[y] > 0)
				Root = 5;
			else
				Root = 7;
		}

		if ((map[x1][y1]) % (Root) != 0)
		{
			map[x1][y1] *= Root;
			return true;
		}
		else
			return false;
	}
	void initialize(int** Omap,int range)
	{
		for (int i = 0; i < range; i++) { 
			Omap[i] = new int[range]; 
			for (int j = 0; j < range; j++) Omap[i][j] = 1; 
		}
	}
public:
	Type *Equation;
	double *y,t=0;
	int N;
	chaos(Type* p)
	{
		Equation = p;
		N =Equation->N ; t = 0;
		y = new double[N];
	}
	~chaos()
	{
		delete[] y;
	}
	void SetY0(double *y0)
	{
		int i;
		for (i = 0; i < N; i++)
			y[i] = y0[i];
	}
	double ARK4FractalSolver(double *t0, double Accuracy, double IniStep, int range)
	{
		double *ys = new double[N];
		int i, enlarge,Ne=0,Neb=0,beta=2; 
		bool Over = false;
		double ts = t0[0],
			pace = IniStep;
		double(*k)[4] = new double[N][4];
		int** map = new int*[range];
		int** mapb = new int*[range / beta];
		this->initialize(map,range);
		this->initialize(mapb, range / beta);
		//ofstream data_file;
		//data_file.open("Adata.dat", ios::app);
		while (ts<t0[1])
		{
			Over = false;
			enlarge = 0;
			SETY ys[i] = y[i];
			SETY k[i][0] = Equation->f(ys, ts, i);
			SETY ys[i] = y[i] + k[i][0] * 0.5*pace;
			SETY k[i][1] =Equation->f(ys, ts + 0.5*pace,i);
			SETY ys[i] = y[i] + k[i][1] * 0.5*pace;
			SETY k[i][2] = Equation->f(ys, ts + 0.5*pace, i);
			SETY ys[i] = y[i] + pace * k[i][2];
			SETY k[i][3] = Equation->f(ys, ts + pace, i);
			SETY
			{
				ys[i] = (pace / 6)*(k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]);
				if (fabs(ys[i])>Accuracy)
				{
					Over = true;
					pace = 0.5*pace;
					break;
				}
				else
					if (fabs(ys[i])<(Accuracy / 4))
						enlarge++;
			}
			if (!Over)
			{
				SETY y[i] += ys[i];
				ts += pace;
				if (this->Out(0, 2, y, range, Accuracy))
				{
					delete[] ys;
					delete[] k;
					this->clearM(map, range);
					this->clearM(mapb, range / beta);
					return(this->Fractal(Ne, Neb, beta));
				}
				/*for (int i = 0; i < N; i++)
				{
					data_file << y[i] << ' ';
				}
				data_file << ts << endl;*/
				
				if (this->check(0, 2, y, map, Accuracy))
				{
					Ne++;
					//data_file << Ne << endl;
					if (this->check(0, 2, y, mapb, Accuracy*beta)) {
						Neb++;
						//data_file << Neb << endl;
					}
				}
			}
			if (enlarge == N) pace *= 2;
			if (pace<0.001*IniStep)
			{
				delete[] ys;
				delete[] k;
				this->clearM(map, range);
				this->clearM(mapb, range / beta);
				return(this->Fractal(Ne, Neb, beta));
			}
		}
		t = ts;
		//data_file.close();
		delete[] ys;
		delete[] k;
		this->clearM(map, range);
		this->clearM(mapb, range / beta);
		return(this->Fractal(Ne, Neb, beta));
	}
};
