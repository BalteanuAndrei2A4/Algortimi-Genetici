#include <iostream>
#include <math.h>
#include <random>
#include <chrono>
using namespace std;

mt19937_64 g_randomGenerator;
float Jong(float* valori, int D)
{
	float s = 0;
	for (int i = 0; i < D; i++)
	{
		s = s + valori[i] * valori[i];
	}
	return s;
}
float decodeDimesion(bool* vc, int begin, double a, double b, int l)
{
	float s = 0;
	float ret = 0;
	int putere = pow(2, l);
	for (int i = begin; i != begin + l; i++)
	{
		s = s * 2;
		s = s + vc[i];
	}
	ret = s / putere;
	ret = ret * (b - a) + a;

	return ret;
}
float* decodeSolution(bool* vc, int begin, double a, double b, int l, int D)
{
	float valori[500];
	for (int i = 0; i < D; i++)
	{
		valori[i] = decodeDimesion(vc, begin, a, b, l);
		begin = begin + l;
	}
	return valori;
}
bool* copiere_vector(bool* sursa, bool* dest, double L)
{
	for (int i = 0; i < L; i++)
	{
		dest[i] = sursa[i];
	}
	return dest;
}
void initializareRandomGenerator()
{
	mt19937_64 initializer;
	initializer.seed(static_cast<long unsigned int>(chrono::high_resolution_clock::now().time_since_epoch().count()));
	initializer.discard(10000);
	g_randomGenerator.seed(initializer());
}
void generateRandomVector(int L, bool vc[])
{
	for (int i = 0; i < L; i++)
	{
		vc[i] = g_randomGenerator() % 2;
	}
}
int generareRandomPozitie(int L)
{
	return g_randomGenerator() % L;
}
int main()
{
	int begin = 0;
	long long int seed = 0;
	double a = -5.12;
	double b = 5.12;
	int D = 5;//aici se inlocuieste cu 10 sau 30 in functie de cate dimensiuni are functia
	long double L = b - a;
	L = L * pow(10, 5);
	L = log2(L);
	L = long long int(L * D) + 1;
	int l = L / D;
	bool vc[1000] = { 0 };
	bool vbest[1000] = { 0 };

	initializareRandomGenerator();
	generateRandomVector(L, vc);
	copiere_vector(vc, vbest, L);
	long double best = Jong(decodeSolution(vc, begin, a, b, l, D), D);


	for (int i = 0; i < 1000; i++)
	{
		generateRandomVector(L, vc);
		float ec = Jong(decodeSolution(vc, begin, a, b, l, D), D);
		bool localBest[1000];
		copiere_vector(vc, localBest, L);
		bool local = false;
		while (local == false)
		{
			vc[0] = !vc[0];
			float eLocalBest = Jong(decodeSolution(vc, begin, a, b, l, D), D);
			vc[0] = !vc[0];
			for (int j = 0; j < L; j++)
			{
				vc[j] = !vc[j];
				float en = Jong(decodeSolution(vc, begin, a, b, l, D), D);
				if (en < eLocalBest)
				{
					eLocalBest = en;
					copiere_vector(vc, localBest, L);
				}
				vc[j] = !vc[j];
			}
			if (eLocalBest < ec)
			{
				ec = eLocalBest;
				copiere_vector(localBest, vc, L);
			}
			else
			{
				local = true;
			}

		}
		cout << i << ":" << ec << endl;
		if (ec < best)
		{
			best = ec;
			copiere_vector(vc, vbest, L);
		}
	}


	return 0;
}