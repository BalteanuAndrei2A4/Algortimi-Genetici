#include <iostream>

using namespace std;
struct variabila {
	int coeficient_variabila;
	int putere_variabila;// aici doar 1 sau 2-
	char operatie;
}parametru[50];
int main()
{
	int minim = INT_MAX;
	int maxim = INT_MIN;
	int numar_variabile,numar_elemente_codomeniu = 0;
	int vector[20] = { 0 };
	cout << "INTRODUCETI NUMARUL DE VARIABILE:";
	cin >> numar_variabile;
	cout << endl;
	for (int i = 0; i < numar_variabile; i++)
	{
		cout << "INTRODUCETI PUTEREA VARIABILEI CU NUMARUL " << i + 1<<":";
		cin >> parametru[i].putere_variabila;
		cout << "INTRODUCETI COEFICIENTUL VARIABILEI CU NUMARUL " << i + 1<<":";
		cin >> parametru[i].coeficient_variabila;
		cout << "SEMNUL DIN FATA VARIABILEI CU NUMARUL " << i + 1 << ":";
		cin >> parametru[i].operatie;
	}
	while (vector[numar_variabile - 1] != 16)
	{
		int valoare_functie=0;
		for (int i = 0; i < numar_variabile; i++)
		{
			if (parametru[i].operatie == '+')
			{
				valoare_functie = valoare_functie + pow(vector[i], parametru[i].putere_variabila) * parametru[i].coeficient_variabila;
			}
			else if (parametru[i].operatie == '-')
			{
				valoare_functie = valoare_functie - pow(vector[i], parametru[i].putere_variabila) * parametru[i].coeficient_variabila;

			}
		}
		if (valoare_functie > maxim)
			maxim = valoare_functie;
		if (valoare_functie < minim)
			minim = valoare_functie;

		vector[0]++;
		for (int i = 0; i < numar_variabile - 1; i++)
		{
			if (vector[i] == 16)
			{
				vector[i] = 0;
				vector[i + 1]++;
			}

		}

		
	}
	cout << endl << "minimul functiei este egal cu: " << minim<<endl;
	cout<<"maximul functiei este egal cu: " << maxim;
	
	


	return 0;

}