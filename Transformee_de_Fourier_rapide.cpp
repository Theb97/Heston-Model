#include "Transformee_de_Fourier_rapide.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>


using namespace std;

const double PI = 3.1415926535;


Transformee_de_Fourier_rapide::Transformee_de_Fourier_rapide(vector<double> const& signal)
{
	

	size_t q(signal.size());
	unsigned int b = 0;
	unsigned int puissance_de_2 = 1;
	while (puissance_de_2 < q)
	{
		b += 1;
		puissance_de_2 *= 2;
	}
	if (puissance_de_2 != q)
	{
		unsigned int compteur(0);
		unsigned int j(0);
		for (j = 0; j < q; j++)
		{
			m_signal_original.push_back(NombreComplexe(signal[j], 0));
		}
		
		m_signal_ralonge = m_signal_original;
		for (compteur = 0; compteur < puissance_de_2 - q; compteur++)
		{
			m_signal_ralonge.push_back(NombreComplexe(0, 0));
		}
	}
	else if (puissance_de_2 == q)
	{
		unsigned int j(0);
		for (j = 0; j < q; j++)
		{
			m_signal_original.push_back(NombreComplexe(signal[j], 0));
		}
		m_signal_ralonge = m_signal_original;
	}


}


Transformee_de_Fourier_rapide::Transformee_de_Fourier_rapide(vector<NombreComplexe> const& signal)
{
	

	size_t q(signal.size());
	unsigned int b = 0;
	unsigned int puissance_de_2 = 1;
	while (puissance_de_2 < q)
	{
		b += 1;
		puissance_de_2 *= 2;
	}
	if (puissance_de_2 != q)
	{
		m_signal_original = signal;
		unsigned int compteur(0);
		m_signal_ralonge = m_signal_original;
		for (compteur = 0; compteur < puissance_de_2 - q; compteur++)
		{
			m_signal_ralonge.push_back(NombreComplexe(0, 0));
		}
	}
	else if (puissance_de_2 == q)
	{
		m_signal_original = signal;
		m_signal_ralonge = m_signal_original;
	}
	

	

}

Transformee_de_Fourier_rapide Transformee_de_Fourier_rapide::compute_fft()
{
	return Transformee_de_Fourier_rapide(_recursive_fft(m_signal_ralonge));
}

Transformee_de_Fourier_rapide Transformee_de_Fourier_rapide::compute_inverse_fft()
{
	size_t n = m_signal_ralonge.size();
	double nn = static_cast<double>(n);
	vector <NombreComplexe> L;
	for (int k = 0; k < n; k++)
	{
		L.push_back(m_signal_ralonge[k].conjuguee()*(1/nn));
	}
	vector<NombreComplexe> inverse = _recursive_fft(L);
	for (int i = 0; i < L.size(); i++)
	{
		inverse[i] = inverse[i].conjuguee();
	}

	return Transformee_de_Fourier_rapide(inverse);
}

 vector<NombreComplexe> Transformee_de_Fourier_rapide::_recursive_fft(vector<NombreComplexe> const& x)
{
	 size_t n(x.size());
	if (n <= 1)
	{
		return x;
	}
	else
	{
		// Séparer les parties paires et impaires
		vector<NombreComplexe> even(n / 2);
		vector<NombreComplexe> odd(n / 2);

		for (int i = 0; i < n / 2; i++) {
			even[i] = x[i * 2];
			odd[i] = x[i * 2 + 1];
		}

		// Appels récursifs
		even = _recursive_fft(even);
		odd = _recursive_fft(odd);

		// Calcul de la transformée FFT
		vector<NombreComplexe> T(n / 2);
		vector<NombreComplexe> resultat(n);

		for (int k = 0; k < n / 2; k++)
		{
			T[k] = expi(-2 * PI * k / n) * odd[k];
			resultat[k] = even[k] + T[k];
			resultat[k + n / 2] = even[k] - T[k];
		}

		return resultat;
	}
	
}
 vector<NombreComplexe> Transformee_de_Fourier_rapide::signal_original() const
 {
	 return m_signal_original;
 }

 vector<NombreComplexe> Transformee_de_Fourier_rapide::signal_ralonge() const
 {
	 return m_signal_ralonge;
 }

void Transformee_de_Fourier_rapide::afficher_signal_original()
{
	unsigned int j(0);
	for (j = 0; j < m_signal_original.size(); j++)
	{
		cout << m_signal_original[j] << endl;
	}
}

void Transformee_de_Fourier_rapide::afficher_signal_ralonge()
{
	unsigned int j(0);
	for (j = 0; j < m_signal_ralonge.size(); j++)
	{
		cout << m_signal_ralonge[j] << endl;
	}
}

Transformee_de_Fourier_rapide Transformee_de_Fourier_rapide::produit_de_convolution(Transformee_de_Fourier_rapide fft_1)
{
	vector<NombreComplexe> x;
	vector<NombreComplexe> y;
	x = _recursive_fft((*this).m_signal_ralonge);
	y = _recursive_fft(fft_1.m_signal_ralonge);
	if (x.size() < y.size())
	{
		for (int k = 0; k < y.size() - x.size(); k++)
		{
			x.push_back(NombreComplexe(0, 0));
		}
	}
	else if (y.size() < x.size())
	{
		for (int k = 0; k < x.size() - y.size(); k++)
		{
			y.push_back(NombreComplexe(0, 0));
		}
	}
	vector<NombreComplexe> prod(x.size(), NombreComplexe(0, 0));
	double u = static_cast<double>(x.size());
	double t(1 / u);
	for (int k = 0; k < y.size(); k++)
	{
		prod[k] = ((x[k] * y[k]).conjuguee())*t;
	}
	vector<NombreComplexe> L = _recursive_fft(prod);
	for (int i = 0; i < L.size(); i++)
	{
		L[i] = L[i].conjuguee();
	}

	return Transformee_de_Fourier_rapide(L);
}

vector<NombreComplexe> Transformee_de_Fourier_rapide::fft_fractionnaire(double a)
{
	vector<NombreComplexe> X = this->m_signal_ralonge;
	size_t n = X.size();
	vector<NombreComplexe> y;
	vector<NombreComplexe> z;
	for (int k = 0; k < n; k++)
	{
		y.push_back(X[k]*expi(-PI*k*k*a));
		z.push_back(expi(PI * k * k * a));
	}
	y.push_back(NombreComplexe(0,0));
	z.push_back(NombreComplexe(0, 0));
	for (int j = n - 1; 0; j--)
	{
		y.push_back(NombreComplexe(0, 0));
		z.push_back(z[j]);
	}
	vector<NombreComplexe> M = y * z; // *= produit de convolution entre deux tableaux de NombreComplexe de taille 2n
	for (int j = 0; j < n; j++)
	{
		M.pop_back();
		M[j] *= expi(-PI * j * j * a);
	}
	return M;
}

vector<NombreComplexe> Transformee_de_Fourier_rapide::transformee_de_Fourier_continue_methode_classique_via_fft(vector<double> const& f, double T)
//f est une liste à n éléments allant de - T / 2 à T / 2 (2k - n) / nT / 2 représentant une fonction f(x) définie sur l'intervalle [-T/2,T/2]
{
	vector<double> M;
	for (size_t k = 1; k < f.size(); k += 2) {
		M.push_back(f[k-1]);
		M.push_back(-f[k]);
	}
	vector<NombreComplexe> N = Transformee_de_Fourier_rapide(M).compute_fft().m_signal_ralonge;
	for (size_t k = 1; k < f.size(); k += 2) {
		M.pop_back();
		M.pop_back();
		N[k] *= -1;
	}
	NombreComplexe y = NombreComplexe(0, -1).puissance_complexe(f.size());
	for (int j = 0; j < f.size(); j++)
	{
		N[j] *= NombreComplexe(T / f.size(), 0) * y;
	}
	return N;
}

vector<NombreComplexe> Transformee_de_Fourier_rapide::transformee_de_Fourier_continue_methode_via_fft_fractionnaire(vector<double> const& f, double a, double b)
//f est une liste à n éléments allant de -a/2 à a/2 (2k-n)/na/2 représentant une fonction f(x) définie sur l'intervalle [-a/2,a/2]
{
	vector<NombreComplexe> M;
	for (size_t k = 1; k < f.size(); k++) {
		M.push_back(f[k]*expi(a*b*k/(2*f.size())));
	}
	vector<NombreComplexe> N = Transformee_de_Fourier_rapide(M).fft_fractionnaire(a*b/(2*PI* f.size()* f.size()));
	for (size_t k = 1; k < f.size(); k++) {
		M.pop_back();
		N[k] *= ((a/ f.size())*expi(a/2*(-b/2+k*b/ f.size())));
	}
	return N;
}

Transformee_de_Fourier_rapide::~Transformee_de_Fourier_rapide()
{

}

Transformee_de_Fourier_rapide operator*(Transformee_de_Fourier_rapide  x, Transformee_de_Fourier_rapide  y)
{
	return x.produit_de_convolution(y);
}

vector<NombreComplexe> operator*(vector<NombreComplexe> x, vector<NombreComplexe>  y)
{
	Transformee_de_Fourier_rapide fft_x = Transformee_de_Fourier_rapide(x);
	Transformee_de_Fourier_rapide fft_y = Transformee_de_Fourier_rapide(y);
	
	return (fft_x.produit_de_convolution(fft_y)).signal_ralonge();
}

vector<NombreComplexe> operator*(vector<double> x, vector<double>  y)
{
	Transformee_de_Fourier_rapide fft_x = Transformee_de_Fourier_rapide(x);
	Transformee_de_Fourier_rapide fft_y = Transformee_de_Fourier_rapide(y);

	return (fft_x.produit_de_convolution(fft_y)).signal_ralonge();
}

