#ifndef TRANSFORMEE_DE_FOURIER_RAPIDE_H_INCLUDED
#define TRANSFORMEE_DE_FOURIER_RAPIDE_H_INCLUDED

#include <string>
#include <vector>
#include "C:/Users/theob/OneDrive/Documents/Métier Quant/C++/Modele de Heston/FonctionsMathematiques/FonctionsMathematiques.h"

class Transformee_de_Fourier_rapide
{
	public:
		Transformee_de_Fourier_rapide(std::vector<double> const& signal);
		Transformee_de_Fourier_rapide(std::vector < NombreComplexe > const& signal);
		Transformee_de_Fourier_rapide compute_fft();
		Transformee_de_Fourier_rapide compute_inverse_fft();
		std::vector<NombreComplexe> _recursive_fft(std::vector<NombreComplexe> const& x);
		std::vector<NombreComplexe> signal_original() const;
		std::vector<NombreComplexe> signal_ralonge() const;
		void afficher_signal_original();
		void afficher_signal_ralonge();
		Transformee_de_Fourier_rapide produit_de_convolution(Transformee_de_Fourier_rapide fft_1);
		std::vector<NombreComplexe> fft_fractionnaire(double a);
		std::vector<NombreComplexe> transformee_de_Fourier_continue_methode_classique_via_fft(std::vector<double> const&f, double T);
//f est une liste à n éléments allant de - T / 2 à T / 2 (2k - n) / nT / 2 représentant une fonction f(x) définie sur l'intervalle [-T/2,T/2]
		std::vector<NombreComplexe> transformee_de_Fourier_continue_methode_via_fft_fractionnaire(std::vector<double> const& f, double a, double b);
//f est une liste à n éléments allant de - a / 2 à a / 2 (2k - n) / n a / 2 représentant une fonction f(x) définie sur l'intervalle [-a/2,a/2]
		~Transformee_de_Fourier_rapide();


	private:
		std::vector<NombreComplexe> m_signal_original;
		std::vector<NombreComplexe> m_signal_ralonge;
};
Transformee_de_Fourier_rapide operator*(Transformee_de_Fourier_rapide x, Transformee_de_Fourier_rapide y); //Extension du produit comme étant le produit de convolution
std::vector<NombreComplexe> operator*(std::vector<NombreComplexe> x, std::vector<NombreComplexe>  y); //Extension du produit comme étant le produit de convolution
std::vector<NombreComplexe> operator*(std::vector<double> x, std::vector<double>  y); //Extension du produit comme étant le produit de convolution


#endif // TRANSFORMEE_DE_FOURIER_RAPIDE_H_INCLUDED
