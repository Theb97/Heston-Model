#ifndef CLASSE_DU_MODELE_DE_HESTON_H_INCLUDED
#define CLASSE_DU_MODELE_DE_HESTON_H_INCLUDED

#include <string>
#include <vector>
#include "C:/Modele de Heston/FonctionsMathematiques/FonctionsMathematiques.h"
#include "C:/Modele de Heston/Transformee_de_Fourier_rapide/Transformee_de_Fourier_rapide.h"

class modele_de_Heston
{
public:
	modele_de_Heston();
	modele_de_Heston(double kappa, double theta, double sigma, double rho, double V_0, double r, double lambda, double T, double t, double S_t);
	NombreComplexe d_2(NombreComplexe const& phi);
	NombreComplexe d_2(double const& phii);
	NombreComplexe g_2(NombreComplexe const& phi);
	NombreComplexe g_2(double const& phii);
	NombreComplexe D_2_tau_phi(NombreComplexe const& phi);
	NombreComplexe D_2_tau_phi(double const& phii);
	NombreComplexe C_2_tau_phi(NombreComplexe const& phi);
	NombreComplexe C_2_tau_phi(double const& phii);
	NombreComplexe fonction_caracteristique_modele_de_Heston_2(NombreComplexe const& phi);
	NombreComplexe fonction_caracteristique_modele_de_Heston_2(double const& phii);
	NombreComplexe psi(double const& phi, double const& alpha); // voir article Carr et Madan pour la definition de la fonction psi
	NombreComplexe zeta(NombreComplexe const& phi); // voir article Carr et Madan pour la definition de la fonction zeta
	NombreComplexe zeta(double const& phii);
	NombreComplexe gamma(double const& phi, double const& alpha); // voir article Carr et Madan pour la definition de la fonction gamma
	~modele_de_Heston();


private:
	double m_kappa;
	double m_theta;
	double m_sigma;
	double m_rho;
	double m_V_0;
	double m_r;
	double m_lambda;
	double m_T;
	double m_t;
	double m_S_t;

};

#endif // CLASSE_DU_MODELE_DE_HESTON_H_INCLUDED
