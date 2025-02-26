// voir article de Carr et Madan (Option valuation using the fast Fourier transform) pour les méthodes en fin de code 


#include "classe_du_modele_de_Heston.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <float.h>


using namespace std;

const double PI = 3.1415926535;

modele_de_Heston::modele_de_Heston() : m_kappa(1), m_theta(0.16), m_sigma(0.5), m_rho(-0.8), m_V_0(0.16), m_r(0), m_lambda(0), m_T(10), m_t(0), m_S_t(1)
{

}
modele_de_Heston::modele_de_Heston(double kappa, double theta, double sigma, double rho, double V_0, double r, double lambda, double T, double t, double S_t) : m_kappa(kappa), m_theta(theta), m_sigma(sigma), m_rho(rho), m_V_0(V_0), m_r(r), m_lambda(lambda), m_T(T), m_t(t), m_S_t(S_t)
{

}

NombreComplexe modele_de_Heston::d_2(NombreComplexe const& phi)
{
	double u_2 = -0.5;
	double b_2 = m_kappa + m_lambda;
	NombreComplexe i = NombreComplexe(0, 1);
	NombreComplexe d_2 = ((m_rho * m_sigma * phi * i - b_2)*(m_rho * m_sigma * phi * i - b_2)-m_sigma*m_sigma*(2*u_2*phi*i-phi*phi)).racine_carre_complexe();
	return d_2;
}
NombreComplexe modele_de_Heston::d_2(double const& phii)
{
	NombreComplexe phi = NombreComplexe(phii, 0);
	return d_2(phi);
}
NombreComplexe modele_de_Heston::g_2(NombreComplexe const& phi)
{
	double u_2 = -0.5;
	double b_2 = m_kappa + m_lambda;
	NombreComplexe i = NombreComplexe(0, 1);
	NombreComplexe d = this->d_2(phi);
	NombreComplexe g_2 = (b_2 - m_rho * m_sigma * phi * i + d) / (b_2 - m_rho * m_sigma * phi * i - d);
	return g_2;
}

NombreComplexe modele_de_Heston::g_2(double const& phii)
{
	NombreComplexe phi = NombreComplexe(phii, 0);
	return g_2(phi);
}

NombreComplexe modele_de_Heston::D_2_tau_phi(NombreComplexe const& phi)
{

	if (phi.module() < DBL_EPSILON)
		return NombreComplexe(0, 0);
	else
	{
		double tau = m_T - m_t;
		double b_2 = m_kappa + m_lambda;
		NombreComplexe i = NombreComplexe(0, 1);
		NombreComplexe d = this->d_2(phi);
		NombreComplexe g = this->g_2(phi);
		NombreComplexe D_tau_phi = NombreComplexe();
		if (d.partie_reelle() * tau >= -log(DBL_EPSILON))
		{
			D_tau_phi = (b_2 - m_rho * m_sigma * phi * i + d) / m_sigma / m_sigma / g;
		}
		else
		{
			D_tau_phi = (b_2 - m_rho * m_sigma * phi * i + d) / m_sigma / m_sigma * (1 - (d * tau).expz()) / (1 - g * ((d * tau).expz()));
		}


		return D_tau_phi;
	}
	
}

NombreComplexe modele_de_Heston::D_2_tau_phi(double const& phii)
{
	NombreComplexe phi = NombreComplexe(phii, 0);
	return D_2_tau_phi(phi);
}

NombreComplexe modele_de_Heston::C_2_tau_phi(NombreComplexe const& phi)
{

	if (phi.module() < DBL_EPSILON)
		return NombreComplexe(0, 0);
	else
	{
		double tau = m_T - m_t;
		double b_2 = m_kappa + m_lambda;
		NombreComplexe i = NombreComplexe(0, 1);
		NombreComplexe d = this->d_2(phi);
		NombreComplexe g = this->g_2(phi);
		double t_g = g.argument();
		NombreComplexe G_D = g - 1;
		int mm = static_cast<int>((t_g + PI) / (2 * PI));
		double m = static_cast<double>(mm);
		NombreComplexe G_N = NombreComplexe();
		NombreComplexe ln_G = NombreComplexe();
		if (d.partie_reelle() * tau >= -log(DBL_EPSILON))
		{
			G_N = g * ((d * tau).expz());
			ln_G = log((g / (g - 1)).module()) + i * (g / (g - 1)).argument() + d * tau;

		}
		else
		{
			G_N = g * ((d * tau).expz()) - 1;
			int nn = static_cast<int>((t_g + (d.partie_imaginaire()) * tau + PI) / (2 * PI));
			double n = static_cast<double>(nn);
			ln_G = log(G_N.module() / G_D.module()) + i * (G_N.argument() - G_D.argument() + 2 * PI * (n - m));
		}
		NombreComplexe C_tau_phi = (m_r)*phi * i * tau + m_kappa * m_theta / m_sigma / m_sigma * ((b_2 - m_rho * m_sigma * phi * i + d) * tau - 2 * ln_G);
		return C_tau_phi;
	}
	
}

NombreComplexe modele_de_Heston::C_2_tau_phi(double const& phii)
{
	NombreComplexe phi = NombreComplexe(phii, 0);
	return C_2_tau_phi(phi);
}

NombreComplexe modele_de_Heston::fonction_caracteristique_modele_de_Heston_2(NombreComplexe const& phi)
{
	double tau = m_T - m_t;
	double x_t = log(m_S_t) + m_r * tau;
	NombreComplexe i(0, 1);
	return (C_2_tau_phi(phi) + D_2_tau_phi(phi) * m_V_0 + i * phi * x_t).expz();

}

NombreComplexe modele_de_Heston::fonction_caracteristique_modele_de_Heston_2(double const& phii)
{
	NombreComplexe phi(phii, 0);
	return fonction_caracteristique_modele_de_Heston_2(phi);
}

// voir article Carr et Madan pour la definition de la fonction psi
NombreComplexe modele_de_Heston::psi(double const& phi, double const& alpha)
{
	double tau = m_T - m_t;
	NombreComplexe i(0, 1);
	return (exp(-m_r*tau)*fonction_caracteristique_modele_de_Heston_2(phi-(alpha+1)*i))/(alpha*alpha+alpha-phi*phi+i*phi*(2*alpha+1));
}
// voir article Carr et Madan pour la definition de la fonction zeta
NombreComplexe modele_de_Heston::zeta(NombreComplexe const& phi)
{
	double tau = m_T - m_t;
	NombreComplexe i(0, 1);
	return ((exp(-m_r * tau)) * (1 / (1 + i * phi) - exp(m_r * tau) / (i * phi) - fonction_caracteristique_modele_de_Heston_2(phi - i) / (phi * phi - i * phi)));
}
NombreComplexe modele_de_Heston::zeta(double const& phii)
{
	NombreComplexe phi(phii, 0);
	return zeta(phi);
}

// voir article Carr et Madan pour la definition de la fonction gamma
NombreComplexe modele_de_Heston::gamma(double const& phi,double const& alpha)
{
	NombreComplexe i(0, 1);
	return ((zeta(phi - i * alpha) - zeta(phi + i * alpha))/2);
}


modele_de_Heston::~modele_de_Heston()
{
	/* Rien à mettre ici car on ne fait pas d'allocation dynamique
	dans la classe modele_de_Heston. Le destructeur est donc inutile mais
	je le mets pour montrer à quoi cela ressemble.
	En temps normal, un destructeur fait souvent des delete et quelques
	autres vérifications si nécessaire avant la destruction de l'objet. */
}




