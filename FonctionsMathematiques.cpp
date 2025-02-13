#include "FonctionsMathematiques.h"
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

int ajouteDeux(int nombreRecu)
{
    int valeur(nombreRecu + 2);

    return valeur;
}

void afficherMessage(const string & message="")
{
    cout << message << std::endl;
}

const double PI = 3.1415926535;

NombreComplexe::NombreComplexe() : m_partie_reelle(1), m_partie_imaginaire(0)
{
    
    //Rien à mettre dans le corps du constructeur, tout a déjà été fait !
}

NombreComplexe::NombreComplexe(double partie_reelle, double partie_imaginaire):m_partie_reelle(partie_reelle), m_partie_imaginaire(partie_imaginaire)
{
    //Rien à mettre dans le corps du constructeur, tout a déjà été fait !
}

double NombreComplexe::module() const
{
    double a = m_partie_reelle;
    double b = m_partie_imaginaire;
    return sqrt(a * a + b * b);
}
double NombreComplexe::partie_reelle() const
{
    return m_partie_reelle;
}
double NombreComplexe::partie_imaginaire() const
{
    return m_partie_imaginaire;
}

void NombreComplexe::afficher(ostream& flux) const
{
    flux << m_partie_reelle << " + " << m_partie_imaginaire << "i" ;
}

NombreComplexe NombreComplexe::conjuguee() const
{
    return NombreComplexe(m_partie_reelle, -m_partie_imaginaire);
}

NombreComplexe NombreComplexe::addition(NombreComplexe const& T) const
{
    double a = m_partie_reelle + T.partie_reelle();
    double b = m_partie_imaginaire + T.partie_imaginaire();
    return NombreComplexe(a, b);
}

NombreComplexe NombreComplexe::soustraction(NombreComplexe const& S) const
{
    double a = m_partie_reelle - S.partie_reelle();
    double b = m_partie_imaginaire - S.partie_imaginaire();
    return NombreComplexe(a, b);
}

NombreComplexe NombreComplexe::lineaire(double t)
{
    return NombreComplexe(t * m_partie_reelle, t * m_partie_imaginaire);
}

NombreComplexe NombreComplexe::produit(NombreComplexe Q) const
{
    double a = m_partie_reelle * Q.partie_reelle() - m_partie_imaginaire * Q.partie_imaginaire();
    double b = m_partie_imaginaire * Q.partie_reelle() + m_partie_reelle * Q.partie_imaginaire();
    return NombreComplexe(a, b);
}

NombreComplexe NombreComplexe::division(NombreComplexe R) const
{
    NombreComplexe Q = NombreComplexe(m_partie_reelle, m_partie_imaginaire);
    Q = Q.produit(R.conjuguee());
    double t = 1 / (R.partie_reelle() * R.partie_reelle() + R.partie_imaginaire() * R.partie_imaginaire());
    return Q.lineaire(t);
}



NombreComplexe NombreComplexe::expz() const
{
    double a = exp(m_partie_reelle);
    return NombreComplexe(a * cos(m_partie_imaginaire), a * sin(m_partie_imaginaire));
}

NombreComplexe NombreComplexe::logz() const
{
    double c = (*this).module();
    double d = 2 * atan(m_partie_imaginaire / (m_partie_reelle + c));
    return NombreComplexe(log(c), d);
}

NombreComplexe NombreComplexe::racine_carre_complexe() const
{
    double a = 0;
    double b = 0;
    if (((*this).m_partie_imaginaire)*((*this).m_partie_imaginaire) < 0.01 && (*this).m_partie_reelle < 0)
    {
        b = sqrt((*this).module());
        a = 0;   
    }
    else if (((*this).m_partie_imaginaire) * ((*this).m_partie_imaginaire) < 0.01 && (*this).m_partie_reelle < 0.01 && (*this).m_partie_reelle > 0 )
    {
        a = 0;
        b = 0;
    }
    else if (((*this).m_partie_imaginaire) * ((*this).m_partie_imaginaire) < 0.01 &&  (*this).m_partie_reelle > 0.01)
    {
        a = sqrt((*this).module());
        b = 0;
    }
    else
    {
        double c = (*this).module();
        double e = sqrt(c);
        double d = atan(m_partie_imaginaire / (m_partie_reelle + c));

        a = e * cos(d);
        b = e * sin(d);
    }
    return NombreComplexe(a,b) ;
}

NombreComplexe NombreComplexe::puissance_complexe(int exposant) const
{
    NombreComplexe resultat = NombreComplexe(1, 0);
    NombreComplexe base = (*this);
    while (exposant > 0)
    {
        if (exposant % 2 == 1)
        {
            resultat *= base;
            exposant -= 1;
        }
        base *= base;
        exposant /= 2;
    }
    return resultat;
}

double NombreComplexe::argument() const
{
    if (m_partie_reelle == 0 && m_partie_imaginaire == 0) {
        throw std::invalid_argument("z ne doit pas être nul pour que l'argument de z soit bien défini.");
    }
    else if (m_partie_reelle < 0 && m_partie_imaginaire == 0) {
        return -PI;
    }
    else {
        double c = module();
        double d = 2 * atan(m_partie_imaginaire / (m_partie_reelle + c));
        return d;
    }
}

bool NombreComplexe::estEgal(NombreComplexe const& z_2) const
{
    //Teste si a.m_heure == b.m_heure etc.  
    if (m_partie_reelle == z_2.m_partie_reelle && m_partie_imaginaire == z_2.m_partie_imaginaire)
        return true;
    else
        return false;
}

NombreComplexe& NombreComplexe::operator+=(NombreComplexe const& T)
{
    m_partie_reelle +=  T.partie_reelle();
    m_partie_imaginaire += T.partie_imaginaire();
    return *this;
}

NombreComplexe& NombreComplexe::operator-=(NombreComplexe const& T)
{
    m_partie_reelle -= T.partie_reelle();
    m_partie_imaginaire -= T.partie_imaginaire();
    return *this;
}
NombreComplexe& NombreComplexe::operator*=(NombreComplexe const& Q)
{
    double a = m_partie_reelle * Q.partie_reelle() + m_partie_imaginaire * Q.partie_imaginaire();
    double b = m_partie_imaginaire * Q.partie_reelle() - m_partie_reelle * Q.partie_imaginaire();
    m_partie_reelle = a;
    m_partie_imaginaire = b;
    return *this;
}

NombreComplexe& NombreComplexe::operator/=(NombreComplexe const& Q)
{
    double a = m_partie_reelle * Q.partie_reelle() - m_partie_imaginaire * Q.partie_imaginaire();
    double b = m_partie_imaginaire * Q.partie_reelle() + m_partie_reelle * Q.partie_imaginaire();
    m_partie_reelle = a / (Q.partie_reelle() * Q.partie_reelle() + Q.partie_imaginaire() * Q.partie_imaginaire());
    m_partie_imaginaire = b / (Q.partie_reelle() * Q.partie_reelle() + Q.partie_imaginaire() * Q.partie_imaginaire());
    return *this;
}

NombreComplexe& NombreComplexe::operator=(NombreComplexe const& NombreComplexeAcopier)
{
    if (this != &NombreComplexeAcopier)
        //On vérifie que l'objet n'est pas le même que celui reçu en argument
    {
        m_partie_reelle = NombreComplexeAcopier.m_partie_reelle; //On copie tous les champs
        m_partie_imaginaire = NombreComplexeAcopier.m_partie_imaginaire;
        
    }
    return *this; //On renvoie l'objet lui-même
}

NombreComplexe::~NombreComplexe()
{
    /* Rien à mettre ici car on ne fait pas d'allocation dynamique
    dans la classe NombreComplexe. Le destructeur est donc inutile mais
    je le mets pour montrer à quoi cela ressemble.
    En temps normal, un destructeur fait souvent des delete et quelques
    autres vérifications si nécessaire avant la destruction de l'objet. */
}

NombreComplexe expi(double theta)
{
    return NombreComplexe(cos(theta), sin(theta));
}

bool operator==(NombreComplexe const& z_1, NombreComplexe const& z_2) 
{
    return z_1.estEgal(z_2);
}

bool operator!=(NombreComplexe const& z_1, NombreComplexe const& z_2)
{
    return not (z_1 == z_2);
}

ostream &operator<<(ostream& flux, NombreComplexe const& z)
{
    z.afficher(flux);
    return flux;
}

NombreComplexe operator+(NombreComplexe const& z_1, NombreComplexe const& z_2)
{
    NombreComplexe resultat(z_1);
    resultat += z_2;
    return resultat ;
}

NombreComplexe operator-(NombreComplexe const& z_1, NombreComplexe const& z_2)
{
    NombreComplexe resultat(z_1);
    resultat -= z_2;
    return resultat;
}

NombreComplexe operator*(NombreComplexe const& z_1, NombreComplexe const& z_2)
{
    
    NombreComplexe resultat = z_1.produit(z_2);
    return resultat;
}

NombreComplexe operator/(NombreComplexe const& z_1, NombreComplexe const& z_2)
{

    NombreComplexe resultat = z_1.division(z_2);
    return resultat;
}



