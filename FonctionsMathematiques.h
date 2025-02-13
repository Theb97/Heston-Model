#ifndef FONCTIONSMATHEMATIQUES_H_INCLUDED
#define FONCTIONSMATHEMATIQUES_H_INCLUDED

#include <string>

int ajouteDeux(int nombreRecu);

void afficherMessage(const std::string & message);

// fonction qui ajoute deux au nombre reçu en entree

class NombreComplexe
{
    // Tout ce qui suit est public (accessible depuis l'extérieur)
public:
    // Méthodes
    NombreComplexe();
    NombreComplexe(double a, double b);
    double module() const ;
    double partie_reelle() const;
    double partie_imaginaire() const ;
    void afficher(std::ostream& flux) const;
    NombreComplexe conjuguee() const;
    NombreComplexe addition(NombreComplexe const& T) const;
    NombreComplexe soustraction(NombreComplexe const& S) const;
    NombreComplexe lineaire(double t);
    NombreComplexe produit(NombreComplexe Q) const;
    NombreComplexe division(NombreComplexe R) const;
    NombreComplexe expz() const;
    NombreComplexe logz() const;
    NombreComplexe racine_carre_complexe() const;
    NombreComplexe puissance_complexe(int exposant) const;
    double argument() const;
    bool estEgal(NombreComplexe const& z_2) const;
    NombreComplexe& operator+=(NombreComplexe const& T);
    NombreComplexe& operator-=(NombreComplexe const& T);
    NombreComplexe& operator*=(NombreComplexe const& Q);
    NombreComplexe& operator/=(NombreComplexe const& Q);
    NombreComplexe& operator=(NombreComplexe const& NombreComplexeAcopier);
    ~NombreComplexe();



    // Tout ce qui suit est prive (inaccessible depuis l'extérieur)
private:
    // Attributs

    double m_partie_reelle;
    double m_partie_imaginaire;
};

NombreComplexe expi(double theta);
bool operator==(NombreComplexe const& z_1, NombreComplexe const& z_2);
bool operator!=(NombreComplexe const& z_1, NombreComplexe const& z_2);
std::ostream& operator<<(std::ostream& flux, NombreComplexe const& z);
NombreComplexe operator+(NombreComplexe const& a, NombreComplexe const& b);
NombreComplexe operator-(NombreComplexe const& z_1, NombreComplexe const& z_2);
NombreComplexe operator*(NombreComplexe const& z_1, NombreComplexe const& z_2);
NombreComplexe operator/(NombreComplexe const& z_1, NombreComplexe const& z_2);


#endif // FONCTIONSMATHEMATIQUES_H_INCLUDED
