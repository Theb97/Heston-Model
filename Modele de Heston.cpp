// Modele de Heston.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//
#include <iostream>
#include "FonctionsMathematiques/FonctionsMathematiques.h"
#include "Transformee_de_Fourier_rapide/Transformee_de_Fourier_rapide.h"
#include "Classe_du_modele_de_Heston/classe_du_modele_de_Heston.h"
#include <fstream>
#include <vector>
#include <string>
#include <float.h>

using namespace std;

const double PI = 3.1415926535;




int main()
{
    // Différents tests pour s'approprier le langage C++ ...      
    NombreComplexe z_1(1,1), z_2(1,1); // Création de deux nombres complexes : z_1 et z_2
    double c(z_1.module());
    modele_de_Heston xxx;
    cout << xxx.D_2_tau_phi(DBL_EPSILON) << endl;
    cout << DBL_EPSILON << endl;
    cout << (z_1 * z_2) / z_1 << endl;
    cout << z_1.logz()  << endl;
    vector<double> signal(5);
    int j(0);
    for (j = 0; j < 5; j++)
    {
        double q = static_cast<double>(j);
        signal[j] = q;
    };
    vector<NombreComplexe> y = (signal * signal);
    for (int k = 0; k < 8; k++)
    {
        cout << y[k] << endl;
    }
    
    int aa(2);
    int b(2);
    cout << "Valeur de a : " << aa << endl;
    cout << "Valeur de b : " << b << endl;
    b = ajouteDeux(aa);                    // Appel de la fonction
    afficherMessage("Je m'appelle Theophile Bourrat");
    cout << "Valeur de a : " << aa << endl;
    cout << "Valeur de b : " << b << endl;
    string const MonFichier("C++.py");
    ofstream MonFlux(MonFichier.c_str(), ios::out | ios::trunc);
    //Déclaration d'un flux permettant d'écrire dans le fichier python

    int n = 8192; // n doit être une puissance de 2 
    vector<NombreComplexe> gamma_k(n);

    vector<NombreComplexe>::iterator ut;
    vector<NombreComplexe>::iterator utt;
    modele_de_Heston z(1.0, 0.0027, 0.2, 0.3, 0.12, 0.0, 0.0, 0.0027, 0.0, 1.0);
    double a =8;
    double eta = a / static_cast<double>(n);
    double bb = 500;
    double lambda = bb / static_cast<double>(n);
    
   

    
    

    for (ut = gamma_k.begin(); ut != gamma_k.end(); ++ut)
    {
        double index = static_cast<double>(distance(gamma_k.begin(), ut));
        NombreComplexe i(0, 1);
        *ut = z.gamma(-a/2+ eta*index,1.1);
    };
    double nn = static_cast<double>(n);

   
    
    Transformee_de_Fourier_rapide gamma_k_bis(gamma_k);
    vector<NombreComplexe> gamma_bis = ((gamma_k_bis.transformee_de_Fourier_continue_methode_via_fft_fractionnaire(a,bb)));
    
    if (MonFlux)  //On teste si tout est OK
    {
        MonFlux << "import numpy as np" << endl;
        MonFlux << "import matplotlib.pyplot as plt" << endl;
        MonFlux << "   " << endl;
        MonFlux << "pi=3.1415926535" << endl;
        MonFlux << "n = " << n << endl;
        MonFlux << "b = " <<bb/2<< endl;
        MonFlux << "x = np.linspace(-b,b,n)" << endl;
        MonFlux << "##included sinh" << endl;
        MonFlux << "Y1 = [";

        vector<NombreComplexe>::iterator it;

        for (it = gamma_bis.begin(); it != gamma_bis.end(); ++it)
        {
            double index = static_cast<double>(distance(gamma_bis.begin(), it));
            MonFlux <<(*it).partie_reelle()/(2*PI);
            if (it != gamma_bis.end()-1)
            {
                MonFlux << "," ;
            }
            else
            {
                MonFlux << "]" << endl;
            }
        };
        MonFlux << "   " << endl;
        vector<NombreComplexe>::iterator utt;
        

        MonFlux << "##without sinh" << endl;
        MonFlux << "Y2 = [";
        vector<NombreComplexe>::iterator itt;
        for (itt = gamma_bis.begin() ; itt != gamma_bis.end(); ++itt)
        {
            double index = static_cast<double>(distance(gamma_bis.begin(), itt));
            if (index != n/2)
            {
                MonFlux << (*itt).partie_reelle() / (2 * PI * sinh(1.1 * (-bb / 2 + lambda * index)));
            }
            else
            {
                MonFlux << 0;
            }
            
            if (itt != gamma_bis.end() - 1)
            {
                MonFlux << ",";
            }
            else
            {
                MonFlux << "]"<< endl;
            }
            
        };
        MonFlux << "plt.plot(x,Y1,color='red')" << endl;
        MonFlux << "plt.xlabel('log strike transform variate')" << endl;
        MonFlux << "plt.ylabel('values including sinh')" << endl;
        MonFlux << "plt.show()" << endl;
        MonFlux << "plt.plot(x,Y2,color='green')" << endl;
        MonFlux << "plt.xlabel('log strike transform variate')" << endl;
        MonFlux << "plt.ylabel('values without sinh')" << endl;
        MonFlux << "plt.show()" << endl;

    }
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
    }
    gamma_k.clear();
    gamma_bis.clear();
    
    const int taille_tableau(100);
    char tableau[taille_tableau];
    tableau[0] = 65;
    tableau[1] = 66;
    for (int i = 2; i < taille_tableau; ++i)
    {
        tableau[i] = 65 + i;
    }
    cout << tableau << endl;
    return 0;
}



// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.
