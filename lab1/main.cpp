#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "gauss_leg.h"
#include "Eigen/Dense"

double phi(int i, double xi){
    switch (i){
        case 0:
            return xi*(xi-1)/2.;
        case 1:
            return -(xi + 1.)*(xi - 1.);
        case 2: 
            return xi*(xi + 1)/2.;
    }
    return 0.0;
}

double x_wspolrzedna(double xi, double x_up, double x_down){
    return (x_down+x_up)/2. + xi*(x_up - x_down)/2.;
}

double Jakobian(double x_up, double x_down){
    return (x_up-x_down)/2.;
}

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


template <typename FunkcjaT>
double pochodna(double x, FunkcjaT funkcja){
    double delta_x = 0.001;
    return (funkcja(x + delta_x) - funkcja(x - delta_x))/(2*delta_x);
}

template <typename FunkcjaT>
double druga_pochodna(double x, FunkcjaT funkcja){
    double delta_x = 0.001;
    return (funkcja(x + delta_x) - 2*funkcja(x) + funkcja(x - delta_x))/(delta_x*delta_x);
}

void zapisz_rozwiazanie_do_pliku(const std::string& nazwa_pliku, 
                                 int M, 
                                 double x_k[], 
                                 const Eigen::MatrixXd& wektory_wlasne) {
    
    std::ofstream plik(nazwa_pliku);
    
    // Ustalmy gęsty krok do rysowania (np. 1000 punktów w całym zakresie)
    double x_start = x_k[0];
    double x_end = x_k[2*M];
    double dx = 0.01; // Gęstość próbkowania

    // Pętla po osi X (zgodnie z instrukcją pkt 5)
    for (double x = x_start; x <= x_end; x += dx) {
        
        plik << x << " "; // Zapisujemy x

        // Dla każdego z 5 stanów własnych (kolumn w macierzy wektory_wlasne)
        for (int mu = 0; mu < 5; mu++) { 
            double u = 0.0;
            
            // Znajdź w którym elemencie 'm' znajduje się obecny 'x'
            // Element m rozciąga się od węzła 2*m do 2*m+2
            for (int m = 0; m < M; m++) {
                double xa = x_k[2 * m];     // Początek elementu
                double xb = x_k[2 * m + 2]; // Koniec elementu
                
                // Sprawdzamy czy x jest w tym elemencie
                if (x >= xa && x <= xb) {
                    double Jm = (xb - xa) / 2.0;       // Jakobian
                    double xi = (x - (xa + xb)/2.0) / Jm; // Przeliczenie na wsp. lokalną [-1, 1]
                    
                    // Sumowanie funkcji kształtu wewnątrz elementu (wzór 1 z PDF)
                    for (int i = 0; i < 3; i++) {
                        int global_idx = 2 * m + i; // Indeks globalny węzła
                        double c_p = wektory_wlasne(global_idx, mu); // Współczynnik wektora własnego
                        
                        u += c_p * phi(i, xi); // phi to twoja funkcja kształtu
                    }
                    break; // Znaleźliśmy element, nie musimy szukać dalej dla tego x
                }
            }
            plik << u << " ";
        }
        plik << "\n";
    }
    plik.close();
}

int main(){
    int M = 5;
    int N = 2*M + 1;

    double x_k[N]; // węzły pomiędzy którymi szukamy rozwiązań
    double x_max = 6., alpha = 0.4;
    //zapis do pliku
    std::ofstream plik ("energie.dat");
    std::ofstream plik2 ("wekt_wlas.dat");

    do{
        // Poprawiona generacja siatki
        for(int k = 0; k < N; k++){
            // t zmienia się liniowo od -1 do 1
            double t = (2.0 * k) / (N - 1.0) - 1.0;
            
            // Wzór z uwzględnieniem potęgi alpha i znaku
            x_k[k] = x_max * std::pow(std::abs(t), alpha) * sgn(t);
        }
        
        // Upewnij się, że skrajne węzły to dokładnie -x_max i x_max
        x_k[0] = -x_max;
        x_k[N-1] = x_max;

        double S[2*M+1][2*M+1];
        double O[2*M+1][2*M+1];
        double tablicaWag[5];
        double tablicaZer[5];

        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                S[i][j] = 0.0;
                O[i][j] = 0.0;
            }
        }

        for(int m=0;m<M;m++){
            for(int i=0;i<3;i++){
                for(int j=0;j<  3;j++){
                    int p=2*m+i;
                    int q=2*m+j;
                    double xa=x_k[2*m];
                    double xb=x_k[2*m+2];
                    double Jm=Jakobian(xb, xa);
                
                    auto funkcja_podcalkowa_S = [=](double x) {
                        return 1/(2.*Jm) * pochodna(x, [i](double t){return phi(i,t);}) * pochodna(x, [j](double t){return phi(j,t);}) + 1/2.*Jm * x_wspolrzedna(x, xb, xa) * x_wspolrzedna(x, xb, xa) * phi(i,x) * phi(j,x);
                    };
                    S[p][q]+= kwad_gauss_leg(4, -1., 1., funkcja_podcalkowa_S, tablicaWag, tablicaZer);

                    auto funkcja_podcalkowa_O = [=](double x) {
                    return Jm * phi(j, x) * phi(i, x);
                    };
                    O[p][q]+= kwad_gauss_leg(4, -1., 1., funkcja_podcalkowa_O, tablicaWag, tablicaZer);
                }
            }
        }

        Eigen::MatrixXd matS(N, N);
        Eigen::MatrixXd matO(N, N);

        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                matS(i, j) = S[i][j]; // Eigen używa nawiasów (wiersz, kolumna)
                matO(i, j) = O[i][j];
            }
        }

        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(matS, matO);

        if(solver.info() != Eigen::Success) {
            std::cerr << "Blad diagonalizacji!" << std::endl;
            return 1;
        }

        // 3. Pobieramy wyniki
        // values() zwraca wektor wartości własnych (energii) posortowany rosnąco
        Eigen::VectorXd energie = solver.eigenvalues();
        
        // eigenvectors() zwraca macierz, gdzie k-ta KOLUMNA to k-ty wektor własny
        Eigen::MatrixXd wektory = solver.eigenvectors();


        zapisz_rozwiazanie_do_pliku("wekt_wlas.dat", M, x_k, wektory);

        plik.precision(5); 
        // plik2.precision(5);
        for(int i = 0; i<N; i++){
            plik << std::fixed << std::setprecision(4) << std::setw(4) <<  energie(i) << std::endl;
            // for(int j = 0; j<N; j++){
            //     plik2 << std::fixed << std::setprecision(4) << std::setw(4) << wektory(i,j) << ' ';
            // }
            // plik2 << std::endl;
        }
        alpha += 0.05;
        plik << std::endl;
        // plik2 << std::endl;
    }while(alpha <= 2);

    return 0;
}