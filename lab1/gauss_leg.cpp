#include <iostream>
#include <cmath>
#include <stdexcept>
#include "gauss_leg.h"


void wiel_leg(int stopien_wielo, double tablicaWag[], double tablicaZer[]){
    double* J_down = new double[stopien_wielo];
    double* J_diag = new double[stopien_wielo];
    double* wagi = new double[stopien_wielo];

    for(int i = 0; i < stopien_wielo; i++){
        J_down[i] = (i+1)/sqrt(4.*(i+1)*(i+1)-1.);
    }

    wagi[0] = 1.;
    J_down[stopien_wielo - 1] = 0.0;

    for(int i = 1; i < stopien_wielo; i++){
        wagi[i] = 0.;
    }

    for(int i = 0; i < stopien_wielo+1; i++){
        J_diag[i] = 0.0;
    }

    for(int index_wart = 0; index_wart<stopien_wielo; index_wart++){
        int IT = 0, IT_MAX = 50;
        do{
            int m = index_wart;

            while(m < stopien_wielo - 1){
                double dd = fabs(J_diag[m]) + fabs(J_diag[m+1]);
                if (fabs(J_down[m]) + dd == dd) break;
                m++;
            }

            if(m == index_wart) break;
            
            if(IT >= IT_MAX){
                delete [] J_diag; delete [] J_down; delete [] wagi;
                throw std::runtime_error("Brak zbieznosci - wart. wlasne");
            }

            double g, r, shift;
            g = (J_diag[index_wart] - J_diag[index_wart+1]) / (J_down[index_wart]*2.);
            r = sqrt(g*g + 1);

            if(g>= 0){
                shift = J_diag[m] - J_diag[index_wart] + J_down[index_wart] / (g + r);  
            }
            else {
                shift = J_diag[m] - J_diag[index_wart] + J_down[index_wart] / (g - r);  
            }

            double s = 1., c = 1., p = 0.;
            double f, b;

            for(int i = m-1; i >= index_wart; i--){
                f = s*J_down[i];
                b = c*J_down[i];
                if(fabs(f) >= fabs(shift)){
                    c = shift / f;
                    r = sqrt(c*c + 1);
                    J_down[i+1] = f * r;
                    s = 1 / r;
                    c = c * s;
                }
                else{
                    s = f / shift;
                    r = sqrt(s*s + 1);
                    J_down[i+1] = shift * r;
                    c = 1 / r;
                    s = s * c;
                }
                shift = J_diag[i+1] - p;
                r = (J_diag[i] - shift) * s + 2.0 * c * b;
                p = s * r;
                J_diag[i+1] = shift + p;
                shift = c * r - b;

                double temp_wagi = wagi[i+1]; 

                wagi[i+1] = s * wagi[i] + c * temp_wagi;
                wagi[i]   = c * wagi[i] - s * temp_wagi;
                
            }

            J_diag[index_wart] = J_diag[index_wart] - p;
            J_down[index_wart] = shift;
            J_down[m] = 0.0;

            IT++;
        }while(IT < IT_MAX);
    }

    for (int i = 0; i < stopien_wielo; i++) {
    tablicaWag[i] = 2.0 * wagi[i] * wagi[i];
    tablicaZer[i] = J_diag[i];
    }

    delete [] J_diag;
    delete [] J_down;
    delete [] wagi;
}


// template <typename FunkcjaT>
// double kwad_gauss_leg(int N, double x_a, double x_b, FunkcjaT funkcja, double wezly[], double wagi[]){
//     double t, wynik = 0.;

//     wiel_leg(N, wagi, wezly);
    
//     for(int i = 0; i<N; i++){
//         t = (x_a + x_b)/2. + wezly[i]*(x_b - x_a)/2.;

//         wynik += (x_b - x_a)/2. * wagi[i] * funkcja(t);
//     }

//     return wynik;
// }
