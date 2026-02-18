#pragma once

void wiel_leg(int stopien_wielo, double tablicaWag[], double tablicaZer[]);


template <typename FunkcjaT>
double kwad_gauss_leg(int N, double x_a, double x_b, FunkcjaT funkcja, double wezly[], double wagi[]){
    double t, wynik = 0.;

    wiel_leg(N, wagi, wezly);
    
    for(int i = 0; i<N; i++){
        t = (x_a + x_b)/2. + wezly[i]*(x_b - x_a)/2.;

        wynik += (x_b - x_a)/2. * wagi[i] * funkcja(t);
    }

    return wynik;
}