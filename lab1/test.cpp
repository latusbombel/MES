#include <iostream>
#include <cmath>
#include "gauss_leg.h"

double kwadratowa(double x){
    return x*x*x;
}

int main(){
    double wezly[4] = {0. , 0., 0., 0.};
    double wagi[4] = {0. , 0., 0., 0.};


    std::cout << "całka z x^2 na przedziale [0;2] to = " << kwad_gauss_leg(4, 0, 2, kwadratowa, wezly, wagi) << std::endl;

    return 0;
}