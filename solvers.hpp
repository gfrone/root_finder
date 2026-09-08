#pragma once

#include <math.h>
#include <iostream>
#include <vector>
#include "dual.hpp"
#include "types.hpp"
#include <functional>

/*
Método da bisseção:
    - Dado um intervalo [a, b], percorrer através do ponto médio do intervalo
    - Esta implementação, ja calcula o numero maximo de iterações com:
        n >= log(b - a) - log(tol)/ log(2)
*/
SolverResult bissecao(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol
);


SolverResult falsaPosicao(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol,
    int max_iter
);

SolverResult newton_raphsen(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double tol,
    int max_iter
);

SolverResult secante(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol,
    int max_iter
);