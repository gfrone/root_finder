#pragma once

#include <math.h>
#include <iostream>
#include <vector>
#include "dual.hpp"
#include "types.hpp"
#include <functional>

// Declaração dos solvers com valores padrão de tolerância e iterações
SolverResult bissecao(const std::function<Dual(const Dual&)>& f, double a, double b, double tol = 1e-6);
SolverResult falsaPosicao(const std::function<Dual(const Dual&)>& f, double a, double b, double tol = 1e-6, int max_iter = 100);
SolverResult newton_raphson(const std::function<Dual(const Dual&)>& f, double a, double tol = 1e-6, int max_iter = 100);
SolverResult secante(const std::function<Dual(const Dual&)>& f, double x0, double x1, double tol = 1e-6, int max_iter = 100);

// Declaração do estimador de ordem de convergência (usado pelo reporter)
double estimate_convergence_order(const SolverResult& result);

/*
    Tratamento de estagnação: Para métodos, exceto a bisseção, podem
    sofrer com este tipo de comportamento. Para isso, façamos uma comparaçao entre
    o x anterior e o proximo ponto. O 8.0 dá uma margem de aproximadamente oito ULPs (units in the last place): pequenas diferenças causadas 
    por arredondamentos ainda são consideradas “sem progresso”.

    Definimos uma escala numerica para comparar:
        - Longe de zero, usa o maior módulo entre o x atual e o do proximo passo
        - Proximo do zero, usa 1 evitando uma tolerancia nula
        - Retorna que são praticamente iguais quando tem diferença minima de 8 × epsilon × escala
    
*/