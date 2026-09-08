#include <iostream>
#include <functional>
#include <string>

#include "dual.hpp"
#include "types.hpp"
#include "solvers.hpp" 
#include "reporter.hpp"

void exercicio_7(int a, int b, int tol) {
    if(a == b || a > b) {
        printf("Intervalos inválidos para Bisseção");
        return;
    }
    
    // Utilizazr método da bisseção para encontrar solução com precisão de 10e-3 para as funções:
    auto f1 = [](const Dual& x) {
        return 4 * cos(x) - exp(2*x);
    };

    // f2(x) = cos(x) - x = 0  (Raiz em x ≈ 0.739085 - Ponto fixo do cosseno)
    auto f2 = [](const Dual& x) {
        return x/2 - tan(x);
    };
    
    // f3(x) = exp(x) - 3*x = 0  (Possui duas raízes: uma em ~0.619 e outra em ~1.512)
    auto f3 = [](const Dual& x) {
        return 1 - x * log(x);
    };

    // f4(x) = x² - 2 = 0  (Raiz em x = sqrt(2) ≈ 1.41421356)
    auto f4 = [](const Dual& x) {
        return x*x*x + x - 10;
    };

    report_bissecao(f1, a, b, tol);
    report_bissecao(f2, a, b, tol);
    report_bissecao(f3, a, b, tol);
    report_bissecao(f4, a, b, tol);


}


int main() {
    int a, b;
    std::cin >> a >> b;
    // exercicio 7 tolerancia é de 10-4
    exercicio_7(a, b, 1e-4);
}