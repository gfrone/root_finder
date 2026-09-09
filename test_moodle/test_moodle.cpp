#include <iostream>
#include <functional>
#include <string>

#include "dual.hpp"
#include "types.hpp"
#include "solvers.hpp" 
#include "reporter.hpp"


void exericicio_7(double tol){
    auto f1 = [](const Dual& x){
        return x - pow(2, x);
    };
    
    auto f2 = [](const Dual& x){
        return exp(x) - x*x + 3*x - 2;
    };
    
    auto f3 = [](const Dual& x){
        return 2 * x * cos(2 * x) - ((x + 1)*(x + 1));
    };

    auto f4 = [](const Dual& x){
        return x * cos(x) - 2*x*x + 3*x - 1;
    };

    report_bissecao(f1, 0, 1, tol, 1);
    report_bissecao(f2, 0, 1, tol, 1);
    report_bissecao(f3, -3, -2, tol, 1);
    report_bissecao(f3, -1, 0, tol, 1);
    report_bissecao(f4, 0.2, 0.3, tol, 1);
    report_bissecao(f4, 1.2 , 1.3, tol, 1);
}

void exercicio_8(double tol){
    auto f1 = [](const Dual& x) {
        return x - pow(25, 1/3);
    };
    // raiz de 25, pra uma aproximação grosseira, está entre 2(8) - 3(27)
    // para o ponto < 0, basta qualuqer negativo, como -1
    report_bissecao(f1, -1, 2, tol, 1);
}

void exercicio_9(double tol) {
    auto f1 = [](const Dual& x) {
        return exp(x) + 2 - x + 2 * cos(x) - 6;
    };

    auto f2 = [](const Dual& x) {
        return log(x - 1) + cos(x - 1);
    };

    auto f3 = [](const Dual& x) {
        return 2 * x * cos(2 * x) - (x - 2) * (x - 2);
    };

    auto f4 = [](const Dual& x) {
        return (x - 2) * (x - 2) - log(x);
    };

    auto f5 = [](const Dual& x) {
        return exp(x) - 3 * x * x;
    };

    report_falsa_posicao(f1, 1, 2, tol);
    report_falsa_posicao(f2, 1.3, 2, tol);
    report_falsa_posicao(f3, 2, 3, tol);
    report_falsa_posicao(f3, 3, 4, tol);
    report_falsa_posicao(f4, 1, 2, tol);
    report_falsa_posicao(f4, 2, 4, tol);
    report_falsa_posicao(f5, 0, 1, tol);
    report_falsa_posicao(f5, 3, 5, tol);
}

void exercicio_10(double tol) {
    auto f1 = [](const Dual& x) {
        return exp(x) + 2 - x + 2 * cos(x) - 6;
    };

    auto f2 = [](const Dual& x) {
        return log(x - 1) + cos(x - 1);
    };

    auto f3 = [](const Dual& x) {
        return 2 * x * cos(2 * x) - (x - 2) * (x - 2);
    };

    auto f4 = [](const Dual& x) {
        return (x - 2) * (x - 2) - log(x);
    };

    auto f5 = [](const Dual& x) {
        return exp(x) - 3 * x * x;
    };

    // Newton usa o ponto médio de cada intervalo como chute inicial.
    report_newton(f1, 1.5, tol);
    report_newton(f2, 1.65, tol);
    report_newton(f3, 2.5, tol);
    report_newton(f3, 3.5, tol);
    report_newton(f4, 1.5, tol);
    report_newton(f4, 3.0, tol);
    report_newton(f5, 0.5, tol);
    report_newton(f5, 4.0, tol);

    // Secantes usa os extremos de cada intervalo como aproximações iniciais.
    report_secante(f1, 1, 2, tol);
    report_secante(f2, 1.3, 2, tol);
    report_secante(f3, 2, 3, tol);
    report_secante(f3, 3, 4, tol);
    report_secante(f4, 1, 2, tol);
    report_secante(f4, 2, 4, tol);
    report_secante(f5, 0, 1, tol);
    report_secante(f5, 3, 5, tol);
}

void exercicio_11(double tol) {
    auto f1 = [](const Dual& x) {
        return x * x * x - 9 * x * x + 12;
    };
    auto df1 = [](const Dual& x) {
        return 3 * x * x - 18 * x;
    };

    auto f2 = [](const Dual& x) {
        return x * x * x * x - 2 * x * x * x - 5 * x * x + 12 * x - 5;
    };
    auto df2 = [](const Dual& x) {
        return 4 * x * x * x - 6 * x * x - 10 * x + 12;
    };

    auto f3 = [](const Dual& x) {
        return x * x / 2 + x * (log(x) - 1);
    };
    auto df3 = [](const Dual& x) {
        return x + log(x);
    };

    std::cout << "\nExercicio 11(a) - zeros:\n";
    report_newton(f1, -1.0, tol);
    report_newton(f1, 1.0, tol);
    report_newton(f1, 9.0, tol);
    std::cout << "\nExercicio 11(a) - pontos criticos:\n";
    report_newton(df1, 1.0, tol);
    report_newton(df1, 7.0, tol);

    std::cout << "\nExercicio 11(b) - zeros:\n";
    report_newton(f2, -2.0, tol);
    report_newton(f2, 0.5, tol);
    report_newton(f2, 1.5, tol);
    report_newton(f2, 2.5, tol);
    std::cout << "\nExercicio 11(b) - pontos criticos:\n";
    report_newton(df2, -2.0, tol);
    report_newton(df2, 0.5, tol);
    report_newton(df2, 2.5, tol);

    // O domínio de f3 é x > 0; há um zero e um ponto crítico nesse domínio.
    std::cout << "\nExercicio 11(c) - zero:\n";
    report_newton(f3, 1.5, tol);
    std::cout << "\nExercicio 11(c) - ponto critico:\n";
    report_newton(df3, 0.5, tol);
}

void exercicio_12(double tol) {
    auto populacao = [](const Dual& lambda) {
        return 1000000 * exp(lambda)
            + (435000 / lambda) * (exp(lambda) - 1)
            - 1564000;
    };

    std::cout << "\nExercicio 12 - lambda:\n";
    report_newton(populacao, 0.1, tol);
}

int main() {
    // exericicio_7(1e-3);
    // exercicio_8(1e-8);
    // exercicio_9(1e-3);
    // exercicio_10(1e-3);
    exercicio_11(1e-4);
    exercicio_12(1e-4);
    return 0;
}
