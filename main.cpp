#include <iostream>
#include <vector>
#include <functional>
#include <string>

#include "dual.hpp"
#include "types.hpp"
#include "solvers.hpp"
#include "reporter.hpp"

// ============================================================
// Definição das Funções de Teste usando Números Duais
// ============================================================

// f1(x) = x³ - x - 2 = 0  (Raiz em x ≈ 1.52138)
Dual f1(const Dual& x) {
    return x * x * x - x - 2.0;
}

// f2(x) = cos(x) - x = 0  (Raiz em x ≈ 0.739085 - Ponto fixo do cosseno)
Dual f2(const Dual& x) {
    return cos(x) - x;
}

// f3(x) = exp(x) - 3*x = 0  (Possui duas raízes: uma em ~0.619 e outra em ~1.512)
Dual f3(const Dual& x) {
    return exp(x) - 3.0 * x;
}

// f4(x) = x² - 2 = 0  (Raiz em x = sqrt(2) ≈ 1.41421356)
Dual f4(const Dual& x) {
    return x * x - 2.0;
}

Dual f5(const Dual& x) {
    return -1564000 + 1000000 * exp(x) + (435000 / x) * (exp(x) - 1);
}

// ============================================================
// Execução dos Relatórios
// ============================================================

void run_detailed_report_example() {
    std::cout << "\n############################################################\n";
    std::cout << "      RELATORIO DETALHADO POR METODO (PASSO A PASSO)       \n";
    std::cout << "############################################################\n";

    double a = 1.0;
    double b = 2.0;
    double tol = 1e-8;
    int max_iter = 100;
    int print_every = 1; // Imprime cada iteração

    std::cout << "\nAnalisando detalhadamente f(x) = x^3 - x - 2 no intervalo:\n";

    // 1. Bisseção
    SolverResult res_bis = bissecao(f1, a, b, tol);
    report_method("Bissecao", res_bis, print_every);

    // 2. Newton-Raphson (chute a = 1.0)
    SolverResult res_newton = newton_raphsen(f1, a, tol, max_iter);
    report_method("Newton-Raphson", res_newton, print_every);

    // 3. Secante (chutes x0 = 1.0, x1 = 2.0)
    SolverResult res_secante = secante(f1, a, b, tol, max_iter);
    report_method("Secante", res_secante, print_every);
}

void run_edge_cases() {
    std::cout << "\n############################################################\n";
    std::cout << "             TESTES DE BORDAS E STATUS CODES                \n";
    std::cout << "############################################################\n";

    double tol = 1e-6;

    // Caso 1: Chute inicial já é a raiz exata em f(x) = x^2 - 4 (raiz x = 2)
    auto f_exact = [](const Dual& x) { return x * x - 4.0; };
    std::cout << "\n[Caso 1] Chute inicial sendo a raiz exata (x = 2.0):\n";
    SolverResult res_exact = newton_raphsen(f_exact, 2.0, tol, 100);
    std::cout << "Newton-Raphson Status: " << res_exact.get_status_message()
              << " | Iteracoes: " << res_exact.iterations_used << '\n';

    // Caso 2: Derivada nula no chute inicial (f(x) = x^2 - 2 em x = 0)
    std::cout << "\n[Caso 2] Derivada nula no ponto inicial (f'(0) = 0 em x = 0):\n";
    SolverResult res_div_zero = newton_raphsen(f4, 0.0, tol, 100);
    std::cout << "Newton-Raphson Status: " << res_div_zero.get_status_message() 
              << " (Esperado: DIVISION_BY_ZERO)\n";

    // Caso 3: Fora de domínio (f(x) = log(x) com intervalo negativo)
    auto f_log = [](const Dual& x) { return log(x); };
    std::cout << "\n[Caso 3] Avaliacao fora de dominio (f(x) = ln(x) em [-2, -1]):\n";
    SolverResult res_domain = bissecao(f_log, -2.0, -1.0, tol);
    std::cout << "Bissecao Status: " << res_domain.get_status_message() 
              << " (Esperado: DOMAIN_INVALID)\n";
}

// ============================================================
// Ponto de Entrada Principal
// ============================================================

int main() {
    // Lista de casos de teste para o benchmark comparativo
    std::vector<TestCase> test_cases = {
        {
            "Polinomio Cubico",
            "f(x) = x^3 - x - 2",
            f1,
            1.0, 2.0, 1e-8, 100
        },
        {
            "Funcao Transcendente",
            "f(x) = cos(x) - x",
            f2,
            0.0, 1.0, 1e-8, 100
        },
        {
            "Exponencial Nao-Linear",
            "f(x) = exp(x) - 3x",
            f3,
            0.0, 1.0, 1e-8, 100
        },
        {
            "Calculo de Raiz Quadrada",
            "f(x) = x^2 - 2",
            f4,
            1.0, 2.0, 1e-12, 50
        },
        {
            "Exercicio 12",
            "f(x) = -1564000 + 1000000*exp(x) e etc ",
            f5,
            1.0, 2.0, 1e-4, 100
        }
    };

    // 1. Executa a tabela comparativa para todos os casos
    run_suite_comparison(test_cases);

    // 2. Executa a exibição passo a passo de um caso específico
    // run_detailed_report_example();

    // 3. Executa testes específicos de códigos de erro
    // run_edge_cases();

    // benchmark_metodos(f5, 0.2, 2.0, 1e-8, 100);

    std::cout << "\nTodos os testes foram executados com sucesso.\n";
    return 0;
}