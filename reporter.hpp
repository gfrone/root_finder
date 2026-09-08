#ifndef REPORTER_HPP
#define REPORTER_HPP

#include <functional>
#include <string>

#include "dual.hpp"
#include "types.hpp"
#include "solvers.hpp"


// ============================================================
// Estrutura para descrever um Caso de Teste
// Sera usado para facilitar os testes
// ============================================================
struct TestCase {
    std::string name;
    std::string expression_str;
    std::function<Dual(const Dual&)> f;
    double a;        // Limite inferior (ou chute inicial / x0)
    double b;        // Limite superior (ou chute x1)
    double tol;      // Tolerância desejada
    int max_iter;    // Limite máximo de iterações
};


// ============================================================
// Utilitários
// ============================================================

std::string status_to_string(StatusCode status);

std::string convergence_to_string(ConvergenceType convergence);

void print_history(
    const SolverResult& result,
    int every
);

// ============================================================
// Estimativa da ordem de convergência
// ============================================================

double estimate_convergence_order(const SolverResult& result);


// ============================================================
// Relatório de um método específico
// ============================================================

void report_method(
    const std::string& method_name,
    const SolverResult& result,
    int every = 1
);


// ============================================================
// Relatório comparativo de todos os métodos
//
// Os parâmetros a e b são usados como:
//   - Bisseção:       [a, b]
//   - Falsa Posição:  [a, b]
//   - Newton:         a como chute inicial
//   - Secante:        a e b como chutes iniciais
// ============================================================

void report_all(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol = 1e-7,
    int max_iter = 100
);

void benchmark_metodos(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol,
    int max_iter
);

void run_suite_comparison(const std::vector<TestCase>& test_cases);

#endif
