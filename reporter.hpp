#pragma once

#include <functional>
#include <string>
#include <vector>

#include "dual.hpp"
#include "solvers.hpp"
#include "types.hpp"

// Caso de teste para baterias de comparação
struct TestCase {
    std::function<Dual(const Dual&)> f;
    double a = 0.0;
    double b = 0.0;
    double tol = 1e-6;
    int max_iter = 100;
    std::string name = "";
    std::string expression_str = "";

    TestCase() = default;

    // Construtor de TestCase
    TestCase(
        std::function<Dual(const Dual&)> func,
        double a_val,
        double b_val,
        double tol_val = 1e-6,
        int max_iter_val = 100,
        std::string test_name = "",
        std::string expr = ""
    ) : f(func), a(a_val), b(b_val), tol(tol_val), max_iter(max_iter_val), name(test_name), expression_str(expr) {}
};

// ============================================================
// Impressão do histórico de iterações
// ============================================================
void print_history(
    const SolverResult& result,
    int every = 1
);

// ============================================================
// Relatório de um método específico a partir do SolverResult
// (O nome é deduzido automaticamente de result.method_name, sem strings)
// ============================================================
void report_method(
    const SolverResult& result,
    int every = 1
);

// Sobrecarga opcional caso queira forçar um nome customizado
void report_method(
    const std::string& method_name,
    const SolverResult& result,
    int every = 1
);

// ============================================================
// Aplicar um método e já imprimir o relatório diretamente
// ============================================================

// Bisseção
SolverResult report_bissecao(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol = 1e-6,
    int every = 1
);

// Falsa Posição
SolverResult report_falsa_posicao(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol = 1e-6,
    int max_iter = 100,
    int every = 1
);

// Newton-Raphson
SolverResult report_newton(
    const std::function<Dual(const Dual&)>& f,
    double x0,
    double tol = 1e-6,
    int max_iter = 100,
    int every = 1
);

// Secante
SolverResult report_secante(
    const std::function<Dual(const Dual&)>& f,
    double x0,
    double x1,
    double tol = 1e-6,
    int max_iter = 100,
    int every = 1
);

// Despacho via enum (sem strings)
SolverResult report_method(
    Method method,
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b = 0.0,
    double tol = 1e-6,
    int max_iter = 100,
    int every = 1
);

// ============================================================
// Relatório comparativo de todos os métodos
// ============================================================
void report_all(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol = 1e-6,
    int max_iter = 100
);

// ============================================================
// Benchmark com ranking por ordem de convergência
// ============================================================
void benchmark_metodos(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol = 1e-6,
    int max_iter = 100
);

// ============================================================
// Bateria de testes
// ============================================================
void run_suite_comparison(const std::vector<TestCase>& test_cases);

// Sobrecarga simples apenas com vetor de funções e parâmetros gerais
void run_suite_comparison(
    const std::vector<std::function<Dual(const Dual&)>>& funcs,
    double a,
    double b,
    double tol = 1e-6,
    int max_iter = 100
);