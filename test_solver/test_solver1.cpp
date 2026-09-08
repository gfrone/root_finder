#include <iostream>
#include <functional>
#include <string>
#include "dual.hpp"
#include "types.hpp"
#include "solvers.hpp" 

// Função auxiliar para converter o enum StatusCode em texto (para printar bonito)
std::string status_to_string(StatusCode s) {
    switch(s) {
        case StatusCode::SUCCESS: return "SUCCESS";
        case StatusCode::DOMAIN_INVALID: return "DOMAIN_INVALID";
        case StatusCode::DIVISION_BY_ZERO: return "DIVISION_BY_ZERO";
        case StatusCode::DIVERGENCE: return "DIVERGENCE";
        case StatusCode::MAX_ITERATIONS_REACHED: return "MAX_ITERATIONS_REACHED";
        default: return "UNKNOWN";
    }
}
#define MAXITER 100
// Macro simples para validar os resultados
#define ASSERT_STATUS(test_name, expected, obtained) \
    if ((expected) == (obtained)) { \
        std::cout << "[PASS] " << test_name << "\n"; \
    } else { \
        std::cout << "[FAIL] " << test_name \
                  << "\n       Esperado: " << status_to_string(expected) \
                  << "\n       Obtido:   " << status_to_string(obtained) << "\n"; \
    }

int main() {
    double tol = 1e-6;

    // Função 1: f(x) = x^2 - 4 (Raiz exata em x = 2 e x = -2)
    auto f_parabola = [](const Dual& x) {
        return x * x - Dual(4.0, 0.0);
    };

    // Função 2: f(x) = ln(x) - 1 (Raiz exata em x = e ~= 2.718)
    auto f_log = [](const Dual& x) {
        return log(x) - Dual(1.0, 0.0);
    };

    std::cout << "--- 1. Testes de Sucesso (Caminho Feliz) ---\n";
    
    SolverResult res_newton = newton_raphsen(f_parabola, 3.0, tol, MAXITER);
    ASSERT_STATUS("Newton-Raphson (Raiz de x^2 - 4)", StatusCode::SUCCESS, res_newton.status);
    
    SolverResult res_bissecao = bissecao(f_parabola, 0.0, 5.0, tol);
    ASSERT_STATUS("Bissecao (Raiz de x^2 - 4)", StatusCode::SUCCESS, res_bissecao.status);


    std::cout << "\n--- 2. Testes de Erros e Comportamentos Problematicos ---\n";

    // Teste de Divisão por Zero no Newton-Raphson
    // Se f(x) = x^2 - 4, f'(x) = 2x. Se começarmos em x=0, a derivada é 0. O método deve falhar.
    SolverResult res_newton_zero = newton_raphsen(f_parabola, 0.0, tol, MAXITER);
    ASSERT_STATUS("Newton-Raphson: Derivada nula (Divisao por zero)", StatusCode::DIVISION_BY_ZERO, res_newton_zero.status);

    // Teste de Intervalo Inválido (Bisseção)
    // O intervalo [3, 5] não contém a raiz de x^2 - 4 (ambos f(a) e f(b) são positivos)
    SolverResult res_bissecao_inv = bissecao(f_parabola, 3.0, 5.0, tol);
    ASSERT_STATUS("Bissecao: Intervalo f(a)*f(b) >= 0", StatusCode::DOMAIN_INVALID, res_bissecao_inv.status);

    // Teste de Erro de Domínio (Secante ou Newton)
    // Se tentarmos avaliar ln(x) em x = -1, o método log() vai lançar std::domain_error
    SolverResult res_secante_dom = secante(f_log, -2.0, -1.0, tol, MAXITER);
    ASSERT_STATUS("Secante: Erro de dominio (ln de numero negativo)", StatusCode::DOMAIN_INVALID, res_secante_dom.status);

    // Teste de Divisão por Zero na Secante
    // Secante divide por (f(x1) - f(x0)). Se chutarmos -2 e 2 na parábola, f(2) = f(-2) = 0.
    // 0 - 0 = 0. O método deve pegar a divisão por zero.
    SolverResult res_secante_zero = secante(f_parabola, -2.0, 2.0, tol, MAXITER);
    ASSERT_STATUS("Secante: Divisao por zero (f(x1) == f(x0))", StatusCode::DIVISION_BY_ZERO, res_secante_zero.status);


    std::cout << "\nTestes de solvers finalizados!\n";
    return 0;
}