/*
    Teste da implementação dos numeros dual, com varias funções
    verificando se está tratando possiveis erros corretamente
*/

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <string>
#include "dual.hpp"

// Função auxiliar para comparar valores de ponto flutuante
bool is_close(double a, double b, double tol = 1e-6) {
    return std::fabs(a - b) < tol;
}

// Macro simples para facilitar a escrita dos testes
#define ASSERT_DUAL(d, expected_val, expected_der, test_name) \
    if (is_close((d).val, (expected_val)) && is_close((d).der, (expected_der))) { \
        std::cout << "[PASS] " << test_name << "\n"; \
    } else { \
        std::cout << "[FAIL] " << test_name \
                  << " -> Esperado: (" << (expected_val) << ", " << (expected_der) \
                  << ") Obtido: (" << (d).val << ", " << (d).der << ")\n"; \
    }

void test_basic_math() {
    std::cout << "--- Testando Operacoes Basicas ---\n";
    Dual x(2.0, 1.0); // f(x) = x, no ponto x = 2
    Dual y(3.0, 0.0); // Constante 3

    ASSERT_DUAL(x + y, 5.0, 1.0, "Soma Dual + Dual");
    // (2 + e) - (3 + 0e) = (-1 + e)
    ASSERT_DUAL(x - y, -1.0, 1.0, "Subtracao Dual - Dual");
    ASSERT_DUAL(x * y, 6.0, 3.0, "Multiplicacao Dual * Dual");
    ASSERT_DUAL(x / Dual(2.0, 0.0), 1.0, 0.5, "Divisao Dual / Dual");
    // Operações com double (testando os bugs mencionados)
    ASSERT_DUAL(5.0 - x, 3.0, -1.0, "Subtracao Double - Dual");
    ASSERT_DUAL(10.0 / x, 5.0, -2.5, "Divisao Double / Dual");
}

void test_trigonometry() {
    std::cout << "\n--- Testando Trigonometria ---\n";
    Dual x(0.0, 1.0); // x = 0
    
    ASSERT_DUAL(sin(x), 0.0, 1.0, "Seno de 0");
    ASSERT_DUAL(cos(x), 1.0, 0.0, "Cosseno de 0");

    Dual half_pi(M_PI / 6.0, 1.0); // x = pi/6 (30 graus)
    ASSERT_DUAL(sin(half_pi), 0.5, std::cos(M_PI/6.0), "Seno de pi/6");
}

void test_exceptions() {
    std::cout << "\n--- Testando Captura de Excecoes ---\n";

    // 1. Divisão por zero
    try {
        Dual a(5.0, 1.0);
        Dual b(0.0, 0.0);
        Dual c = a / b;
        std::cout << "[FAIL] Divisao por zero nao capturada.\n";
    } catch (const std::domain_error& e) {
        std::cout << "[PASS] Divisao por zero capturada: " << e.what() << "\n";
    }

    // 2. Logaritmo de número negativo
    try {
        Dual a(-2.0, 1.0);
        Dual b = log(a);
        std::cout << "[FAIL] Logaritmo de negativo nao capturado.\n";
    } catch (const std::domain_error& e) {
        std::cout << "[PASS] Log de negativo capturado: " << e.what() << "\n";
    }

    // 3. Raiz quadrada de número negativo
    try {
        Dual a(-4.0, 1.0);
        Dual b = sqrt(a);
        std::cout << "[FAIL] Raiz quadrada de negativo nao capturada.\n";
    } catch (const std::domain_error& e) {
        std::cout << "[PASS] Raiz de negativo capturada: " << e.what() << "\n";
    }

    // 4. Arco seno fora do domínio (>= 1)
    try {
        Dual a(1.5, 1.0);
        Dual b = asin(a);
        std::cout << "[FAIL] asin fora do dominio nao capturado.\n";
    } catch (const std::domain_error& e) {
        std::cout << "[PASS] asin fora do dominio capturado: " << e.what() << "\n";
    }
}

int main() {
    test_basic_math();
    test_trigonometry();
    test_exceptions();
    
    std::cout << "\nTestes finalizados!\n";
    return 0;
}