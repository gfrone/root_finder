#include "reporter.hpp"
#include "solvers.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>


// ============================================================
// Impressão do histórico
// ============================================================

void print_history(
    const SolverResult& result,
    int every
) {
    if (result.history.empty()) {
        std::cout << "Nenhuma iteracao registrada.\n";
        return;
    }
    
    if (every <= 0) {
        every = 1;
    }

    std::cout
        << std::left
        << std::setw(10) << "Iteracao"
        << std::setw(20) << "x"
        << std::setw(20) << "f(x)"
        << std::setw(20) << "Erro X"
        << std::setw(20) << "Erro Y"
        << '\n';

    std::cout << std::string(90, '-') << '\n';

    const std::size_t last = result.history.size() - 1;

    for (std::size_t i = 0; i < result.history.size(); ++i) {
        const auto& data = result.history[i];

        // mostra sempre a primeira e a final, depois a cada k iterações
        bool show =
            (i == 0) ||
            (data.it % every == 0) ||
            (i == last);

        if (!show) {
            continue;
        }

        std::cout
            << std::left
            << std::setw(10) << data.it
            << std::scientific << std::setprecision(10)
            << std::setw(20) << data.x_k
            << std::setw(20) << data.f_x
            << std::setw(20) << data.error_x
            << std::setw(20) << data.error_y
            << '\n';
    }

    std::cout << std::defaultfloat; // Restaura o formato padrão de ponto flutuante
}


// ============================================================
// Relatório de um método específico usado por funções 
// ============================================================
void report_method(
    const std::string& method_name,
    const SolverResult& result,
    int every
) {
    std::cout << '\n';
    std::cout << "============================================================\n";
    std::cout << "                    " << method_name << '\n';
    std::cout << "============================================================\n";

    std::cout
        << std::fixed
        << std::setprecision(10);

    std::cout << "Status:              "
              << result.get_status_message()
              << '\n';

    std::cout << "Iteracoes usadas:    "
              << result.iterations_used
              << '\n';

    std::cout << "Convergencia:        "
              << result.get_convergence_message()
              << '\n';

    std::cout << "Raiz aproximada:     "
              << result.root
              << '\n';

    if (result.status == StatusCode::SUCCESS) {
        double order = estimate_convergence_order(result);

        if (std::isfinite(order)) {
            std::cout << "Ordem estimada:      "
                      << std::fixed << std::setprecision(6)
                      << order
                      << '\n';

            if (!result.history.empty()) {
                const auto& last = result.history.back();

                std::cout << std::scientific << std::setprecision(10);
                std::cout << "Erro final em X:     "
                          << last.error_x
                          << '\n';
                std::cout << "Erro final em Y:     "
                          << last.error_y
                          << '\n';
            }
        } else {
            std::cout << "Ordem estimada:      N/A\n";
        }
    }

    std::cout << '\n';
    std::cout << "Historico:\n";
    std::cout << std::string(90, '-') << '\n';

    print_history(result, every);

    std::cout << "============================================================\n";
    std::cout << std::defaultfloat;
}


// ============================================================
// Relatório resumido de todos os métodos
// ============================================================

void report_all(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol,
    int max_iter
) {
    SolverResult result_bissecao =
        bissecao(f, a, b, tol);

    SolverResult result_falsa_posicao =
        falsaPosicao(f, a, b, tol, max_iter);

    SolverResult result_newton =
        newton_raphsen(f, a, tol, max_iter);

    SolverResult result_secante =
        secante(f, a, b, tol, max_iter);

    std::cout << '\n';
    std::cout << "========================================================================================\n";
    std::cout << "                                COMPARACAO DOS METODOS\n";
    std::cout << "========================================================================================\n";

    std::cout
        << std::left
        << std::setw(20) << "Metodo"
        << std::setw(12) << "Iteracoes"
        << std::setw(26) << "Status"
        << std::setw(18) << "Convergencia"
        << std::setw(14) << "Ordem"
        << '\n';

    std::cout << std::string(88, '-') << '\n';

    auto print_summary = [](
        const std::string& name,
        const SolverResult& result
    ) {
        std::cout
            << std::left
            << std::setw(20) << name
            << std::setw(12) << result.iterations_used
            << std::setw(26) << result.get_status_message()
            << std::setw(18) << result.get_convergence_message();

        if (result.status == StatusCode::SUCCESS) {
            double order = estimate_convergence_order(result);

            if (std::isfinite(order)) {
                std::cout
                    << std::setw(14)
                    << std::fixed
                    << std::setprecision(4)
                    << order;
            } else {
                std::cout
                    << std::setw(14)
                    << "N/A";
            }
        } else {
            std::cout
                << std::setw(14)
                << "N/A";
        }

        std::cout << '\n';
    };

    print_summary("Bissecao", result_bissecao);
    print_summary("Falsa Posicao", result_falsa_posicao);
    print_summary("Newton-Raphson", result_newton);
    print_summary("Secante", result_secante);

    std::cout << "========================================================================================\n";

    std::cout << '\n';
    std::cout << "Erros finais:\n";
    std::cout << std::string(88, '-') << '\n';

    auto print_errors = [](
        const std::string& name,
        const SolverResult& result
    ) {
        if (result.status != StatusCode::SUCCESS || result.history.empty()) {
            return;
        }

        const auto& last = result.history.back();

        std::cout
            << std::left
            << std::setw(20) << name
            << "Erro X = "
            << std::scientific
            << std::setprecision(6)
            << last.error_x
            << "    Erro Y = "
            << last.error_y
            << '\n';
    };

    print_errors("Bissecao", result_bissecao);
    print_errors("Falsa Posicao", result_falsa_posicao);
    print_errors("Newton-Raphson", result_newton);
    print_errors("Secante", result_secante);

    std::cout << "========================================================================================\n";
    std::cout << std::defaultfloat;
}

 
// ============================================================
// Benchmark: ranking dos metodos pela ordem de convergencia
// ============================================================
 
namespace {
 
    // Entrada auxiliar do ranking: nome do metodo, resultado bruto
    // e a ordem de convergencia estimada (NaN se nao aplicavel).
    struct BenchmarkEntry {
        std::string name;
        SolverResult result;
        double order;
    };
 
    // Um método só é "valido" para fins de ranking se convergiu E se
    // foi possível estimar uma ordem de convergência finita a partir
    // do seu histórico (precisa de pelo menos 3 erros úteis).
    bool is_ranked_valid(const BenchmarkEntry& entry) {
        return entry.result.status == StatusCode::SUCCESS &&
               std::isfinite(entry.order);
    }
 
} // namespace
 
void benchmark_metodos(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol,
    int max_iter
) {
 
    /*
        Executa os quatro métodos, exatamente como em report_all.
 
        Bisseção e Falsa Posição usam o intervalo [a, b];
        Newton usa a como chute inicial;
        Secante usa a e b como os dois chutes iniciais.
    */
 
    std::vector<BenchmarkEntry> entries = {
        { "Bissecao",       bissecao(f, a, b, tol),                  0.0 },
        { "Falsa Posicao",  falsaPosicao(f, a, b, tol, max_iter),    0.0 },
        { "Newton-Raphson", newton_raphsen(f, a, tol, max_iter),     0.0 },
        { "Secante",        secante(f, a, b, tol, max_iter),         0.0 },
    };
 
    // Estima a ordem de convergência apenas para quem convergiu de fato.
    for (auto& entry : entries) {
        if (entry.result.status == StatusCode::SUCCESS) {
            entry.order = estimate_convergence_order(entry.result);
        } else {
            entry.order = std::numeric_limits<double>::quiet_NaN();
        }
    }
 
    /*
        Ranking:
          1) Métodos que convergiram e têm ordem estimável, ordenados
             da maior para a menor ordem (convergência mais rápida primeiro).
          2) Métodos que convergiram mas não têm ordem estimável
             (histórico curto demais, denominador ~0, etc).
          3) Métodos que não convergiram (MAX_ITERATIONS_REACHED,
             DOMAIN_INVALID, DIVISION_BY_ZERO, DIVERGENCE...).
 
        std::stable_sort preserva a ordem original (Bissecao, Falsa Posicao,
        Newton, Secante) como desempate dentro de cada grupo.
    */
 
    std::stable_sort(
        entries.begin(),
        entries.end(),
        [](const BenchmarkEntry& x, const BenchmarkEntry& y) {
 
            bool x_valid = is_ranked_valid(x);
            bool y_valid = is_ranked_valid(y);
 
            if (x_valid != y_valid) {
                return x_valid; // válidos vêm antes dos inválidos
            }
 
            if (x_valid && y_valid) {
                return x.order > y.order; // maior ordem primeiro
            }
 
            // Nenhum dos dois é "válido" para ranking por ordem: entre
            // não convergidos, prioriza quem pelo menos convergiu
            // (só não deu pra estimar ordem) antes de quem falhou de vez.
            bool x_convergiu = x.result.status == StatusCode::SUCCESS;
            bool y_convergiu = y.result.status == StatusCode::SUCCESS;
 
            return x_convergiu && !y_convergiu;
        }
    );
 
    std::cout << '\n';
    std::cout << "============================================================\n";
    std::cout << "     BENCHMARK: RANKING POR ORDEM DE CONVERGENCIA\n";
    std::cout << "============================================================\n";
 
    std::cout
        << std::left
        << std::setw(6)  << "Rank"
        << std::setw(18) << "Metodo"
        << std::setw(14) << "Ordem"
        << std::setw(12) << "Iteracoes"
        << std::setw(18) << "Erro" // - Pode ser x ou em y
        << std::setw(14) << "Status"
        << '\n';
 
    std::cout << std::string(82, '-') << '\n';
 
    int rank = 1;
 
    for (const auto& entry : entries) {
 
        std::cout
            << std::left
            << std::setw(6)  << rank++
            << std::setw(18) << entry.name;
 
        // Ordem estimada (ou N/A se o método não convergiu ou a ordem
        // não pôde ser calculada a partir do histórico).
        if (entry.result.status == StatusCode::SUCCESS && std::isfinite(entry.order)) {
            std::cout
                << std::setw(14)
                << std::fixed << std::setprecision(4)
                << entry.order;
        } else {
            std::cout << std::setw(14) << "N/A";
        }
 
        std::cout << std::setw(12) << entry.result.iterations_used;
 
        // abs(erro_y ou y) final: só existe se houve pelo menos uma iteração
        // registrada no histórico.
        if (!entry.result.history.empty()) {
            double error = (entry.result.convergence == ConvergenceType::BY_X)
            ? std::fabs(entry.result.history.back().error_x)
            : std::fabs(entry.result.history.back().error_y);

            std::cout
                << std::setw(18)
                << std::scientific << std::setprecision(6)
                << error;
        } else {
            std::cout << std::setw(18) << "N/A";
        }
        
        std::string status_s = entry.result.get_status_message();
        std::cout
            << std::setw(14) << status_s
            << '\n';
    }
 
    std::cout << "============================================================\n";
 
    // Um aviso explícito se nenhum método convergiu, pra deixar claro
    // que o "ranking" acima não reflete velocidade de convergência real.
    bool algum_convergiu = std::any_of(
        entries.begin(),
        entries.end(),
        [](const BenchmarkEntry& e) {
            return e.result.status == StatusCode::SUCCESS;
        }
    );
 
    if (!algum_convergiu) {
        std::cout << "Aviso: nenhum metodo convergiu para esta funcao/intervalo/chute.\n";
        std::cout << "O ranking acima reflete apenas a ordem de execucao, nao velocidade de convergencia.\n";
        std::cout << "============================================================\n";
    }
 
    std::cout << std::defaultfloat;
}

void run_suite_comparison(const std::vector<TestCase>& test_cases) {
    std::cout << "\n############################################################\n";
    std::cout << "          BATERIA DE TESTES: COMPARACAO GERAL              \n";
    std::cout << "############################################################\n";

    for (const auto& tc : test_cases) {
        std::cout << "\n>>> " << tc.name << ": " << tc.expression_str << '\n';
        std::cout << "    Parametros: a = " << tc.a 
                  << ", b = " << tc.b 
                  << ", tol = " << tc.tol 
                  << ", max_iter = " << tc.max_iter << '\n';

        report_all(tc.f, tc.a, tc.b, tc.tol, tc.max_iter);
    }
}


