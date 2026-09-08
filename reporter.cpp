#include "reporter.hpp"
#include "solvers.hpp"
#include "types.hpp"
#include <algorithm>
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

    std::cout << std::defaultfloat;
}

// ============================================================
// Relatório de um método específico
// Função suporte para tratar o nome do método
// ============================================================

void report_method(
    const SolverResult& result,
    int every
) {
    std::string name = result.method_name.empty() ? "Metodo" : result.method_name;
    report_method(name, result, every);
}


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
// Aplicar um método e imprimir o relatório diretamente
// ============================================================

SolverResult report_bissecao(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol,
    int every
) {
    SolverResult result = bissecao(f, a, b, tol);
    report_method(result, every);
    return result;
}

SolverResult report_falsa_posicao(
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol,
    int max_iter,
    int every
) {
    SolverResult result = falsaPosicao(f, a, b, tol, max_iter);
    report_method(result, every);
    return result;
}

SolverResult report_newton(
    const std::function<Dual(const Dual&)>& f,
    double x0,
    double tol,
    int max_iter,
    int every
) {
    SolverResult result = newton_raphson(f, x0, tol, max_iter);
    report_method(result, every);
    return result;
}

SolverResult report_secante(
    const std::function<Dual(const Dual&)>& f,
    double x0,
    double x1,
    double tol,
    int max_iter,
    int every
) {
    SolverResult result = secante(f, x0, x1, tol, max_iter);
    report_method(result, every);
    return result;
}

SolverResult report_method(
    Method method,
    const std::function<Dual(const Dual&)>& f,
    double a,
    double b,
    double tol,
    int max_iter,
    int every
) {
    switch (method) {
        case Method::BISSECAO:
            return report_bissecao(f, a, b, tol, every);
        case Method::FALSA_POSICAO:
            return report_falsa_posicao(f, a, b, tol, max_iter, every);
        case Method::NEWTON:
            return report_newton(f, a, tol, max_iter, every);
        case Method::SECANTE:
            return report_secante(f, a, b, tol, max_iter, every);
    }
    return SolverResult();
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
        newton_raphson(f, a, tol, max_iter);

    SolverResult result_secante =
        secante(f, a, b, tol, max_iter);

    std::cout << '\n';
    std::cout << "===============================================================================================\n";
    std::cout << "                                   COMPARACAO DOS METODOS\n";
    std::cout << "===============================================================================================\n";

    std::cout
        << std::left
        << std::setw(18) << "Metodo"
        << std::setw(12) << "Iteracoes"
        << std::setw(26) << "Status"
        << std::setw(24) << "Convergencia"
        << std::setw(14) << "Ordem"
        << '\n';

    std::cout << std::string(94, '-') << '\n';

    auto print_summary = [](
        const std::string& name,
        const SolverResult& result
    ) {
        std::cout
            << std::left
            << std::setw(18) << name
            << std::setw(12) << result.iterations_used
            << std::setw(26) << result.get_status_message()
            << std::setw(24) << result.get_convergence_message();

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

    std::cout << "===============================================================================================\n";

    std::cout << '\n';
    std::cout << "Erros finais:\n";
    std::cout << std::string(94, '-') << '\n';

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
            << std::setw(18) << name
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

    std::cout << "===============================================================================================\n";
    std::cout << std::defaultfloat;
}

// ============================================================
// Benchmark: ranking dos métodos pela ordem de convergência
// ============================================================

namespace {

    struct BenchmarkEntry {
        std::string name;
        SolverResult result;
        double order;
    };

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
    std::vector<BenchmarkEntry> entries = {
        { "Bissecao",       bissecao(f, a, b, tol),                  0.0 },
        { "Falsa Posicao",  falsaPosicao(f, a, b, tol, max_iter),    0.0 },
        { "Newton-Raphson", newton_raphson(f, a, tol, max_iter),     0.0 },
        { "Secante",        secante(f, a, b, tol, max_iter),         0.0 },
    };

    for (auto& entry : entries) {
        if (entry.result.status == StatusCode::SUCCESS) {
            entry.order = estimate_convergence_order(entry.result);
        } else {
            entry.order = std::numeric_limits<double>::quiet_NaN();
        }
    }

    std::stable_sort(
        entries.begin(),
        entries.end(),
        [](const BenchmarkEntry& x, const BenchmarkEntry& y) {
            bool x_valid = is_ranked_valid(x);
            bool y_valid = is_ranked_valid(y);

            if (x_valid != y_valid) {
                return x_valid;
            }

            if (x_valid && y_valid) {
                return x.order > y.order;
            }

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
        << std::setw(18) << "Erro"
        << std::setw(14) << "Status"
        << '\n';

    std::cout << std::string(82, '-') << '\n';

    int rank = 1;

    for (const auto& entry : entries) {
        std::cout
            << std::left
            << std::setw(6)  << rank++
            << std::setw(18) << entry.name;

        if (entry.result.status == StatusCode::SUCCESS && std::isfinite(entry.order)) {
            std::cout
                << std::setw(14)
                << std::fixed << std::setprecision(4)
                << entry.order;
        } else {
            std::cout << std::setw(14) << "N/A";
        }

        std::cout << std::setw(12) << entry.result.iterations_used;

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

// ============================================================
// Bateria de testes
// ============================================================

void run_suite_comparison(const std::vector<TestCase>& test_cases) {
    std::cout << "\n############################################################\n";
    std::cout << "          BATERIA DE TESTES: COMPARACAO GERAL              \n";
    std::cout << "############################################################\n";

    int index = 1;
    for (const auto& tc : test_cases) {
        std::string title = tc.name.empty() ? ("Caso de Teste #" + std::to_string(index)) : tc.name;
        std::cout << "\n>>> " << title;
        if (!tc.expression_str.empty()) {
            std::cout << ": " << tc.expression_str;
        }
        std::cout << '\n';
        std::cout << "    Parametros: a = " << tc.a 
                  << ", b = " << tc.b 
                  << ", tol = " << tc.tol 
                  << ", max_iter = " << tc.max_iter << '\n';

        report_all(tc.f, tc.a, tc.b, tc.tol, tc.max_iter);
        index++;
    }
}

void run_suite_comparison(
    const std::vector<std::function<Dual(const Dual&)>>& funcs,
    double a,
    double b,
    double tol,
    int max_iter
) {
    std::vector<TestCase> test_cases;
    test_cases.reserve(funcs.size());
    int i = 1;
    for (const auto& f : funcs) {
        test_cases.emplace_back(f, a, b, tol, max_iter, "Funcao #" + std::to_string(i++), "");
    }
    run_suite_comparison(test_cases);
}