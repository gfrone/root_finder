#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>
#include <string>
#include <cmath>

enum class StatusCode {
    SUCCESS = 0,
    MAX_ITERATIONS_REACHED,
    DOMAIN_INVALID,
    DIVISION_BY_ZERO,
    DIVERGENCE,
    OTHERS  // - inclui a == b, ou outros
};

enum class ConvergenceType {
    NONE,
    BY_X,
    BY_Y
};

struct IterationData {
    int it;           // Iteração atual
    double x_k;       // Valor de x na iteração atual
    double f_x;       // Valor da função em x_k
    double error_x;     // Erro estimado (ex: |x_k - x_{k-1}| ou |f(x)|)
    double error_y; 
};

enum class Method {
    BISSECAO,
    FALSA_POSICAO,
    NEWTON,
    SECANTE,
};

struct SolverResult {
    StatusCode status = StatusCode::OTHERS;
    double root = 0.0;
    int iterations_used = 0;
    std::vector<IterationData> history;
    ConvergenceType convergence = ConvergenceType::NONE;
    std::string method_name = ""; 

    SolverResult() = default;

    SolverResult(
        StatusCode stat, 
        double root_val = 0.0, 
        int iter = 0, 
        const std::vector<IterationData>& hist = {}, 
        const std::string& name = "",
        ConvergenceType conv = ConvergenceType::NONE
    ) : status(stat), 
        root(root_val), 
        iterations_used(iter), 
        history(hist), 
        convergence(conv),
        method_name(name) {}

    bool is_sucess() const {
        return status == StatusCode::SUCCESS;
    }

    std::string get_status_message() const {
        switch (status)
        {
        case StatusCode::SUCCESS: 
            return "Raiz encontrada";
        case StatusCode::MAX_ITERATIONS_REACHED: 
            return "Aviso: Número máximo de iterações atingido antes da convergência. ";
        case StatusCode::DOMAIN_INVALID: 
            return "Erro Fatal: Avaliação fora do domínio matemático. ";
        case StatusCode::DIVISION_BY_ZERO: 
            return "Erro Fatal: Divisão por zero detectada (possível derivada nula). ";
        case StatusCode::DIVERGENCE: 
            return "Erro Fatal: O método divergiu para o infinito (Inf ou NaN). ";
        default: 
            return "Erro Fatal: Ocorrência não mapeada. ";
        }
    }

    std::string get_convergence_message() const {
        switch (convergence)
        {
        case ConvergenceType::BY_X:
            return "Convergência por x. ";

        case ConvergenceType::BY_Y:
            return "Convergência por y. ";

        default:
            return "Convergência não atingida. ";
        }
    }
};


#endif