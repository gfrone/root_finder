#include <math.h>
#include <iostream>
#include <vector>
#include "types.hpp"
#include "solvers.hpp"

bool safe_eval(
    const std::function<Dual(const Dual&)>& f,
    double x,
    Dual& result,
    double der = 0.0
) {
    try {
        result = f(Dual(x, der));
        return true;
    }
    catch (const std::domain_error&) {
        return false;
    }
}

SolverResult bissecao(const std::function<Dual(const Dual&)>& f, double a, double b, double tol) {
    SolverResult result;
    result.method_name = "Bissecao";
    result.iterations_used=0;
    double fa, fb;

    if(a == b) {
        result.root = -1;
        result.convergence = ConvergenceType::NONE;
        result.status = StatusCode::OTHERS;
        return result;
    }

    try{
        fa = f(Dual(a)).val;
        fb = f(Dual(b)).val;
    } catch (const std::domain_error& e){
        result.status = StatusCode::DOMAIN_INVALID;
        result.root = a;
        result.convergence = ConvergenceType::NONE;
        return result;
    }

    if(fa == 0 || fb == 0) {
        double root;
        root = (fa == 0) ? a : b;
        result.root = root;
        result.convergence = ConvergenceType::BY_Y;
        result.status = StatusCode::SUCCESS;
        return result;
    }
    
    if (fa * fb > 0) {
        result.status = StatusCode::DOMAIN_INVALID;
        result.root = a;
        result.convergence = ConvergenceType::NONE;
        return result;
    } 

    double intervalo = std::fabs(b - a);

    // 4. Cálculo seguro do número máximo de iterações
    const int max_iter = 100;
    int num_iters = max_iter;

    if(intervalo > tol) {
        double n = std::log(intervalo / tol) / std::log(2.0);
        if(std::isfinite(n) && n > 0) {
            num_iters = (n < max_iter) ? n : max_iter;
        }
    }
    else {
        // se intervalo ja é maior que  tol, basta somente avaliar no ponto medio
        num_iters = 1;
    }
    

    double meio = (a + b) / 2;
    double fm;

    for(int i = 0; i < num_iters; i++) {
        meio = (a + b) / 2; 

        try{
            fm = f(Dual(meio)).val;
        } catch (const std::domain_error& e){
            result.status = StatusCode::DOMAIN_INVALID;
            result.root = meio;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }

        // Verifica se a avaliação resultou em um valor infinito ou NaN
        if (std::isnan(fm) || std::isinf(fm)) {
            result.status = StatusCode::DIVERGENCE;
            result.root = meio;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }

        // Erro em função do eixo x e y (|f(meio)|)
        double error_x = (b - a) / 2;
        double error_y = std::fabs(fm);

        result.history.push_back({i + 1, meio, fm, error_x, error_y});

        // Verifica qual dos dois critérios converge primeiro. Como os dois são
        // checados na mesma iteração, priorizo error_y (raiz "exata" em y) sobre
        // error_x — ajuste a ordem se preferir o contrário.
        if (error_y < tol) {
            result.status = StatusCode::SUCCESS;
            result.root = meio;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::BY_Y;
            return result;
        }
        if (error_x < tol) {
            result.status = StatusCode::SUCCESS;
            result.root = meio;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::BY_X;
            return result;
        }

        if (fa * fm < 0.0) {
            b = meio;
            fb = fm;
        } else {
            a = meio;
            fa = fm;
        }
    }
    
    result.root = meio;
    result.iterations_used = num_iters;
    
    if(!result.history.empty() && result.history.back().error_x < tol) {
        // Se chegou aqui, o número de iterações calculado por num_iters garante,
        // por construção, erro em x < tol — então a convergência é por x.
        result.status = StatusCode::SUCCESS;
        result.convergence = ConvergenceType::BY_X;
    }
    
    else {
        result.status = StatusCode::MAX_ITERATIONS_REACHED;
        result.convergence = ConvergenceType::NONE;
    }

    return result;
}

SolverResult falsaPosicao(const std::function<Dual(const Dual&)>& f, double a, double b, double tol, int max_iter) {
    SolverResult result;
    result.method_name= "Falsa Posicao";
    result.iterations_used=0;
    double fa, fb;

    if(a == b) {
        result.root = -1;
        result.convergence = ConvergenceType::NONE;
        result.status = StatusCode::OTHERS;
        return result;
    }

    try{
        fa = f(Dual(a)).val;
        fb = f(Dual(b)).val;
    } catch (const std::domain_error& e){
        result.status = StatusCode::DOMAIN_INVALID;
        result.root = a;
        result.convergence = ConvergenceType::NONE;
        return result;
    }
    
    if(fa == 0 || fb == 0) {
        double root;
        root = (fa == 0) ? a : b;
        result.root = root;
        result.convergence = ConvergenceType::BY_Y;
        result.status = StatusCode::SUCCESS;
        return result;
    }

    if (fa * fb > 0) {
        result.status = StatusCode::DOMAIN_INVALID;
        result.root = a;
        result.convergence = ConvergenceType::NONE;
        return result;
    } 

    bool has_prev = false;
    double prev_x = 0.0;
    double x = a;
    double error_x = 0.0;
    double error_y = 0.0;

    for(int i = 0; i < max_iter; i++) {
        double denom = fb - fa;

        if(denom == 0) {
            result.status = StatusCode::DIVISION_BY_ZERO;
            result.root = has_prev ? prev_x : a;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }
        
        x = a - fa * (b - a) / denom;
        double fx;
        try {
            fx = f(Dual(x)).val;
        } catch(const std::domain_error& e) {
            result.status = StatusCode::DOMAIN_INVALID;
            result.root = x;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        } 

        if (std::isnan(fx) || std::isinf(fx)) {
            result.status = StatusCode::DIVERGENCE;
            result.root = x;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }

        error_x = has_prev ? std::fabs(x - prev_x) : std::fabs(b - a);
        error_y = std::fabs(fx);

        result.history.push_back({i + 1, x, fx, error_x, error_y});

        if (error_y < tol) {
            result.status = StatusCode::SUCCESS;
            result.root = x;
            result.convergence = ConvergenceType::BY_Y;
            result.iterations_used = i + 1;
            return result;
        }
        if (error_x < tol) {
            result.status = StatusCode::SUCCESS;
            result.root = x;
            result.convergence = ConvergenceType::BY_X;
            result.iterations_used = i + 1;
            return result;
        }

        if (fa * fx< 0.0) {
            b = x;
            fb = fx;
        } else {
            a = x;
            fa = fx;
        }

        prev_x = x;
        has_prev= true;
    }

    result.status = StatusCode::MAX_ITERATIONS_REACHED;
    result.root = x;
    result.convergence = ConvergenceType::NONE;
    result.iterations_used = max_iter;
    return result;
}

SolverResult newton_raphson(const std::function<Dual(const Dual&)>& f, double a, double tol, int max_iter) {
    SolverResult result;
    result.method_name= "Newton-Raphson";
    result.iterations_used = 0;
    double current_x = a;

    Dual fd(0.0, 1.0);

    if(!safe_eval(f, current_x, fd, 1.0)) {
        result.status = StatusCode::DOMAIN_INVALID;
        result.convergence = ConvergenceType::NONE;
        result.root = current_x;
        result.iterations_used = 0;
        return result;
    }
    // chute inicial e a raiz
    if(std::fabs(fd.val) < tol) {
        result.status = StatusCode::SUCCESS;
        result.convergence = ConvergenceType::BY_Y;
        result.root = current_x;
        result.iterations_used = 0;
        return result;
    }

    result.history.push_back({0, current_x, fd.val, current_x, abs(fd.val)});

    for(int i = 1; i < max_iter; i++) {
        double fx = fd.val;
        double dfx = fd.der;

        if(std::isinf(fx) || std::isnan(fx)) {
            result.status = StatusCode::DIVERGENCE;
            result.convergence = ConvergenceType::NONE;
            result.root = current_x;
            result.iterations_used = i;
            return result;
        }

        if(std::isinf(dfx) || std::isnan(dfx)) {
            result.status = StatusCode::DIVERGENCE;
            result.convergence = ConvergenceType::NONE;
            result.root = current_x;
            result.iterations_used = i;
            return result;
        }

        // derivada mt proxima de zero, n deve prosseguir
        if(abs(dfx) <= 1e-15) {
            result.status = StatusCode::DIVISION_BY_ZERO;
            result.convergence = ConvergenceType::NONE;
            result.root = current_x;
            result.iterations_used = i;
            return result;
        }

        double delta = fx / dfx;
        double next_x = current_x - delta;

        // inicializa com zero por enquanto para avaliarmos somente fx
        Dual next_fd(0.0, 0.0);

        // derivada ocmo 1 porque preciso da derivada
        if(!safe_eval(f, next_x, next_fd, 1.0)) {
            result.status = StatusCode::DOMAIN_INVALID;
            result.convergence = ConvergenceType::NONE;
            result.root = next_x;
            result.iterations_used = i;
            return result;
        }

        double next_fx = next_fd.val;

        if(std::isinf(next_fx) || std::isnan(next_fx)) {
            result.status = StatusCode::DIVERGENCE;
            result.convergence = ConvergenceType::NONE;
            result.root = next_x;
            result.iterations_used = i;
            return result;
        }

        double error_x = std::fabs(delta);
        double error_y = std::fabs(next_fx);

        result.history.push_back({i, next_x, next_fx, error_x, error_y});

        if(error_y < tol) {
            result.status = StatusCode::SUCCESS;
            result.convergence = ConvergenceType::BY_Y;
            result.root = next_x;
            result.iterations_used = i;
            return result;
        }
        if(error_x < tol) {
            result.status = StatusCode::SUCCESS;
            result.convergence = ConvergenceType::BY_X;
            result.root = next_x;
            result.iterations_used = i;
            return result;
        }

        current_x = next_x;
        fd = next_fd;
    }

    result.status = StatusCode::MAX_ITERATIONS_REACHED;
    result.convergence = ConvergenceType::NONE;
    result.root = current_x;
    result.iterations_used = max_iter;
    
    return result;
}

SolverResult secante(const std::function<Dual(const Dual&)>& f, double x0, double x1, double tol, int max_iter) {
    SolverResult result;
    result.method_name = "Secante";
    result.iterations_used = 0;

    double f0, f1;

    try {
        f0 = f(Dual(x0)).val;
        f1 = f(Dual(x1)).val;
    } catch (const std::domain_error& e) {
        result.status = StatusCode::DOMAIN_INVALID;
        result.root = x0;
        result.convergence = ConvergenceType::NONE;
        return result;
    }

    // Checagem de NaN/Inf nas avaliações iniciais
    if (std::isnan(f0) || std::isinf(f0) || std::isnan(f1) || std::isinf(f1)) {
        result.status = StatusCode::DIVERGENCE;
        result.root = x1;
        result.convergence = ConvergenceType::NONE;
        return result;
    }

    if (std::fabs(f1) < tol) {
        result.status = StatusCode::SUCCESS;
        result.root = x1;
        result.convergence = ConvergenceType::BY_Y;
        result.iterations_used = 0;
        return result;
    }
    if (std::fabs(f0) < tol) {
        result.status = StatusCode::SUCCESS;
        result.root = x0;
        result.convergence = ConvergenceType::BY_Y;
        result.iterations_used = 0;
        return result;
    }

    for (int i = 0; i < max_iter; i++) {
        double denom = f1 - f0;

        if (denom == 0) {
            result.status = StatusCode::DIVISION_BY_ZERO;
            result.root = x1;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }

        double x2 = x1 - f1 * (x1 - x0) / denom;
        double f2;

        try {
            f2 = f(Dual(x2)).val;
        } catch (const std::domain_error& e) {
            result.status = StatusCode::DOMAIN_INVALID;
            result.root = x2;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }

        if (std::isnan(f2) || std::isinf(f2)) {
            result.status = StatusCode::DIVERGENCE;
            result.root = x2;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }

        double error_x = std::fabs(x2 - x1);
        double error_y = std::fabs(f2);

        result.history.push_back({i + 1, x2, f2, error_x, error_y});

        if (error_y < tol) {
            result.status = StatusCode::SUCCESS;
            result.root = x2;
            result.convergence = ConvergenceType::BY_Y;
            result.iterations_used = i + 1;
            return result;
        }
        if (error_x < tol) {
            result.status = StatusCode::SUCCESS;
            result.root = x2;
            result.convergence = ConvergenceType::BY_X;
            result.iterations_used = i + 1;
            return result;
        }

        // Desloca a janela: (x0, x1) -> (x1, x2)
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
    }

    result.status = StatusCode::MAX_ITERATIONS_REACHED;
    result.root = x1;
    result.convergence = ConvergenceType::NONE;
    result.iterations_used = max_iter;
    return result;
}

// ============================================================
// Estimativa da ordem de convergência
// ============================================================

double estimate_convergence_order(const SolverResult& result) {
    // Para estimar p precisamos de pelo menos 3 iterações registradas
    if (result.history.size() < 3) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> errors;
    errors.reserve(result.history.size());

    for (const auto& iteration : result.history) {
        double error = std::fabs(iteration.error_x);
        if (std::isfinite(error) && error > 0.0) {
            errors.push_back(error);
        }
    }

    if (errors.size() < 3) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    /*
        p ≈ log(e[k+1] / e[k]) / log(e[k] / e[k-1])
        Usamos os últimos três erros disponíveis.
    */
    double e0 = errors[errors.size() - 3];
    double e1 = errors[errors.size() - 2];
    double e2 = errors[errors.size() - 1];

    if (e0 <= 0.0 || e1 <= 0.0 || e2 <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double denominator = std::log(e1 / e0);
    if (std::fabs(denominator) < 1e-15) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double p = std::log(e2 / e1) / denominator;

    // Descarta ordens não-finitas ou espúrias decorrentes de ruído numérico
    if (!std::isfinite(p) || p <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return p;
}