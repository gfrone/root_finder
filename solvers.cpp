#include <math.h>
#include <iostream>
#include <limits>
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

namespace {
// Esse trecho implementa uma verificação de estagnação numérica 
//nos métodos iterativos. 
//A função compara duas aproximações consecutivas: oito vezes a precisão relativa do tipo \texttt{double}.

constexpr double STAGNATION_ULPS = 8.0;

bool has_no_numeric_progress(double current_x, double previous_x) {
    const double scale = std::fmax(
        1.0,
        std::fmax(std::fabs(current_x), std::fabs(previous_x))
    );

    return std::fabs(current_x - previous_x) <=
           STAGNATION_ULPS * std::numeric_limits<double>::epsilon() * scale;
}

} // namespace

SolverResult bissecao(const std::function<Dual(const Dual&)>& f, double a, double b, double tol) {
    SolverResult result;
    result.method_name = "Bissecao";
    result.iterations_used=0;
    double fa, fb;

    if (tol <= 0.0 || !std::isfinite(tol) || !std::isfinite(a) || !std::isfinite(b)) {
        result.root = a;
        result.convergence = ConvergenceType::NONE;
        result.status = StatusCode::DOMAIN_INVALID;
        return result;
    }

    if (a == b) {
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

    if (!std::isfinite(fa) || !std::isfinite(fb)) {
        result.status = StatusCode::DIVERGENCE;
        result.root = a;
        result.convergence = ConvergenceType::NONE;
        return result;
    }

    if (fa == 0.0 || fb == 0.0) {
        double root;
        root = (fa == 0) ? a : b;
        result.root = root;
        result.convergence = ConvergenceType::BY_Y;
        result.status = StatusCode::SUCCESS;
        return result;
    }

    if (std::signbit(fa) == std::signbit(fb)) {
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
            // se numero n de iterações é finito e maior que zero, pega ele, garantido que sempre sera
            // um inteiro maior, caso der numero fracionario o ceil trata acima
            num_iters = (n < max_iter) ? std::ceil(n) : max_iter;
        }
    }
    else {
        // se intervalo ja é maior que  tol, basta somente avaliar no ponto medio
        num_iters = 1;
    }
    

    double meio = a / 2.0 + b / 2.0;
    double previous_meio = 0.0;
    double fm;

    for(int i = 0; i < num_iters; i++) {
        meio = a / 2.0 + b / 2.0;

        // Se o ponto médio for igual a um dos extremos, chegamos ao limite do hardware (double)
        //O tipo double (padrão IEEE 754) possui 53 bits de mantissa 
        // (aproximadamente 15 a 17 dígitos significativos).
    
        if (meio == a || meio == b) {
            const double error_x = std::fabs(b - a);
            result.root = meio;
            result.iterations_used = i + 1;
            
            if(error_x <= tol) {
                result.status = StatusCode::SUCCESS;
                result.convergence = ConvergenceType::BY_X;
            }
            else {
                result.status = StatusCode::STAGNATION;
                result.convergence = ConvergenceType::NONE;
            }

            return result;
        }

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
        double error_x = std::fabs(b - a) / 2;
        double error_y = std::fabs(fm);
        // Step para calculo da convergencia
        double step_x = (i == 0) ? 0.0 : std::fabs(meio - previous_meio);

        result.history.push_back({i + 1, meio, fm, error_x, error_y, step_x});
        previous_meio = meio;


        /* Atenção: Como o numero de iterações garante, por definição, a covergencia pelo eixo x
        Façamos uma checagem em relação ao eixo y, para verificar se convergiu mais rapidamente.                            */
        if (error_y <= tol) {
            result.status = StatusCode::SUCCESS;
            result.root = meio;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::BY_Y;
            return result;
        }

        // Ainda assim, decidi manter o calculo do erro_x
        if (error_x <= tol) {
            result.status = StatusCode::SUCCESS;
            result.root = meio;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::BY_X;
            return result;
        }

        if (fm == 0.0 || std::signbit(fa) != std::signbit(fm)) {
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

    if (tol <= 0.0 || !std::isfinite(tol) || !std::isfinite(a) || !std::isfinite(b)) {
        result.root = a;
        result.convergence = ConvergenceType::NONE;
        result.status = StatusCode::DOMAIN_INVALID;
        return result;
    }

    // intervalos iguais já são tratados antes de chamar as funções. É preciso de um intervalo [a,b] tal que b > a
    if (a == b) {
        result.root = a;
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
    
    if (!std::isfinite(fa) || !std::isfinite(fb)) {
        result.status = StatusCode::DIVERGENCE;
        result.root = a;
        result.convergence = ConvergenceType::NONE;
        return result;
    }

    if (std::fabs(fa) <= tol || std::fabs(fb) <= tol) {
        double root;
        root = (std::fabs(fa) <= tol) ? a : b;
        result.root = root;
        result.convergence = ConvergenceType::BY_Y;
        result.status = StatusCode::SUCCESS;
        return result;
    }

    if (std::signbit(fa) == std::signbit(fb)) {
        result.status = StatusCode::DOMAIN_INVALID;
        result.root = a;
        result.convergence = ConvergenceType::NONE;
        return result;
    } 

    int max_it = (max_iter <= 0) ? 100 : max_iter;
    int stagnation_count = 0;
    bool has_prev = false;
    double prev_x = 0.0;
    double x = a;
    double error_x = 0.0;
    double error_y = 0.0;

    for(int i = 0; i < max_it; i++) {
        double denom = fb - fa;

        // Proteção contra denominadores quase nulos
        if (denom == 0.0) {
            result.status = StatusCode::DIVISION_BY_ZERO;
            result.root = has_prev ? prev_x : a;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }
        
        x = a - (fa / denom) * (b - a);

        if (!std::isfinite(x)) {
            result.status = StatusCode::DIVERGENCE;
            result.root = has_prev ? prev_x : a;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }

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
        double step_x = has_prev ? std::fabs(x - prev_x) : 0.0;

        result.history.push_back({i + 1, x, fx, error_x, error_y, step_x});

        if (error_y <= tol) {
            result.status = StatusCode::SUCCESS;
            result.root = x;
            result.convergence = ConvergenceType::BY_Y;
            result.iterations_used = i + 1;
            return result;
        }

        // Um extremo repetido significa que a interpolacao nao consegue mais
        // produzir um novo double. Como o residuo nao atingiu a tolerancia,
        // o metodo estagnou.
        if (x == a || x == b || (has_prev && x == prev_x)) {
            result.status = StatusCode::STAGNATION;
            result.root = x;
            result.convergence = ConvergenceType::NONE;
            result.iterations_used = i + 1;
            return result;
        }

        // Proximidade de poucas ULPs pode ser apenas ruido de arredondamento;
        // exigimos repeticao para nao encerrar por uma unica iteracao.
        if (has_prev && has_no_numeric_progress(x, prev_x)) {
            stagnation_count++;
            if (stagnation_count >= 3) {
                result.status = StatusCode::STAGNATION;
                result.root = x;
                result.convergence = ConvergenceType::NONE;
                result.iterations_used = i + 1;
                return result;
            }
        } else {
            stagnation_count = 0;
        }

        if (fx == 0.0 || std::signbit(fa) != std::signbit(fx)) {
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
    result.iterations_used = max_it;
    return result;
}

SolverResult newton_raphson(const std::function<Dual(const Dual&)>& f, double a, double tol, int max_iter) {
    SolverResult result;
    result.method_name= "Newton-Raphson";
    int stagnation_count = 0;
    result.iterations_used = 0;
    double current_x = a;

    if (tol <= 0.0 || !std::isfinite(tol) || !std::isfinite(a)) {
        result.status = StatusCode::DOMAIN_INVALID;
        result.convergence = ConvergenceType::NONE;
        result.root = a;
        result.iterations_used = 0;
        return result;
    }

    // limite max iter
    if (max_iter <= 0) {
        max_iter = 100;
    }

    Dual fd(0.0, 1.0);

    if(!safe_eval(f, current_x, fd, 1.0)) {
        result.status = StatusCode::DOMAIN_INVALID;
        result.convergence = ConvergenceType::NONE;
        result.root = current_x;
        result.iterations_used = 0;
        return result;
    }

    // chute inicial é a raiz
    if(std::fabs(fd.val) <= tol || fd.val == 0.0) {
        result.status = StatusCode::SUCCESS;
        result.convergence = ConvergenceType::BY_Y;
        result.root = current_x;
        result.iterations_used = 0;
        return result;
    }

    // Na iteracao inicial nao existe x_{k-1}; portanto, nao ha erro/step em X.
    result.history.push_back({0, current_x, fd.val, 0.0, std::fabs(fd.val), 0.0});

    for(int i = 1; i <= max_iter; i++) {
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

        // derivada igual a zero, n deve prosseguir
        if (dfx == 0.0) {
            result.status = StatusCode::DIVISION_BY_ZERO;
            result.convergence = ConvergenceType::NONE;
            result.root = current_x;
            result.iterations_used = i;
            return result;
        }

        double delta = fx / dfx;

        // Porteção contra double_max
        if (!std::isfinite(delta)) {
            result.status = StatusCode::DIVERGENCE;
            result.convergence = ConvergenceType::NONE;
            result.root = current_x;
            result.iterations_used = i;
            return result;
        }


        double next_x = current_x - delta;

        // Problema percebido: Em algumas situações, newton pode crescer indefinidamente
        // Problemas que mostram a divergencia de newton = x/2 - tan(x)
        // Para isso, vamos colocar um limite na divergencia deste método

        if (!std::isfinite(next_x)) {
            result.status = StatusCode::DIVERGENCE;
            result.root = current_x;
            result.iterations_used = i;
            return result;
        }

        // inicializa com zero por enquanto para avaliarmos somente fx
        Dual next_fd(0.0, 0.0);

        // derivada como 1 porque preciso da derivada
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

        double step_x = std::fabs(next_x - current_x);
        double error_x = step_x;
        double error_y = std::fabs(next_fx);

        result.history.push_back({i, next_x, next_fx, error_x, error_y, step_x});

        // Primeiro confirmamos que o novo ponto realmente e uma raiz.
        if (error_y <= tol) {
            result.status = StatusCode::SUCCESS;
            result.convergence = ConvergenceType::BY_Y;
            result.root = next_x;
            result.iterations_used = i;
            return result;
        }

        // Um ponto repetido nao pode gerar novo progresso em double.
        if (next_x == current_x) {
            result.status = StatusCode::STAGNATION;
            result.root = next_x;
            result.convergence = ConvergenceType::NONE;
            result.iterations_used = i;
            return result;
        }

        // Proximidade de poucas ULPs exige recorrencia para caracterizar
        // estagnacao, pois um unico passo pequeno pode ser transitorio.
        if (has_no_numeric_progress(next_x, current_x)) {
            stagnation_count++;

            // Um unico passo pequeno nao e suficiente para concluir estagnacao.
            if (stagnation_count >= 3) {
                result.status = StatusCode::STAGNATION;
                result.root = next_x;
                result.convergence = ConvergenceType::NONE;
                result.iterations_used = i;
                return result;
            }
        } else {
            stagnation_count = 0;
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
    int stagnation_count = 0; // critério para definir quando um método estagnou

    if (tol <= 0.0 || !std::isfinite(tol) || !std::isfinite(x0) || !std::isfinite(x1)) {
        result.status = StatusCode::DOMAIN_INVALID;
        result.convergence = ConvergenceType::NONE;
        result.root=x1;
        return result;
    }
    if (max_iter <= 0) {
        max_iter = 100;
    }

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

    if (std::fabs(f1) <= tol) {
        result.status = StatusCode::SUCCESS;
        result.root = x1;
        result.convergence = ConvergenceType::BY_Y;
        result.iterations_used = 0;
        return result;
    }
    if (std::fabs(f0) <= tol) {
        result.status = StatusCode::SUCCESS;
        result.root = x0;
        result.convergence = ConvergenceType::BY_Y;
        result.iterations_used = 0;
        return result;
    }

    for (int i = 0; i < max_iter; i++) {
        double denom = f1 - f0;

        // Proteção contra denominador nulo
        if (denom == 0.0) {
            result.status = StatusCode::DIVISION_BY_ZERO;
            result.root = x1;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }

        double x2 = x1 - f1 * (x1 - x0) / denom;
        if(!std::isfinite(x2)) {
            result.status = StatusCode::DIVERGENCE;
            result.root = x1;
            result.iterations_used = i + 1;
            result.convergence = ConvergenceType::NONE;
            return result;
        }

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

        double step_x = std::fabs(x2 - x1);
        double error_y = std::fabs(f2);

        // Na Secante, error_x e mantido apenas para o relatorio; o deslocamento
        // entre iteracoes e representado explicitamente por step_x.
        result.history.push_back({i + 1, x2, f2, step_x, error_y, step_x});

        if (error_y <= tol) {
            result.status = StatusCode::SUCCESS;
            result.root = x2;
            result.convergence = ConvergenceType::BY_Y;
            result.iterations_used = i + 1;
            return result;
        }
        if (x2 == x1) {
            result.status = StatusCode::STAGNATION;
            result.root = x2;
            result.convergence = ConvergenceType::NONE;
            result.iterations_used = i + 1;
            return result;
        }

        if (has_no_numeric_progress(x2, x1)) {
            // Proximidade de poucas ULPs so caracteriza estagnacao quando se
            // repete em iteracoes consecutivas.
            stagnation_count++;
            if (stagnation_count >= 3) {
                result.status = StatusCode::STAGNATION;
                result.root = x2;
                result.convergence = ConvergenceType::NONE;
                result.iterations_used = i + 1;
                return result;
            }
        } else {
            stagnation_count = 0;
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
        errors.push_back(error);
    }

    if (errors.size() < 3) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    /*
        p ≈ log(e[k+1] / e[k]) / log(e[k] / e[k-1])
        Usamos os últimos três erros disponíveis.
    */
    const auto& h0 = result.history[result.history.size() - 3];
    const auto& h1 = result.history[result.history.size() - 2];
    const auto& h2 = result.history[result.history.size() - 1];

    if(h1.it != h0.it + 1 || h2.it != h1.it + 1){
        return std::numeric_limits<double>::quiet_NaN();
    }

    double e0 = std::fabs(h0.error_x);
    double e1 = std::fabs(h1.error_x);
    double e2 = std::fabs(h2.error_x);

    // Rejeitar estimativa se algum erro ou zero ou não finito 
    if (!std::isfinite(e0) || !std::isfinite(e1) || !std::isfinite(e2) ||
        e0 <= 0.0 || e1 <= 0.0 || e2 <= 0.0 ||
        e1 >= e0 || e2 >= e1) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double denom = std::log(e1 / e0);
    if (std::fabs(denom) < 1e-15) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return  std::log(e2 / e1) / denom;
}
