#include "dual.hpp"
#include <math.h>
#include <stdexcept>

using namespace std;

// operações entre Dual, ex: Dual + Dual
// Falta implementar POW, tipo 2^x
Dual::Dual(double val, double der)
    : val(val), der(der) {}

Dual Dual::operator+(const Dual& other) const {
    // (a + be) + (c + de) = (a + c) + e(b + d)
    return Dual(
        val + other.val,
        der + other.der
    );
}

Dual Dual::operator-(const Dual& other) const {
    // (a + be) - (c + de) = (a - c) + e(b - d)
    return Dual(
        val - other.val,
        der - other.der
    );
}

Dual Dual::operator*(const Dual& other) const {
    // (a + be) * (c + de) = ac + ade + cbe = (ac) + e(ad + cb). Derivada é ad + cb, igual regra da multiplicação
    return Dual(
        val * other.val,
        der * other.val + val * other.der
    );
}

Dual Dual::operator/(const Dual& other) const {
    if (other.val == 0.0) {
        throw std::domain_error("Divisao por zero no operador / (Dual / Dual)");
    }

    return Dual(
        val / other.val,
        (der * other.val - val * other.der)
            / (other.val * other.val)
    );
}

// -Dual = -val, -der
Dual Dual::operator-() const {
    return Dual(-val, -der);
}
// Dual + float - so altera neste caso o val. Derivada de x - c em x continua sendo x'
Dual Dual::operator+(double other) const {
    return Dual(val + other, der);
}

// Dual - float
Dual Dual::operator-(double other) const {
    return Dual(val - other, der);
}

// Dual * float -> Derivada de C*x = x' * C -> der * C 
Dual Dual::operator*(double other) const {
    return Dual(val * other, der * other);
}

// Dual / float -> Derivada de x/C = 1/C * x' = der / other
Dual Dual::operator/(double other) const {
    if (other == 0.0) {
        throw std::domain_error("Divisao por zero no operador / (Dual / double)");
    }
    return Dual(val / other, der / other);
}

// Dual Dual::operator-(const Dual& x) const {
//     return Dual(-val, -der);
// }

Dual exp(const Dual& x) {
    double v = std::exp(x.val);

    //  - Derivada de e^x é = d/dx * (e^x)'
    
    return Dual(v, x.der * v);
}

Dual log(const Dual& x) {
    if(x.val <= 0.0) {
        throw std::domain_error("log: argumento deve ser maior que zero");
    }

    double v = std::log(x.val);
    //  - Derivada de ln(x) é = d/dx * 1/x
    return Dual(v, x.der / x.val);
}

Dual log(const Dual& x, double base) {

    if(x.val <= 0.0) {
        throw std::domain_error("log: argumento deve ser maior que zero");
    }

    if (base <= 0 || base == 1){
        throw std::domain_error("base do logarítimo deve ser > 0 e diferente de 1");
    }

    double ln_base = std::log(base);

    return Dual(
        std::log(x.val) / ln_base,
        x.der / (x.val * ln_base)
    );
}

Dual sqrt(const Dual& x) {
    if(x.val <= 0.0) {
        throw std::domain_error("domínio deve ser >= 0");
    }

    if(x.val == 0.0) {
        throw std::domain_error("A derivada não existe neste ponto - divisão por zero");
    }

    double sqrt_val = std::sqrt(x.val);

    return Dual(sqrt_val, x.der * (1 / (2 * sqrt_val)));
}

// Trigonométricas

Dual sin(const Dual &x){
    return Dual(std::sin(x.val), x.der * std::cos(x.val));
}

Dual cos(const Dual &x){
    return Dual(std::cos(x.val), -x.der * std::sin(x.val));
}

Dual sec(const Dual &x){
    double cos = std::cos(x.val);

    if (std::fabs(cos) < 1e-9) {
        throw std::domain_error("Sec proximo à assintota");
    }

    double sec_v = 1 / cos;

    return Dual(sec_v, x.der * (std::tan(x.val) * sec_v));
}

Dual csc(const Dual &x) {
    // nao implementada
    return Dual(0.0, 0.0);
}

Dual tan(const Dual &x){
    double cos = std::cos(x.val);

    if (std::fabs(cos) < 1e-9) {
        throw std::domain_error("tan: Valor proximo à assintota");
    }
    // Assíntota em pi/2 + kpi. Atualmente não implementa algo sobre a isso
    double sec_val = 1 / cos;
    return Dual(std::tan(x.val), x.der * (sec_val * sec_val));
}



// Inversas Trigonométricas
Dual asin(const Dual& x){
    // intervalo asin é de [-1 a 1]
    // O intervalo da derivada é (-1, 1), por isso o operador >=
    if (abs(x.val) >= 1) {
        throw std::domain_error("para asin, deve ser no intervalo de [-1, 1]");
    }
    // derivada é d/dx * (asin(x)) = 1/sqrt(1+x²)'
    double asin_value = std::asin(x.val);
    double x_squared = std::pow(x.val, 2);

    return Dual(asin_value, x.der * (1.0 / std::sqrt(1.0 - x_squared)));
}

Dual acos(const Dual& x){
    // intervalo acos é de [-1 a 1]
    if (abs(x.val) >= 1) {
        throw std::domain_error("para asin, deve ser no intervalo de (-1, 1)");
    }
    // derivada é inverso do asin'
    double acos_value = std::acos(x.val);
    double x_squared = std::pow(x.val, 2);

    return Dual(acos_value, - x.der * (1.0 / std::sqrt(1 + x_squared)));
}

Dual atan(const Dual& x){
    
    // derivada é d/dx * 1/ 1+x²'
    double atan_value = std::atan(x.val);
    double x_squared = std::pow(x.val, 2);

    return Dual(atan_value,  x.der  / (1 + x_squared));
}

// Operações com Dual e Float

Dual operator+(double lhs, const Dual& x){
    return Dual(x.val + lhs, x.der);
}

Dual operator-(double lhs, const Dual& x){
    return Dual(lhs - x.val, -x.der);
}

Dual operator*(double lhs, const Dual& x){
    return Dual(lhs * x.val, lhs * x.der);
}

Dual operator/(double lhs, const Dual& x){
    if (x.val == 0.0) {
        throw std::domain_error("Divisao por zero no operador / (Dual / Dual)");
    }
    // derivada de c/ x é  (-c * x'/ x²)
    return Dual(lhs / x.val, -lhs * x.der / (x.val * x.val));
}
