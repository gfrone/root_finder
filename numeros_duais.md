# Números Duais e Diferenciação Automática

## 1. Motivação

Quando implementamos métodos como Newton-Raphson, precisamos da derivada `f'(x)` em cada iteração. Existem três formas de obter isso:

1. **Derivada analítica manual** — precisa, mas trabalhosa e propensa a erro se `f` mudar.
2. **Diferenças finitas** (`(f(x+h) - f(x-h)) / 2h`) — fácil de implementar, mas é uma *aproximação*, sensível à escolha de `h` (erro de truncamento vs. erro de arredondamento).
3. **Números duais** — dão a derivada **exata** (até a precisão de ponto flutuante), sem aproximação nenhuma, calculando `f(x)` e `f'(x)` ao mesmo tempo, na mesma passada.

A classe `Dual` do projeto é exatamente essa terceira abordagem, também chamada de **diferenciação automática em modo direto** (*forward-mode automatic differentiation*).

## 2. Definição algébrica

Um número dual é um par ordenado:

```
a + b·ε
```

onde `a` e `b` são números reais, e `ε` (épsilon) é um novo tipo de número que obedece a uma única regra:

```
ε² = 0        com        ε ≠ 0
```

Isso é estranho à primeira vista — nos reais, só `0² = 0`. `ε` é chamado de **infinitesimal nilpotente**: não é zero, mas seu quadrado é. Pense nele como "um número tão pequeno que qualquer coisa multiplicada por `ε²` é desprezível/zero".

Na struct do projeto, isso vira:

```cpp
class Dual {
public:
    double val;   // parte "real" -> a
    double der;   // parte "infinitesimal" -> b (coeficiente de ε)
};
```

- `val` guarda o **valor** da função naquele ponto.
- `der` guarda a **derivada** da função naquele ponto.

## 3. Operações algébricas (de onde vêm as fórmulas do código)

Toda operação com números duais nasce de aplicar álgebra comum e depois usar `ε² = 0` para simplificar. Tome dois duais `x = a + b·ε` e `y = c + d·ε`.

### Soma e subtração

```
x + y = (a + b·ε) + (c + d·ε) = (a + c) + (b + d)·ε
```

Direto: soma as partes reais, soma as partes derivadas. É exatamente `operator+` no código.

### Multiplicação

```
x·y = (a + b·ε)(c + d·ε)
    = ac + ad·ε + bc·ε + bd·ε²
    = ac + (ad + bc)·ε + bd·(0)      // porque ε² = 0
    = ac + (ad + bc)·ε
```

Parte real: `a·c` (o produto dos valores).
Parte derivada: `a·d + b·c` — que é **exatamente a regra do produto** da derivada: `(uv)' = u'v + uv'`.

Isso não é coincidência — é o ponto central dos números duais: a álgebra de `ε² = 0` **reproduz automaticamente as regras de derivação** que normalmente decoramos separadamente (regra do produto, do quociente, da cadeia).

### Divisão

Buscamos `x/y = p + q·ε` tal que `(p + q·ε)(c + d·ε) = a + b·ε`. Expandindo o lado esquerdo:

```
pc + (pd + qc)·ε = a + b·ε
```

Igualando parte real e parte de `ε`:

```
p = a/c
pd + qc = b   →   q = (b - pd)/c = (bc - ad)/c²
```

Ou seja:

```
x/y = (a/c) + [(b·c - a·d) / c²]·ε
```

Compare com a regra do quociente `(u/v)' = (u'v - uv') / v²` — é a mesma fórmula, com `a↔u`, `b↔u'`, `c↔v`, `d↔v'`. É exatamente o que está implementado em `Dual::operator/`.

## 4. A dedução "de verdade": série de Taylor

As contas acima mostram *que* as fórmulas batem com as regras de derivação, mas o motivo *por que* isso funciona em geral — para qualquer função `f` (não só `+`, `-`, `*`, `/`), inclusive `sin`, `exp`, `log`, etc — vem da série de Taylor.

### 4.1 Expansão de Taylor em torno de um ponto `a`

Para uma função suave `f`, a expansão de Taylor ao redor de `a`, avaliada em `a + h`, é:

```
f(a + h) = f(a) + f'(a)·h + f''(a)/2!·h² + f'''(a)/3!·h³ + ...
```

Essa série é exata (supondo convergência); ela expressa `f(a+h)` como uma soma infinita de termos, cada um envolvendo uma derivada de ordem cada vez mais alta.

### 4.2 Substituindo `h` por `ε`

Agora, em vez de um número real pequeno `h`, usamos o infinitesimal `ε`, com a propriedade `ε² = 0` (e, por extensão, `ε³ = ε²·ε = 0`, `ε⁴ = 0`, etc — toda potência de `ε` a partir da segunda é zero):

```
f(a + ε) = f(a) + f'(a)·ε + f''(a)/2!·ε² + f'''(a)/3!·ε³ + ...
         = f(a) + f'(a)·ε + f''(a)/2!·(0) + f'''(a)/3!·(0) + ...
         = f(a) + f'(a)·ε
```

**Todos os termos de ordem 2 em diante desaparecem sozinhos**, porque cada um carrega uma potência de `ε` maior ou igual a 2. Sobra só:

```
f(a + ε) = f(a) + f'(a)·ε
```

Essa é a fórmula fundamental por trás de cada função em `dual.hpp`. Ela diz: **se você avalia `f` num número dual `a + 1·ε` (isto é, `val = a`, `der = 1`), o resultado é um número dual cuja parte real é `f(a)` e cuja parte infinitesimal é `f'(a)`** — a derivada exata, sem nenhuma aproximação numérica, sem escolher nenhum `h`.

### 4.3 Verificando com um exemplo do código: `exp`

```cpp
Dual exp(const Dual& x) {
    double v = std::exp(x.val);
    return Dual(v, x.der * v);
}
```

Pela fórmula de Taylor: `exp(a + ε) = exp(a) + exp'(a)·ε = exp(a) + exp(a)·ε` (pois a derivada de `exp` é ela mesma). Isso bate exatamente com o código: `val = exp(x.val)` e `der = x.der * exp(x.val)` — o `x.der` aqui é a regra da cadeia entrando em ação (ver seção 4.4).

### 4.4 Regra da cadeia "de graça"

Quando `x` já é o resultado de outra expressão (não a variável independente pura), `x = g(t) = c + d·ε` onde `d = g'(t)` (não necessariamente `1`). Aplicando a mesma expansão de Taylor de `f` em torno de `c`:

```
f(x) = f(c + d·ε) = f(c) + f'(c)·(d·ε) = f(c) + [f'(c)·d]·ε
```

A parte infinitesimal do resultado é `f'(c)·d = f'(g(t))·g'(t)` — **exatamente a regra da cadeia**. É por isso que toda função do tipo `Dual f(const Dual& x)` no código multiplica a derivada interna (`x.der`) pela derivada externa avaliada em `x.val`:

```cpp
Dual sin(const Dual& x) {
    return Dual(std::sin(x.val), x.der * std::cos(x.val));
    //                            ^^^^^^   ^^^^^^^^^^^^^^^
    //                            g'(t)     f'(g(t))
}
```

Isso significa que, ao encadear várias operações (`sin(x*x + 1)`, por exemplo), a regra da cadeia se aplica automaticamente, camada por camada, sem que você precise derivá-la manualmente — é a própria mecânica dos números duais fazendo isso.

## 5. Propriedades algébricas resumidas

| Propriedade | Números reais | Números duais |
|---|---|---|
| Estrutura | Corpo (todo não-zero tem inverso) | Anel comutativo (**não** é corpo) |
| Elemento problemático | — | `ε` não tem inverso multiplicativo |
| Nilpotência | Só `0` satisfaz `xⁿ = 0` para algum `n` | `ε` é nilpotente: `ε² = 0`, mas `ε ≠ 0` |
| Divisores de zero | Não existem (`ab=0 → a=0 ou b=0`) | `ε` é um divisor de zero: `ε·ε = 0` sem que `ε = 0` |
| Divisão por `x = a + b·ε` | — | Só é definida quando `a ≠ 0` (a parte real não pode ser zero) — coerente com o código, que lança exceção quando `other.val == 0.0` |

O fato de os duais formarem um **anel** (e não um corpo) é justamente por causa de `ε`: ele quebra a existência de inverso multiplicativo para todo elemento não-nulo. Isso é consistente com a checagem de domínio em `operator/`: a divisão só é bem definida quando a parte real do divisor é diferente de zero.

## 6. Resumo

- Um número dual é `a + b·ε`, com a regra algébrica `ε² = 0`.
- Expandindo `f(a + ε)` em série de Taylor, todos os termos de ordem `≥ 2` desaparecem por conta de `ε² = 0`, sobrando `f(a) + f'(a)·ε`.
- Isso significa: **avaliar uma função num número dual calcula, na mesma operação, o valor da função (parte real) e sua derivada exata (parte de `ε`)**.
- As regras de derivação que normalmente aplicamos manualmente (produto, quociente, cadeia) emergem automaticamente da álgebra de `ε² = 0` — é por isso que `Dual::operator*`, `operator/` e cada função em `dual.hpp` têm exatamente a forma das regras de derivação clássicas.
- Isso é a base do que se chama **diferenciação automática em modo direto (forward-mode)**, usada no projeto para fornecer `f'(x)` exato ao método de Newton-Raphson, sem aproximação numérica.
