# Root Finder

Projeto em C++17 para encontrar raízes de funções numéricas com os métodos de bisseção, falsa posição, Newton-Raphson e secante. O projeto usa números duais para obter automaticamente a derivada necessária pelo método de Newton.

## Requisitos e execução

São necessários `g++` com suporte a C++17 e `make`.

```bash
# aplicação principal
make
./main_app

# testes de números duais
make test
./run_test_dual

# exercícios do Moodle
make moodle
./moodle_test
```

O programa principal está em `main.cpp`; os exercícios estão em `test_moodle/test_moodle.cpp`. No `main` deste último arquivo, descomente a chamada do exercício desejado e informe a tolerância, por exemplo `exercicio_11(1e-4)`.

Para remover apenas os executáveis gerados:

```bash
make clean
```

## Resultado dos solvers

Cada método retorna um `SolverResult` com:

| Campo | Descrição |
| --- | --- |
| `root` | Melhor aproximação de raiz encontrada. |
| `status` | Situação da execução, como `SUCCESS`, `DOMAIN_INVALID` ou `DIVERGENCE`. |
| `convergence` | Critério de término: `BY_Y`, `BY_X` ou `NONE`. |
| `iterations_used` | Número de iterações realizadas. |
| `history` | Histórico de `IterationData` usado nos relatórios. |
| `method_name` | Nome do método executado. |

Cada item de `history` contém `it` (iteração), `x_k` (aproximação), `f_x` (valor da função), `error_x` (erro em `x`), `error_y` (`|f(x_k)|`) e `step_x` (`|x_k - x_(k-1)|`).

`error_y` é o resíduo: indica se a equação está sendo satisfeita. `step_x` indica o quanto a aproximação ainda mudou, mas, isoladamente, não comprova que uma raiz foi encontrada.

### Validação e falhas numéricas

Tolerâncias não positivas, valores não finitos e pontos iniciais inválidos retornam `DOMAIN_INVALID`. Uma `std::domain_error` lançada pela função — por exemplo, `log(x)` fora do domínio ou divisão por zero — também produz esse status.

Valores ou aproximações não finitos produzem `DIVERGENCE`. Divisores exatamente nulos produzem `DIVISION_BY_ZERO`: no denominador da falsa posição e secante, ou na derivada de Newton. Um denominador pequeno, mas não nulo, ainda pode ser válido; se causar uma aproximação não finita, será classificado como divergência.

## Métodos implementados

### Bisseção

Recebe um intervalo `[a, b]` que precisa conter mudança de sinal. O código compara os sinais com `std::signbit(fa) == std::signbit(fb)`, em vez de multiplicar os valores; assim evita overflow em `fa * fb`. Se um extremo já zera a função, ele é retornado imediatamente.

O intervalo é dividido ao meio e a metade que preserva a mudança de sinal é mantida. Há sucesso quando `|f(meio)| <= tol` (`BY_Y`) ou quando metade da largura do intervalo atende a tolerância (`BY_X`). Neste método, `error_x` é metade da largura do intervalo e `step_x` é a distância entre pontos médios consecutivos. Há limite de 100 iterações.

### Falsa posição

Também exige mudança de sinal em `[a, b]`, mas usa a interseção da secante entre os extremos:

```text
x = a - (f(a) / (f(b) - f(a))) * (b - a)
```

Após calcular `f(x)`, substitui o extremo que tem o mesmo sinal de `f(x)`, preservando o intervalo que contém a raiz. Na primeira iteração, `error_x` é a largura inicial; nas demais, `error_x` e `step_x` são `|x - x_anterior|`. O sucesso é confirmado pelo resíduo. Se `x` repetir um extremo ou a aproximação anterior sem que o resíduo seja aceitável, o resultado é `STAGNATION`.

### Newton-Raphson

Usa a derivada obtida com números duais:

```text
x_(k+1) = x_k - f(x_k) / f'(x_k)
```

O ponto inicial é avaliado com derivada `1.0`. Se seu resíduo já atende a tolerância, o método termina sem iterações; caso contrário, ele entra no histórico como iteração 0. `error_x` e `step_x` são o módulo da correção de Newton. Derivada exatamente zero retorna `DIVISION_BY_ZERO`; o sucesso depende sempre de `|f(x)| <= tol`.

### Secante

Não precisa de derivada e usa os dois pontos anteriores:

```text
x_(k+1) = x_k - f(x_k) * (x_k - x_(k-1)) / (f(x_k) - f(x_(k-1)))
```

Os dois chutes são validados antes das iterações. Denominador zero retorna `DIVISION_BY_ZERO`; valores não finitos retornam `DIVERGENCE`. `error_x` e `step_x` são `|x_(k+1) - x_k|`, e a confirmação de sucesso depende do resíduo.

### Estagnação e ordem de convergência

Falsa posição, Newton e secante testam estagnação depois de verificar o sucesso por resíduo. Uma repetição exata da próxima aproximação causa `STAGNATION`. Também há uma verificação de progresso na escala de `double`: passos de até oito ULPs (`8 * epsilon * max(1, |x_atual|, |x_anterior|)`) em três iterações consecutivas são considerados falta de progresso.

A ordem local é estimada pelos três últimos erros válidos em `x`:

```text
p ~= log(e_(k+1) / e_k) / log(e_k / e_(k-1))
```

Sem três erros adequados, ou se o cálculo não for finito, a ordem é exibida como `N/A`.

## Números duais e diferenciação automática

Um dual tem a forma `a + b*ε`, em que `ε² = 0`, mas `ε != 0`. Na classe `Dual`, `val` armazena `a` e `der` armazena `b`.

Para `x = a + b*ε` e `y = c + d*ε`, as operações reproduzem as regras de derivação:

```text
x + y = (a + c) + (b + d)*ε
x * y = ac + (ad + bc)*ε
x / y = (a/c) + ((bc - ad)/c²)*ε
```

A razão está na série de Taylor. Como todas as potências de `ε` a partir da segunda são zero:

```text
f(a + ε) = f(a) + f'(a)*ε
```

Assim, ao avaliar uma função com `Dual(a, 1.0)`, a parte real do resultado é `f(a)` e a parte derivada é `f'(a)`, até a precisão de ponto flutuante. Para uma expressão intermediária `g(t)`, a parte derivada se torna `f'(g(t)) * g'(t)`: a regra da cadeia é aplicada automaticamente. Esse mecanismo é a diferenciação automática em modo direto e permite a Newton obter valor e derivada numa única avaliação.

Os duais formam um anel, e não um corpo: a divisão por `a + b*ε` só existe quando `a != 0`. Isso corresponde à checagem de domínio de `Dual::operator/`.

## Mini guia de `reporter.cpp`

Inclua `reporter.hpp` para usar as funções abaixo. O parâmetro `every` controla a periodicidade das linhas do histórico; valores menores ou iguais a zero são tratados como `1`.

| Função | Uso |
| --- | --- |
| `print_history(result, every)` | Imprime a tabela de iterações de um `SolverResult`; sempre mostra a primeira e a última linha. |
| `report_method(result, every)` | Imprime o relatório completo de um resultado e usa `result.method_name`. |
| `report_method(nome, result, every)` | Igual à anterior, mas permite informar um nome personalizado. |
| `report_bissecao(f, a, b, tol, every)` | Executa a bisseção, imprime o relatório e devolve o resultado. |
| `report_falsa_posicao(f, a, b, tol, max_iter, every)` | Executa e relata falsa posição. |
| `report_newton(f, x0, tol, max_iter, every)` | Executa e relata Newton-Raphson. |
| `report_secante(f, x0, x1, tol, max_iter, every)` | Executa e relata secante. |
| `report_method(method, f, a, b, tol, max_iter, every)` | Despacha para um método a partir do enum `Method`. Para Newton, `a` é o chute inicial; para secante, `a` e `b` são os dois chutes. |
| `report_all(f, a, b, tol, max_iter)` | Executa todos os métodos, imprime tabela comparativa e os erros finais dos que convergiram. |
| `benchmark_metodos(f, a, b, tol, max_iter)` | Executa todos os métodos e os ordena pela ordem de convergência estimada; métodos sem ordem válida ficam como `N/A`. |
| `run_suite_comparison(casos)` | Executa `report_all` para cada `TestCase`, exibindo nome, expressão e parâmetros. |
| `run_suite_comparison(funcoes, a, b, tol, max_iter)` | Sobrecarga que cria casos simples para uma lista de funções. |

Exemplo de uso:

```cpp
auto f = [](const Dual& x) {
    return x * x - 2.0;
};

report_newton(f, 1.0, 1e-8);
report_all(f, 0.0, 2.0, 1e-8);
```

Os relatórios individuais mostram status, número de iterações, critério de convergência, raiz aproximada, ordem estimada e erros finais. A ordem e os erros finais são detalhados quando a execução tem sucesso.
