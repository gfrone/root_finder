# Implementação dos métodos numéricos

Este documento descreve como os métodos de busca de raízes foram implementados em `solvers.cpp`. Todos recebem uma função `f` que opera sobre `Dual`, mas usam principalmente `f(x).val` como o valor numérico da função. No Newton-Raphson, `f(x).der` também fornece a derivada automaticamente.

## Resultado de um método

Cada método retorna um `SolverResult`, contendo:

- `root`: melhor aproximação de raiz encontrada;
- `status`: resultado da execução (`SUCCESS`, `DOMAIN_INVALID`, `DIVERGENCE`, etc.);
- `convergence`: critério que permitiu encerrar com sucesso;
- `iterations_used`: número de iterações realizadas;
- `history`: dados de cada iteração para relatório e análise;
- `method_name`: nome do método.

Os possíveis critérios de convergência são:

- `BY_Y`: o resíduo foi suficientemente pequeno, isto é, `|f(x)| <= tol`;
- `BY_X`: a incerteza ou erro estimado em `x` ficou menor que a tolerância;
- `NONE`: não houve convergência confirmada.

## Dados registrados em cada iteração

Cada item de `history` é um `IterationData`:

| Campo | Significado |
| --- | --- |
| `it` | Número da iteração. |
| `x_k` | Aproximação atual da raiz. |
| `f_x` | Valor de `f(x_k)`. |
| `error_x` | Estimativa de erro no eixo `x`, específica de cada método. |
| `error_y` | `|f(x_k)|`, também chamado de resíduo. |
| `step_x` | Distância percorrida desde a aproximação anterior: `|x_k - x_(k-1)|`. |

`error_y` responde à pergunta “o ponto encontrado realmente satisfaz a equação?”. Já `step_x` responde “o método ainda está mudando sua aproximação?”. Um passo pequeno, sozinho, não prova que existe uma raiz: o método pode ter parado perto de um ponto que não zera a função. Por isso, os métodos iterativos verificam primeiro `error_y` e somente depois investigam estagnação.

## Tratamento de erros numéricos e de domínio

Antes de iniciar, os métodos validam tolerância e pontos iniciais. Tolerância não positiva, `NaN`, `Inf` ou entradas não finitas resultam em `DOMAIN_INVALID`.

As funções de números duais podem lançar `std::domain_error`, por exemplo em `log(x)` para `x <= 0`, divisão por zero ou tangente próxima de uma assíntota. Quando isso ocorre durante a avaliação de `f`, o solver encerra com `DOMAIN_INVALID`.

Mesmo sem lançar exceção, uma função pode retornar `NaN` ou `Inf`. Toda avaliação relevante é verificada com `std::isfinite`; se falhar, o método retorna `DIVERGENCE`. Também são verificadas as aproximações calculadas (`x`, `x2` ou `next_x`) antes de usá-las.

Divisão por zero é tratada explicitamente nos métodos que usam denominador:

- Falsa posição e secante retornam `DIVISION_BY_ZERO` quando a diferença entre os valores de função no denominador é exatamente zero.
- Newton-Raphson retorna `DIVISION_BY_ZERO` quando a derivada calculada é exatamente zero.

O teste é exato de propósito: um denominador pequeno pode ser matematicamente válido. Caso ele cause uma aproximação não finita, o método retorna `DIVERGENCE` no teste seguinte.

## Bisseção

A bisseção exige um intervalo `[a, b]` com mudança de sinal. Essa verificação não é feita com a expressão usual `f(a) * f(b) > 0`. Em vez disso, o código compara os sinais:

```cpp
std::signbit(fa) == std::signbit(fb)
```

Em `double`, o sinal é armazenado em um bit próprio da representação IEEE 754: `0` para número não negativo e `1` para número negativo. `std::signbit` consulta esse bit de sinal. Assim:

- sinais iguais significam que não há mudança de sinal no intervalo;
- sinais diferentes significam que há mudança de sinal;
- os extremos que valem exatamente zero são tratados antes dessa comparação como raízes imediatas.

Além de expressar diretamente a intenção, esse teste evita multiplicar valores muito grandes. Por exemplo, mesmo que `fa` e `fb` sejam finitos, o produto `fa * fb` pode transbordar para `Inf`; comparar os bits de sinal não realiza nenhuma multiplicação e, portanto, não sofre esse problema.

Fluxo principal:

1. Rejeita extremos iguais e valida as avaliações nos extremos.
2. Se `f(a)` ou `f(b)` for exatamente zero, retorna aquele extremo como raiz.
3. Confirma que os valores nos extremos possuem sinais opostos.
4. Avalia o ponto médio e mantém a metade do intervalo que preserva a mudança de sinal.
5. Para quando `|f(meio)| <= tol` (`BY_Y`) ou quando metade da largura do intervalo atende a tolerância (`BY_X`).

Na bisseção, `error_x` é metade da largura atual do intervalo, pois existe uma garantia teórica de que a raiz está nesse intervalo. `step_x` é a diferença entre dois pontos médios consecutivos.

O número de iterações é limitado a 100 e é calculado usando a redução do intervalo pela metade. Se o ponto médio se tornar igual a um dos extremos, o tipo `double` não consegue mais representar um ponto interno diferente.

## Falsa posição

A falsa posição também requer mudança de sinal em `[a, b]`. A validação inicial usa a mesma comparação de sinais da bisseção:

```cpp
std::signbit(fa) == std::signbit(fb)
```

Se os sinais forem iguais, o intervalo não é válido para o método. Se forem diferentes, a raiz permanece delimitada. Durante as iterações, o extremo substituído também é escolhido comparando sinais:

```cpp
std::signbit(fa) != std::signbit(fx)
```

Quando os sinais são diferentes, `x` passa a ser o novo extremo `b`; caso contrário, ele substitui `a`. Tal como na bisseção, valores exatamente zero são tratados antes, e a comparação por `std::signbit` evita risco de overflow em produtos como `fa * fb` ou `fa * fx`.

Em vez de escolher o ponto médio, ela usa a interseção da reta que liga `(a, f(a))` a `(b, f(b))` com o eixo `x`:

```text
x = a - (f(a) / (f(b) - f(a))) * (b - a)
```

Depois de avaliar `f(x)`, o extremo com o mesmo sinal de `f(x)` é substituído por `x`, preservando o intervalo que contém a raiz.

Neste método:

- na primeira iteração, `error_x` é a largura inicial do intervalo;
- nas demais, `error_x` e `step_x` são `|x - x_anterior|`;
- o sucesso é confirmado por `|f(x)| <= tol`.

A falsa posição pode ficar presa porque um extremo do intervalo muda muito pouco. Se o novo `x` for idêntico a `a`, `b` ou ao `x` anterior, e o resíduo ainda não for pequeno, o método retorna `STAGNATION` imediatamente.

## Newton-Raphson

Newton-Raphson usa números duais para obter valor e derivada numa única avaliação:

```text
x_(k+1) = x_k - f(x_k) / f'(x_k)
```

O ponto inicial é avaliado com derivada inicial igual a `1.0`. Caso já satisfaça `|f(x_0)| <= tol`, o método termina sem iterações. Caso contrário, o ponto inicial é registrado no histórico como iteração `0`.

Em cada iteração, o método verifica se valor, derivada, passo e nova avaliação são finitos. Derivada exatamente zero impede a divisão e produz `DIVISION_BY_ZERO`. O sucesso ocorre apenas quando o resíduo da nova aproximação atende `tol`.

Para Newton, `error_x` e `step_x` são ambos `|x_(k+1) - x_k|`. Eles informam o tamanho da correção de Newton, mas não são usados isoladamente como prova de raiz.

## Secante

A secante não precisa da derivada. Ela usa os dois últimos pontos:

```text
x_(k+1) = x_k - f(x_k) * (x_k - x_(k-1)) / (f(x_k) - f(x_(k-1)))
```

Antes das iterações, são testadas as duas aproximações iniciais. Se uma delas já satisfaz a tolerância em `y`, o método retorna sucesso. Durante a execução, denominador exatamente zero resulta em `DIVISION_BY_ZERO`; valores não finitos produzem `DIVERGENCE`.

Para a secante, `error_x` e `step_x` são `|x_(k+1) - x_k|`, e a confirmação de raiz também depende de `|f(x_(k+1))| <= tol`.

## Critério de estagnação

Newton-Raphson, secante e falsa posição usam dois níveis de detecção, sempre depois do teste de sucesso por resíduo.

### Repetição exata

Se a próxima aproximação é exatamente igual à anterior, não existe progresso representável em `double`:

- Newton: `next_x == current_x`;
- secante: `x2 == x1`;
- falsa posição: `x == a`, `x == b` ou, a partir da segunda iteração, `x == prev_x`.

Nesses casos, o método retorna `STAGNATION` imediatamente, desde que `|f(x)| > tol`.

### Proximidade na escala do `double`

Também existe a função auxiliar `has_no_numeric_progress`:

```cpp
constexpr double STAGNATION_ULPS = 8.0;

bool has_no_numeric_progress(double current_x, double previous_x) {
    const double scale = std::fmax(
        1.0,
        std::fmax(std::fabs(current_x), std::fabs(previous_x))
    );

    return std::fabs(current_x - previous_x) <=
           STAGNATION_ULPS *
           std::numeric_limits<double>::epsilon() * scale;
}
```

`epsilon` é aproximadamente `2.22e-16`, a precisão relativa do tipo `double`. A escala cresce com a magnitude de `x`, pois números muito grandes possuem espaçamento maior entre valores representáveis. O piso `1.0` cria uma tolerância absoluta mínima perto de zero. A margem de oito ULPs reduz falsos positivos devido aos arredondamentos das operações intermediárias.

Uma única aproximação muito próxima pode ser transitória. Portanto, esse teste precisa ocorrer em três iterações consecutivas para retornar `STAGNATION`. Se houver um passo numericamente significativo, o contador é zerado.

Esse critério é diferente de `step_x <= tol`: `tol` é uma exigência do problema, enquanto `epsilon * scale` descreve o limite de representação do computador. Separar os dois evita declarar sucesso ou estagnação apenas porque um passo foi pequeno em relação à tolerância escolhida.

## Ordem de convergência

Para resultados com pelo menos três erros válidos em `x`, `estimate_convergence_order` estima a ordem local usando os últimos três valores:

```text
p ~= log(e_(k+1) / e_k) / log(e_k / e_(k-1))
```

O valor é usado pelo `reporter.cpp` para exibir relatórios e rankings. Caso não existam dados suficientes, haja divisão/logaritmo inválido ou o resultado não seja finito, a ordem é apresentada como `N/A`.
