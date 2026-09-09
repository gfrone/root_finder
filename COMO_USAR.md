# Como usar o projeto

## Requisitos

É necessário ter `g++` com suporte a C++17 e `make` instalados.

## Compilar e executar a aplicação principal

Na raiz do projeto, execute:

```bash
make
./main_app
```

O código da aplicação principal fica em `main.cpp`.

## Compilar e executar os testes de números duais

```bash
make test
./run_test_dual
```

## Compilar e executar os exercícios do Moodle

```bash
make moodle
./moodle_test
```

Os exercícios estão em `test_moodle/test_moodle.cpp`.

No fim desse arquivo, no `main`, descomente a chamada do exercício que deseja executar. Por exemplo:

```cpp
int main() {
    exercicio_11(1e-4);
    exercicio_12(1e-4);
    return 0;
}
```

As tolerâncias usadas nos exercícios são passadas como argumento, como `1e-3` ou `1e-4`.

## Limpar executáveis compilados

```bash
make clean
```

Esse comando remove `main_app`, `run_test_dual` e `moodle_test`. Os arquivos-fonte não são removidos.
