# Compilador e flags de compilação
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -I.

# Arquivos base que são compartilhados entre o main e os testes
CORE_SRC = dual.cpp solvers.cpp reporter.cpp

# Nomes dos executáveis que serão gerados
APP = main_app
TEST_APP = run_test_dual
TEST_MOODLE = moodle_test

# Alvo padrão (roda quando você digita apenas 'make')
all: $(APP)

# Regra para compilar o programa principal
$(APP): main.cpp $(CORE_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Regra para compilar os testes da pasta test_solver
test: test_dual/test_dual.cpp $(CORE_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $(TEST_APP)
	@echo "Testes compilados! Execute com: ./$(TEST_APP)"

moodle: test_moodle/test_moodle.cpp $(CORE_SRC)
	$(CXX) $(CXXFLAGS) $^ -o $(TEST_MOODLE)
	@echo "Testes compilados! Execute com: ./$(TEST_MOODLE)"

# Regra para limpar os executáveis gerados
clean:
	rm -f $(APP) $(TEST_APP) $(TEST_MOODLE)