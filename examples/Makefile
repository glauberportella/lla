CC      ?= gcc
CFLAGS  ?= -Og -g -Wall -Wextra -Werror -I../src
LDLIBS  ?= -lm

EXAMPLES = 01_vetores 02_matrizes 03_determinante_inversa \
           04_decomposicao_lu 05_sistema_linear 06_minimos_quadrados 07_svd

LIB_SRC = ../src/vector.c ../src/matrix.c ../src/lu.c ../src/svd.c ../src/lsystem.c

.PHONY: all run clean $(EXAMPLES)

all: $(EXAMPLES)

01_vetores: 01_vetores.c $(LIB_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

02_matrizes: 02_matrizes.c $(LIB_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

03_determinante_inversa: 03_determinante_inversa.c $(LIB_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

04_decomposicao_lu: 04_decomposicao_lu.c $(LIB_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

05_sistema_linear: 05_sistema_linear.c $(LIB_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

06_minimos_quadrados: 06_minimos_quadrados.c $(LIB_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

07_svd: 07_svd.c $(LIB_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

run: all
	@for ex in $(EXAMPLES); do \
	  echo "### executando $$ex ###"; \
	  ./$$ex || true; \
	  echo; \
	done

clean:
	rm -f $(EXAMPLES)
