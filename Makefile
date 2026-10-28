CC      ?= gcc
CFLAGS  ?= -O2 -g -Wall -Wextra -Werror -Isrc
LDLIBS  ?= -lm

# make OPENMP=1: paraleliza o produto matricial (e as operações por coluna)
# com OpenMP. Não muda o comportamento do build padrão.
ifeq ($(OPENMP),1)
  CFLAGS += -fopenmp -DLLA_OPENMP
  LDLIBS += -lm -lpthread
endif

LIB_SRC = src/vector.c src/matrix.c src/lu.c src/svd.c src/lsystem.c
LIB_NAME = lla

OBJ_DIR  = obj
LIB_OBJS = $(patsubst src/%.c,$(OBJ_DIR)/%.o,$(LIB_SRC))

.PHONY: all build libs static shared bench lim run bench-run lim-run clean

all: build libs bench lim

# objetos com -fPIC, usados pelos dois formatos de biblioteca
obj:
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: src/%.c | obj
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

build: $(LIB_SRC) src/main.c
	$(CC) $(CFLAGS) $^ -o al $(LDLIBS)

libs: static shared
static: $(LIB_OBJS)
	ar rcs lib$(LIB_NAME).a $(LIB_OBJS)

shared: $(LIB_OBJS)
	$(CC) -shared $(LIB_OBJS) -o lib$(LIB_NAME).so $(LDLIBS)

bench: $(LIB_SRC) src/bech.c
	$(CC) $(CFLAGS) $^ -o bech $(LDLIBS)

lim: $(LIB_SRC) src/lim.c
	$(CC) $(CFLAGS) $^ -o lim $(LDLIBS)

run: build
	./al

bench-run: bench
	./bech

lim-run: lim
	./lim

clean:
	rm -f al bech lim lib$(LIB_NAME).a lib$(LIB_NAME).so
	rm -rf al.dSYM bech.dSYM lim.dSYM $(OBJ_DIR)
