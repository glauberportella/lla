CC      ?= gcc
CFLAGS  ?= -g -Wall -Wextra -Werror -Isrc
LDLIBS  ?= -lm

LIB_SRC = src/vector.c src/matrix.c src/lu.c src/svd.c

.PHONY: all build bench lim run bench-run lim-run clean

all: build bench lim

build: $(LIB_SRC) src/main.c
	$(CC) $(CFLAGS) $^ -o al $(LDLIBS)

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
	rm -f al bech lim
	rm -rf al.dSYM bech.dSYM lim.dSYM
