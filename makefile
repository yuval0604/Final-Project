CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors

all: symnmf

symnmf: main.o symnmf.o
	$(CC) $(CFLAGS) main.o symnmf.o -o symnmf -lm

main.o: main.c symnmf.h
	$(CC) $(CFLAGS) -c main.c

symnmf.o: symnmf.c symnmf.h
	$(CC) $(CFLAGS) -c symnmf.c

clean:
	rm -f *.o symnmf
