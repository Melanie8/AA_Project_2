CC=g++
CFLAGS=-Werror
EXEC=fact
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(EXEC)

fact: $(OBJ)
	    @$(CC) -lgmp -o $@ $^

search: $(OBJ)
	    @$(CC) -lgmp -o $@ $^

%.o: %.c
	@$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	@rm -rf *.o

mrproper: clean
	@rm -rf $(EXEC)