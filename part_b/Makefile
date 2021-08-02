CC = cc
BINARIES = ./red_blue_movement

THREADS = 2
N = 24
T = 12
C = 40
MAXITERS = 100

all: main run

main: red_blue_movement.c
	$(CC) -Wall -Werror red_blue_movement.c -lpthread -o red_blue_movement

run:
	$(BINARIES) $(THREADS) $(N) $(T) $(C) $(MAXITERS)

clean:
	rm -f $(BINARIES) *.o