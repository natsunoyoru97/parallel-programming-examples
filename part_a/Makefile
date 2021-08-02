CC = mpicc
RUN = mpirun
BINARIES = ./pairwise_inner_product_parallel

NUMPROCS = 25
N = 100
M = 30

all: main run

main: pairwise_inner_product_parallel.c
	$(CC) -Wall -Werror pairwise_inner_product_parallel.c -o pairwise_inner_product_parallel

run:
	$(RUN) -np $(NUMPROCS) $(BINARIES) $(N) $(M)

clean:
	rm -f $(BINARIES) *.o