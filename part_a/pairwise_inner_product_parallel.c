#include "parallel_header.h"

#define pp(process) ((process > myid) ? (process - 1) : (process))

int main(int argc, char** argv) {
    assert(argc == 3);

    int numprocs, myid;
    int previd, nextid;
    MPI_Status* statuses = NULL;
    MPI_Request* requests = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int N = (int)strtol(argv[1], (char**)NULL, 10);
    int M = (int)strtol(argv[2], (char**)NULL, 10);
    printf("N: %d, M:%d, numprocs:%d\n", N, M, numprocs);

    assert(N > 0 && M > 0);
    // numprocs must be odd to acheive best load balancing
    assert(numprocs % 2 == 1);
    assert(N % numprocs == 0);

    int is_proc_zero, is_parallel;
    int iter_curr, iter_cnt;
    int BLOCKROWS, BLOCKCOLS;
    int i, j;
    int cnt, steps, total;

    double** full_matrix = NULL;
    double **block_curr = NULL, **block_to_comp = NULL;
    double** block_tmp = NULL;
    double *p_arr = NULL, *s_arr = NULL;

    is_proc_zero = (myid == 0);
    is_parallel = (numprocs > 1);
    previd = is_proc_zero ? (numprocs - 1) : (myid - 1);
    nextid = (myid + 1) % numprocs;

    BLOCKROWS = N / numprocs;
    iter_cnt = (numprocs + 1) / 2;
    steps = (N % 2 == 0) ? ((N + 2) / 2) : ((N + 1) / 2);

    total = (N * 2 - N + 1) * N / 2;

    statuses = (MPI_Status*)malloc((numprocs) * sizeof(MPI_Status));
    requests = (MPI_Request*)malloc((numprocs) * sizeof(MPI_Request));

    for (i = 0; i < numprocs; i++) {
        requests[i] = MPI_REQUEST_NULL;
    }

    block_curr = init_matrix(BLOCKROWS, M);
    block_to_comp = init_matrix(BLOCKROWS, M);

    p_arr = init_res_arr(total);

    int* index_arr = init_index_arr(N);

    /*
     * Process 0 allocate blocks to other processes
     * and process 0 also PERFORMS COMPUTATION
     */
    if (is_proc_zero) {
        full_matrix = init_matrix(N, M);
        fill_matrix(full_matrix, N, M);
        print_matrix(full_matrix, N, M);

        if (is_parallel) {
            // Note: result = block_curr * block_to_compute (persudocode)

            // Send allcated blocks to each process
            memcpy(block_curr[0], full_matrix[0],
                   BLOCKROWS * M * sizeof(double));
            for (i = 0; i < numprocs; i++) {
                MPI_Isend(full_matrix[i * BLOCKROWS], BLOCKROWS * M, MPI_DOUBLE,
                          i, 0, MPI_COMM_WORLD, &requests[pp(i)]);
            }
            MPI_Waitall(numprocs, requests, statuses);

            // Send blocks to compute to each process
            memcpy(block_curr[0], full_matrix[BLOCKROWS],
                   BLOCKROWS * M * sizeof(double));
            for (i = 0; i < numprocs; i++) {
                j = i % numprocs;
                MPI_Isend(full_matrix[j * BLOCKROWS], BLOCKROWS * M, MPI_DOUBLE,
                          i, 0, MPI_COMM_WORLD, &requests[pp(i)]);
            }
        }
    }

    if (is_parallel) {
        // Other processes receive allocated blocks including the process 0
        MPI_Recv(block_curr[0], BLOCKROWS * M, MPI_DOUBLE, 0, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &statuses[pp(0)]);
        MPI_Recv(block_to_comp[0], BLOCKROWS * M, MPI_DOUBLE, 0, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &statuses[pp(0)]);

        // Perform vectors computation
        //printf("Process %d is performing parallel computation.\n", myid);

        BLOCKCOLS = iter_cnt * BLOCKROWS;
        block_tmp = init_matrix(BLOCKROWS, steps);
        fill_matrix_with_zero(block_tmp, BLOCKROWS, steps);

        comp_matrix(block_curr, BLOCKROWS, M, block_tmp);

        int row_curr;

        iter_curr = 0;
        cnt = 0;

        for (;;) {
            row_curr = BLOCKROWS;
            for (i = 0; i < BLOCKROWS; i++, row_curr--) {
                // At this time block_curr == block_comp, to avoid duplicate
                // computation
                if (iter_curr == 0) break;

                for (j = 0; j < BLOCKROWS; j++) {
                    int bound = row_curr + (iter_curr - 1) * BLOCKROWS;
                    if (bound + j < steps) {
                        block_tmp[i][bound + j] =
                            inner_product(block_curr[i], block_to_comp[j], M);
                    }
                }
            }

            iter_curr++;
            if (iter_curr > iter_cnt) break;

            // Shift right block
            MPI_Wait(&requests[pp(previd)], &statuses[pp(previd)]);
            MPI_Isend(block_to_comp[0], BLOCKROWS * M, MPI_DOUBLE, previd, 0,
                      MPI_COMM_WORLD, &requests[pp(previd)]);
            MPI_Recv(block_to_comp[0], BLOCKROWS * M, MPI_DOUBLE, nextid,
                     MPI_ANY_TAG, MPI_COMM_WORLD, &statuses[pp(nextid)]);
        }

        // Merge the compuation result
        if (is_proc_zero) {
            iter_curr = 0;

            int row_curr = N;
            int index_curr = 0;

            for (;;) {
                for (i = 0; i < BLOCKROWS; i++) {
                    for (j = 0; j < steps; j++) {
                        if (j < row_curr) {
                            p_arr[j + index_curr] = block_tmp[i][j];
                        } else if (j == row_curr) {
                            p_arr[iter_curr * BLOCKROWS + i] = block_tmp[i][j];
                        } else if (j > row_curr) {
                            p_arr[index_arr[j - row_curr - 1] - (j - row_curr) +
                                  iter_curr * BLOCKROWS + i] = block_tmp[i][j];
                        }
                    }
                    row_curr--;
                    index_curr = index_arr[iter_curr * BLOCKROWS + i];
                }

                iter_curr++;

                if (iter_curr >= numprocs) break;

                MPI_Recv(block_tmp[0], BLOCKROWS * steps, MPI_DOUBLE, iter_curr,
                         MPI_ANY_TAG, MPI_COMM_WORLD, &statuses[pp(i)]);
            }
        }
        // Other processes send results to the process 0
        else {
            MPI_Wait(&requests[pp(0)], &statuses[pp(0)]);
            MPI_Send(block_tmp[0], BLOCKROWS * steps, MPI_DOUBLE, 0, 0,
                     MPI_COMM_WORLD);

            //printf("Send blocks from proc %d\n", myid);
        }
    }

    if (is_proc_zero) {
        // Compute results of sequential computing
        s_arr = init_res_arr(total);
        comp_arr(full_matrix, N, M, s_arr, index_arr);
        printf("The sequential computation result is:\n");
        print_upper_tri_matrix(s_arr, index_arr, total, N);

        // Compare results of parallel computing and results of sequential
        // computing
        if (is_parallel) {
            printf("\nThe parallel computation result is:\n");
            print_upper_tri_matrix(p_arr, index_arr, total, N);

            int diff;
            diff = 0;
            for (i = 0; i < total; i++) {
                if (p_arr[i] != s_arr[i]) {
                    diff++;
                    printf("p_arr[%d]: %f is different from s_arr[%d]: %f\n", i,
                           p_arr[i], i, s_arr[i]);
                }
            }

            if (diff == 0) {
                printf("p_arr is same with s_arr.\n");
            }
        } else {
            printf(
                "There is only 1 process, so sequential computing has "
                "performed.\n");
        }
    }
    
    free(block_curr);
    free(block_to_comp);
    free(block_tmp);

    free(p_arr);
    free(s_arr);

    MPI_Finalize();
    return 0;
}