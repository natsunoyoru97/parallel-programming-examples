#include "red_blue_movement.h"

int **full_grid;  // the full grid to allocate/compute
int **s_grid;     // the grid runs sequential computation
int **p_grid;     // the grid runs parallel computation

int NUM_THREADS;
int n, t, c, MAX_ITERS;
int rows_each_block;
int finished = 0;

mylib_barrier_t barrier;

int rc, task, len;

void* Compute(void *threadarg);

int main(int argc, char** argv) {
    NUM_THREADS = (int)strtol(argv[1], (char **)NULL, 10); // number of threads
    n = (int)strtol(argv[2], (char **)NULL, 10); // cell grid size
    t = (int)strtol(argv[3], (char **)NULL,
                        10); // tile grid size, size of group
    c = (int)strtol(argv[4], (char **)NULL,
                        10); // terminating threshold, e.g. 40 stands for 40%
    MAX_ITERS = (int)strtol(argv[5], (char **)NULL,
                                10); // maximum number of iterations

    assert(n > 0 && t > 0);
    assert(n % t == 0); // the board should be prefectly overlaid with tiles
    assert(NUM_THREADS > 0 && NUM_THREADS <= t);
    assert(c > 0 && MAX_ITERS > 0);

    rows_each_block = n / t;
    int d = (NUM_THREADS > 1) ? (t / NUM_THREADS) : 0;
    int r = (NUM_THREADS > 1) ? (t % NUM_THREADS) : 0;
    int d_ex = d + 1;

    int n_itrs;
    int is_parallel = (NUM_THREADS > 1);

    struct thread_data thread_data_arr[NUM_THREADS];
    pthread_t threads[NUM_THREADS];

    len = 0;

    // Initialize grids
    s_grid = init_block(n, n);
    cells_init(s_grid, n, n * n / 3);
    p_grid = init_block(n, n);
    memcpy(&p_grid[0][0], &s_grid[0][0], n * n * sizeof(int));

    board_print(p_grid, n, n);

    if (is_parallel) {
        int cnt, index, tile_index;
        cnt = r;
        index = 0;
        tile_index = 0;

        mylib_barrier_init(&barrier);

        for (task = 0; task < NUM_THREADS; task++) {
            printf("In main: creating thread %d\n", task);
            thread_data_arr[task].thread_id = task;
            thread_data_arr[task].row_upper_bound = index;
            thread_data_arr[task].cols = n;

            // Assign attributes for each thread
            if (cnt > 0) {
                thread_data_arr[task].t_start = tile_index;
                tile_index += d_ex;
                thread_data_arr[task].t_end = tile_index;

                thread_data_arr[task].rows = d_ex * rows_each_block;
                index += rows_each_block * d_ex;

                cnt--;
            }
            else {
                thread_data_arr[task].t_start = tile_index;
                tile_index += d;
                thread_data_arr[task].t_end = tile_index;

                thread_data_arr[task].rows = d * rows_each_block;
                index += rows_each_block * d;
            }
        }

        for (task = 0; task < NUM_THREADS; task++) {
            // Create threads
            rc = pthread_create(&threads[task], NULL, Compute, (void *)&thread_data_arr[task]);
            if (rc) {
                printf("ERROR; return code from pthread_create() is %d\n", rc);
                exit(-1);
            }
        }
        
        for (task = 0; task < NUM_THREADS; task++) {
            // Specify the sequence to execute the tasks
            rc = pthread_join(threads[task], NULL);
            if (rc) {
                printf("ERROR; return code from pthread_join() is %d\n", rc);
                exit(-1);
            }
        }

        board_print(p_grid, n, n);
    }

    // Perform sequential computing
    finished = false;
    for (n_itrs = 0; !finished && n_itrs < MAX_ITERS; n_itrs++) {
        red_cell_movement(s_grid, 0, n, n);
        blue_cell_movement(s_grid, 0, n, n, n);

        finished = find_colored_tiles(s_grid, 0, t, 0, t, rows_each_block, c);
    }

    printf("s_grid:\n");
    board_print(s_grid, n, n);

    if (is_parallel) {
        int i, j;
        int cnt;
        cnt = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (p_grid[i][j] != s_grid[i][j]) {
                    printf(
                        "p_grid[%d][%d]: %d is different from s_grid[%d][%d]: "
                        "%d\n",
                        i, j, p_grid[i][j], i, j, s_grid[i][j]);
                    cnt++;
                }
            }
        }
        if (cnt == 0) printf("p_grid is same with s_grid.\n");
        mylib_barrier_destroy(&barrier);
    }
    
    pthread_exit(NULL);
}

void *Compute(void *threadarg) {
    int row_upper_bound;
    int rows, cols;
    int t_start, t_end;
    int finished_p;
    struct thread_data *my_data;

    my_data = (struct thread_data *)threadarg;
    row_upper_bound = my_data->row_upper_bound;
    rows = my_data->rows;
    cols = my_data->cols;
    t_start = my_data->t_start;
    t_end = my_data->t_end;

    finished_p = false;

    // Perform computation
    int n_itrs;
    for (n_itrs = 0; !finished && n_itrs < MAX_ITERS; n_itrs++) {
        red_cell_movement(p_grid, row_upper_bound, rows, cols);
        mylib_barrier_wait(&barrier, NUM_THREADS);
        blue_cell_movement(p_grid, row_upper_bound, rows, cols, n);
        mylib_barrier_wait(&barrier, NUM_THREADS);
        finished_p =
            find_colored_tiles(p_grid, t_start, t_end, t_start, t_end,
                               rows_each_block, c);
        if (finished_p) {
            finished = finished_p;
        }
        mylib_barrier_wait(&barrier, NUM_THREADS);
    }

    pthread_exit(NULL);
}

// Initialize all other functions
int **init_block(int row, int col) {
    int *data = (int *)malloc(row * col * sizeof(int));
    int **block = (int **)malloc(row * sizeof(int *));

    int i, j;
    for (i = 0; i < row; i++) {
        block[i] = &(data[col * i]);
    }

    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            block[i][j] = 0;
        }
    }

    return block;
}

void cells_init(int **grid, int n, int init_num) {
    /* init red cells */
    int i;
    for (i = 0; i < init_num; i++) {
        int num = rand() % (n * n) - 1;
        int col = num % n, row = num / n;
        while (grid[row][col] != 0) {
            num = rand() % (n * n) - 1;
            col = num % n, row = num / n;
        }
        grid[row][col] = red;
    }

    /* init blue cells */
    for (i = 0; i < init_num; i++) {
        int num = rand() % (n * n) - 1;
        int col = num % n, row = num / n;
        while (grid[row][col] != 0) {
            num = rand() % (n * n) - 1;
            col = num % n, row = num / n;
        }
        grid[row][col] = blue;
    }
}

void board_print(int **grid, int row, int col) {
    int i, j;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%d ", grid[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/* movement functions */
void red_cell_movement(int **grid, int row_upper_bound, int rows, int cols) {
    int i, j;
    for (i = row_upper_bound; i < rows + row_upper_bound; i++) {
        if (grid[i][0] == red && grid[i][1] == white) {
            grid[i][0] = out;
            grid[i][1] = in;
        }
        for (j = 1; j < cols; j++) {
            if (grid[i][j] == red && grid[i][(j + 1) % cols] == white) {
                grid[i][j] = white;
                grid[i][(j + 1) % cols] = in;
            } else if (grid[i][j] == in) {
                grid[i][j] = red;
            }
        }
        if (grid[i][0] == in) {
            grid[i][0] = red;
        } else if (grid[i][0] == out) {
            grid[i][0] = white;
        }
    }
}

void blue_cell_movement(int **grid, int row_upper_bound, int rows, int cols,
                        int n) {
    int i, j;
    for (j = row_upper_bound; j < row_upper_bound + rows; j++) {
        if (grid[0][j] == blue && grid[1][j] == white) {
            grid[0][j] = out;
            grid[1][j] = in;
        }
        for (i = 1; i < n; i++) {
            if (grid[i][j] == blue && grid[(i + 1) % n][j] == white) {
                grid[i][j] = white;
                grid[(i + 1) % n][j] = in;
            } else if (grid[i][j] == in) {
                grid[i][j] = blue;
            } 
        }
        // border check
        if (grid[0][j] == in) {
            grid[0][j] = blue;
        } else if (grid[0][j] == out) {
            grid[0][j] = white;
        }
    }
}

int cal_num_colored(int **grid, int t_row_itr, int t_col_itr, int len,
                        int color) {
    int cnt = 0;

    if (len == 1 && grid[1][0] == color) {
        cnt += 1;
    }

    int i, j;
    for (i = t_row_itr; i < len; i++) {
        for (j = t_col_itr; j < len; j++) {
            if (grid[i][j] == color) {
                cnt += 1;
            }
        }
    }

    return cnt;
}

int find_colored_tiles(int **grid, int row_start, int row_end, int col_start,
                       int col_end, int tile_len, int c) {
    int finished = false;
    int max_cells = tile_len * tile_len * (float)(c * 0.01);
    int i, j;

    for (i = row_start; i < row_end; i += tile_len) {
        for (j = col_start; j < col_end; j += tile_len) {
            int redcount = cal_num_colored(grid, i, j, tile_len, red);
            int bluecount = cal_num_colored(grid, i, j, tile_len, blue);

            if (redcount > max_cells || bluecount > max_cells) {
                finished = true;
                printf(
                    "The tile from row %d, col %d has more than %d%% one color, the iteration is terminated.\n",
                    i, j, c);
                return finished;
            }
        }
    }

    return finished;
}