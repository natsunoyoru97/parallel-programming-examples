#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum bool { false = 0, true = 1 };
enum cells { white = 0, red = 1, blue = 2, in = 3, out = 4};

struct thread_data {
    int thread_id;
    int row_upper_bound;
    int rows, cols;
    int t_start, t_end;
};

typedef struct {
    pthread_mutex_t count_lock;
    pthread_cond_t ok_to_proceed;
    int count;
} mylib_barrier_t;

/* function declaration */
// barrier function declaration
void mylib_barrier_init(mylib_barrier_t *b);
void mylib_barrier_destroy(mylib_barrier_t *b);
void mylib_barrier (mylib_barrier_t *b, int num_threads);

int **init_block(int row, int col);
void cells_init(int **grid, int n, int init_num);
void board_print(int **grid, int row, int col);

void red_cell_movement(int **grid, int row_upper_bound, int rows, int cols);
void blue_cell_movement(int **grid, int row_upper_bound, int rows, int cols, int n);

int cal_num_colored(int **grid, int t_row_itr, int t_col_itr, int len,
                        int color);
int find_colored_tiles(int **grid, int row_start, int row_end, int col_start,
                           int col_end, int tile_len, int c);

/* function initialzation */
void mylib_barrier_init(mylib_barrier_t *b) {
    b->count = 0;
    pthread_mutex_init(&(b->count_lock), NULL);
    pthread_cond_init(&(b->ok_to_proceed), NULL);
}

void mylib_barrier_destroy(mylib_barrier_t *b) {
    pthread_mutex_destroy(&(b->count_lock));
    pthread_cond_destroy(&(b->ok_to_proceed));
}

void mylib_barrier_wait(mylib_barrier_t *b, int num_threads) {
    pthread_mutex_lock(&(b->count_lock));
    b->count += 1;
    
    if ((b->count) == num_threads) {
        b->count = 0;
        pthread_cond_broadcast(&(b->ok_to_proceed));
    }
    else {
        pthread_cond_wait(&(b->ok_to_proceed), &(b->count_lock));
    }
    pthread_mutex_unlock(&(b->count_lock));
}