#include <png.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <pthread.h>

#define RECINIT 0.1f
#define RMFIRST 200
#define ITLIMIT 250
#define ABS(X) ((X) >= 0.0 ? (X) : -(X))
#define MOD1(X) ((X) - floor(X))
#define A_RES 10.0

#define N_THREADS 8
/* A coloured pixel. */

typedef struct
{
    uint8_t red;
    uint8_t green;
    uint8_t blue;
}
pixel_t;

/* A picture. */
    
typedef struct
{
    pixel_t *pixels;
    size_t width;
    size_t height;
}
bitmap_t;
    
typedef struct
{
    unsigned *data;
    size_t rows; // Range of Omega
    size_t cols; // Range of K
    double o_low;
    double k_low;
    double o_stepsz;
    double k_stepsz;
    double eps;
} plot_t;

typedef struct
{
    plot_t* plot;
    size_t offset;
    size_t n;
} args_t;

static unsigned * point_at (plot_t * plot, size_t r, size_t c) 
{
    return plot->data + plot->rows * c + r;
}    


/* Given "value" and "max", the maximum value which we expect "value"
   to take, this returns an integer between 0 and 255 proportional to
   "value" divided by "max". */

static int pix (int value, int max)
{
    if (value < 0) {
        return 0;
    }
    return (int) (256.0 *((double) (value)/(double) max));
}

static inline double circle_map(double O, double K, double theta) {
    return MOD1(theta + O - K / (2*M_PI) * sin(2*M_PI * theta));
}

unsigned get_recn(double O, double K, double eps, double init) {
    /*Ignore first few iterations to get to convergence */
    double theta_0 = init;
    for (int i = 0; i < RMFIRST; i++)
        theta_0 = circle_map(O, K, theta_0);
    /*calculate recurrence number*/
    double theta_1 = circle_map(O,K, theta_0);
    unsigned n = 0;
    for (; n < ITLIMIT && ABS(theta_1 - theta_0) >= eps; n++)
        theta_1 = circle_map(O, K, theta_1);
    return n;
}

unsigned get_avg_recn(double O, double K, double eps) {
    /*unsigned acc = 0;
    double stepsz = 1.0/A_RES;
    for (double init = 0.0; init < 1.0; init += stepsz) {
        acc += get_recn(O, K, eps, init);
    }
    return acc / A_RES;*/
    //unsigned acc = 0;
    return get_recn(O, K, eps, 0.1);
    /*acc += get_recn(O, K, eps, 0.9);
    acc += get_recn(O, K, eps, 0.5);*/
    //return acc / 1; 
}

pixel_t color_from_recn(int n) {
    pixel_t colo;
    if (n < 10) {

    } else if (n < 50) {
        colo.blue = 255;
    } else if (n < 140) {
        colo.green = 255;
    } else {
        colo.red = 255;
    }
    return colo;
}

void * color_segment(void* args) {
    args_t *a = (args_t *) args;
    plot_t plot = *(a->plot);
    size_t offset = a->offset;
    size_t n = a->n;
    
    double o = plot.o_low + offset * plot.o_stepsz;
    for (size_t oidx = offset; oidx < offset + n; oidx++, o += plot.o_stepsz) {
        double k = plot.k_low;
        for (size_t kidx = 0; kidx < plot.rows; kidx++, k += plot.k_stepsz) {
            unsigned n = get_recn(o, k, plot.eps, 0.5);
            unsigned *pt = point_at (&plot, kidx, oidx);
            *pt = n;
        }
    }
    return NULL;
}

void save_to_csv(plot_t plot, const char *path) {
    FILE * fp = fopen (path, "w");
    if (!fp) {
        printf("Error: Couldn't open \"%s\"", path);
    }
    for (size_t r = 0; r < plot.rows; r++) {
        size_t c;
        for (c = 0; c < plot.cols - 1; c++) {
            fprintf(fp,"%d,", plot.data[c * plot.rows + r]);
        }
        fprintf(fp, "%d\n", plot.data[c * plot.rows + r]);
    }
    fclose(fp);
}

void mk_plot(size_t k_steps, size_t o_steps, double k_low, double k_high, double o_low, double o_high, double eps) {
    /* Set up bitmap and calculate parameter steps*/
    plot_t plot;
    plot.cols = o_steps;
    plot.rows = k_steps;
    plot.data = calloc(o_steps * k_steps, sizeof(unsigned));
    plot.o_low = o_low;
    plot.k_low = k_low;
    plot.k_stepsz = (k_high - k_low) / k_steps;
    plot.o_stepsz = (o_high - o_low) / o_steps;
    plot.eps = eps;
    if (!plot.data)
        exit(1);
    
    /* Parallelize Workload */
    pthread_t tids[N_THREADS];
    args_t argss[N_THREADS];
    size_t batchsz = o_steps / N_THREADS;
    size_t final_n = o_steps - (N_THREADS - 1) * batchsz;
    for (size_t i = 0; i < N_THREADS; i++) {
        argss[i].plot = &plot;
        argss[i].offset = i * batchsz;
        argss[i].n = batchsz;
    }
    argss[N_THREADS-1].n = final_n;
    printf("Dispatching to %d threads, batch size: %zu\n", N_THREADS, batchsz);
    for (size_t i = 0; i < N_THREADS; i++) {
        pthread_create(tids + i, NULL, color_segment, (void *)(argss + i));
    }
    for (size_t i = 0; i < N_THREADS; i++) {
        pthread_join(tids[i], NULL);
    }

    /* Save File */
    char fname[100];
    snprintf(fname, 100, "output/arnold_%d.csv", (unsigned)time(NULL));
    save_to_csv(plot, fname);
    free(plot.data);
}

int main ()
{
    mk_plot(3000, 1500, 0.0f, 4*M_PI, 0.0f, 1.0f, 0.0001f);
  
}