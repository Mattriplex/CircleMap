#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <pthread.h>
#include <gmp.h>
#include <mpfr.h>

#define RECINIT 0.1f
#define RMFIRST 200
#define ITLIMIT 10000
#define ABS(X) ((X) >= 0.0 ? (X) : -(X))
#define MOD1(X) ((X) - floor(X))
#define A_RES 10.0

#define DEF_OL 0.0f
#define DEF_OH 0.5f
#define DEF_KL 0.0f
#define DEF_KH (4*M_PI)

#define N_THREADS 8

#define MULTIPREC
#ifdef MULTIPREC
#define PRECISION 500
mpfr_t TAU;
mpfr_t EPS;
#endif

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
    size_t convcnt;
    size_t n;
} args_t;

static unsigned * point_at (plot_t * plot, size_t r, size_t c) 
{
    return plot->data + plot->rows * c + r;
}    


/* theta_n+1 = theta_n + O - K/(2*PI)*sin(2*pi*theta_n)*/
static inline double circle_map(double O, double K, double theta) {
    return MOD1(theta + O - K / (2*M_PI) * sin(2*M_PI * theta));
}

#ifdef MULTIPREC
static inline void circle_map_mp(double o, mpfr_t kot, mpfr_t theta, mpfr_t tmp) {
    mpfr_add_d(theta, theta, o, MPFR_RNDN);
    mpfr_mul(tmp, theta, TAU, MPFR_RNDN);
    mpfr_sin(tmp, tmp, MPFR_RNDN);
    mpfr_mul(tmp, tmp, kot, MPFR_RNDN);
    mpfr_sub(theta, theta, tmp, MPFR_RNDN);
    mpfr_floor(tmp, theta);
    mpfr_sub(theta, theta, tmp, MPFR_RNDN);
}
#endif

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

#ifdef MULTIPREC
unsigned get_recn_mp(double o, mpfr_t kot, mpfr_t theta_0, mpfr_t theta_1, mpfr_t tmp, double init) {
    /*Ignore first few iterations to get to convergence */
    mpfr_set_d(theta_0, init, MPFR_RNDN);
    for (int i = 0; i < RMFIRST; i++)
        circle_map_mp(o, kot, theta_0, tmp);
    /*calculate recurrence number*/
    mpfr_set(theta_1, theta_0, MPFR_RNDN);
    circle_map_mp(o, kot, theta_1, tmp);
    for (unsigned n = 0; n < ITLIMIT; n++) {
        mpfr_sub(tmp, theta_0, theta_1, MPFR_RNDN);
        mpfr_abs(tmp, tmp, MPFR_RNDN);
        if (mpfr_cmp(tmp, EPS) < 0)
            return n;
        circle_map_mp(o, kot, theta_1, tmp);
    }
    return ITLIMIT;
}
#endif

unsigned get_avg_recn(double O, double K, double eps) {
    /*unsigned acc = 0;
    double stepsz = 1.0/A_RES;
    for (double init = 0.0; init < 1.0; init += stepsz) {
        acc += get_recn(O, K, eps, init);
    }
    return acc / A_RES;*/
    return get_recn(O,K,eps, 0.5);
}

#ifdef MULTIPREC
void * color_segment_mp(void* args) {
    args_t *a = (args_t *) args;
    plot_t plot = *(a->plot);
    size_t offset = a->offset;
    size_t n = a->n;
    size_t pct = n / 10;
    size_t pctcnt = 0;
    mpfr_t theta_0, theta_1, kmp, tmp, delta_k_ot, k_low_ot;
    mpfr_inits2(PRECISION, theta_0, theta_1, kmp, tmp, delta_k_ot, k_low_ot, NULL);
    mpfr_d_div(delta_k_ot, plot.k_stepsz, TAU, MPFR_RNDN);
    mpfr_d_div(k_low_ot, plot.k_low, TAU, MPFR_RNDN);
    double o = plot.o_low + offset * plot.o_stepsz;
    for (size_t oidx = offset; oidx < offset + n; oidx++) {
        mpfr_set(kmp, k_low_ot, MPFR_RNDN);
        for (size_t kidx = 0; kidx < plot.rows; kidx++) {
            unsigned n = get_recn_mp(o, kmp, theta_0, theta_1, tmp, 0.5);
            if (n != ITLIMIT)
                a->convcnt++;
            unsigned *pt = point_at (&plot, kidx, oidx);
            *pt = n;
            mpfr_add(kmp, kmp, delta_k_ot, MPFR_RNDN);

        }
        o += plot.o_stepsz;
        if ((oidx - offset) * 100 / n > pctcnt) {
            pctcnt = (oidx - offset) * 100 / n;
            printf("%zu%% completed\n", pctcnt);
        }
    }
    printf("segment complete.\n");
    mpfr_clears(theta_0, theta_1, kmp, tmp, delta_k_ot, k_low_ot, NULL);
    return NULL;
}
#endif

void * color_segment(void* args) {
    args_t *a = (args_t *) args;
    plot_t plot = *(a->plot);
    size_t offset = a->offset;
    size_t n = a->n;
    
    double o = plot.o_low + offset * plot.o_stepsz;
    for (size_t oidx = offset; oidx < offset + n; oidx++) {

        double k = plot.k_low;
        for (size_t kidx = 0; kidx < plot.rows; kidx++) {
            unsigned n = get_avg_recn(o, k, plot.eps); 
            if (n != ITLIMIT)
                a->convcnt++;
            unsigned *pt = point_at (&plot, kidx, oidx);
            *pt = n;
            k += plot.k_stepsz;
        }
        o += plot.o_stepsz;
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
    
#ifdef MULTIPREC
    /*Calculate 2*PI*/
    mpfr_inits2(PRECISION, TAU, EPS, NULL);
    mpfr_set_si(TAU, -1, MPFR_RNDN);
    mpfr_acos(TAU, TAU, MPFR_RNDN);
    mpfr_add(TAU, TAU, TAU, MPFR_RNDN);
    mpfr_set_d(EPS, eps, MPFR_RNDN);
#endif
    /* Parallelize Workload */
    pthread_t tids[N_THREADS];
    args_t argss[N_THREADS];
    size_t batchsz = o_steps / N_THREADS;
    size_t final_n = o_steps - (N_THREADS - 1) * batchsz;
    for (size_t i = 0; i < N_THREADS; i++) {
        argss[i].plot = &plot;
        argss[i].offset = i * batchsz;
        argss[i].n = batchsz;
        argss[i].convcnt = 0;
    }
    argss[N_THREADS-1].n = final_n;
    printf("Dispatching to %d threads, batch size: %zu\n", N_THREADS, batchsz);
    for (size_t i = 0; i < N_THREADS; i++) {
#ifdef MULTIPREC
        pthread_create(tids + i, NULL, color_segment_mp, (void *)(argss + i));
#else
        pthread_create(tids + i, NULL, color_segment, (void *)(argss + i));
#endif
    }
    size_t convcnt = 0;
    for (size_t i = 0; i < N_THREADS; i++) {
        pthread_join(tids[i], NULL);
        convcnt += argss[i].convcnt;
    }
    /* Save File */
    printf("Convergence rate: %3f\n", (float)convcnt/ (float)(o_steps * k_steps));
    char fname[100];
    snprintf(fname, 100, "output/arnold_%d.csv", (unsigned)time(NULL));
    save_to_csv(plot, fname);
    free(plot.data);
}

int main (int argc, char **argv)
{   

    printf("bpl: %d\n", mp_bits_per_limb);
    size_t h = 3000, w = 1500;
    double ol = DEF_OL, oh = DEF_OH, kl = DEF_KL, kh = DEF_KH, eps = 0.0003f;
    if (argc > 1) {
        w = atol(argv[1]);
        h = atol(argv[2]);
    }
 
    printf("Generating data on a %zu * %zu grid.\n", w, h);
    mk_plot(h, w, kl, kh, ol, oh, eps);
  
}