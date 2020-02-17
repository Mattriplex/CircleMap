#include <png.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <pthread.h>

#define RECINIT 0.1f
#define RMFIRST 200
#define ITLIMIT 300
#define ABS(X) ((X) >= 0.0 ? (X) : -(X))
#define MOD1(X) ((X) - floor(X))

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
    bitmap_t bitmap;
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

/* Given "bitmap", this returns the pixel of bitmap at the point 
   ("x", "y"). */

static pixel_t * pixel_at (bitmap_t * bitmap, int x, int y)
{
    return bitmap->pixels + bitmap->width * y + x;
}
    
/* Write "bitmap" to a PNG file specified by "path"; returns 0 on
   success, non-zero on error. */

static int save_png_to_file (bitmap_t *bitmap, const char *path)
{
    FILE * fp;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    size_t x, y;
    png_byte ** row_pointers = NULL;
    /* "status" contains the return value of this function. At first
       it is set to a value which means 'failure'. When the routine
       has finished its work, it is set to a value which means
       'success'. */
    int status = -1;
    /* The following number is set by trial and error only. I cannot
       see where it it is documented in the libpng manual.
    */
    int pixel_size = 3;
    int depth = 8;
    
    fp = fopen (path, "wb");
    if (! fp) {
        goto fopen_failed;
    }

    png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        goto png_create_write_struct_failed;
    }
    
    info_ptr = png_create_info_struct (png_ptr);
    if (info_ptr == NULL) {
        goto png_create_info_struct_failed;
    }
    
    /* Set up error handling. */

    if (setjmp (png_jmpbuf (png_ptr))) {
        goto png_failure;
    }
    
    /* Set image attributes. */

    png_set_IHDR (png_ptr,
                  info_ptr,
                  bitmap->width,
                  bitmap->height,
                  depth,
                  PNG_COLOR_TYPE_RGB,
                  PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_DEFAULT,
                  PNG_FILTER_TYPE_DEFAULT);
    
    /* Initialize rows of PNG. */

    row_pointers = png_malloc (png_ptr, bitmap->height * sizeof (png_byte *));
    for (y = 0; y < bitmap->height; y++) {
        png_byte *row = 
            png_malloc (png_ptr, sizeof (uint8_t) * bitmap->width * pixel_size);
        row_pointers[y] = row;
        for (x = 0; x < bitmap->width; x++) {
            pixel_t * pixel = pixel_at (bitmap, x, y);
            *row++ = pixel->red;
            *row++ = pixel->green;
            *row++ = pixel->blue;
        }
    }
    
    /* Write the image data to "fp". */

    png_init_io (png_ptr, fp);
    png_set_rows (png_ptr, info_ptr, row_pointers);
    png_write_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    /* The routine has successfully written the file, so we set
       "status" to a value which indicates success. */

    status = 0;
    
    for (y = 0; y < bitmap->height; y++) {
        png_free (png_ptr, row_pointers[y]);
    }
    png_free (png_ptr, row_pointers);
    
 png_failure:
 png_create_info_struct_failed:
    png_destroy_write_struct (&png_ptr, &info_ptr);
 png_create_write_struct_failed:
    fclose (fp);
 fopen_failed:
    return status;
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

int get_recurrence(double O, double K, double eps, double init) {
    /*Ignore first few iterations to get to convergence */
    double theta_0 = init;
    for (int i = 0; i < RMFIRST; i++)
        theta_0 = circle_map(O, K, theta_0);
    /*calculate recurrence number*/
    double theta_1 = circle_map(O,K, theta_0);
    int n = 0;
    for (; n < ITLIMIT && ABS(theta_1 - theta_0) >= eps; n++)
        theta_1 = circle_map(O, K, theta_1);
    return n;
}

pixel_t color_from_recn(int n) {
    pixel_t colo;
    if (n < 6) {

    } else if (n < 50) {
        colo.blue = 255;
    } else if (n < 200) {
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
        for (size_t kidx = 0; kidx < plot.bitmap.height; kidx++, k += plot.k_stepsz) {
            int n = get_recurrence(o, k, plot.eps, RECINIT);
            pixel_t colo = color_from_recn(n);
            pixel_t *px = pixel_at (&plot.bitmap, oidx, kidx);
            px->red = colo.red;
            px->green = colo.green;
            px->blue = colo.blue;
        }
    }
    return NULL;
}

void mk_plot(size_t k_steps, size_t o_steps, double k_low, double k_high, double o_low, double o_high, double eps) {
    /* Set up bitmap and calculate parameter steps*/
    plot_t plot;
    plot.bitmap.width = o_steps;
    plot.bitmap.height = k_steps;
    plot.bitmap.pixels = calloc(o_steps * k_steps, sizeof(pixel_t));
    plot.k_stepsz = (k_high - k_low) / k_steps;
    plot.o_stepsz = (o_high - o_low) / o_steps;
    if (!plot.bitmap.pixels)
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
    snprintf(fname, 100, "arnold_%d.png", (unsigned)time(NULL));
    if (save_png_to_file (& plot.bitmap, fname)) {
	fprintf (stderr, "Error writing file.\n");
    }
    free(plot.bitmap.pixels);
}

int main ()
{
    /* Write the image to a file 'fruit.png'. */
    mk_plot(3000, 1500, 0.0f, 4*M_PI, 0.0f, 1.0f, 0.0001f);
  
}