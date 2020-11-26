/*******************************************************************************
 * Includes
 ******************************************************************************/
#include <stdlib.h>
#include <math.h>

/*******************************************************************************
 * Definitions
 ******************************************************************************/
#define FLOAT_T double

/*******************************************************************************
 * Structures
 ******************************************************************************/

typedef struct bm_parameters_s {
    int num_sections;
    FLOAT_T x_lo;
    FLOAT_T x_hi;
    FLOAT_T damping;
} bm_parameters_s;


typedef struct ihc_parameters_s {
    FLOAT_T hpf_cutoff;
    FLOAT_T tau_in;
    FLOAT_T tau_out;
    FLOAT_T tau_ihc;
} ihc_parameters_s;

typedef struct ohc_parameters_s {
    FLOAT_T scale;
    FLOAT_T offset;
    FLOAT_T * b;

} ohc_parameters_s;

typedef struct sai_parameters_s {
    FLOAT_T trig_win_t;         //length of trigger window in seconds
    FLOAT_T adv_t;          // how far to advance the base time between frames
    int num_trig_win;       // number of trigger windows in a frame
    int num_sections;       //
    int num_samples;        //
} sai_parameters_s;

typedef struct sai_s {
    int num_frames;
    int num_sections;
    int frame_len_n;
    FLOAT_T * images;
    FLOAT_T * delay_t;
    FLOAT_T * times;
} sai_s;

typedef struct car_state_s {
    int num_sections;
    int block_size;

    // derived parameters from BM settings
    FLOAT_T * f;
    FLOAT_T * a0;
    FLOAT_T * c0;
    FLOAT_T * r;
    FLOAT_T * h;
    FLOAT_T * g;

    // derived parameters from IHC settings
    FLOAT_T q;          //hpf coefficient
    FLOAT_T c_in;
    FLOAT_T c_out;
    FLOAT_T c_ihc;

    // state parameters
    FLOAT_T * bm;
    FLOAT_T * bm_hpf;
    FLOAT_T * ihc_out;
    FLOAT_T * ihc_state;
    FLOAT_T * w0;
    FLOAT_T * w1;
    FLOAT_T * trans;

    // various dummy variables (so processors doesn't have to keep allocating)
    FLOAT_T * ihc_new;
    FLOAT_T * z;
    FLOAT_T * v_mem;
    FLOAT_T prev;
    FLOAT_T w0_new;
    

} car_state_s;


typedef struct carfac_state_s {
    int num_sections;
    int block_size;

    // derived parameters from BM settings
    FLOAT_T * f;
    FLOAT_T * a0;
    FLOAT_T * c0;
    FLOAT_T * r;
    FLOAT_T * r1;
    FLOAT_T * h;
    FLOAT_T * g;

    // derived parameters from IHC settings
    FLOAT_T q;          //hpf coefficient
    FLOAT_T c_in;
    FLOAT_T c_out;
    FLOAT_T c_ihc;

    // derived parameters from OHC settings
    FLOAT_T scale;
    FLOAT_T offset;
    FLOAT_T * b;
    FLOAT_T * d_rz;

    // state parameters
    FLOAT_T * bm;
    FLOAT_T * bm_hpf;
    FLOAT_T * ihc_out;
    FLOAT_T * ihc_state;
    FLOAT_T * w0;
    FLOAT_T * w1;
    FLOAT_T * w1_old;
    FLOAT_T * trans;

    // various dummy variables (so processors doesn't have to keep allocating)
    FLOAT_T * ihc_new;
    FLOAT_T * z;
    FLOAT_T * v_mem;
    FLOAT_T * v_ohc;
    FLOAT_T * sqrd;
    FLOAT_T * nlf;
    FLOAT_T prev;
    FLOAT_T w0_new;
    

} carfac_state_s;


typedef struct carfacagc_state_s {
    int num_sections;
    int block_size;

    // derived parameters from BM settings
    FLOAT_T * f;
    FLOAT_T * a0;
    FLOAT_T * c0;
    FLOAT_T * r;
    FLOAT_T * r1;
    FLOAT_T * h;
    FLOAT_T * g;

    // derived parameters from IHC settings
    FLOAT_T q;          //hpf coefficient
    FLOAT_T c_in;
    FLOAT_T c_out;
    FLOAT_T c_ihc;

    // derived parameters from OHC settings
    FLOAT_T scale;
    FLOAT_T offset;
    FLOAT_T * b;
    FLOAT_T * d_rz;

    // AGC parameters
    FLOAT_T * c_agc;
    FLOAT_T * sa;
    FLOAT_T * sb;
    FLOAT_T * sc;

    // state parameters
    FLOAT_T * bm;
    FLOAT_T * bm_hpf;
    FLOAT_T * ihc_out;
    FLOAT_T * ihc_state;
    FLOAT_T * w0;
    FLOAT_T * w1;
    FLOAT_T * w1_old;
    FLOAT_T * trans;
    FLOAT_T * acc8;
    FLOAT_T * acc16;
    FLOAT_T * acc32;
    FLOAT_T * acc64;
    FLOAT_T * agc;
    FLOAT_T * agc0;
    FLOAT_T * agc1;
    FLOAT_T * agc2;
    FLOAT_T * agc3;

    // various other internal variables (so processors doesn't have to keep allocating)
    FLOAT_T * ihc_new;
    FLOAT_T * z;
    FLOAT_T * v_mem;
    FLOAT_T * v_ohc;
    FLOAT_T * sqrd;
    FLOAT_T * nlf;
    FLOAT_T prev;
    FLOAT_T w0_new;
    

} carfacagc_state_s;


/*******************************************************************************
 * Function Declarations
 *******************************************************************************/

/*
 * -----------------
 * Utility Functions
 * ------------------
 */

FLOAT_T * linspace(FLOAT_T lo, FLOAT_T hi, int num_points);

FLOAT_T * zeros(int num);

FLOAT_T * ones(int num);

FLOAT_T * pmul(FLOAT_T *v1, FLOAT_T *v2, int len);

int smul(FLOAT_T *v1, FLOAT_T g, int len);

int argmax(FLOAT_T *v, int len);

FLOAT_T * greenwood(FLOAT_T *points, int num_points);

FLOAT_T clip_lo(FLOAT_T val, FLOAT_T thresh); 


/* 
 * -------------
 * SAI functions
 * -------------
 */

sai_s * sai_init(int num_frames, int num_sections, int frame_len_n);

int sai_free(sai_s * sai);

int sai_free_except_images(sai_s * sai);

sai_s * sai_generate(FLOAT_T * nap, FLOAT_T fs, sai_parameters_s * sai_params);



/* 
 * ----------------
 * CARFAC functions
 * ----------------
 */

void spatial_filter(FLOAT_T * states, FLOAT_T sa, FLOAT_T sb, FLOAT_T sc, int num_states);

car_state_s * car_init(bm_parameters_s *bmp, ihc_parameters_s *ihcp,
                       int block_size, FLOAT_T fs);

int car_free(car_state_s * cs);

int car_process_block(car_state_s * cs, FLOAT_T * sig);

carfac_state_s * carfac_init(bm_parameters_s *bmp, ihc_parameters_s *ihcp,
                             ohc_parameters_s *ohcp, int block_size, FLOAT_T fs);

int carfac_free(carfac_state_s * cs);

int carfac_process_block(carfac_state_s * cs, FLOAT_T * sig);

carfacagc_state_s * carfacagc_init(bm_parameters_s *bmp, ihc_parameters_s *ihcp,
                             ohc_parameters_s *ohcp, int block_size, FLOAT_T fs);

int carfacagc_free(carfacagc_state_s * cs);

int carfacagc_free_except_signal(carfacagc_state_s * cs);

int carfacagc_process_block(carfacagc_state_s * cs, FLOAT_T * sig);
