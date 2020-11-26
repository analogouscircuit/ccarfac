#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sndfile.h>

#include "siggen.h"
#include "carfac.h"
#include "ctopy.h"

/* Note that FLOAT_T is defined in carfac.h */

int main( int argc, char *argv[] )
{
    FLOAT_T dur = 1.0;
    FLOAT_T fs;
    FLOAT_T * sig;

    if( argc == 2) {
        /* Read a wav file and prepare for processing */
        SNDFILE * sf;
        SF_INFO info;
        int num_channels;
        int num, num_frames;

        info.format = 0;
        sf = sf_open(argv[1], SFM_READ, &info);
        if (sf == NULL) {
            printf("Error opening file %s! Exiting...\n", argv[1]);
            exit(0);
        }
        fs = (FLOAT_T) info.samplerate;
        num_channels = info.channels;
        num_frames = info.frames;
        if(num_channels != 1) {
            printf("Can only read wav files with one channel! Exiting...\n");
            exit(0);
        }
        dur = num_frames/fs;
        sig = malloc(sizeof(FLOAT_T)*num_frames);
        num = sf_read_double(sf, sig, num_frames);
        sf_close(sf);

    } else {
        /* If no file given, generate a test signal */
        FLOAT_T sig_gain = 1.0;
        FLOAT_T f_lo = 100.0;
        FLOAT_T f_hi = 1000.0;
        dur = 1.0;
        fs = 32000.0;

        sig = swept_sine(dur, fs, f_lo, f_hi, 0.0, sig_gain, LOGARITHMIC);
        //sig = sinusoid_f(dur, fs, 500.0, 0.0, sig_gain);
    }

    /* General processing parameters */
    int num_points = fs*dur;
    int block_size = num_points;
    int num_blocks = (int) num_points/block_size;
    int num_sections = 50;

    /* Establish parameters for the basilar membrane (bm) */
    bm_parameters_s bm_params = { .num_sections = num_sections,
                                  .x_lo = 0.1, // 0.1 about 100 Hz
                                  .x_hi = 0.8, // 0.8 about 7750 Hz
                                  .damping = 0.2};

    /* Establish parameters for the inner hair cells (ihc) */
    ihc_parameters_s ihc_params = { .hpf_cutoff = 20.0,
                                    .tau_in = 10.0e-3,
                                    .tau_out = 0.5e-3,
                                    .tau_ihc = 80.0e-6 }; 

    /* Establish parameters for the outer hair cells (ohc) */
    FLOAT_T b_init = 1.0;
    FLOAT_T * b = malloc(sizeof(FLOAT_T)*num_sections);
    for(int k=0; k < num_sections; k++) {
        b[k] = b_init;
    }

    ohc_parameters_s ohc_params = { .scale = 0.1,
                                    .offset = 0.04,
                                    .b = b };

    /* initialize a CARFAC structure */
    carfacagc_state_s * carfac = carfacagc_init(&bm_params,
                                                &ihc_params, 
                                                &ohc_params, 
                                                block_size, 
                                                fs);


    /* Process through CARFAC model */
    FLOAT_T * bm_out = malloc(sizeof(FLOAT_T)*num_blocks*block_size*bm_params.num_sections);
    FLOAT_T * ihc_out = malloc(sizeof(FLOAT_T)*num_blocks*block_size*bm_params.num_sections);
    for(int k = 0; k < num_blocks; k++) {
        carfacagc_process_block(carfac, (sig+k*block_size)); 
        /* copy the state into a record */
        for(int t=0; t < block_size; t++) {
            for(int s=0; s < bm_params.num_sections; s++) {
                bm_out[(k*block_size) + s*(block_size*num_blocks) + t] =
                                                        carfac->bm_hpf[s*block_size+t];
                ihc_out[(k*block_size) + s*(block_size*num_blocks) + t] =
                                                        carfac->ihc_out[s*block_size+t];
            }
        }
    }
    

    /* Now generate SAI images */
    sai_parameters_s sai_params = { .trig_win_t = 0.010,
                                    .adv_t = 0.005,
                                    .num_trig_win = 5,
                                    .num_sections = num_sections,
                                    .num_samples = num_points };

    sai_s * sai = sai_generate(ihc_out, fs, &sai_params);


    /* Write data in Python readable format */
    int dims1[2] = {num_sections, num_blocks*block_size};
    int dims2[1] = {1};
    int dims3 [1] = {num_points};
    int dims4 [3] = { sai->num_frames, sai->num_sections, sai->frame_len_n};
    int dims5 [1] = { sai->num_sections };
    char data_type[7] = "double";

    pydata_list * data_list = pydata_list_init();
    pydata_list_push(bm_out, data_type, 2, dims1, "bm_out", data_list);
    pydata_list_push(ihc_out, data_type, 2, dims1, "ihc_out", data_list);
    pydata_list_push(&fs, data_type, 1, dims2, "fs", data_list);
    pydata_list_push(sig, data_type, 1, dims3, "sig", data_list);
    pydata_list_push(sai->images, data_type, 3, dims4, "sai_images", data_list);
    pydata_list_push(carfac->f, data_type, 1, dims5, "f_vals", data_list); 
    pydata_write(data_list, "carfac_test_data");
    pydata_list_free(data_list);
                     

    /* Free heap-allocated data */
    carfacagc_free(carfac);
    free(bm_out);
    free(ihc_out);
    free(sig);
    free(b);
    sai_free(sai);

    return 0;
}
