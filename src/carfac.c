#include "carfac.h"


/*
 * Various convenenience functions
 */

FLOAT_T *linspace(FLOAT_T lo, FLOAT_T hi, int num_points)
{
    FLOAT_T delta = (hi-lo)/(num_points-1);
    FLOAT_T *vals = malloc(sizeof(FLOAT_T)*num_points);

    vals[0] = lo;
    for(int k=1; k<num_points; k++) {
        vals[k] = vals[k-1] + delta;    
    }
    return vals;
}

FLOAT_T *zeros(int num)
{
    FLOAT_T *vec = calloc(sizeof(FLOAT_T), num);

    for(int k=0; k < num; k++) {
        vec[k] = 0.0;
    }
    return vec;
}

FLOAT_T *ones(int num)
{
    FLOAT_T *vec = malloc(sizeof(FLOAT_T)*num);

    for(int k=0; k<num; k++) {
        vec[k] = 1.0;
    }
    return vec;
}

/*
 * Function: pmul
 * --------------
 *
 *  Computes pointwise multiplication of two vectors (arrays of FLOAT_Ts)
 *
 */

FLOAT_T * pmul(FLOAT_T *v1, FLOAT_T *v2, int len)
{
    FLOAT_T *out = malloc(sizeof(FLOAT_T)*len);

    for(int k=0; k<len; k++) {
        out[k] = v1[k]*v2[k];
    }
    return out;
}

/*
 * Function: smul
 * --------------
 *
 *  Multiplies a vector (array of FLOAT_Ts) by a scalar (FLOAT_T). 
 *
 */

int smul(FLOAT_T *v1, FLOAT_T g, int len)
{
    for(int k=0; k<len; k++) {
        v1[k] *= g;
    }
    return 0;
}


/*
 * Function: argmax
 * ----------------
 *
 *  Just what it sounds like (returns index of first instance of maximum value)
 *
 */

int argmax(FLOAT_T *v, int len) 
{
    int idx = 0;
    FLOAT_T max_val = v[0];

    for(int k=1; k < len; k++) {
        if(v[k] > max_val) {
            max_val = v[k];
            idx = k;
        }
    }
    return idx;
}


FLOAT_T * greenwood(FLOAT_T *points, int num_points)
{
    FLOAT_T * vals = malloc(sizeof(FLOAT_T)*num_points);

    for(int k=0; k<num_points; k++) {
        vals[k] = 165.4*(pow(10.,(2.1*points[k]))-1.0);
    }
    return vals;
}

FLOAT_T inline clip_lo(FLOAT_T val, FLOAT_T thresh)
{
    return val < thresh ? thresh : val;
}


/*
 * ----------------
 * CARFAC-AGC Model
 * ----------------
 */
carfacagc_state_s *carfacagc_init(bm_parameters_s *bmp, ihc_parameters_s *ihcp,
                             ohc_parameters_s *ohcp, int block_size, FLOAT_T fs)
{
    carfacagc_state_s *state = malloc(sizeof(carfacagc_state_s));
    
    // basic dimensions
    state->block_size = block_size;
    state->num_sections = bmp->num_sections;

    // init BM-derived parameters
    FLOAT_T * x = linspace(bmp->x_hi, bmp->x_lo, state->num_sections);
    state->f = greenwood(x, state->num_sections);
    free(x);
    state->a0 = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->c0 = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->h = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->r = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->r1 = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->d_rz = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->g = malloc(sizeof(FLOAT_T)*state->num_sections);
    for(int k=0; k<state->num_sections; k++) {
        state->a0[k] = cos(2.0*M_PI*state->f[k]/fs);
        state->c0[k] = sin(2.0*M_PI*state->f[k]/fs);
        state->h[k] = state->c0[k];
        state->r[k] = 1.0 - bmp->damping*2.0*M_PI*state->f[k]/fs;  // actual
        state->r1[k] = 1.0 - bmp->damping*2.0*M_PI*state->f[k]/fs; // minimum set point
        state->d_rz[k] = 0.7*(1.0 - state->r1[k]);
        state->g[k] = (1.0-2.0*state->a0[k]*state->r[k]+state->r[k]*state->r[k]) / 
                      (1.0 - (2.0*state->a0[k]-state->h[k]*state->c0[k])*state->r[k] 
                       + state->r[k]*state->r[k]);
    }

    // init IHC-derived parameters
    state->q = 1.0/(1.0 + (2.0*M_PI*ihcp->hpf_cutoff/fs));
    state->c_in = 1.0/(fs*ihcp->tau_in);
    state->c_out = 1.0/(fs*ihcp->tau_out); 
    state->c_ihc = 1.0/(fs*ihcp->tau_ihc);

    // init OHC-derived parameters
    state->scale = ohcp->scale;
    state->offset = ohcp->offset;
    state->b = malloc(sizeof(FLOAT_T)*state->num_sections);
    for(int s=0; s < state->num_sections; s++) {
        state->b[s] = ohcp->b[s];
    }
    
    // AGC parameters
    FLOAT_T tau_agc[4] = {0.002, 0.008, 0.032, 0.128};
    FLOAT_T shift_agc[4];
    FLOAT_T spread_sq_agc[4];
    state->c_agc = malloc(sizeof(FLOAT_T)*4);
    state->sa = malloc(sizeof(FLOAT_T)*4);
    state->sb = malloc(sizeof(FLOAT_T)*4);
    state->sc = malloc(sizeof(FLOAT_T)*4);
    for(int k=0; k<4; k++) {
        state->c_agc[k] = 8.0 * pow(2.0, k) / (fs*tau_agc[k]);
        shift_agc[k] = state->c_agc[k] * 0.65 * pow(sqrt(2.0), k);
        spread_sq_agc[k] = state->c_agc[k] * (pow(1.65,2)+1.0) * pow(2.0, k);
        state->sa[k] = (spread_sq_agc[k] + pow(shift_agc[k],2.0) - shift_agc[k])/2.0;
        state->sb[k] = (spread_sq_agc[k] + pow(shift_agc[k],2.0) + shift_agc[k])/2.0;
        state->sc[k] = 1.0 - state->sa[k] - state->sb[k];
    }


    // init model state variables
    int dim2d = state->block_size*state->num_sections;
    state->bm = zeros(dim2d);
    state->bm_hpf = zeros(dim2d);
    state->ihc_out = zeros(dim2d);
    state->ihc_state = zeros(dim2d);
    state->w0 = zeros(state->num_sections);
    state->w1 = zeros(state->num_sections);
    state->w1_old = zeros(state->num_sections);
    state->trans = ones(state->num_sections);
    state->acc8 = zeros(state->num_sections);
    state->acc16 = zeros(state->num_sections);
    state->acc32 = zeros(state->num_sections);
    state->acc64 = zeros(state->num_sections);
    state->agc0 = zeros(state->num_sections);
    state->agc1 = zeros(state->num_sections);
    state->agc2 = zeros(state->num_sections);
    state->agc3 = zeros(state->num_sections);
    state->agc = zeros(state->num_sections*state->block_size);

    // init other internal variables
    state->ihc_new = zeros(state->num_sections);
    state->z = zeros(state->num_sections);
    state->v_mem = zeros(state->num_sections);
    state->v_ohc = zeros(state->num_sections);
    state->sqrd = zeros(state->num_sections);
    state->nlf = zeros(state->num_sections);
    state->prev = 0;
    state->w0_new = 0;

    return state;   
}

int carfacagc_free_except_signal(carfacagc_state_s * cs)
{

    // free parameters
    free(cs->f);
    free(cs->a0);
    free(cs->c0);
    free(cs->r);
    free(cs->r1);
    free(cs->d_rz);
    free(cs->b);
    free(cs->h);
    free(cs->g);
    free(cs->c_agc);
    free(cs->sa);
    free(cs->sb);
    free(cs->sc);

    // free state parameters
    free(cs->bm);
    free(cs->bm_hpf);
    free(cs->ihc_state);
    free(cs->w0);
    free(cs->w1);
    free(cs->w1_old);
    free(cs->trans);
    free(cs->acc8);
    free(cs->acc16);
    free(cs->acc32);
    free(cs->acc64);
    free(cs->agc);
    free(cs->agc0);
    free(cs->agc1);
    free(cs->agc2);
    free(cs->agc3);

    // free dummy variables
    free(cs->ihc_new);
    free(cs->z);
    free(cs->v_mem);
    free(cs->v_ohc);
    free(cs->sqrd);
    free(cs->nlf);

    // free car_state itself
    free(cs);

    return 0;
}

int carfacagc_free(carfacagc_state_s * cs)
{

    // free parameters
    free(cs->f);
    free(cs->a0);
    free(cs->c0);
    free(cs->r);
    free(cs->r1);
    free(cs->d_rz);
    free(cs->b);
    free(cs->h);
    free(cs->g);
    free(cs->c_agc);
    free(cs->sa);
    free(cs->sb);
    free(cs->sc);

    // free state parameters
    free(cs->bm);
    free(cs->bm_hpf);
    free(cs->ihc_out);
    free(cs->ihc_state);
    free(cs->w0);
    free(cs->w1);
    free(cs->w1_old);
    free(cs->trans);
    free(cs->acc8);
    free(cs->acc16);
    free(cs->acc32);
    free(cs->acc64);
    free(cs->agc);
    free(cs->agc0);
    free(cs->agc1);
    free(cs->agc2);
    free(cs->agc3);

    // free dummy variables
    free(cs->ihc_new);
    free(cs->z);
    free(cs->v_mem);
    free(cs->v_ohc);
    free(cs->sqrd);
    free(cs->nlf);

    // free car_state itself
    free(cs);

    return 0;
}


void inline spatial_filter(FLOAT_T * states, FLOAT_T sa, FLOAT_T sb, FLOAT_T sc, int num_states)
{
    /* Rigid boundary conditions */
    states[0] = (sa+sc)*states[0] + sb*states[1];   

    /* Circular boundary conditions */
    //states[0] = sa*states[num_states-1] + sc*states[0] + sb*states[1];

    /* Main loop */
    for(int s=1; s < num_states-1; s++) {
        states[s] = sa*states[s-1] + sc*states[s] + sb*states[s+1];
    }

    /* Rigid boundary conditions */
    states[num_states-1] = sa*states[num_states-2] + (sb+sc)*states[num_states-1];

    /* Circular boundary conditions */
    //states[num_states-1] = sa*states[num_states-2]+sc*states[num_states-1] + sb*states[0];
}


int carfacagc_process_block(carfacagc_state_s * cs, FLOAT_T * sig)
{
    int tm1;    // (t-1) index
    int offset;
    for(int t=0; t < cs->block_size; t++) {
        // set up index for t-1 (circularly)
        tm1 = (t==0) ? cs->block_size-1 : (t-1);
        for(int s=0; s < cs->num_sections; s++) {
            // indexing chores
            cs->prev = (s==0) ? sig[t] : cs->bm[(s-1)*cs->block_size + t];
            offset = s*cs->block_size;              // avoid repeated calculation

            // cascading of signal down BM
            cs->w0_new = cs->prev + cs->r[s]*(cs->a0[s]*cs->w0[s] - cs->c0[s]*cs->w1[s]);
            cs->w1[s] = cs->r[s] * (cs->a0[s]*cs->w1[s] + cs->c0[s]*cs->w0[s]);
            cs->w0[s] = cs->w0_new;
            cs->bm[offset + t] = cs->g[s]*(cs->prev + cs->h[s]*cs->w1[s]);

            // local filtering and haircells
            cs->bm_hpf[offset+t] = cs->q*(cs->bm_hpf[offset+tm1] 
                                            + cs->bm[offset+t]
                                            - cs->bm[offset+tm1]);
            cs->z[s] = clip_lo((cs->bm_hpf[offset+t]+0.175),0);
            cs->v_mem[s] = pow(cs->z[s],3)/(pow(cs->z[s],3)+pow(cs->z[s],2)+0.1);
            cs->ihc_new[s] = cs->v_mem[s]*cs->trans[s];
            cs->trans[s] += cs->c_in*(1.0-cs->trans[s]) - cs->c_out*cs->ihc_new[s];
            cs->ihc_state[offset+t] = (1.0-cs->c_ihc)*cs->ihc_state[offset+tm1] 
                                       + cs->c_ihc*cs->ihc_new[s];
            cs->ihc_out[offset+t] = (1.0-cs->c_ihc)*cs->ihc_out[offset+tm1] 
                                       + cs->c_ihc*cs->ihc_state[offset+t];
            cs->v_ohc[s] = cs->w1[s] - cs->w1_old[s];
            cs->w1_old[s] = cs->w1[s];
            cs->sqrd[s] = pow((cs->v_ohc[s]*cs->scale + cs->offset),2);
            cs->nlf[s] = 1.0/(1.0 + pow((cs->scale*cs->v_ohc[s] + cs->offset),2));
            cs->acc8[s] += cs->ihc_out[offset+t]/8.0;
            }

        // AGC calculations
        if(t%64 == 0) {
            //LPF in time domain
            for(int s=0; s < cs->num_sections; s++) {
                cs->agc3[s] = (1.0 - cs->c_agc[3])*cs->agc3[s] + cs->c_agc[3]*cs->acc64[s]; 
                cs->acc64[s] *= 0;
            }
            //LPF in spatial domain -- must be done after time filtering
            spatial_filter(cs->agc3, cs->sa[3], cs->sb[3], cs->sc[3], cs->num_sections);
        }
        if(t%32 == 0) {
            //LPF in channel (time)
            for(int s=0; s < cs->num_sections; s++) {
                cs->agc2[s] = (1.0 - cs->c_agc[2])*cs->agc2[s] 
                            + cs->c_agc[2]*(cs->acc32[s] + 2.0*cs->agc3[s]);
                cs->acc64[s] += cs->acc32[s];
                cs->acc32[s] *= 0;
            }
            //LPF spatially
            spatial_filter(cs->agc2, cs->sa[2], cs->sb[2], cs->sc[2], cs->num_sections);
        }
        if(t%16 == 0) {
            //LPF in channel (time)
            for(int s=0; s < cs->num_sections; s++) {
                cs->agc1[s] = (1.0 - cs->c_agc[1])*cs->agc1[s]
                            + cs->c_agc[1]*(cs->acc16[s] + 2.0*cs->agc2[s]);
                cs->acc32[s] += cs->acc16[s];
                cs->acc16[s] *= 0;
            }
            //LPF spatially
            spatial_filter(cs->agc1, cs->sa[1], cs->sb[1], cs->sc[1], cs->num_sections);
        }
        if(t%8 == 0) {
            //LPF in channel (time)
            for(int s=0; s < cs->num_sections; s++) {
                cs->agc0[s] = (1.0 - cs->c_agc[0])*cs->agc0[s]
                            + cs->c_agc[0]*(cs->acc8[s] + 2.0*cs->agc1[s]);
                cs->acc16[s] += cs->acc8[s];
                cs->acc8[s] *= 0;
            }
            //LPF spatially and gain adjustment
            spatial_filter(cs->agc0, cs->sa[0], cs->sb[0], cs->sc[0], cs->num_sections);
            for(int s=0; s < cs->num_sections; s++) {
                cs->agc[s*cs->block_size + t] = cs->agc0[s];
                cs->b[s] = cs->agc0[s];
                cs->r[s] = cs->r1[s] + cs->d_rz[s]*(1.0 - cs->b[s])*cs->nlf[s];
                cs->g[s] = (1.0 - 2.0*cs->a0[s]*cs->r[s] + cs->r[s]*cs->r[s]) /
                           (1.0 - (2.0 * cs->a0[s] - cs->h[s]*cs->c0[s])*cs->r[s]
                                                                + cs->r[s]*cs->r[s]);
            }
        }
        else {
            for(int s=0; s < cs->num_sections; s++) {
                cs->agc[s*cs->block_size + t] = cs->agc[s*cs->block_size + tm1];
            }
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
//CARFAC Model
////////////////////////////////////////////////////////////////////////////////
carfac_state_s * carfac_init(bm_parameters_s *bmp, ihc_parameters_s *ihcp,
                        ohc_parameters_s *ohcp, int block_size, FLOAT_T fs)
{
    carfac_state_s * state = malloc(sizeof(carfac_state_s));
    
    // basic dimensions
    state->block_size = block_size;
    state->num_sections = bmp->num_sections;

    // init BM-derived parameters
    FLOAT_T * x = linspace(bmp->x_hi, bmp->x_lo, state->num_sections);
    state->f = greenwood(x, state->num_sections);
    free(x);
    state->a0 = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->c0 = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->h = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->r = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->r1 = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->d_rz = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->g = malloc(sizeof(FLOAT_T)*state->num_sections);
    for(int k=0; k<state->num_sections; k++) {
        state->a0[k] = cos(2.0*M_PI*state->f[k]/fs);
        state->c0[k] = sin(2.0*M_PI*state->f[k]/fs);
        state->h[k] = state->c0[k];
        state->r[k] = 1.0 - bmp->damping*2.0*M_PI*state->f[k]/fs;  // actual
        state->r1[k] = 1.0 - bmp->damping*2.0*M_PI*state->f[k]/fs; // minimum set point
        state->d_rz[k] = 0.7*(1.0 - state->r1[k]);
        state->g[k] = (1.0-2.0*state->a0[k]*state->r[k]+state->r[k]*state->r[k]) / 
                      (1.0 - (2.0*state->a0[k]-state->h[k]*state->c0[k])*state->r[k] 
                       + state->r[k]*state->r[k]);
    }

    // init IHC-derived parameters
    state->q = 1.0/(1.0 + (2.0*M_PI*ihcp->hpf_cutoff/fs));
    state->c_in = 1.0/(fs*ihcp->tau_in);
    state->c_out = 1.0/(fs*ihcp->tau_out); 
    state->c_ihc = 1.0/(fs*ihcp->tau_ihc);

    // init OHC-derived parameters
    state->scale = ohcp->scale;
    state->offset = ohcp->offset;
    state->b = malloc(sizeof(FLOAT_T)*state->num_sections);
    for(int k=0; k< state->num_sections; k++) {
        state->b[k] = ohcp->b[k];
    }
    //d_rz set above


    // init model state variables
    int dim2d = state->block_size*state->num_sections;
    state->bm = zeros(dim2d);
    state->bm_hpf = zeros(dim2d);
    state->ihc_out = zeros(dim2d);
    state->ihc_state = zeros(dim2d);
    state->w0 = zeros(state->num_sections);
    state->w1 = zeros(state->num_sections);
    state->w1_old = zeros(state->num_sections);
    state->trans = ones(state->num_sections);

    // init dummy variables
    state->ihc_new = zeros(state->num_sections);
    state->z = zeros(state->num_sections);
    state->v_mem = zeros(state->num_sections);
    state->v_ohc = zeros(state->num_sections);
    state->sqrd = zeros(state->num_sections);
    state->nlf = zeros(state->num_sections);
    state->prev = 0;
    state->w0_new = 0;

    return state;   
}


int carfac_free(carfac_state_s * cs)
{
    // free malloc-ed members of car_state
    free(cs->f);
    free(cs->a0);
    free(cs->c0);
    free(cs->r);
    free(cs->r1);
    free(cs->d_rz);
    free(cs->h);
    free(cs->g);

    // free state parameters
    free(cs->bm);
    free(cs->bm_hpf);
    free(cs->ihc_out);
    free(cs->ihc_state);
    free(cs->w0);
    free(cs->w1);
    free(cs->w1_old);
    free(cs->trans);
    free(cs->b);

    // free dummy variables
    free(cs->ihc_new);
    free(cs->z);
    free(cs->v_mem);
    free(cs->v_ohc);
    free(cs->sqrd);
    free(cs->nlf);

    // free car_state itself
    free(cs);

    return 0;
}

int carfac_process_block(carfac_state_s * cs, FLOAT_T * sig)
{
    int tm1;    // (t-1) index
    int offset;
    for(int t=0; t < cs->block_size; t++) {
        // set up index for t-1 (circularly)
        tm1 = (t==0) ? cs->block_size-1 : (t-1);
        for(int s=0; s < cs->num_sections; s++) {
            // indexing chores
            cs->prev = (s==0) ? sig[t] : cs->bm[(s-1)*cs->block_size + t];
            offset = s*cs->block_size;              // avoid repeated calculation

            // cascading of signal down BM
            cs->w0_new = cs->prev + cs->r[s]*(cs->a0[s]*cs->w0[s] - cs->c0[s]*cs->w1[s]);
            cs->w1[s] = cs->r[s] * (cs->a0[s]*cs->w1[s] + cs->c0[s]*cs->w0[s]);
            cs->w0[s] = cs->w0_new;
            cs->bm[offset + t] = cs->g[s]*(cs->prev + cs->h[s]*cs->w1[s]);

            // local filtering and haircells
            cs->bm_hpf[offset+t] = cs->q*(cs->bm_hpf[offset+tm1] 
                                            + cs->bm[offset+t]
                                            - cs->bm[offset+tm1]);
            cs->z[s] = clip_lo((cs->bm_hpf[offset+t]+0.175),0);
            cs->v_mem[s] = pow(cs->z[s],3)/(pow(cs->z[s],3)+pow(cs->z[s],2)+0.1);
            cs->ihc_new[s] = cs->v_mem[s]*cs->trans[s];
            cs->trans[s] += cs->c_in*(1.0-cs->trans[s]) - cs->c_out*cs->ihc_new[s];
            cs->ihc_state[offset+t] = (1.0-cs->c_ihc)*cs->ihc_state[offset+tm1] 
                                       + cs->c_ihc*cs->ihc_new[s];
            cs->ihc_out[offset+t] = (1.0-cs->c_ihc)*cs->ihc_out[offset+tm1] 
                                       + cs->c_ihc*cs->ihc_state[offset+t];
            cs->v_ohc[s] = cs->w1[s] - cs->w1_old[s];
            cs->w1_old[s] = cs->w1[s];
            cs->sqrd[s] = pow((cs->v_ohc[s]*cs->scale + cs->offset),2);
            cs->nlf[s] = 1.0/(1.0 + pow((cs->scale*cs->v_ohc[s] + cs->offset),2));
            cs->r[s] = cs->r1[s] + cs->d_rz[s]*(1.0-cs->b[s])*cs->nlf[s];
            cs->g[s] = (1.0 - 2.0*cs->a0[s]*cs->r[s] + cs->r[s]*cs->r[s]) / 
                            (1.0 - (2.0*cs->a0[s] - cs->h[s]*cs->c0[s])*cs->r[s] + 
                                    cs->r[s]*cs->r[s]);
        }
    }
    return 0;
}



////////////////////////////////////////////////////////////////////////////////
//CAR Model
////////////////////////////////////////////////////////////////////////////////

car_state_s * car_init(bm_parameters_s *bmp, ihc_parameters_s *ihcp, int block_size, FLOAT_T fs)
{
    car_state_s * state = malloc(sizeof(car_state_s));
    
    // basic dimensions
    state->block_size = block_size;
    state->num_sections = bmp->num_sections;

    // init BM-derived parameters
    FLOAT_T * x = linspace(bmp->x_hi, bmp->x_lo, state->num_sections);
    state->f = greenwood(x, state->num_sections);
    free(x);
    state->a0 = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->c0 = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->h = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->r = malloc(sizeof(FLOAT_T)*state->num_sections);
    state->g = malloc(sizeof(FLOAT_T)*state->num_sections);
    for(int k=0; k<state->num_sections; k++) {
        state->a0[k] = cos(2.0*M_PI*state->f[k]/fs);
        state->c0[k] = sin(2.0*M_PI*state->f[k]/fs);
        state->h[k] = state->c0[k];
        state->r[k] = 1.0 - bmp->damping*2.0*M_PI*state->f[k]/fs;
        state->g[k] = (1.0-2.0*state->a0[k]*state->r[k]+state->r[k]*state->r[k]) / 
                      (1.0 - (2.0*state->a0[k]-state->h[k]*state->c0[k])*state->r[k] 
                       + state->r[k]*state->r[k]);
    }

    // init IHC-derived parameters
    state->q = 1.0/(1.0 + (2.0*M_PI*ihcp->hpf_cutoff/fs));
    state->c_in = 1.0/(fs*ihcp->tau_in);
    state->c_out = 1.0/(fs*ihcp->tau_out); 
    state->c_ihc = 1.0/(fs*ihcp->tau_ihc);

    // init model state variables
    int dim2d = state->block_size*state->num_sections;
    state->bm = zeros(dim2d);
    state->bm_hpf = zeros(dim2d);
    state->ihc_out = zeros(dim2d);
    state->ihc_state = zeros(dim2d);
    state->w0 = zeros(state->num_sections);
    state->w1 = zeros(state->num_sections);
    state->trans = ones(state->num_sections);

    // init dummy variables
    state->ihc_new = zeros(state->num_sections);
    state->z = zeros(state->num_sections);
    state->v_mem = zeros(state->num_sections);
    state->prev = 0;
    state->w0_new = 0;

    return state;   
}


int car_free(car_state_s * cs)
{
    // free malloc-ed members of car_state
    free(cs->f);
    free(cs->a0);
    free(cs->c0);
    free(cs->r);
    free(cs->h);
    free(cs->g);

    // free state parameters
    free(cs->bm);
    free(cs->bm_hpf);
    free(cs->ihc_out);
    free(cs->ihc_state);
    free(cs->w0);
    free(cs->w1);
    free(cs->trans);

    // free dummy variables
    free(cs->ihc_new);
    free(cs->z);
    free(cs->v_mem);

    // free car_state itself
    free(cs);

    return 0;
}


int car_process_block(car_state_s * cs, FLOAT_T * sig)
{
    int tm1;    // (t-1) index
    int offset;
    for(int t=0; t < cs->block_size; t++) {
        // set up index for t-1 (circularly)
        tm1 = (t==0) ? cs->block_size-1 : (t-1);
        for(int s=0; s < cs->num_sections; s++) {
            // indexing chores
            cs->prev = (s==0) ? sig[t] : cs->bm[(s-1)*cs->block_size + t];
            offset = s*cs->block_size;              // avoid repeated calculation

            // cascading of signal down BM
            cs->w0_new = cs->prev + cs->r[s]*(cs->a0[s]*cs->w0[s] - cs->c0[s]*cs->w1[s]);
            cs->w1[s] = cs->r[s] * (cs->a0[s]*cs->w1[s] + cs->c0[s]*cs->w0[s]);
            cs->w0[s] = cs->w0_new;
            cs->bm[offset + t] = cs->g[s]*(cs->prev + cs->h[s]*cs->w1[s]);

            // local filtering and haircells
            cs->bm_hpf[offset+t] = cs->q*(cs->bm_hpf[offset+tm1] 
                                            + cs->bm[offset+t]
                                            - cs->bm[offset+tm1]);
            cs->z[s] = clip_lo((cs->bm_hpf[offset+t]+0.175),0);
            cs->v_mem[s] = pow(cs->z[s],3)/(pow(cs->z[s],3)+pow(cs->z[s],2)+0.1);
            cs->ihc_new[s] = cs->v_mem[s]*cs->trans[s];
            cs->trans[s] += cs->c_in*(1.0-cs->trans[s]) - cs->c_out*cs->ihc_new[s];
            cs->ihc_state[offset+t] = (1.0-cs->c_ihc)*cs->ihc_state[offset+tm1] 
                                       + cs->c_ihc*cs->ihc_new[s];
            cs->ihc_out[offset+t] = (1.0-cs->c_ihc)*cs->ihc_out[offset+tm1] 
                                       + cs->c_ihc*cs->ihc_state[offset+t];
        }
    }
    return 0;
}


/*
 * --------------
 * SAI generation
 * --------------
 */

sai_s * sai_init(int num_frames, int num_sections, int frame_len_n) 
{
    sai_s * sai = malloc(sizeof(sai_s));

    sai->num_frames = num_frames;
    sai->num_sections = num_sections;
    sai->frame_len_n = frame_len_n;

    return sai;
}


int sai_free(sai_s * sai) 
{
    free(sai->images);
    free(sai->delay_t);
    free(sai->times);
    free(sai);

    return 0;
}

int sai_free_except_images(sai_s * sai) 
{
    free(sai->delay_t);
    free(sai->times);
    free(sai);

    return 0;
}


sai_s * sai_generate(FLOAT_T * nap, FLOAT_T fs, sai_parameters_s * sai_params)
{
    int trig_win_n = sai_params->trig_win_t * fs;   
    int frame_len_n = sai_params->num_trig_win * trig_win_n;
    int num_sections = sai_params->num_sections;
    FLOAT_T frame_len_t = sai_params->num_trig_win * sai_params->trig_win_t;
    
    int delay_n = frame_len_n;
    FLOAT_T * delay_t = malloc(sizeof(FLOAT_T)*frame_len_n);
    FLOAT_T dt = 1.0/fs;
    for(int k = 0; k < frame_len_n; k++) {
        delay_t[k] = ((FLOAT_T) k) * dt;
    }
    int adv_n = (int) (sai_params->adv_t * fs);
    int num_frames = (int) ((sai_params->num_samples - frame_len_n)/adv_n);
    
    FLOAT_T * images = malloc(sizeof(FLOAT_T) * num_frames 
                                            * num_sections
                                            * frame_len_n);
    int frame_n = 0;
    int trig_n;
    FLOAT_T peak_val, alpha, beta;
    FLOAT_T * f_in = zeros(frame_len_n);
    for(int f = 0; f < num_frames; f++) {
        for(int s = 0; s < num_sections; s++) {
            for(int w = 0; w < sai_params->num_trig_win; w++) {

                // find the trigger point within the trigger window
                trig_n = argmax( &nap[s*sai_params->num_samples
                                    + frame_n
                                    + w*trig_win_n], trig_win_n);   
                trig_n += frame_n;  // offset from base
                trig_n = (trig_n > 0) ? trig_n : 0;
                peak_val = nap[s*sai_params->num_samples + trig_n];
                
                // calculate the coefficients for smoothing/resetting
                alpha = peak_val / (1.0 + peak_val);
                beta = 1.0 - alpha;

                // set up the data for the triggered/reset portion
                if(trig_n - delay_n < 0) {
                    smul(f_in, 0.0, frame_len_n);
                    for(int k=0; k < trig_n; k++) {
                        f_in[frame_len_n-trig_n+k] = nap[s*sai_params->num_samples + k];
                    }
                } else {
                    for(int k=0; k < frame_len_n; k++) {
                        f_in[k] = nap[s*sai_params->num_samples + trig_n - frame_len_n + k];
                    }
                }

                // finally set the image
                for(int k=0; k < frame_len_n; k++) {
                    images[f*(num_sections*frame_len_n) + s*(frame_len_n) + k] = 
                        beta * images[f*(num_sections*frame_len_n) + s*(frame_len_n) + k] +
                        alpha * f_in[k];
                }
            }
        }
        frame_n += adv_n;
    }

    FLOAT_T * times = malloc(sizeof(FLOAT_T)*num_frames);
    FLOAT_T coef = ((FLOAT_T) adv_n)/fs;
    for(int k = 0; k < num_frames; k++) {
        times[k] = k*coef + frame_len_t;
    }

    // set up SAI structure for return
    sai_s * sai = sai_init(num_frames, num_sections, frame_len_n);
    sai->times = times;
    sai->images = images;
    sai->delay_t = delay_t;

    // free heap-allocated memory that isn't part of struct
    free(f_in);

    return sai;
}

