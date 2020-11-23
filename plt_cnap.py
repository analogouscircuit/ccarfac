import array
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import scipy.fftpack as fft
from ctopy import ctopy_read

class MyNormalize(mcolors.Normalize):
    def __call__(self, value, clip=None):
        f = lambda x,a: (2*x)**a*(2*x<1)/2. +(2-(2*(1-1*x))**a)*(2*x>=1)/2.
        return np.ma.masked_array(f(value,0.8))

def main(BM_out, IHC_out, sig, fs):
    num_sections, num_points = BM_out.shape
    print(num_sections, num_points)

    ## staggered plot of IHC output
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    t_max = num_points/fs
    num_points = int(t_max*fs)
    t = np.arange(num_points)/fs
    p = 0.01
    for k in range(num_sections):
        ax1.fill_between(t, 0, IHC_out[k,:num_points]+p*(num_sections-k), facecolor='w',
                edgecolor='k', linewidth=0.6)

    ax2 = fig.add_subplot(1,2,2)
    ax2.imshow(IHC_out, aspect='auto', cmap=cm.binary,
            norm=MyNormalize(vmin=0, vmax=np.max(IHC_out)))

    plt.show()

## Main script
if __name__=="__main__":
    # import data
    imported_data_dict = ctopy_read("carfac_test_data")
    BM_out = imported_data_dict['bm_out']
    IHC_out = imported_data_dict['ihc_out']
    fs = imported_data_dict['fs']
    sig = imported_data_dict['sig']
    f_vals = imported_data_dict['f_vals']
    print("Channel CFs: ", f_vals)

    # run script
    main(BM_out, IHC_out, sig, fs)
