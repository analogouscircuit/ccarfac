import array
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.animation import FuncAnimation
import numpy as np
import scipy.fftpack as fft
import importlib.util
import os
from ctopy import ctopy_read

class MyNormalize(mcolors.Normalize):
    def __call__(self, value, clip=None):
        f = lambda x,a: (2*x)**a*(2*x<1)/2. +(2-(2*(1-1*x))**a)*(2*x>=1)/2.
        return np.ma.masked_array(f(value,0.2))

def animate_SAI(I, fs, f, times, colormap = cm.binary, adv_time = 50):
    num_frames = I.shape[0]
    num_sect = I.shape[1]
    frame_len_n = I.shape[2]
    delays = np.flip(np.arange(frame_len_n)/fs,0)

    # Set up the figure
    fig = plt.figure()
    im = plt.imshow(np.flip(I[0,:,:],0),cmap=colormap,
            norm=MyNormalize(vmin=-0.3, vmax=1.2), origin='lower', aspect='auto')
    # text_labels = ["%.1f" % (np.flip(f,0)[k]) for k in range(num_sect)]
    # ticks = []
    # labels = []
    # for k in range(num_sect):
    #     if k % 8 == 0:
    #         ticks.append(k)
    #         labels.append(text_labels[k])
    # plt.yticks(ticks,labels)
    # plt.ylabel("Channel CF")
    adv = frame_len_n//4
    # ticks = np.arange(4)*adv
    # labels = ["%.0f" % (delays[tick]*1000) for tick in ticks]
    # plt.xticks(ticks, labels)
    # plt.xlabel('delay t (in ms)')

    # Animation update function
    def update(frame_num):
        plt.title("t=%.3f" % (times[frame_num]) )
        im.set_array(np.flip(I[frame_num,:,:],0))
        return im

    # Let 'er rip
    anim = FuncAnimation(fig, update, frames=range(1,num_frames), interval=adv_time, repeat=True)
    plt.show()

## Basic parameters (need to import correctly)
if __name__=="__main__":
    data_dict = ctopy_read("carfac_test_data")
    fs = data_dict["fs"]
    I = data_dict["sai_images"]
    num_frames, num_sections, frame_len_n = I.shape
    tau = np.arange(1,frame_len_n+1)/fs
    f = 1./np.flip(tau)

    ## import signal
    I /= np.max(I)
    animate_SAI(I, fs, f, np.arange(num_frames), colormap = cm.viridis, adv_time = 100)

