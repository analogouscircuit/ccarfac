# ccarfac
Very simple pure C implementation of [Richard Lyon's CARFAC model](https://books.google.com/books?id=ENmiDgAAQBAJ&dq=human+and+machine+hearing+lyon&lr=&source=gbs_navlinks_s).  (See also pycarfac and JCARFAC for Python and Julia wrappers.)

This implementation is a straightforward C adaptation of [Andre van Schaik's Python version](https://github.com/vschaik/CARFAC).
I have found it to be stable, but bear in mind this was written by a non-professional programmer for scientific research.  Use at your own risk.

A basic makefile is included with two build targets: `clib` (the default) and `test`.
The former creates a linkable library, used, for example, by the Julia wrapper; the latter compiles a
simple executable which either takes no argument or a (mono) wav file.  It processes a test signal or the given file through the
CARFAC model and writes the results to disk.  To plot the results, use the plt_cnap.py script. This 
has been tested on Ubuntu 18.04 and an Arch Linux installation, both x86.

The interface is quite straightfoward. 
There structures to hold parameter values for the basilar membrane (`bm_parameters_s`), inner hair cells 
(`ihc_parameters_s`), and outer hair cells (`ohc_parameters_s`).  See Lyon's book or test program
for reasonable parameter values.  These are used, together with a block size and sample frequency,
to initialize the CARFAC model state, `carfacagc_state_s`, using the `carfacagc_init` function. A 
mono signal can then be processed in blocks using `carfacagc_process_block`.  The output is contained
in the CARFAC state structure (e.g. `bm_hpf` or `ihc_out` for basilar membrane or hair cell data respectively).

There is also a function to produce a stabilized auditory image (SAI),  namely `sai_generate`.
