This is the code which takes care of setting up an HMM analysis for MATLAB to run, and processing the results.

`run_hmm.py` is the main piece of code.
It is invoked as: `./run_hmm.py --par J1234-1234.par --tim J1234-1234.tim --ini J1234-1234.ini`

The `ini` file contains the parameters of the analysis, as well as paths pointing to the MATLAB code,
and the desired results location. Note that for now all paths should be absolute rather than relative.

`do_matlab_analysis.m` is a piece of wrapper MATLAB code which is invoked by `run_hmm.py`.

A small demo analysis on the UTMOST data release for PSR J1452-6036 is in `demo/`.

This is still very much a work in progress.
