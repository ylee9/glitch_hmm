#!/usr/bin/env python

import pulsar_hmm.HMM as HMM
import glob
import libstempo
import numpy as np
import scipy.special as special
import sys
import os.path
import os
import shutil
import subprocess

def chunk(list, n):
    return [list[i::n] for i in range(n)]

def setup_hmm(par, tim, mjd_range=None):
    with open(par,'r') as f:
        psr_name = None
        red_idx = None
        red_amp = None
        efac = 1
        equad = -100
        for line in f.readlines():
            split = line.split()
            if split[0].strip() == 'PSRJ':
                psr_name = split[-1].strip()
                print(psr_name)
            elif split[0].strip() == 'TNRedAmp':
                red_amp = 10**float(split[-1].strip())
            elif split[0].strip() == 'TNRedGam':
                red_idx = float(split[-1].strip())
            elif split[0].strip().lower() == 'tnglobalef':
                efac = float(split[-1].strip())
            elif split[0].strip().lower() == 'tnglobaleq':
                equad = float(split[-1].strip())
            elif split[0].strip() == 'BINARY': # We can't do binaries at this point
                return (None, psr_name)
        if psr_name == 'J0835-4510':
            return (None, psr_name)
        sigma = 1e-21
        #if red_amp is not None and red_idx is not None:
        if True:
            psr = libstempo.tempopulsar(parfile=par, timfile=tim)

            freqs = np.arange(-3e-7, 3e-7, 1e-9)
            #freqs = np.linspace(-3e-7, 2.5e-5, 1501)
            mean_z = (np.mean(np.diff(sorted(psr.toas()))))*86400
            #sigma = max(np.mean(np.diff(freqs))/(mean_z)**1.5, sigma)
            #fdots = np.linspace(-5*sigma*mean_z**0.5, 5*sigma*mean_z**0.5, 11)
            fdot_range = np.min([1e-15, -0.1*psr['F1'].val])
            fdots = np.linspace(-fdot_range, fdot_range, 11)

            #phase_cov = lambda z: psr['F0'].val**2*red_amp**2/(12*np.pi**2)*(np.pi*1e7)**(4-red_idx)/(red_idx)*z**(red_idx - 1)
            #noise_cov = lambda z: [[phase_cov(z)/z**2/red_idx**2, phase_cov(z)/z**3/(red_idx-1)**2], [phase_cov(z)/z**3/(red_idx-1)**2, phase_cov(z)/z**4/(red_idx-2)**2]]

            #if (np.sqrt(phase_cov(mean_z))/mean_z**2/(red_idx-2)/np.mean(np.diff(fdots))) < 1:
            #    sigma = np.mean(np.diff(fdots))/mean_z**0.5
            #    print(f"Adjusting sigma for {psr_name}: inferred sigma was {np.sqrt(phase_cov(1))/mean_z**2/(red_idx-2)}, setting to {sigma}")
            #    sys.stdout.flush()
            #    noise_cov = lambda z: [[sigma**2*z**3/3, sigma**2*z**2/2],[sigma**2*z**2/2, sigma**2*z]]


            sigma = max(np.mean(np.diff(fdots))/(mean_z)**0.5, sigma)
            noise_cov = lambda z: [[sigma**2*z**3/3, sigma**2*z**2/2],[sigma**2*z**2/2, sigma**2*z]]
            print(f"{psr_name}: amp = {red_amp}, idx = {red_idx}, sigma={sigma}")
            sys.stdout.flush()
            hmm = HMM.HMM.from_tempo2(par, tim, freqs, fdots, noise_cov,[], min_toa_gap = 0, mjd_range=mjd_range)
            return hmm, psr_name
        else:
            return (None, psr_name)

            psr = libstempo.tempopulsar(parfile=par, timfile=tim)
            freqs = np.linspace(-5e-7, 5e-7, 1001)
            mean_z = (np.mean(np.diff(sorted(psr.toas()))))*86400
            sigma = max(np.mean(np.diff(freqs))/(mean_z)**1.5, sigma)
            fdots = np.linspace(-2*sigma*mean_z**0.5, 2*sigma*mean_z**0.5, 11)
            noise_cov = lambda z: [[sigma**2*z**3/3, sigma**2*z**2/2],[sigma**2*z**2/2, sigma**2*z]]
            hmm = HMM.HMM.from_tempo2(par, tim, freqs, fdots, noise_cov,[], efac=efac, equad=10**equad*1e6)
            return hmm, psr_name

def do_psr(hmm):
    evs = []
    for n in range(len(hmm.zs)):
        hmm.glitches = [n]
        hmm.gen_all_obs_loglikes()
        hmm.forward()
        evs.append(hmm.evidence[-1])
    print(evs)
    if max(evs - evs[0]) > 1:
        print(f"Glitch at idx {np.argmax(evs)} with log BF {max(evs-evs[0])}")

    #return {"PSRJ": psr_name, "idx": np.argmax(evs), "log_BF": max(evs-evs[0])}
    return np.argmax(evs), max(evs-evs[0])
    #return hmm

def do_psr_matlab(hmm, psr_name, analysis_file):
    np.savetxt(f'test_kappas.dat', hmm.kappas)
    np.savetxt(f'test_zs.dat', hmm.zs)
    np.savetxt(f'test_fdots.dat', hmm.fdots)
    np.savetxt(f'test_freqs.dat', hmm.freqs)
    np.savetxt(f'test_freq_fiducial.dat', [hmm.f_fiducial])
    np.savetxt(f'test_fdot_fiducial.dat', [hmm.fd_fiducial])
    np.savetxt(f'test_fddot_fiducial.dat', [hmm.fdd_fiducial])
    cmd = f"matlab -nosplash -nodesktop -r \"do_matlab_analysis; save {analysis_file}; exit\""
    subprocess.run(cmd, shell=True)


#for par in glob.glob("*.par"):
#for par in [sys.argv[1]]:
par = sys.argv[1]

print(par)
tim = par[:-3] + 'tim'
#if not os.path.exists(par[:-4]):

#if not os.path.exists(f"{par[:-4]}/analysis.mat"):
#shutil.copyfile(par, f"{par[:-4]}/{par}")
#shutil.copyfile(tim, f"{par[:-4]}/{tim}")
#shutil.copyfile("do_analysis.py", f"{par[:-4]}/do_analysis.py")
#shutil.copyfile("do_matlab_analysis.m", f"{par[:-4]}/do_matlab_analysis.m")
mjd_range = None
if len(sys.argv) == 4:
    mjd_range = (float(sys.argv[2]), float(sys.argv[3]))

hmm, psr_name = setup_hmm(par,tim, mjd_range = mjd_range)
if mjd_range:
    do_psr_matlab(hmm, psr_name, analysis_file=f"analysis_{mjd_range[0]}_{mjd_range[1]}.mat")
else:
    do_psr_matlab(hmm, psr_name, analysis_file=f"analysis.mat")
