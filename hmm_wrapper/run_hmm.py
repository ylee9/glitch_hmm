#!/usr/bin/env python

import pulsar_hmm.HMM as HMM
import libstempo
import numpy as np
import configparser
import shutil
import sys
import argparse
import tempfile
import subprocess
import matplotlib.pyplot as plt
import matplotlib
plt.rcParams['figure.figsize'] = [12.0, 9.0]
from matplotlib import rc
#matplotlib.rcParams["text.latex.preamble"] += r'\usepackage[dvips]{graphicx}\usepackage{amsmath}\usepackage{amssymb}'
#rc('text', usetex=True)
rc('font', size=24.0)
rc('font',**{'family':'serif'})

def setup_hmm(par, tim, config):
    psr = libstempo.tempopulsar(parfile=par, timfile=tim)

    freq_min = float(config['doi']['freq_min'])
    freq_max = float(config['doi']['freq_max'])
    dfreq = float(config['doi']['dfreq'])
    fdot_min = float(config['doi']['fdot_min'])
    fdot_max = float(config['doi']['fdot_max'])
    dfdot = float(config['doi']['dfdot'])

    freqs = np.arange(freq_min, freq_max, dfreq)
    fdots = np.arange(fdot_min, fdot_max, dfdot)
    if len(fdots) % 2 == 0:
        fdots = np.linspace(fdot_min, fdot_max, len(fdots)+1)
        dfdot = np.diff(fdots)[0]
        print(f'Odd number of fdot bins required, adjusting dfdot. New dfdot = {dfdot}')

    min_toa_gap = 0
    mjd_range = None
    if 'toas' in config:
        if 'min_toa_gap' in config['toas']:
            min_toa_gap = float(config['toas']['min_toa_gap'])
        if 'mjd_min' in config['toas'] and 'mjd_max' in config['toas']:
            mjd_range = [float(config['toas']['mjd_min']), float(config['toas']['mjd_max'])]

    psr_name = None
    red_idx = None
    red_amp = None
    red_comp = None
    efac = 1
    equad = -100

    with open(par,'r') as f:
        for line in f.readlines():
            split = line.split()
            if split[0].strip() == 'PSRJ':
                psr_name = split[-1].strip()
                print(psr_name)
            elif split[0].strip() == 'TNRedAmp':
                red_amp = 10**float(split[-1].strip())
            elif split[0].strip() == 'TNRedGam':
                red_idx = float(split[-1].strip())
            elif split[0].strip() == 'TNRedC':
                red_comp = int(split[-1].strip())
            elif split[0].strip().lower() == 'tnglobalef':
                efac = float(split[-1].strip())
            elif split[0].strip().lower() == 'tnglobaleq':
                equad = float(split[-1].strip())

    tn_params = {'red_amp': red_amp, 'red_idx': red_idx, 'red_comp': red_comp, 'efac': efac, 'equad': equad}

    hmm = HMM.HMM.from_tempo2(par, tim, freqs, fdots, None, [], min_toa_gap = min_toa_gap, mjd_range=mjd_range)

    if 'tn' in config:
        if 'sigma' in config['tn']:
            print('Overriding default value of sigma')
            sigma = float(config['tn']['sigma'])
    else:
        mean_z = np.mean(hmm.zs)
        sigma = max(dfdot/np.sqrt(mean_z), 1e-21)

    return hmm, psr, sigma, tn_params

def save_hmm_files(hmm, config):
    working_prefix = config['matlab']['working_prefix'] if 'working_prefix' in config['matlab'] else None
    if not working_prefix:
        tmp_dir = tempfile.mkdtemp()
        working_prefix = tmp_dir + "/"

    np.savetxt(f"{working_prefix}kappas.dat", hmm.kappas)
    np.savetxt(f"{working_prefix}freqs.dat", hmm.freqs)
    np.savetxt(f"{working_prefix}fdots.dat", hmm.fdots)
    np.savetxt(f"{working_prefix}zs.dat", hmm.zs)

    return working_prefix

def do_psr(hmm, sigma, config, extra_matlab_cmd=None, ul=False):
    matlab_wrapper = config['matlab']['matlab_wrapper'] if not ul else config['matlab']['matlab_ul_wrapper']
    matlab_path = config['matlab']['matlab_path']
    out_prefix = config['out']['out_prefix']
    working_prefix = config['matlab']['working_prefix'] if 'working_prefix' in config['matlab'] else None


    matlab_cmd = "global f_fiducial fd_fiducial fdd_fiducial sigma kappa_per_toa;"
    matlab_cmd += f"f_fiducial = {hmm.f_fiducial};"
    matlab_cmd += f"fd_fiducial = {-hmm.fd_fiducial};"
    matlab_cmd += f"fdd_fiducial = {hmm.fdd_fiducial};"
    matlab_cmd += f"sigma = {sigma};"
    matlab_cmd += f"out_prefix = '{out_prefix}';"
    matlab_cmd += f"working_prefix = '{working_prefix}';"
    matlab_cmd += f"matlab_path = '{matlab_path}';"
    if extra_matlab_cmd:
        matlab_cmd += extra_matlab_cmd
    matlab_cmd += f"run('{matlab_wrapper}');"
    print(matlab_cmd)
    cmd = f"matlab -nosplash -nodesktop -r \"{matlab_cmd} exit\""
    subprocess.run(cmd, shell=True)
    return

def make_plots(hmm, out_prefix):
    freqs = hmm.freqs
    fdots = hmm.fdots
    zs = hmm.zs
    f_posterior = np.loadtxt(f"{out_prefix}f_posterior.dat")
    freq_ticklabels = ["{:.2E}".format(freqs[int(idx)] if np.abs(freqs[int(idx)]) > 1e-25 else 0) for idx in np.floor(np.linspace(0, len(freqs)-1, 11))]
    plt.yticks(np.floor(np.linspace(0, len(freqs)-1, 11)), freq_ticklabels)
    plt.imshow(f_posterior, aspect='auto')
    cbar = plt.colorbar()
    cbar.set_label(r"$\ln[\gamma_f(t_n)]$")
    plt.xlabel('ToA index')
    plt.ylabel('Frequency bin')
    plt.gca().invert_yaxis()
    plt.savefig(f"{out_prefix}f_posterior.pdf", bbox_inches='tight')
    plt.clf()

    fdot_posterior = np.loadtxt(f"{out_prefix}fdot_posterior.dat")
    plt.imshow(np.exp(fdot_posterior), aspect='auto')
    cbar = plt.colorbar()
    cbar.set_label(r"$\gamma_{\dot{f}}(t_n)$")
    fdot_ticklabels = ["{:.2E}".format(fdots[int(idx)] if np.abs(fdots[int(idx)]) > 1e-25 else 0) for idx in np.floor(np.linspace(0, len(fdots)-1, 5))]
    plt.yticks(np.floor(np.linspace(0, len(fdots)-1, 5)), fdot_ticklabels)
    plt.xlabel('ToA index')
    plt.ylabel('Frequency derivative bin')
    plt.gca().invert_yaxis()
    plt.savefig(f"{out_prefix}fdot_posterior.pdf", bbox_inches='tight')
    plt.clf()

    f_path = np.loadtxt(f"{out_prefix}f_path.dat", dtype=np.float)
    plt.plot(np.cumsum(zs)/86400, f_path)
    plt.ylim(min(freqs), max(freqs))
    plt.xlabel('Days since first ToA')
    plt.ylabel('Frequency (Hz)')
    plt.savefig(f"{out_prefix}f_path.pdf", bbox_inches='tight')
    plt.clf()

    fdot_path = np.loadtxt(f"{out_prefix}fdot_path.dat", dtype=np.float)
    print(fdot_path)
    plt.plot(np.cumsum(zs)/86400, fdot_path)
    plt.ylim(min(fdots), max(fdots))
    plt.xlabel('Days since first ToA')
    plt.ylabel('Frequency derivative (Hz/s)')
    plt.savefig(f"{out_prefix}fdot_path.pdf", bbox_inches='tight')
    plt.clf()

    bfs = np.atleast_2d(np.loadtxt(f"{out_prefix}bfs.dat"))
    #num_models = int(len(bfs)/len(zs))
    #bfs = np.transpose(np.reshape(bfs, (len(zs), num_models)))
    print(bfs.shape)
    for i in range(0, bfs.shape[0]):
        plt.plot(np.cumsum(zs)/86400, bfs[i, :])
        plt.xlabel('Days since first ToA')
        plt.ylabel('ln Bayes factor')
        plt.tight_layout()
        plt.savefig(f"{out_prefix}bfs_{i}.pdf", bbox_inches='tight')
        plt.clf()

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--par', '-P', type=str, help='The .par file to be used', required=True)
    parser.add_argument('--tim', '-T', type=str, help='The .tim file to be used', required=True)
    parser.add_argument('--ini', '-I', type=str, help='The .ini file containing HMM parameters to be used', required=True)

    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.ini)
    hmm, psr, sigma, tn_params = setup_hmm(args.par, args.tim, config)
    if not 'working_prefix' in config['matlab']:
        delete_working = True
    else:
        delete_working = False
    working_prefix = save_hmm_files(hmm, config)
    config['matlab']['working_prefix'] = working_prefix
    do_psr(hmm, sigma, config)

    make_plots(hmm, config['out']['out_prefix'])

    if delete_working:
        shutil.rmtree(config['matlab']['working_prefix'])
