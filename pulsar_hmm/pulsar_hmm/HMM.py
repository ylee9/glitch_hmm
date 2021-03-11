import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import scipy.misc as misc
import scipy.special as special
from skimage.util.shape import view_as_windows
import itertools

import libstempo

class HMM:

    def __init__(self, toas, toa_errs, freqs, fdots, noise_cov, glitches, f_fiducial=0, fd_fiducial=0, fdd_fiducial=0):
        self.zs = np.diff(sorted(toas))*86400
        self.noise_cov = noise_cov
        self.glitches = glitches
        self.num_timesteps = len(self.zs)
        self.freqs = freqs
        self.df = np.diff(freqs)[0]
        self.fdots = fdots
        self.dfd = np.diff(fdots)[0]
        self.f_fiducial = f_fiducial
        self.fd_fiducial = fd_fiducial
        self.fdd_fiducial = fdd_fiducial

        err_df = np.diff(freqs)[0]*self.zs
        err_dfd = 0.5*np.diff(fdots)[0]*self.zs**2
        # err_sigma_toas = 5*np.sqrt((f_fiducial + np.mean(freqs))**2*(toa_errs[:-1]**2 + toa_errs[1:]**2))
        err_sigma_toas = np.sqrt((f_fiducial + np.mean(freqs))**2*(toa_errs[:-1]**2 + toa_errs[1:]**2))

        #err_df = 2*np.diff(freqs)[0]*self.zs
        #err_dfd = 2*0.5*np.diff(fdots)[0]*self.zs**2
        #err_sigma_toas = 5*np.sqrt((f_fiducial + np.mean(freqs))**2*(toa_errs[:-1]**2 + toa_errs[1:]**2))
        self.kappas = 1/4/np.pi**2/(err_df**2 + err_dfd**2 + err_sigma_toas**2)

    @staticmethod
    def from_tempo2(parfile, timfile, freqs, fdots, noise_cov, glitches, efac=1, equad=0, min_toa_gap=0, mjd_range=None):
        psr = libstempo.tempopulsar(parfile=parfile, timfile=timfile)
        psr_toas = sorted(psr.toas())

        try:
            phase_jumps = psr.flagvals('phaseJ').astype(np.float)
            psr_toas += phase_jumps/86400
        except Exception as e:
            print(e)

        psr_toaerrs = [x for _,x in sorted(zip(psr_toas, psr.toaerrs))]
        toas = np.array([psr_toas[0]])
        toaerrs = np.array([np.sqrt(efac**2*psr_toaerrs[0]**2 + equad**2)])

        for i in range(1, len(psr.toaerrs)):
            if psr_toas[i] - toas[-1] >= min_toa_gap:
                toas = np.append(toas, psr_toas[i])
                toaerrs = np.append(toaerrs, np.sqrt(efac**2*psr_toaerrs[i]**2 + equad**2))

        #toas = psr.toas()
        #toaerrs = np.sqrt(efac**2*psr.toaerrs**2 + equad**2)
        print(toaerrs)
        if mjd_range is not None:
            filt = [toa >= mjd_range[0] and toa <= mjd_range[1] for toa in toas]
            toas = list(itertools.compress(toas,filt))
            toaerrs = np.array(list(itertools.compress(toaerrs,filt)))
        pepoch = psr['PEPOCH'].val
        f_fiducial = psr['F0'].val + psr['F1'].val * (min(toas) - pepoch)*86400 + 0.5*psr['F2'].val * ((min(toas) - pepoch)*86400)**2
        fd_fiducial = psr['F1'].val + psr['F2'].val * (min(toas) - pepoch)*86400

        return HMM(toas, toaerrs*1e-6, freqs, fdots, noise_cov, glitches, f_fiducial, fd_fiducial, psr['F2'].val)


    def fokker_planck_pdf(self, freqs, fdots, z):
        cov_matrix = self.noise_cov(z)
        cov_matrix[0][0] /= self.df**2
        cov_matrix[1][0] /= self.df*self.dfd
        cov_matrix[0][1] /= self.df*self.dfd
        cov_matrix[1][1] /= self.dfd**2
        freq_size = len(freqs)
        fdot_size = len(fdots)
        grid_f, grid_fd = np.meshgrid(np.linspace(-(freq_size-1)/2, (freq_size-1)/2, freq_size), np.linspace(-(fdot_size-1)/2, (fdot_size-1)/2, fdot_size))
        pos = np.empty(grid_f.shape + (2,))
        pos[:,:,0] = grid_f
        pos[:,:,1] = grid_fd
        try:
            rand_var = stats.multivariate_normal(mean=[0,0], cov=cov_matrix)
            pdf = rand_var.logpdf(pos)
            return pdf - special.logsumexp(pdf)
        except Exception as e:
            pdf = np.ones(grid_f.shape)*-np.inf
            pdf[int((fdot_size-1)/2)][int((freq_size-1)/2)] = 0
            return pdf

    def gen_trans_matrix_block(self, z):
        freq_size = int(np.max([3, np.min([len(self.freqs), (2*3 + 1)*np.sqrt(self.noise_cov(z)[0][0])/self.df])]))
        if freq_size % 2 == 0: # Pad freq_size to be odd so that we have a good center
            freq_size += 1

        fdot_size = len(self.fdots)

        fs = np.linspace(-self.df*(freq_size-1)/2, self.df*(freq_size-1)/2, freq_size)
        fds = np.linspace(-self.dfd*(fdot_size-1), self.dfd*(fdot_size-1), fdot_size*2-1)
        return self.fokker_planck_pdf(fs[:, None], fds[:, None], z)

    def step(self, prev_loglikes, z, glitch, direction='fwd'):
        new_loglikes = np.zeros(np.shape(prev_loglikes))
        if not glitch:
            trans_matrix_block = np.flipud(self.gen_trans_matrix_block(z))
            window_shape = (len(self.fdots), np.shape(trans_matrix_block)[1])
            new_loglikes_unsummed = np.zeros((prev_loglikes.shape[0]*trans_matrix_block.shape[1], prev_loglikes.shape[0], prev_loglikes.shape[1]))
            for fdot_idx in range(np.shape(prev_loglikes)[0]):
                trans_matrix_block_fdot_lower = -fdot_idx + len(self.fdots) - 1
                trans_matrix_block_fdot_upper = -fdot_idx + 2*len(self.fdots) - 1
                trans_matrix_block_selected_fdots = trans_matrix_block[trans_matrix_block_fdot_lower:trans_matrix_block_fdot_upper, :]
                trans_matrix_block_replicated = np.repeat(trans_matrix_block_selected_fdots[:,:, np.newaxis], np.shape(prev_loglikes)[1], axis=2)

                if direction == 'fwd':
                    freq_idx_offset = int(np.rint(self.fdots[fdot_idx]*z/self.df))
                elif direction == 'bwd':
                    freq_idx_offset = int(np.rint(self.fdots[fdot_idx]*z/self.df)) * -1
                else:
                    raise ValueError('Invalid value for direction')

                rolled_loglikes = np.roll(prev_loglikes, freq_idx_offset, axis=1)
                #if freq_idx_offset < 0:
                #    rolled_loglikes[:, freq_idx_offset:] = -np.inf
                #elif freq_idx_offset > 0:
                #    rolled_loglikes[:, :(freq_idx_offset)] = -np.inf
                padding_width = int((window_shape[1] - 1)/2)
                rolled_loglikes = np.pad(rolled_loglikes, ((0,0), (padding_width, padding_width)), constant_values=(-np.inf, -np.inf), mode='constant')

                rolled_loglikes_win = np.moveaxis(view_as_windows(rolled_loglikes, window_shape).squeeze(), [0,1,2], [2,0,1])
                #summed = special.logsumexp(rolled_loglikes_win + trans_matrix_block_replicated, axis=(0,1))
                summed = rolled_loglikes_win + trans_matrix_block_replicated
                new_loglikes_unsummed[:, fdot_idx, :] = summed.reshape((trans_matrix_block_selected_fdots.size, prev_loglikes.shape[1]))
                #print(summed.shape)
                #summed = special.logsumexp(summed, axis=0)
                #print(summed.shape)

            new_loglikes = np.logaddexp.reduce(new_loglikes_unsummed, axis=0)

            for fdot_idx in range(np.shape(prev_loglikes)[0]):
                if direction == 'fwd':
                    freq_idx_offset = int(np.rint(self.fdots[fdot_idx]*z/self.df))
                elif direction == 'bwd':
                    freq_idx_offset = int(np.rint(self.fdots[fdot_idx]*z/self.df)) * -1
                else:
                    raise ValueError('Invalid value for direction')

                if freq_idx_offset < 0:
                    new_loglikes[fdot_idx, freq_idx_offset:] = -np.inf
                elif freq_idx_offset > 0:
                    new_loglikes[fdot_idx, :(freq_idx_offset)] = -np.inf

        else:
            all_glitch_loglikes = np.ones((prev_loglikes.shape[0], prev_loglikes.shape[1], prev_loglikes.shape[0]*prev_loglikes.shape[1] + 1))*-np.inf
            all_glitch_loglikes[:,:,0] = prev_loglikes
            for i in range(prev_loglikes.size):
                fdot, f = np.unravel_index(i, prev_loglikes.shape)
                if direction == 'fwd':
                    freq_idx_offset = int(np.rint(self.fdots[fdot]*z/self.df))
                    freq_min = np.clip(f + freq_idx_offset, 1, prev_loglikes.shape[1])
                    all_glitch_loglikes[:, freq_min:, i] = prev_loglikes[fdot, f] - np.log((prev_loglikes.shape[1]-freq_min)*prev_loglikes.shape[0])
                elif direction == 'bwd':
                    freq_idx_offset = int(np.rint(self.fdots[fdot]*z/self.df)) * -1
                    freq_max = np.clip(f + freq_idx_offset, 1, prev_loglikes.shape[1])
                    all_glitch_loglikes[:, :freq_max, i] = prev_loglikes[fdot, f] #- np.log((freq_max)*prev_loglikes.shape[0])
                else:
                    raise ValueError('Invalid value for direction')


            new_loglikes = special.logsumexp(all_glitch_loglikes, axis=2)

            #for fdot_idx in range(prev_loglikes.shape[0]):
            #    for f_idx in range(1, prev_loglikes.shape[1]):
            #
            #        for fdot_prev_idx in range(prev_loglikes.shape[0]):
            #            if direction == 'fwd':
            #                freq_idx_offset = int(np.rint(self.fdots[fdot_idx]*z/self.df))
            #            elif direction == 'bwd':
            #                freq_idx_offset = int(np.rint(self.fdots[fdot_idx]*z/self.df)) * -1
            #            else:
            #                raise ValueError('Invalid value for direction')
            #
            #            freq_max = np.clip(f_idx + freq_idx_offset, 1, prev_loglikes.shape[1])
            #            new_loglikes[fdot_idx][f_idx] = special.logsumexp(prev_loglikes[:, :f_idx].flatten()) #- np.log((prev_loglikes.shape[0] + 1)*freq_max)
        return new_loglikes

    def obs_loglikes(self, z, freqs, fdots, kappa, running_f, running_fd):
        obs_loglikes = np.zeros((len(self.fdots), len(self.freqs)))
        fs, fds = np.meshgrid(freqs, fdots)
        phase_fiducial = 2*np.pi*(z*running_f - 0.5*running_fd*z**2)
        phases = 2*np.pi*(z*fs - 0.5*z**2*fds) + phase_fiducial
        log_besseli = np.log(special.i0(kappa.astype(np.float64)))
        if np.isinf(log_besseli): # For large argument, I_0(x) ~ exp(x)/sqrt(2pi*x)
            log_besseli = kappa - 0.5*np.log(2*np.pi*kappa)

        obs_loglikes = kappa*np.cos(phases) - log_besseli - np.log(2*np.pi)
        return obs_loglikes

    def gen_all_obs_loglikes(self):
        self.all_obs_loglikes = np.zeros((len(self.zs), len(self.fdots), len(self.freqs)))
        running_f = self.f_fiducial
        running_fd = self.fd_fiducial
        for n in range(len(self.zs)):
            running_f += self.fd_fiducial*self.zs[n]
            self.all_obs_loglikes[n, :, :] = self.obs_loglikes(self.zs[n], self.freqs, self.fdots, self.kappas[n], running_f, running_fd)

    def forward(self):
        self.forward_loglikes = np.zeros((len(self.zs)+1, len(self.fdots), len(self.freqs)))
        self.evidence = np.zeros(len(self.zs))
        running_f = self.f_fiducial
        running_fd = self.fd_fiducial
        for n in range(1, len(self.zs)):
            self.forward_loglikes[n, :, :] = self.step(self.forward_loglikes[n-1,:,:], self.zs[n], n in self.glitches, direction='fwd')
            self.forward_loglikes[n,:,:] += self.all_obs_loglikes[n, :, :]
            self.evidence[n] = special.logsumexp(self.forward_loglikes[n,:,:].flatten())


    def backward(self):
        self.backward_loglikes = np.zeros((len(self.zs)+1, len(self.fdots), len(self.freqs)))

        for n in range(len(self.zs)-1, 0, -1):
            new_loglikes = self.all_obs_loglikes[n,:,:] + self.backward_loglikes[n+1,:,:]
            self.backward_loglikes[n, :, :] = self.step(new_loglikes, self.zs[n], n in self.glitches, direction='bwd')

    def fw_bw(self):
        self.gen_all_obs_loglikes()
        self.forward()
        self.backward()

        self.combined_loglikes = self.forward_loglikes[:,:,:] + self.backward_loglikes[:,:,:]
        for n in range(len(self.zs)):
            self.combined_loglikes[n,:,:] -= special.logsumexp(self.combined_loglikes[n,:,:].flatten())

        self.gen_path()
        g = np.zeros((len(self.zs), len(self.freqs)))
        for n in range(len(self.zs)):
            g[n, :] = self.forward_loglikes[n, 0, :]

        np.savetxt('g.dat', g)

    def gen_path(self):
        path = []
        for n in range(1, len(self.zs)):
            (fd, f) = np.unravel_index(self.combined_loglikes[n,:,:].argmax(), self.combined_loglikes[n,:,:].shape)
            path.append((fd,f))
        self.path = path

    def get_residuals(self):
        running_f = self.f_fiducial
        running_fd = self.fd_fiducial
        residuals = []
        for n in range(0,len(self.zs)):
            running_f += running_fd*self.zs[n]
            phase_fiducial = running_f*self.zs[n] - 0.5*running_fd*self.zs[n]**2
            phase = self.freqs[self.path[n][1]]*self.zs[n] - 0.5*self.fdots[self.path[n][0]]*self.zs[n]**2 + phase_fiducial

            residuals.append(phase - np.round(phase))

        return residuals
