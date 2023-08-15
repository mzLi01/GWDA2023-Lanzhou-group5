# %%
import os
import sys
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

import pyswarms as ps
from bilby.gw.utils import noise_weighted_inner_product, optimal_snr_squared
from bilby.gw.detector import PowerSpectralDensity, InterferometerStrainData as StrainData
from bilby.core.likelihood import Likelihood
from bilby.core.utils import nfft

if len(sys.argv)==1:
    outdir='.'
else:
    outdir=sys.argv[1]
if not os.path.exists(outdir):
    os.mkdir(outdir)

noise = np.loadtxt('../../MDC/traindata.txt')
noise_fs = 1024
noise_psd_f, noise_psd_value = signal.welch(noise, noise_fs, nperseg=128)

plt.loglog(noise_psd_f, noise_psd_value)
plt.xlabel('f')
plt.ylabel('PSD')
plt.title('Noise PSD')
plt.savefig(os.path.join(outdir, 'noise_psd.png'), bbox_inches='tight')

# %%

def quadratic_chirp(t, A, a1, a2, a3, scale=1):
    phi = sum([ai*t**(i+1) for i, ai in enumerate([a1, a2, a3])])
    s = A*scale * np.sin(2*np.pi*phi)
    return s


class NoisedQuadracticChirp(Likelihood):
    def __init__(self, strain: StrainData, noise_psd: PowerSpectralDensity, scale=1, parameters=None):
        super().__init__(parameters)
        self.strain = strain
        self.noise_psd = noise_psd
        self.scale = scale

    @property
    def noise_psd_array(self):
        return self.noise_psd.get_power_spectral_density_array(
            self.strain.frequency_array) * self.scale**2 * self.strain.window_factor

    def frequency_domain_signal(self, parameters):
        time_domain = quadratic_chirp(
            self.strain.time_array, scale=self.scale, **parameters)
        freq_domain, _ = nfft(time_domain, self.strain.sampling_frequency)
        return freq_domain

    def inject_signal(self, parameters):
        signal = self.frequency_domain_signal(parameters)
        self.strain.frequency_domain_strain += signal

    def noise_log_likelihood(self):
        log_l = -noise_weighted_inner_product(
            self.strain.frequency_domain_strain, self.strain.frequency_domain_strain,
            self.noise_psd_array, self.strain.duration)/2
        return float(np.real(log_l))

    def log_likelihood_ratio(self):
        signal = self.frequency_domain_signal(self.parameters)
        inner = noise_weighted_inner_product(
            self.strain.frequency_domain_strain, signal,
            self.noise_psd_array, self.strain.duration)
        optimal_snr2 = optimal_snr_squared(
            signal, self.noise_psd_array, self.strain.duration)
        return np.real(inner-optimal_snr2/2)

    def log_likelihood(self):
        return self.noise_log_likelihood()+self.log_likelihood_ratio()


# %%
data = np.loadtxt('../../MDC/analysisData.txt')
fs = 1024
duration = data.shape[0]/fs

inject_para = {'A': 1e-22, 'a1': 50, 'a2': 10, 'a3': 3}

noise_psd = PowerSpectralDensity(
    frequency_array=noise_psd_f, psd_array=noise_psd_value)

strain = StrainData()
# strain.set_from_time_domain_strain(data*1e20, sampling_frequency=fs, duration=duration)
strain.set_from_power_spectral_density(
    noise_psd, sampling_frequency=fs, duration=duration)

likelihood = NoisedQuadracticChirp(strain, noise_psd)
likelihood.inject_signal(inject_para)

likelihood.parameters = inject_para
print(f'inject llr: {likelihood.log_likelihood_ratio():.2f}')

# %%
fig, ax = plt.subplots(1, 3, figsize=(20, 5))
fig.subplots_adjust(wspace=0.1)
for axi, x, title in zip(ax, [likelihood.strain.time_domain_strain, noise, data], ['Inject Signal', 'Noise', 'Data']):
    f, t, Sxx = signal.spectrogram(x, likelihood.strain.sampling_frequency)
    color = axi.pcolormesh(t, f, np.log(Sxx), shading='gouraud')
    axi.set_ylabel('Frequency [Hz]')
    axi.set_xlabel('Time [sec]')
    fig.colorbar(color, ax=axi)
    axi.set_title(title)

plt.savefig(os.path.join(outdir, 'inject_specgram.png'), bbox_inches='tight')

# %%

swarm_size = 40
dim = 4
epsilon = 1.0
options = {'c1': 2, 'c2': 2, 'w': 0.9}

bounds = np.array([(0, 1e-21), (40, 100), (1, 50), (1, 15)]).T


def opt_func(x):
    llr = np.zeros(x.shape[0])
    for i, xi in enumerate(x):
        likelihood.parameters = {
            key: value for key, value in zip(['A', 'a1', 'a2', 'a3'], xi)}
        llr[i] = likelihood.log_likelihood_ratio()
    return -llr


optimizer = ps.single.GlobalBestPSO(
    n_particles=swarm_size, dimensions=dim, bounds=bounds, options=options)

cost, joint_vars = optimizer.optimize(opt_func, iters=2000)



fit_parameters = {key: value for key, value in zip(
    ['A', 'a1', 'a2', 'a3'], joint_vars)}
likelihood.parameters = fit_parameters

signal_fit = likelihood.frequency_domain_signal(fit_parameters)
snr2 = np.real(optimal_snr_squared(
    signal_fit, likelihood.noise_psd_array, likelihood.strain.duration))

print(f'PSO estimated parameters: A={joint_vars[0]:.2e}, [a1,a2,a3]={joint_vars[1:]}\n'
      f'Log likelihood ratio={likelihood.log_likelihood_ratio():.2f},SNR^2={snr2:.2e}')

data = likelihood.strain.time_domain_strain.copy()

signal_time = quadratic_chirp(likelihood.strain.time_array, **fit_parameters)
likelihood.strain.frequency_domain_strain -= likelihood.frequency_domain_signal(
    fit_parameters)

fig, ax = plt.subplots(1, 3, figsize=(20, 5), sharey=True)
fig.subplots_adjust(wspace=0.1)
for axi, x, title in zip(ax, [data, signal_time, likelihood.strain.time_domain_strain], ['Data', 'Fit Signal', 'Residual']):
    f, t, Sxx = signal.spectrogram(x, likelihood.strain.sampling_frequency)

    real_freq = sum([inject_para[ai]*t**(i+1)
                     for i, ai in enumerate(['a1', 'a2', 'a3'])])
    axi.plot(t, real_freq, 'r', color='r')

    color = axi.pcolormesh(t, f, np.log(Sxx), shading='gouraud')
    axi.set_ylabel('Frequency [Hz]')
    axi.set_xlabel('Time [sec]')
    fig.colorbar(color, ax=axi)
    axi.set_title(title)

plt.savefig(os.path.join(outdir, 'fit_result.png'), bbox_inches='tight')

# %%
cost_history = optimizer.cost_history
plt.figure()
plt.plot(-np.array(cost_history))
plt.xlabel('Iteration')
plt.ylabel('Log Liklihood Ratio')

plt.savefig(os.path.join(outdir, 'pso_convergence.png'), bbox_inches='tight')
# %%
