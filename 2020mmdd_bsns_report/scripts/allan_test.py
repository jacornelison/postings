import allantools as at
import numpy as np
import matplotlib.pyplot as plt


def gennoise(beta):
    """
    creates noise with a certain power spectrum.
    input f the slope of the spectra
    for technical definitin of noise colors, see wikipedia article
    https://en.wikipedia.org/wiki/Colors_of_noise
    This function was originally gennoise.c http://paulbourke.net/fractals/noise/

    """
    np.random.seed(42)
    beta = int(beta)
    N = int(1e5)
    N2 = int(N/2+1)

    x = np.arange(0, N2)
    mag = x ** (beta / 2.) * np.random.randn(N2)
    pha = 2 * np.pi * np.random.rand(N2)
    real = mag * np.cos(pha);
    imag = mag * np.sin(pha);
    real[0] = 0;
    imag[0] = 0;
    # real = concatenate((real,real[::-1]))
    # imag = concatenate((imag,-imag[::-1]))
    # set real and imag to 0 for f = 0
    # imag[N/2] = 0;

    c = real + 1j * imag
    # plot(abs(c**2))
    # xscale('log')
    # yscale('log')
    tod = np.fft.irfft(c)

    return tod

p = gennoise(-3)

[tau, ad, ade, ns] = at.adev(p,2,'freq','decade')

plt.figure()
ax = plt.gca()
plt.errorbar(tau,np.array(ad),np.array(ade),fmt='--o')
ax.set_yscale('log')
ax.set_xscale('log')



a = at.Dataset(p,2,'freq','decade')
a.compute('adev')
b = at.Plot()
b.plot(a,errorbars=True,grid=True)


plt.figure()
plt.plot(p)

plt.show()




