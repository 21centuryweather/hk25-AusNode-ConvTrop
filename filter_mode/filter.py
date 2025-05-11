import numpy as np
import pandas as pd
import xarray as xr

# For all functions defined below, TLL means 3D shaped (time, lat, lon) and TLLL means 4D shaped (time, lev, lat, lon)

def calcfft2dfTLL(x, alpha=0.1, norm="forward"):    # Calculate forward 2D FFT (i.e., from temporal space to spectral space)
    # alpha is the proportion of data to be tapered (alpha = 0.1 means 5% at the beginning and 5% at the end of data is tapered to zero)
    from scipy.signal import detrend
    from scipy.signal.windows import tukey

    x = detrend(x, axis=0)
    x = (x.T * tukey(x.shape[0], alpha=alpha)).T

    cf = np.fft.rfft2(x, axes=(2, 0), norm=norm)

    return cf

def calcfft2dbTLL(x, ntime, norm="forward"):    # Calculate backward 2D FFT (i.e., from spectral space back to temporal space)
    y = np.real(np.fft.irfft2(x, s=[x.shape[2], ntime], axes=(2, 0), norm=norm))

    return y

# Two functions below are basically the same as the two functions above but for 4D shaped array

def calcfft2dfTLLL(x, alpha=0.1, norm="forward"):
    from scipy.signal import detrend
    from scipy.signal.windows import tukey

    x = detrend(x, axis=0)
    x = (x.T * tukey(x.shape[0], alpha=alpha)).T

    cf = np.fft.rfft2(x, axes=(3, 0), norm=norm)

    return cf

def calcfft2dbTLLL(x, ntime, norm="forward"):
    y = np.real(np.fft.irfft2(x, s=[x.shape[3], ntime], axes=(3, 0), norm=norm))

    return y

def calc_frq_Nmode(k, Nmode, alpha=0.9, c=50):  # At zonal wavenumber k, calculate the frequency corresponding to the Nmode value (phase speed = frequency / wavenumber)
    a = 6.371e+06

    frq = (24 * 3600) / (2 * np.pi * a) * np.sqrt(Nmode * alpha * (1 - alpha)) * c * k

    return frq

def filter_mode_TLL(x, mode, ntime, spd, frqmin=1/120, kmax=20):    # Retain wavenumber and frequency belonging to a defined mode and set the rest to zero, by default only synoptic to planetary scale signals are considered (frequency > 1/120 cycles per day and wavenumber < 20)
    cf = calcfft2dfTLL(x.fillna(0).data)
    frq = np.fft.rfftfreq(x.data.shape[0], d=1/spd)

    cf[frq < frqmin, :, :] = 0
    cf[:, :, kmax+1:-kmax] = 0
    cf[:, :, 0] = 0

    if mode == "Moisture Mode":
        for k in range(-kmax, kmax+1):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[frq >= frq_top, :, k] = 0
    elif mode == "Mixed System":
        for k in range(-kmax, kmax+1):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[(frq >= frq_top) | (frq < frq_bottom), :, k] = 0
    elif mode == "IG Wave":
        for k in range(-kmax, kmax+1):
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            cf[frq < frq_bottom, :, k] = 0
    elif mode == "Westward Moisture Mode":
        for k in range(-kmax, 0):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[frq >= frq_top, :, k] = 0
    elif mode == "Eastward Moisture Mode":
        for k in range(1, kmax+1):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[frq >= frq_top, :, k] = 0
    elif mode == "Westward Mixed System":
        for k in range(-kmax, 0):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[(frq >= frq_top) | (frq < frq_bottom), :, k] = 0
    elif mode == "Eastward Mixed System":
        for k in range(1, kmax+1):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[(frq >= frq_top) | (frq < frq_bottom), :, k] = 0
    elif mode == "Westward IG Wave":
        for k in range(-kmax, 0):
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            cf[frq < frq_bottom, :, k] = 0
    elif mode == "Eastward IG Wave":
        for k in range(1, kmax+1):
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            cf[frq < frq_bottom, :, k] = 0

    y = calcfft2dbTLL(cf, ntime)

    y = xr.DataArray(y, dims=x.dims, coords=x.coords, attrs=dict(long_name=mode+" in "+x.attrs["long_name"], units=x.attrs["units"]))

    return y

def filter_mode_TLLL(x, mode, ntime, spd, frqmin=1/120, kmax=20):
    cf = calcfft2dfTLLL(x.fillna(0).data)
    frq = np.fft.rfftfreq(x.data.shape[0], d=1/spd)

    cf[frq < frqmin, :, :, :] = 0
    cf[:, :, :, kmax+1:-kmax] = 0
    cf[:, :, :, 0] = 0

    if mode == "Moisture Mode":
        for k in range(-kmax, kmax+1):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[frq >= frq_top, :, :, k] = 0
    elif mode == "Mixed System":
        for k in range(-kmax, kmax+1):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[(frq >= frq_top) | (frq < frq_bottom), :, :, k] = 0
    elif mode == "IG Wave":
        for k in range(-kmax, kmax+1):
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            cf[frq < frq_bottom, :, :, k] = 0
    elif mode == "Westward Moisture Mode":
        for k in range(-kmax, 0):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[frq >= frq_top, :, :, k] = 0
    elif mode == "Eastward Moisture Mode":
        for k in range(1, kmax+1):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[frq >= frq_top, :, :, k] = 0
    elif mode == "Westward Mixed System":
        for k in range(-kmax, 0):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[(frq >= frq_top) | (frq < frq_bottom), :, :, k] = 0
    elif mode == "Eastward Mixed System":
        for k in range(1, kmax+1):
            frq_top = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** -0.5)
            cf[(frq >= frq_top) | (frq < frq_bottom), :, :, k] = 0
    elif mode == "Westward IG Wave":
        for k in range(-kmax, 0):
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            cf[frq < frq_bottom, :, :, k] = 0
    elif mode == "Eastward IG Wave":
        for k in range(1, kmax+1):
            frq_bottom = calc_frq_Nmode(np.abs(k), 10 ** 0.5)
            cf[frq < frq_bottom, :, :, k] = 0

    y = calcfft2dbTLLL(cf, ntime)

    y = xr.DataArray(y, dims=x.dims, coords=x.coords, attrs=dict(long_name=mode+" in "+x.attrs["long_name"], units=x.attrs["units"]))

    return y
