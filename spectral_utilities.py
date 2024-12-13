#!/usr/bin/python

import numpy as np
from scipy import ndimage
from pstats import *
from numpy import trapz,linspace,logspace,log,log10
from mstats import *
import os, datetime
import scipy.signal

def get_power_spectral_density(timeseries, rate=1., nfft=4096, noverlap=None):
    """
    Returns the power spectral density of the given timeseries.

    Parameters
    ----------
    timeseries : ndarray
        One-dimensional timeseries for which we want the PSD.
    rate : float, optional
        Sampling rate of the data, in Hz. (or 1/whatever consistent with below)
    nfft : int, optional
        Window size for sampling power spectral density, in number of samples.
    noverlap : int, optional
        Number of points to overlap between segments. Defaults to nfft/2.

    Returns
    -------
    f : ndarray
        The frequency axis corresponding to Pxx.
    Pxx : ndarray
        The power spectral density of the given timeseries.
    """
    f, Pxx = scipy.signal.welch(
        timeseries,
        fs=rate,
        nperseg=nfft,
        noverlap=noverlap,
        window=np.hanning(nfft),
        detrend=lambda x: x,  # no detrending. detrend=False should also work
        scaling='density')
    return f, Pxx


def apply_high_pass_filter(timeseries, rate=1., cutoff_frequency=0.05):
    """
    Takes in a timeseries, and returns a timeseries with a 4th-order
    butterworth high-pass filter of the given cutoff frequency (in Hz) applied.
    Rate is the frequency of the incoming timeseries, in Hz.
    OR rate and cutoff_freq can be (1/whatever) as long as they are the same
    since they always appear as a ratio
    """
    b, a = scipy.signal.butter(N=4, Wn=cutoff_frequency/rate, btype='high')
    high_passed = scipy.signal.filtfilt(b=b, a=a, x=timeseries)
    initialization_error_cutoff = round(rate/cutoff_frequency)
    return high_passed[initialization_error_cutoff:]