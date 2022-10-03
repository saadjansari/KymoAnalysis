#!/usr/bin/env python

import os, pdb
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from .Kymograph import *
import shutil

'''
Name: breakBipolar.py
Description: Plots the pole separation of a bipolar file
'''

folds = ['wild type']
savepath = Path('/Users/saadjansari/Documents/Projects/ImageAnalysis/KymoAnalysis/Analysis/result_smoothing')
parent_path = Path('/Users/saadjansari/Documents/Projects/ImageAnalysis/KymoAnalysis/data/temp')

def smooth_data(arr, span):
    re = np.convolve(arr, np.ones(span * 2 + 1) / (span * 2 + 1), mode="same")

    # The "my_average" part: shrinks the averaging window on the side that
    # reaches beyond the data, keeps the other side the same size as given
    # by "span"
    re[0] = np.average(arr[:span])
    for i in range(1, span + 1):
        re[i] = np.average(arr[:i + span])
        re[-i] = np.average(arr[-i - span:])
    return re

if not Path.exists( savepath):
    os.mkdir( savepath)

for jj,jfold in enumerate(folds):
    print('Data: {0}'.format(jfold))
    mainpath = parent_path / jfold
    # txt files
    files2break = mainpath.glob('*txt')

    # Pole separation vs time
    print('Calculating pole separations...')
    for jj, fil in enumerate(files2break):
        # print(fil)
        kymo = Kymograph(fname=str(fil))

        if len(kymo.poles) == 2:
            fig, ax = plt.subplots()
            time = np.array( sorted( np.hstack( (kymo.poles[0].time, kymo.poles[1].time) ) )[1::10] )
            time = np.linspace(time[0], time[-1], int(np.ceil(time[-1]-time[0])))
            spindleLength = np.array( np.absolute( kymo.poles[1].ifunc(time)- kymo.poles[0].ifunc(time)) )
            slope_windows = [5,25,50]
            for slope in slope_windows:
                spindleLength_cnv = smooth_data(spindleLength, slope)
                ax.plot(time, spindleLength_cnv, label='Window = {0}'.format(slope))

            ax.plot(time, spindleLength, 'k:', lw=2,label='Original')
            ax.legend()
            ax.set(xlabel='Time (s)', ylabel=r'Pole separation ($\mu m$)')

            plt.tight_layout()
            plt.savefig(savepath / 'smoothing_{0}_{1}.pdf'.format(mainpath.stem, jj))
            plt.close(fig)

