#!/usr/bin/env python

import os, pdb
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from .Kymograph import *
import shutil
from random import sample
import seaborn as sns
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
import pickle

'''
Name: breakBipolar.py
Description: Plots the pole separation of a bipolar file
'''

parent_path = Path('/Users/saadjansari/Documents/Projects/ImageAnalysis/KymoAnalysis/data/temp')

# Strain folders
folds = ['wild type']
# folds = ['cut7-989TD,pkl1D,klp2D']
# folds = ['wild type','cut7-989TD,pkl1D,klp2D']

# savepath = Path('/Users/saadjansari/Documents/Projects/ImageAnalysis/KymoAnalysis/Analysis/result_wt')
# savepath = Path('/Users/saadjansari/Documents/Projects/ImageAnalysis/KymoAnalysis/Analysis/result_mutant')
savepath = Path('/Users/saadjansari/Documents/Projects/ImageAnalysis/KymoAnalysis/Analysis/blahblah')

# slope_window = 25
# nsamples = 10
# thresh = [0.008, 0.005]


# get_data_from_files {{{
def get_data_from_files(parent_path, folds):

    # lists
    times = []
    length = []
    velocity = []
    acceleration = []
    strain_tag = []
    file_tag = []

    # folds is list containing different strain names
    for jj,jfold in enumerate(folds):

        # txt files
        mainpath = parent_path / jfold
        files2break = mainpath.glob('*txt')

        for jfile, fil in enumerate(files2break):
            kymo = Kymograph(fname=str(fil))

            # Only do stuff if its bipolar 
            if len(kymo.poles) == 2:

                # Times
                time = np.array( sorted( np.hstack( (kymo.poles[0].time, kymo.poles[1].time) ) )[1::10] )
                time = np.linspace(time[0], time[-1], int(np.ceil(time[-1]-time[0])))
                # times.append(time)

                # Calculate spindle length, velocity, acceleration
                clen = np.absolute( kymo.poles[1].ifunc( time)- kymo.poles[0].ifunc(time))
                cvel = list( (clen[1:]-clen[:-1]) / (time[1:] - time[:-1]) )
                cvel.insert(0, cvel[0])
                cvel = np.array(cvel)
                cacc = list( (cvel[1:]-cvel[:-1]) / (time[1:] - time[:-1]) )
                cacc.insert(0, cacc[0])
                cacc = np.array(cacc)
                for jt in range(len(time)):
                    times.append(time[jt])
                    length.append(clen[jt])
                    velocity.append(cvel[jt])
                    acceleration.append(cacc[jt])
                    strain_tag.append( jfold)
                    file_tag.append( os.path.basename(kymo.label) )
    df = pd.DataFrame({'strain':strain_tag,
        'index':file_tag,
        'time':times,
        'length':length, 
        'velocity':velocity,
        'acceleration':acceleration,
        })
    return df
# }}}

# pair_plot {{{
def pair_plot(datframe, vars_compare, label=None, savePath=None, title=''):
    
    fig,ax = plt.subplots(figsize=(12,9))
    if label is None:
        sns.pairplot(datframe, vars=vars_compare,
                     plot_kws=dict(marker="+", s=50,linewidth=3, alpha=0.1),
                     diag_kind='kde',
                     palette='Dark2', height=3)
    else:
        sns.pairplot(datframe, vars=vars_compare, hue=label,
                     plot_kws=dict(marker="+", s=50,linewidth=3, alpha=0.1),
                     diag_kind='kde',
                     palette='Dark2', height=3)
    plt.tight_layout()
    plt.title(title)
    if savePath is not None:
        plt.savefig(savePath)
    plt.close()
# }}}

# Kmeans {{{
def do_KMeans(df, vars_compare, n_clusters=2, display=True, savePath=None):

    data = df[vars_compare].to_numpy()
    # scaler = StandardScaler()
    # X_std = scaler.fit_transform(data)
    X_std = data

    print('KMeans clustering: N_clusters = {}'.format(n_clusters))
    kmeans = KMeans(n_clusters=n_clusters, 
            init='k-means++', 
            max_iter=300, 
            n_init=10, 
            random_state=10)
    model = kmeans.fit(X_std)
    labels = model.predict(X_std)
    sil_score = silhouette_score(X_std,labels)
    print('Silhouette Score = {0:.3f}'.format(sil_score))
    df['label'] = labels
    df = labels_ordered(df,'length')

    if display and savePath is not None:
        pair_plot(df, vars_compare, label="label",title='kmeans',
                  savePath=savePath)
    return df, model
# }}}

# GMM {{{
def do_GaussianMixtureModel(df, vars_compare, n_clusters=2, display=True, savePath=None):

    data=df[vars_compare].to_numpy()
    scaler = StandardScaler()
    X_std = scaler.fit_transform(data)

    # define the model
    print('Gaussian Mixture Model: N_components = {}'.format(n_clusters))
    model = GaussianMixture(n_components=n_clusters).fit(X_std)
    labels = model.predict(X_std)
    sil_score = silhouette_score(X_std,labels)
    print('Silhouette Score = {0:.3f}'.format(sil_score))
    df['label'] = labels
    df = labels_ordered(df,'length')

    if display and savePath is not None:
        pair_plot(df, vars_compare, label="label",title='kmeans',
                  savePath=savePath)
    return df, model
# }}}

# LabelsOrdered {{{
def labels_ordered( df, ref_name):

    label_list_new = []

    # Get unique labels and the mean values of the reference variable
    labels = sorted(df.label.unique())
    mu = np.zeros( len(labels) )
    for jlab in range(len(labels)):
        mu[jlab] = df[df.label == labels[jlab]][ref_name].mean()
    
    # Create mapping from old_label to new
    labels_new = [x for _,x in sorted( list(zip(mu,labels)), key=lambda x:x[0])]
    # pdb.set_trace()
    mapping = {k:v for k,v in zip(labels,labels_new)}
    print(mapping)
    
    for jlab in df.label:
        label_list_new.append( mapping[jlab])
    df.label = label_list_new
    return df
# }}}

# plotClassifiedTracks {{{
def plotClassifiedTracks(df, saveParent=None, nSamples=50, model=None):

    # for each unique strain, make a plot
    # extract tracks
    strains = df.strain.unique().tolist()
    for strain in strains:
        fig,(ax0,ax1,ax2) = plt.subplots(1,3,figsize=(18,4.5), sharey=True)
        
        indices = df[df.strain == strain]['index'].unique()
        # Pick nSamples indices at random
        # indices2plot = sample(list(indices), nSamples)
        indices2plot = indices[:nSamples]
        
        # plot each track
        for ind in indices2plot:

            # get track to plot
            track = df[ (df['strain'] == strain) & (df['index'] == ind)]
            time = np.array(track.time)
            length = np.array(track.length)
            label = np.array(track.label)

            # Plot Axis 0
            ax0.plot(time, length, alpha=0.5, color='k', lw=2)

            # Plot Axis 1
            # len_group0 and len_group1 
            len_0 = length.copy()
            len_1 = length.copy()
            idx0 = np.where(label == 0)[0]
            idx1 = np.where(label == 1)[0]
            len_0[idx1] = np.nan 
            len_1[idx0] = np.nan 
            ax1.plot(time, len_0, alpha=0.5, lw=2, color='green')
            ax1.plot(time, len_1, alpha=0.5, lw=2,color='purple')

            # Plot Axis 2
            label_new = np.array(ForceLabelsOneWay( SmoothClassifiedLabels(label, span=100) ) )
            len_0 = length.copy()
            len_1 = length.copy()
            idx0 = np.where(label_new == 0)[0]
            idx1 = np.where(label_new == 1)[0]
            len_0[idx1] = np.nan 
            len_1[idx0] = np.nan 
            # pdb.set_trace()
            ax2.plot(time, len_0, alpha=0.5, lw=2,color='green')
            ax2.plot(time, len_1, alpha=0.5, lw=2,color='purple')

        # Labels/Legend Axis 0
        ax0.set(ylabel=r'Spindle Length $(\mu m)$', xlabel='Time (s)')

        # Labels/Legend Axis 1
        ax1.plot([],[], alpha=0.7, color='green', label='Group 0')
        ax1.plot([],[], alpha=0.7, color='purple', label='Group 1')
        ax1.legend()
        ax1.set(xlabel='Time (s)')
        
        # Labels/Legend Axis 2
        ax2.plot([],[], alpha=0.7, color='green', label='Group 0')
        ax2.plot([],[], alpha=0.7, color='purple', label='Group 1')
        ax1.legend()
        ax2.set(xlabel='Time (s)')

        plt.suptitle(strain)
        plt.tight_layout()

        if saveParent is not None:
            if model is None:
                plt.savefig( saveParent / 'tracks_{0}.pdf'.format(strain))
            else:
                plt.savefig( saveParent / 'tracks_{0}_{1}.pdf'.format(model,strain))
        plt.close()
# }}}

# SmoothClassifiedLabels {{{
def SmoothClassifiedLabels(label, span=100): 

    # smooth_data {{{
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
    # }}}

    # Smoothed Labels
    label_new = np.where(np.array( smooth_data( label,min([span, int(len(label)/2)]))) >= 0.5, 1, 0)

    # Once 1, always 1
    # label_perm = [max(label_new[:1+jj]) for jj in range(len(label_new))]
    return label_new
# }}}

# ForceLabelsOneWay {{{
def ForceLabelsOneWay( label):
    labels = [np.max(label[:1+idx]) for idx in range(len(label))]
    return np.array(labels)
# }}}

if not Path.exists( savepath):
    os.mkdir( savepath)

# Load data into dataframe
df = get_data_from_files(parent_path, folds)
names = ['velocity']

# Display (pre clustering)
pair_plot(df, names, savePath=savepath/'features_grid_raw.png')

# Kmeans
df_kmean, model_kmean = do_KMeans(df.copy(),names, savePath=savepath/'features_grid_kmeans.png')
print(df_kmean.groupby('label').mean() )

plotClassifiedTracks(df_kmean, model='kmeans',saveParent=savepath)

# Save model
with open(parent_path / 'kmeans.pickle', 'wb') as f:
    pickle.dump(model_kmean, f)

# GMM
# df_gmm, model_gmm = do_GaussianMixtureModel(df.copy(),names, savePath=savepath/'features_grid_gmm.png')
# print(df_gmm.groupby('label').mean() )
# plotClassifiedTracks(df_gmm, model='gmm',saveParent=savepath)


        # if vel_thresh[0]==1:
            # anaphase_time = 'Always'
        # elif vel_thresh[0]==0 and vel_thresh[-1]==1:
            # anaphase_time = timelist[ np.where(np.array(vel_thresh)>0.5)[0][0] ]
        # elif vel_thresh[-1]==0:
            # anaphase_time = 'Never'
    
        # # anaphase_time = timelist[ np.where(np.array(vel_thresh)>0.5)[0][0] ]
        # print( '{0} --> Anaphase B Transition = {1} sec'.format( files2break[idx].stem,anaphase_time))


