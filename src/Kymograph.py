#!/usr/bin/env python

import os, pdb
import math
import numpy as np
from scipy import interpolate
from .node_graph import Graph
import matplotlib.pyplot as plt

from .Features import *
from .ReadFiles import *

'''
Name: Kymograph.py
Description: Parses general, poles and feature information for a single kymograph and stores the data accordingly 
'''

class Kymograph:
    def __init__(self, fname='example.txt'):
        
        # get label without '.txt'
        self.label = fname[:-9] 

        # Read file information 
        self.general,self.poles,self.tracks = ReadTxt( fname)

        # Remove Bad tracks
        self.RemoveBadTracks()

        # Merge tracks whose ends are close enough
        self.MergeTracks(self.tracks)

        # Order track poles
        for track in self.tracks:
            track.OrderPoles()

        # Trim tracks based on kmeans label 
        # self.TrimTracksKmeansLabel()

    def RemoveBadTracks(self):
        # Remove tracks that go backward in time

        # Find bad tracks
        bad_tracks = []
        for track in self.tracks:
            if not track.CheckViability():
                bad_tracks += [track]

        if len( bad_tracks) != 0:
            print('Found some bad tracks')

        # Remove bad tracks
        for track in bad_tracks:
            self.tracks.remove( track)

    def TrimBasedOnTime(self, time_keep=[-1,-1]):
        
        # Trim poles 
        poles_new = []
        for pole in self.poles:
            trimmed = pole.TrimBasedOnTime(time_keep)
            if trimmed is not np.nan and trimmed is not None:
                poles_new.append(trimmed)
        # print(poles_new)
        # if self.label == '/Users/saadjansari/Documents/Projects/ImageAnalysis/KymoAnalysis/data/bipolar/wild type/MAX_1032_100msR_50msG_7Z_004_cell A_KYMOGRAPH':
            # pdb.set_trace()
        self.poles= poles_new 

        # Trim tracks
        tracks_new = []
        for track in self.tracks:
            trimmed = track.TrimBasedOnTime(time_keep)
            if trimmed is not np.nan and trimmed is not None:
                # trimmed.poles = poles_new
                tracks_new.append(trimmed)
        # if self.label == '/Users/saadjansari/Documents/Projects/ImageAnalysis/KymoAnalysis/data/bipolar/wild type/MAX_1032_100msR_50msG_7Z_004_cell A_KYMOGRAPH':
            # pdb.set_trace()
        # print(tracks_new)
        self.tracks = tracks_new

        
    def MergeTracks(self, tracks):
    # Merge linear tracks into a single bidirectional track
    # Represent tracks as nodes in a directional graph
        
        box_half_width = 0.15
        box_height = 2*self.general['time_step']

        g = Graph( len( tracks) )
        matches = [[] for i in range( len(tracks) )]
        dist = [[] for i in range( len(tracks) )]

        # For each node, find prospective matches 
        for v, trackv in enumerate( tracks):

            # Find all possible matches
            for w, trackw in enumerate( tracks):

                # if tracks are close together
                if ( trackv.position[-1]-box_half_width < trackw.position[0] < trackv.position[-1]+box_half_width ) and ( trackv.time[-1] < trackw.time[0] < trackv.time[-1]+box_height ): 

                    # Add as a possible match
                    matches[v].append(w)

                    # find distance of match 
                    t1 = [ trackv.position[-1], trackv.time[-1]]
                    t2 = [ trackw.position[0], trackw.time[0]]
                    dist[v].append( math.sqrt( ((t1[0]-t2[0])**2)+((t1[1]-t2[1])**2) ) )

        # Find the best match
        for v, trackv in enumerate( tracks):

            if len( matches[v]) == 0:
                continue

            # Find match with lowest distance
            w = matches[v][dist[v].index( min( dist[v]) )]

            # Add edge between v and w
            g.addEdge(v,w)

        # Find connected components
        cc = g.connectedComponents()

        # Merge the tracks in time order
        tracks_merged = []
        for comp in cc:
            time = None 
            position = None 
            
            if len( comp) == 1:
                line_type = tracks[comp[0]].line_type
                direction = tracks[comp[0]].direction
            else:
                line_type = 'Curve'
                direction = 'Ambiguous'

            for v in comp:
                if time is None: 
                    time = tracks[v].time
                else:
                    time = np.concatenate( (time, tracks[v].time) )

                if position is None:
                    position = tracks[v].position
                else:
                    position = np.concatenate( (position, tracks[v].position) )

            tracks_merged += [Track(time, position, self.general['image'], self.poles, direction, line_type, time_step = self.general['time_step'], pos_step=self.tracks[0].pos_step)]
            
        return tracks_merged 

    def PlotTracks( self, tracks, poles=[], figName='tracks.pdf'):
        # Plot the given tracks in a figure

        # Number of tracks
        nt = len(tracks)
        # Number of poles
        np = len(poles)
        # Number of plots
        nn = nt+np

        # Colormap
        cm = plt.get_cmap('gist_rainbow')

        # Generate figure and axes and set colors
        fig = plt.figure( figsize=(12,8) )
        ax = fig.add_subplot(111)
        ax.set_prop_cycle(color=[cm( 1.*i/nn) for i in range(nt)])

        for idx,pole in enumerate(poles):
            ax.plot( pole.position, pole.time, linewidth=3, label = 'Pole {}'.format(1+idx))

        for idx,track in enumerate(tracks):
            ax.plot( track.position, track.time, linewidth=2, label = 'Track {}'.format(1+idx))
        plt.legend()

        # Set axes limits
        time_max = max( [max(trk.time) for trk in tracks] + [max(pol.time) for pol in poles] )
        time_min = min( [min(trk.time) for trk in tracks] + [min(pol.time) for pol in poles] )
        x_max = max( [max(trk.position) for trk in tracks] + [max(pol.position) for pol in poles] ) + 0.5
        x_min = min( [min(trk.position) for trk in tracks] + [min(pol.position) for pol in poles] ) - 0.5

        axes = plt.gca()
        axes.set_xlim([x_min, x_max])
        axes.set_ylim([time_min,time_max])
        fig.savefig( figName )

    def FindIntensityAlongSpindle(self, lrange=[0, 10]):

        if len( self.poles) != 2:
            return None

        # pdb.set_trace()
        dimT = np.shape( self.general['image'])[0]
        dimX = np.shape( self.general['image'])[1]

        # interpolation function for image
        try:
            f = interpolate.interp2d( self.tracks[0].pos_step*np.arange(0,dimX), self.tracks[0].time_step*np.arange(0,dimT), self.general['image'])
        except:
            pdb.set_trace()
            print('1')

        # Get times to find pole position
        tStart = max( self.poles[0].time[0], self.poles[1].time[0])
        tEnd = min( self.poles[0].time[-1], self.poles[1].time[-1])
        tVec = np.linspace( tStart, tEnd, math.ceil( (tEnd-tStart)/self.tracks[0].time_step) )

        # Get pole position
        pos0 = self.poles[0].ifunc( tVec)
        pos1 = self.poles[1].ifunc( tVec)

        # pdb.set_trace()
        # Trim to be within range
        pos0c = [i for i,j in zip(pos0,pos1) if np.abs(i-j) > lrange[0] and np.abs(i-j) < lrange[1]]
        pos1c = [j for i,j in zip(pos0,pos1) if np.abs(i-j) > lrange[0] and np.abs(i-j) < lrange[1]]
        tVecc = [k for i,j,k in zip(pos0,pos1,tVec) if np.abs(i-j) > lrange[0] and np.abs(i-j) < lrange[1]]
        if len(pos0c) == 0:
            return None

        # Find intensity between poles for each time value
        intense = np.zeros( (len(tVecc),100) )
        for i, tt in enumerate(tVecc):
            pVec =  np.linspace( pos0c[i], pos1c[i],100)
            ttVec = tt*np.ones((100,))
            intense[i,:] = f( pVec, ttVec)[0,:]
        return intense

    # # Trim tracks based on kmeans label 
    # def TrimTracksKmeansLabel(self, label=-1):

    def DisplayTracks(self, ax=None):

        if ax is None:
            fig, ax = plt.subplots(figsize=(6,6))

        # Display kymograph image
        ax.imshow( self.tracks[0].image)

        # Plot tracks
        for track in self.tracks:
            ax.plot( track.position/track.pos_step, track.time/track.time_step, color='red', linewidth=3)
        plt.show()


    def Print(self):
    # Print information about poles and tracks

        print(' ')
        print(' path: {}'.format(self.general['path_tiff'][0:-1]))
        print(' name: {}'.format(self.label ))
        print(' n_poles_exp: {}'.format(self.general['n_poles']))
        print(' n_poles_found: {}'.format(len(self.poles)))
        print(' n_tracks_exp: {}'.format(self.general['n_tracks']))
        print(' n_tracks_found: {}'.format( len(self.tracks)))
        print(' ')

        for feat in self.poles+self.tracks:
            feat.Print()


##########################################
if __name__ == "__main__":
    print('No default run method')


