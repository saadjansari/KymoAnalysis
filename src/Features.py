#!/usr/bin/env python

import os, pdb
import numpy as np
from scipy import interpolate, signal
from .node_graph import Graph
import matplotlib.pyplot as plt
import math
import uuid

# Superclass for Poles and tracks that stores positional and intensity information
class Feature:
    def __init__(self, time, position, intensity, strain='xxx', time_step=1):
        self.time = np.array( time )
        self.position = np.array( position )
        self.intensity = np.array( intensity )
        self.id = uuid.uuid1()
        self.strain = strain
        self.time_step = time_step
        self.pixel_time = self.time / self.time_step

        # Resample data
        self.ResampleData()

    def ResampleData( self, sample_factor=3):
        # resample data based on time pixels

        # Define an interpolation function for positions
        ifunc_pos = interpolate.interp1d( self.time, self.position, kind='linear')

        # Define a grid of resampled time points        
        self.time = np.linspace( self.time[0], self.time[-1], max([ 2, sample_factor*(self.time[-1]-self.time[0])]) )
        if len(self.time) == 1:
            pdb.set_trace()
            print('oops')
        self.position = ifunc_pos( self.time) 

# Class for a Pole
class Pole(Feature):
    def __init__(self, time, position, intensity=[], time_step=1, strain='xxx'):
        Feature.__init__(self, time, position, intensity, strain=strain, time_step=time_step)

        # Define an interpolation/extrapolation function
        self.ifunc = interpolate.interp1d(self.time, self.position, kind='linear', fill_value='extrapolate')

    def Print(self):
        print('Pole :')
        print(' ID : {}'.format(self.id))
        print(' Time : {}'.format( self.time))
        print(' Position : {}'.format( self.position))
        print(' Intensity : {}'.format( self.intensity))
        print('--------------------------------- ')

# Class for a Track: additionally stores associated poles and track direction
class Track(Feature):
    def __init__(self, time, position, intensity, poles, direction, line_type, time_step=1, strain='xxx'):
        Feature.__init__(self, time, position, intensity, time_step=time_step, strain=strain)

        if time_step == 1:
            pdb.set_trace()
            print('woah')

        self.poles = poles
        self.direction = direction
        self.line_type = line_type
        self.polePosition = []
        self.data = {
                'pos_pole' : np.zeros( (2, np.size(self.position) ) ),
                'pos_track_rel' : [],
                'velocity' : {  'P' : [], 'AP' : [],'I': []},
                'runlength' : { 'P' : [], 'AP' : [],'I': []},
                'lifetime' : {  'P' : [], 'AP' : [],'I': []},
                'lifetime_total' : [],
                'velocity_mean' : [],
                'switch_count' : [],
                'switch_total' : [],
                }
        self.bad = 0

        # Order poles with 1st pole being main pole(closest at start)
        # self.OrderPoles()

        # Calcualte spindle length
        self.CalcSpindleLength()

        if self.line_type == 'Curve' and self.direction != 'Ambiguous':
            self.direction = 'Ambiguous'
            # pdb.set_trace()
            # print('1')

    def Analyze(self, ipole=0):
    # Run useful analysis methods

        self.CalcPositionTrackRelativeToPole()

        # Split the track and save analysis 
        tracks_mini, switches = self.SplitTrack( ipole=ipole)

        for track in tracks_mini:

            if track.direction == 'Poleward':
                label = 'P'
            elif track.direction == 'Antipoleward':
                label = 'AP'
            elif track.direction == 'Inactive':
                label = 'I'
            else:
                pdb.set_trace()
                raise ValueError('line direction is neither poleward nor antipoleward nor inactive')

            # Calculate and append data of the mini track
            # Velocity
            self.data['velocity'][label] += [track.CalcVelocityLinear(ipole=ipole)]
            # Run length
            self.data['runlength'][label] += [track.CalcRunLength(ipole=ipole)]
            # Lifetime
            self.data['lifetime'][label] += [track.CalcLifetime()]

        # Combine data from the split tracks
        self.data['lifetime_total'] = self.CalcLifetime()
        self.data['velocity_mean'] = self.CalcVelocityMean(ipole=ipole)
        self.data['switch_count'] = switches
        # pdb.set_trace()
        # print('woah')

    def OrderPoles(self):
        # Order the poles with the first one being the closest one to the start of the track

        if len(self.poles) != 2:
            return

        pos = self.CalcPositionTrackRelativeToPole()
        if np.absolute( pos[1,0] ) < np.absolute( pos[0,0]):
            self.poles = [self.poles[1], self.poles[0]]
        
    def CalcSpindleLength(self):
        # Calculate the spindle length

        if len(self.poles) != 2:
            return

        # Find the distance between the poles for the extent of this track
        self.spindleLength = np.absolute( self.poles[0].ifunc( self.time) - self.poles[1].ifunc( self.time) )

    def CheckViability(self):
        # Check if the track's time points are always increasing

        self.bad = 0
        # Check track time is always increasing
        if np.any( np.diff( self.time) <= 0 ):
            self.bad = 1 

        return self.bad

    def CheckLinearLifetime( self, min_lt = 0.5):
        # Check lifetime is above a min threshold

        self.bad = 0
        if self.line_type == 'Line' and self.CalcLifetime() < min_lt:
            self.bad = 1

        return self.bad


    def CalcPositionPoleCurrent(self):
    # Get pole position at the current time (i.e at the times of the track) by using the interpolation/extrapolation function of the pole

        for idx, pole in enumerate( self.poles) :
            pos_pole = np.array(pole.ifunc( self.time) )
            self.data['pos_pole'][idx,:] = pos_pole
        return self.data['pos_pole'] 

    def CalcPositionTrackRelativeToPole(self):
    # Calculate track position relative to the pole 

        pos_pole = self.CalcPositionPoleCurrent()

        # If bipolar spindle
        if len( self.poles) == 2:
            pos_track_rel = np.zeros( np.shape(pos_pole))
            for idx,ele in enumerate( pos_pole):
                pos_track_rel[idx,:] = np.array( self.position - ele)

        # If monopolar spindle
        else:
            pos_track_rel = np.array( self.position - pos_pole )
            # pos_track_rel = pos_track_rel[0,:]

        self.data['pos_track_rel'] = pos_track_rel
        return pos_track_rel

    def CalcVelocityLinear(self, ipole=0):
    # Calculate the velocity of this linear track
        
        if self.direction == 'Ambiguous':
            raise Exception('Track.CalcVelocityLinear() is only defined for tracks with a single direction')

        # Calc relative positions if not done already
        if len( self.data['pos_track_rel']) == 0 or not self.data['pos_track_rel'].any():
            pos_track_rel = self.CalcPositionTrackRelativeToPole()
        else:
            pos_track_rel = self.data['pos_track_rel']

        # Check
        if len(self.time) <= 1:
            pdb.set_trace()
            print('oops')

        # Find Velocity
        vel = np.average( np.absolute( np.divide( np.diff( pos_track_rel[ipole,:]) , np.diff( self.time) ) ), weights = np.diff(self.time) )

        # Check
        if np.size( vel) > 1:
            pdb.set_trace()
            print('1')

        return vel
    
    def CalcRunLength(self, ipole=0):
    # Calculate the run length of this track
        
        # Calc relative positions if not done already
        if len( self.data['pos_track_rel']) == 0 or not self.data['pos_track_rel'].any():
            pos_track_rel = self.CalcPositionTrackRelativeToPole()
        else:
            pos_track_rel = self.data['pos_track_rel']

        # Find Run length 
        run_length = np.absolute( pos_track_rel[ipole,-1] - pos_track_rel[ipole,0] )
            
        self.data['run_length'] = run_length

        # Check
        if np.size( run_length) > 1:
            pdb.set_trace()
            print('1')

        return run_length

    def CalcLifetime(self):
    # Calculate the lifetime of this track

        lifetime = self.time[-1] - self.time[0]
        self.data['lifetime'] = lifetime
        return lifetime

    def CalcVelocityMean(self,ipole=0):
    # Calculate the mean velocity of this track

        if self.line_type == 'Curve' and not self.data['velocity']:
            # Split the track 
            tracks_mini = self.SplitTrack()
            for track in tracks_mini:
                vv = track.CalcVelocityLinear()
                if track.direc == 'Poleward':
                    self.data['velocity']['P'] += [ vv[ipole]]
                elif self.direc == 'Antipoleward':
                    self.data['velocity']['AP'] += [ vv[ipole]]
        
        vel_mu = np.mean( np.concatenate( (self.data['velocity']['P'], self.data['velocity']['AP']) ) )
        return vel_mu

    def CalcSwitchingCount(self):
    # Calculate the mean velocity of this track

        if self.line_type == 'Curve' and not self.data['velocity']:
            # Split the track 
            tracks_mini = self.SplitTrack()
            for track in tracks_mini:
                if track.direc == 'Poleward':
                    self.data['velocity']['P'] += [track.CalcVelocityLinear()]
                elif self.direc == 'Antipoleward':
                    self.data['velocity']['AP'] += [track.CalcVelocityLinear()]
        
        vel_mu = np.mean( np.concat( self.data['velocity']['P'], self.data['velocity']['AP']) )
        return vel_mu

    def CalcIntensityMean(self):
    # Calculate the mean intensity of this track
        self.data['intensity_mean'] = np.mean( self.intensity) 
        
    def SplitTrack(self, ipole=0):
    # Spit curved track into multiple mini unidirectional tracks

        switches = {
                'P' : { 'P' : 0, 'AP': 0, 'I' : 0,},
                'AP' : { 'P' : 0, 'AP': 0, 'I' : 0,},
                'I' : { 'P' : 0, 'AP': 0, 'I' : 0,},
                }

        if self.direction != 'Ambiguous':
            return [self], switches

        if self.line_type == 'Line' and self.direction == 'Ambiguous':
            position = np.absolute( self.CalcPositionTrackRelativeToPole() )
            vel = np.mean( np.divide( np.diff( position) , np.diff(self.time) ) )
            if abs( vel) < 0.005:
                self.direction = 'Inactive'
            elif vel > 0:
                self.direction = 'Antipoleward'
            elif vel  < 0:
                self.direction = 'Poleward'
            return [self], switches

        # Find track position relative to the pole
        if len( self.data['pos_track_rel']) == 0 or not self.data['pos_track_rel'].any():
            position = self.CalcPositionTrackRelativeToPole()
        else:
            position = self.data['pos_track_rel'][ipole,:]

        position = np.absolute( position)
        states = []

        # Smoothing window:
        # Use a time-sliding window to find the average velocity, and use that to figure out state
        def FindStates_RollingWindow( positions, times, t_window, v_cutoff=1):

            dt = np.mean( np.diff( times) ) 
            n_hwindow = int( np.ceil( t_window / (2*dt)) )
            states = []

            for i, t in enumerate( times): 
                i_min = max( [ 0, i-n_hwindow])
                i_max = min( [ len(times), i+n_hwindow])

                vel = np.mean( np.divide( np.diff( positions[i_min:i_max] ) , np.diff( times[i_min:i_max] ) ) )
                # pdb.set_trace()

                # Assign labels based on value of vel 
                if abs( vel) < v_cutoff:
                    states += ['I']
                elif vel > 0:
                    states += ['AP']
                elif vel  < 0:
                    states += ['P']

            return states
        
        states = FindStates_RollingWindow(position,self.time,5,v_cutoff=0.005)

        # Remove singly occuring states
        for cnt, st in enumerate(states):

            if cnt > 1 and cnt < len(states)-1:
                if st != states[cnt-1] and st != states[cnt+1]:
                    states[cnt] = states[cnt-1]

            # set first state to second state
            if cnt == 0:
                states[cnt] = states[cnt+1]

            # set last state to second last state
            if cnt == len(states)-1:
                states[cnt] = states[cnt-1]

        # Count switches and get track indices
        p_state = 'XXX' 
        track = { 'pos': [], 'time': [], 'dir':[] }
        idx = [0 , 0]
        for cnt, st in enumerate(states):
            
            if cnt == 0:
                p_state = st
                idx[0] = 0
                continue

            if st == p_state:
                idx[1] += 1 
            
            if st != p_state:

                # store old stuff
                pos = self.position[ idx[0]: idx[1]+2]
                # pos.tolist()
                time = self.time[ idx[0]: idx[1]+2]
                # time.tolist()
                track['pos'] += [pos]
                track['time'] += [time]
                track['dir'] += [p_state]
                p_state = st

                # begin new
                idx[0] = cnt
                idx[1] = cnt

            # Store the last info
            if cnt == len(states)-1:

                pos = self.position[ idx[0]: idx[1]+1]
                # pos.tolist()
                time = self.time[ idx[0]: idx[1]+1]
                # time.tolist()
                track['pos'] += [pos]
                track['time'] += [time]
                track['dir'] += [p_state]

        # record switches
        for cnt, dd in enumerate( track['dir']):
            if cnt == 0:
                continue
            switches[ track['dir'][cnt-1]][track['dir'][cnt]] += 1
            
        # Create track objects from the information
        mini_tracks = []
        for time, pos, direc in zip( track['time'], track['pos'], track['dir']):
            if direc is 'P':
                direction = 'Poleward'
            elif direc is 'AP':
                direction = 'Antipoleward'
            elif direc is 'I':
                direction = 'Inactive'
            pos = pos.tolist()
            time = time.tolist()
            if len(pos) == 1:
                pdb.set_trace()
                print('oops')

            mini_tracks += [Track( time, pos, self.intensity, self.poles, direction, 'Line', time_step=self.time_step, strain=self.strain)]

        # if self.strain == 'B PA-GFP' and mini_tracks[0].direction == 'Inactive':
            # pdb.set_trace()
            # print('1')
        
        for t in mini_tracks:
            if len( t.position) < 2:
                pdb.set_trace()
                print('oops')

        return mini_tracks, switches

    def PlotCurveWithStates(self, figname='curved_track.pdf'):
        # Plot a curve with states( inactive, poleward and antipoleward) in different colors.

        cols = {
                'Inactive'     : 'blue',
                'Poleward'     : 'red',
                'Antipoleward' : 'green',
                }

        minis, switches = self.SplitTrack()

        # Generate figure and axes and set colors
        fig = plt.figure( figsize=(6,4) )
        ax = fig.add_subplot(111)

        pos_pole = self.CalcPositionPoleCurrent()
        for idx,pole in enumerate(self.poles):
            ax.plot( pos_pole[idx,:], self.time, linewidth=3 )

        for trk in minis:
            ax.plot( trk.position, trk.time, linewidth=2, color=cols[trk.direction] )

        plt.text(1,1, 'Poleward', color='red', transform=ax.transAxes, ha='right', va='top')
        plt.text(1,0.95, 'AntiPoleward', color='green', transform=ax.transAxes, ha='right', va='top')
        plt.text(1,0.9, 'Inactive', color='blue', transform=ax.transAxes, ha='right', va='top')
        plt.text(1,0.85, 'MainPole', color='skyblue', transform=ax.transAxes, ha='right', va='top')
        plt.text(1,0.8, 'SecondaryPole', color='orange', transform=ax.transAxes, ha='right', va='top')
        # Set axes limits
        axes = plt.gca()
        x_min = min([ min(self.position), min([min( pol.position) for pol in self.poles]) ]) -0.5
        x_max = max([ max( self.position), max([max( pol.position) for pol in self.poles]) ]) +0.5
        axes.set_xlim([ x_min, x_max])
        axes.set_xlabel('Position')
        axes.set_ylim([ min(self.time)-5,max(self.time)+5])
        axes.set_ylabel('Time')

        fig.savefig( figname)

    def Trim(self,lrange=[0,100]):
        # Trim the track to be inside the range specified
        
        if len( self.poles) == 1:
            return self

        # Get indices of times when spindle length is between the given range values
        lens = self.spindleLength
        idx = np.argwhere( (lens > lrange[0]) & (lens < lrange[1]) ).T[0].tolist()
        if len(idx) == 0:
            return None

        idx = range( idx[0], idx[-1]+1) 

        # Create the new trimmed track
        tracknew = Track( self.time[idx], self.position[idx], self.intensity, self.poles, self.direction, self.line_type, time_step=self.time_step, strain=self.strain)

        return tracknew


    def Print(self):
        print('Feature :')
        print(' ID : {}'.format(self.id))
        print(' Direction : {}'.format( self.direction))
        print(' Line type : {}'.format( self.line_type))
        print(' Time : {}'.format( self.time))
        print(' Position : {}'.format( self.position))
        print(' Intensity : {}'.format( self.intensity))
        print('--------------------------------- ')


