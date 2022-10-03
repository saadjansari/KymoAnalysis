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
    def __init__(self, time, position, image, time_step=1):
        self.time = np.array( time )
        self.position = np.array( position )
        self.id = uuid.uuid1()
        self.time_step = time_step
        self.pixel_time = self.time / self.time_step
        self.image = image

        # Resample data
        self.ResampleData()

    def ResampleData( self, sample_factor=3):
        # resample data based on time pixels

        # Define an interpolation function for positions
        ifunc_pos = interpolate.interp1d( self.time, self.position, kind='linear')

        # Define a grid of resampled time points        
        self.time = np.linspace( self.time[0], self.time[-1], int( np.floor(max([ 2, sample_factor*(self.time[-1]-self.time[0])])) ))
        if len(self.time) == 1:
            pdb.set_trace()
            print('oops')
        self.position = ifunc_pos( self.time) 

# Class for a Pole
class Pole(Feature):
    def __init__(self, time, position, image=[], time_step=1):
        Feature.__init__(self, time, position, image, time_step=time_step)

        # Define an interpolation/extrapolation function
        # self.ifunc = interpolate.interp1d(self.time, self.position, kind='linear', fill_value='extrapolate')
        self.ifunc = interpolate.interp1d(self.time, self.position, kind='linear', fill_value=(self.position[0], self.position[-1]), bounds_error=False)

    def Print(self):
        print('Pole :')
        print(' ID : {}'.format(self.id))
        print(' Time : {}'.format( self.time))
        print(' Position : {}'.format( self.position))
        print('--------------------------------- ')

    def TrimBasedOnTime(self, time_keep):
        # Trim the pole to be inside the time range specified
        
        if np.all(time_keep == -1):
            return np.nan

        # Check if track exists between those times
        start_before = (self.time[0] < time_keep[0])
        start_after = (self.time[0] > time_keep[1])
        end_before = (self.time[-1] < time_keep[0])
        end_after = (self.time[-1] > time_keep[1])
        if start_before and end_before:
            return np.nan
        elif start_after and end_after:
            return np.nan

        # Get indices of times 
        idx = np.argwhere( (self.time > time_keep[0]) & (self.time < time_keep[1]) ).T[0].tolist()
        if len(idx) < 3:
            return None 
        idx = range( idx[0], idx[-1]+1) 

        # Create the new trimmed pole 
        polenew = Pole( self.time[idx], self.position[idx], self.image, time_step=self.time_step)
        # print(time_keep)
        # print(polenew.time[0])
        # print(polenew.time[-1])
        return polenew

# Class for a Track: additionally stores associated poles and track direction
class Track(Feature):
    def __init__(self, time, position, image, poles, direction, line_type, time_step=1, pos_step=1, kymo_file=None):
        Feature.__init__(self, time, position, image, time_step=time_step)

        self.poles = poles
        self.direction = direction
        self.line_type = line_type
        self.pos_step = pos_step
        self.kymo_file = kymo_file

    def CalcPositionPoleCurrent(self):
    # Get pole position at the current time (i.e at the times of the track) by using the interpolation/extrapolation function of the pole
        pos = np.zeros( (len(self.poles), np.size(self.position) ) )
        for idx, pole in enumerate( self.poles) :
            pos[idx,:] = np.array( pole.ifunc( self.time) )
        return pos 

    def CalcPositionRelative(self):
    # Calculate track position relative to the pole 

        pole = self.CalcPositionPoleCurrent()
        pos = np.zeros( np.shape(pole) )
        for idx, ele in enumerate( pole):
            pos[idx,:] = np.abs( np.array( self.position - ele) )
        return pos

    def CalcVelocity(self):
    # Calculate the velocity of this linear track
        
        pos = self.CalcPositionRelative()
        # Find Velocity
        vel = np.zeros( (len(self.poles)) )
        for idx in range( len(self.poles)):
            vel[idx] = np.average( np.absolute( np.divide( np.diff( pos[idx,:]) , np.diff( self.time) ) ), weights = np.diff(self.time) )

        return vel

    def CalcSpindleLength(self):
        # Calculate the spindle length

        if len(self.poles) != 2: 
            return
        # Find the distance between the poles for the extent of this track
        leng = np.absolute( self.poles[0].ifunc( self.time) - self.poles[1].ifunc( self.time) )
        return leng

    def CalcIntensity( self):
        # Interpolate to find the mean intensity of the track

        dimT = np.shape( self.image)[0]
        dimX = np.shape( self.image)[1]
        f = interpolate.interp2d( self.pos_step*np.arange(0,dimX), self.time_step*np.arange(0,dimT), self.image)
        intense = f(self.position, self.time)
        return np.mean(intense)

    def CheckViability(self):
        # Check track time is always increasing
        if np.any( np.diff( self.time) <= 0 ):
            return 0 
        return 1 

    def OrderPoles(self):
        # Order the poles with the first one being the closest one to the start of the track
        if len(self.poles) != 2:
            return

        pos = self.CalcPositionRelative()
        if np.absolute( pos[1,0] ) < np.absolute( pos[0,0]):
            self.poles = [self.poles[1], self.poles[0]]

    def Trim(self, lrange):
        # Trim the track to be inside the range specified
        
        if len( self.poles) == 1:
            return self
        if lrange is None:
            return self

        # Get indices of times when spindle length is between the given range values
        lens = self.CalcSpindleLength()
        idx = np.argwhere( (lens > lrange[0]) & (lens < lrange[1]) ).T[0].tolist()
        if len(idx) < 3:
            return None 
        idx = range( idx[0], idx[-1]+1) 

        # Create the new trimmed track
        tracknew = Track( self.time[idx], self.position[idx], self.image, self.poles, self.direction, self.line_type, time_step=self.time_step, pos_step=self.pos_step)
        return tracknew

    def TrimBasedOnTime(self, time_keep):
        # Trim the track to be inside the time range specified
        
        if np.all(time_keep == -1):
            return np.nan

        # Check if track exists between those times
        start_before = (self.time[0] < time_keep[0])
        start_after = (self.time[0] > time_keep[1])
        end_before = (self.time[-1] < time_keep[0])
        end_after = (self.time[-1] > time_keep[1])
        if start_before and end_before:
            return np.nan
        elif start_after and end_after:
            return np.nan

        # Get indices of times 
        idx = np.argwhere( (self.time > time_keep[0]) & (self.time < time_keep[1]) ).T[0].tolist()
        if len(idx) < 3:
            return None 
        idx = range( idx[0], idx[-1]+1) 

        # Create the new trimmed track
        tracknew = Track( self.time[idx], self.position[idx], self.image, self.poles, self.direction, self.line_type, time_step=self.time_step, pos_step=self.pos_step)
        if tracknew is None:
            pdb.set_trace()
            print('b')
        return tracknew

    def SplitTrack(self, ipole=0, cutoff=0.003):
    # Spit curved track into multiple mini unidirectional segments 
    # cutoff : units micron/sec

        switches = {
                'P' : { 'P' : 0, 'AP': 0, 'I' : 0,},
                'AP' : { 'P' : 0, 'AP': 0, 'I' : 0,},
                'I' : { 'P' : 0, 'AP': 0, 'I' : 0,},
                }

        # If linear directional track, cant split, so exit
        if self.direction != 'Ambiguous':
            return [self], switches

        # If linear ambiguous track, figure out direction, then exit 
        if self.line_type == 'Line' and self.direction == 'Ambiguous':
            if len(self.CalcPositionRelative()) == 0:
                pdb.set_trace()
                print('a')
            position = np.absolute( self.CalcPositionRelative()[ipole,:] )
            vel = np.mean( np.divide( np.diff( position) , np.diff(self.time) ) )
            if abs( vel) < cutoff:
                self.direction = 'Inactive'
            elif vel > 0:
                self.direction = 'Antipoleward'
            elif vel  < 0:
                self.direction = 'Poleward'
            return [self], switches

        # Get track position relative to the pole
        position = np.absolute( self.CalcPositionRelative()[ipole,:] )

        # Use a rolling window to find velocities
        vel = FindGradientRollingWindow( position, self.time, window=16)

        # Assign states based on value of velocity at each timestep
        states = []
        for v in vel:
            if abs( v) < cutoff:
                states += ['I']
            elif v > 0:
                states += ['AP']
            elif v  < 0:
                states += ['P']
        # set first state to second state. last state to second last state
        states[0] = states[1]
        states[-1] = states[-2]
        # Remove singly occuring states
        for i, state in enumerate(states):
            if i>0 and i< len(states)-1:
                if state != states[i-1] and state != states[i+1]:
                    states[i] = states[i-1]

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
        segments = []
        for time, pos, direc in zip( track['time'], track['pos'], track['dir']):
            if direc is 'P':
                direction = 'Poleward'
            elif direc is 'AP':
                direction = 'Antipoleward'
            elif direc is 'I':
                direction = 'Inactive'
            pos = pos.tolist()
            time = time.tolist()
            segments += [Track( time, pos, self.image, self.poles, direction, 'Line', time_step=self.time_step, pos_step=self.pos_step, kymo_file=self.kymo_file)]
        
        return segments, switches

    def DisplayTrack(self, ax=None):

        if ax is None:
            fig, ax = plt.subplots(figsize=(6,6))

        # Display kymograph image
        ax.imshow( self.image)
        ax.plot( self.position/self.pos_step, self.time/self.time_step, color='red')

    def Print(self):
        print('Feature :')
        print(' ID : {}'.format(self.id))
        print(' Direction : {}'.format( self.direction))
        print(' Line type : {}'.format( self.line_type))
        print(' Time : {}'.format( self.time))
        print(' Position : {}'.format( self.position))
        print('--------------------------------- ')

def CountSwitches( states, switches):
    # Given a list of 

    dt = np.mean( np.diff( t) ) 
    nHalfWindow = int( np.ceil( t_window / (2*dt)) )

    for i  in range(len(t)):
        # get upper lower indices of window
        i_lb = max( [ 0, i-nHalfWindow])
        i_ub = min( [ len(t), i+nHalfWindow])

        # Find gradient
        diff = lambda xx : np.diff( xx[i_lb:i_ub])
        grad = np.mean( np.divide( diff(x), diff(t) ) )
    return grad 

def FindGradientRollingWindow( x, t, window=6):

    dt = np.mean( np.diff( t) ) 
    nHalfWindow = int( np.ceil( window / (2*dt)) )
    grads = []
    for i  in range(len(t)):
        # get upper lower indices of window
        i_lb = max( [ 0, i-nHalfWindow])
        i_ub = min( [ len(t), i+nHalfWindow])

        # Find gradient
        diff = lambda xx : np.diff( xx[i_lb:i_ub])
        grads += [np.mean( np.divide( diff(x), diff(t) ) )]
    return grads

if __name__ == "__main__":
    print('Not implemented')
