#!/usr/bin/env python

import os, pdb
import math
import numpy as np
from .Track import *
import matplotlib.pyplot as plt

def ReadTxt( fname, verbose=0):
    # Read data from files and parse into general, poles and feature information

    if verbose:
        PrintFile( fname)

    # Initialize lists
    geninfo = []
    polesinfo = []
    featureinfo = []

    # Add General Information
    with open(fname) as fp:
        
        addLine = None 
        for cnt, line in enumerate(fp):
            
            if line.find( 'General Information') > -1:
                addLine = 'G' 
            if line.find( 'Poles Information') > -1:
                addLine = 'P' 
            if line.find( 'Feature Information') > -1:
                addLine = 'F' 

            # Add General Information
            if addLine == 'G':
                geninfo.append( line)

            # Add Poles Information
            elif addLine == 'P':
                polesinfo.append( line)

            # Add Feature Information
            elif addLine == 'F':
                featureinfo.append( line)
        
        # Parse information
        general = ParseGeneralInfo( fname, geninfo)
        poles = ParsePolesInfo( polesinfo, general)
        tracks = ParseTracksInfo( featureinfo, poles, general)
        if polesinfo == []:
            pdb.set_trace()
            hi = 1

        return general, poles, tracks

def ParseGeneralInfo( fname, geninfo):
    # Parse information about general information
    
    general = {
            'path_tiff' : [],
            'type' : [],
            'time_start' : [],
            'time_end' : [],
            'time_step' : [],
            'n_poles' : [],
            'n_tracks': [],
            'image': [],
            }


    for line in geninfo:

        # Tiff Path
        path_tiff = FindSingleSubstring( line, 'Tiff path : ') 
        if path_tiff is not None:
            general['path_tiff'] = path_tiff
        # Spindle Type
        typ = FindSingleSubstring( line, 'Spindle type : ') 
        if typ is not None:
            general['type'] = typ 
        # Time Start 
        time_start = FindNumbers( line, 'Start time (s) : ') 
        if time_start is not None:
            general['time_start'] = time_start 
        # Time End 
        time_end = FindNumbers( line, 'End time (s) : ') 
        if time_end is not None:
            general['time_end'] = time_end 
        # Time Step 
        time_step = FindNumbers( line, 'Time step (s) : ') 
        if time_step is not None:
            general['time_step'] = time_step[0]
        # Num Poles 
        npoles = FindNumbers( line, 'Num poles : ') 
        if npoles is not None:
            general['n_poles'] = int( npoles[0]) 
        # Num Tracks 
        ntracks = FindNumbers( line, 'Num tracks : ') 
        if ntracks is not None:
            general['n_tracks'] = int( ntracks[0])

    general['image'] = LoadTiff( fname[:-9]+'.tif')
    return general 

def ParsePolesInfo( polesinfo, general):
    # Parse information about poles 

    if not polesinfo or len(polesinfo) == 0:
        print('No poles information here')
        return

    # Determine number of poles and split information
    polelist = []
    idxPole = None
    nPoles = 0
    for line in polesinfo:

        # Look for the next pole
        if line.find( 'Pole number : {}'.format( nPoles+1)) > -1:
            nPoles += 1 
        if nPoles == 0:
            continue

        if nPoles != len(polelist):
            polelist += [[line]]
        else:
            polelist[ nPoles-1] += [line]

    # print('Found {} poles'.format( nPoles) )

    # for each split pole, get useful information and initialize a Pole object
    poles = []
    for pole in polelist:

        for line in pole:

            # Time pixels
            if FindNumbers( line, 'Time pixel : ') is not None:
                time = FindNumbers( line, 'Time pixel : ')
                time = [x * general['time_step'] for x in time]

            # # Times 
            # if FindNumbers( line, 'Time (s) : ') is not None:
                # time = FindNumbers( line, 'Time (s) : ') 
            
            # Position 
            if FindNumbers( line, 'Position (um) : ') is not None:
                position = FindNumbers( line, 'Position (um) : ') 

            # Intensity 
            if FindNumbers( line, 'Intensity : ') is not None:
                intensity = FindNumbers( line, 'Intensity : ') 
            
        poles += [Pole( time, position, general['image'], time_step=general['time_step']) ]

    return poles


def ParseTracksInfo( featureinfo, poles, general):
    # Parse information about tracks 

    if not featureinfo or len(featureinfo) == 0:
        print('No tracks information here')
        return

    # Determine number of tracks and split information
    tracklist = []
    idxTrack = None
    nTracks = 0
    for line in featureinfo:

        # Look for the next track 
        if line.find( 'Feature number : {}'.format( nTracks+1)) > -1:
            nTracks += 1 
        if nTracks == 0:
            continue

        if nTracks != len(tracklist):
            tracklist += [[line]]
        else:
            tracklist[ nTracks-1] += [line]

    # print('Found {} tracks'.format( nTracks) )

    # for each split track, get useful information and initialize a Track object
    tracks = []
    for trck in tracklist:

        for line in trck:

            # Time pixels
            if FindNumbers( line, 'Time pixel : ') is not None:
                time = FindNumbers( line, 'Time pixel : ')
                timePix = time
                time = [x * general['time_step'] for x in time]

            # # Time 
            # if FindNumbers( line, 'Time (s) : ') is not None:
                # time = FindNumbers( line, 'Time (s) : ') 
            
            # Position 
            if FindNumbers( line, 'Position pixel : ') is not None:
                positionPix = FindNumbers( line, 'Position pixel : ') 

            # Position 
            if FindNumbers( line, 'Position (um) : ') is not None:
                position = FindNumbers( line, 'Position (um) : ') 

            # Intensity 
            if FindNumbers( line, 'Intensity : ') is not None:
                intensity = FindNumbers( line, 'Intensity : ') 
            
            # Direction
            if FindSingleSubstring( line, 'Feature direction : ') is not None:
                direction = FindSingleSubstring( line, 'Feature direction : ') 
                direction = direction[0:-1]

            # Line type 
            if FindSingleSubstring( line, 'Feature type : ') is not None:
                line_type = FindSingleSubstring( line, 'Feature type : ') 
                line_type = line_type[0:-1]

        tracks += [Track( time, position, general['image'], poles, 'Ambiguous', line_type, time_step=general['time_step'], pos_step=0.1067) ]
    return tracks

def LoadTiff( fname):
    # load tiff file
    arr = plt.imread( fname)
    if len( arr.shape) == 3:
        arr = np.mean(arr,axis=2)
    return arr

def FindSingleSubstring(strSearch, strLabel):
    # Find a single substring that contains strLabel. We delete the strLabel. 
   
    if strSearch.find( strLabel) > -1:
        strMatch = strSearch.replace( strLabel, '')
        return strMatch
    return None

def FindNumbers(strSearch, strLabel):
    # Find numbers from a string that starts with strLabel
   
    if strSearch.find( strLabel) > -1:
        strMatch = strSearch.replace( strLabel, '')
        strList = strMatch.split(',')
        nums = [float(i) for i in strList]
        return nums 
    return None

def PrintFile(fname):
    # Print all the information from a file to screen 
    with open( fname) as f:
        print( f.read() )

##########################################
if __name__ == "__main__":
    print("no default implementation")
    
