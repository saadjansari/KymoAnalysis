#!/usr/bin/env python

import os, pdb
import yaml
import glob

'''
Name: Load.py
Description: loads and splits the tracks saved by the trackBuilder (kyman.mlapp) into general, poles and feature sections to be parsed by Kymograph.py 
'''


# Class to load data from files 
class Load:
    def __init__(self, verbose=0):

        file_name = 'track_files.yaml'

        with open(file_name) as infile:
            self.data = yaml.load(infile)

        self.verbose = verbose 

        self.GetFilenames()
        self.ReadFromFiles()

    def GetFilenames(self):
        # Expand filenames in the case of special characters

        for strain, dat in self.data['strain'].items():
            for idx,fpath in enumerate( dat['path']):
                files = []
                for fname in dat['files'][idx]:
                    
                    temp = glob.glob( os.path.join(fpath,fname) )                      
                    for fil in temp:
                        head_tail = os.path.split(fil)
                        files += [ head_tail[1] ]
                self.data['strain'][strain]['files'][idx] = files
   
    def ReadFromFiles(self):
        # Read information from all files given yaml data 
        
        for strain, dat in self.data['strain'].items():
            for idx,fpath in enumerate( dat['path']):

                self.data['strain'][strain]['geninfo'] = []
                self.data['strain'][strain]['polesinfo'] = []
                self.data['strain'][strain]['featureinfo'] = []

                for fname in dat['files'][idx]:

                    gen, poles, feats = self.ReadFromFile( fpath, fname)
                    self.data['strain'][strain]['geninfo'] += [gen]
                    self.data['strain'][strain]['polesinfo'] += [poles]
                    self.data['strain'][strain]['featureinfo'] += [feats]


    def ReadFromFile(self, fpath, fname):
        # Read data from files and parse into general, poles and feature information

        # Initialize lists
        geninfo = []
        polesinfo = []
        featureinfo = []

        if self.verbose:
            self.PrintFile( fname)

        # Add General Information
        with open(fpath + fname) as fp:
            
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

        return geninfo, polesinfo, featureinfo

    def PrintFile(self, fname):
        # Print all the information from a file to screen 

        fp = open( self.fpath + fname)
        fc = fp.read()
        print(fc)
        fp.close()

##########################################
if __name__ == "__main__":
    
    x = Load(verbose=1)

