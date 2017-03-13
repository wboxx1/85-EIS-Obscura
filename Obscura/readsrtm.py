# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 09:43:08 2017

@author: Boxx
"""

import numpy as np
import pandas as pd
import readhgt
import collections
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#@profile
def run(coord, res_km, init_elev_km, filenamekml, filenamecsv):
    '''    
    function [azimuths,elevations] = readsrtm(coord,res_km,init_elev_km,filenamekml,filenamecsv)
    %READSRTM       Create 360 degree terrain obscura using srtm data
    %   [AZIMUTHS,ELEVATIONS] = READSRTM(coord,res_km,init_elev_km,fnkml,fncsv)  
    %   This function is used to create a terrain obscura using a NASA srtm 
    %   data. The srtm data will be downloaded and requires internet access.  
    %   The proper srtm file will be located using the center coordinate in 
    %   coord. The function will create a kml file of lat lon points at every 
    %   line given in coord that maximize obscura elevation.  It also creates a
    %   csv file of all the obscura elevations for 360 degrees of azimuth.  
    %   coord is the array of coordinates created by using azimuthLines.m.
    %   res_km is the resolution or distance between points in the line and
    %   must be in kilometers.  init_elev_km is the elevation of the source
    %   point and must be in kilometers.  fnkml is the filename for the kml
    %   file.  fncsv is the filename for the csv file.  This function requires
    %   readhgt to run
    %
    %   Example:  Generate terrain obscura for coordinates in previously
    %               generated vector coord.  Use a resolution of 10 meters and
    %               a source point elevation of 10 meters.  Use kml filename as
    %               testkml and csv filename of testcsv.
    %
    %       [azimuths,elevations] = ...
    %       readsrtm(coord,0.01,0.01,'testkml','testcsv')
    %
    %   Dependency - This function requires that the readhgt program bundle be
    %   in the same directory.
    '''
    
    #read data block of srtm data using mu coordinate
    X = readhgt.run(lat=coord[0,0,0],lon=coord[0,1,0],options=['srtm1']);
    
    #X = pd.read_pickle(path='.\X.pkl') #used for debugging
        
    feet = 1
    meters = 1
    
    lats = X['lat'][0]
    lats = lats[::-1]
    lons = X['lon'][0]
    els = X['z'][0]
    
    [h,w,d] = coord.shape

    mylist = coord.values
               
    elevations = np.zeros([h,d])
    ea = np.zeros(h)
    elAngles = np.zeros(d)
    res_m = res_km*1000.0*meters    #resolution - set in azimuthLines.m
    indexes = np.zeros(d)
    '''
    ##############################
    X,Y = np.meshgrid(lons,lats)
    Z = els

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X,Y,Z)
    plt.show()
    ##############################
    '''
    #mueHeightFT = 30.5*feet;
    init_elev_m = init_elev_km*1000.0*meters  #mueHeightFT * 0.3048*meters;
    
    for i in range(d):
        j = 0;
        lat = mylist[j,0,i]
        lon = mylist[j,1,i]
        
        a = np.abs(lats - lat)
        b = np.abs(lons - lon)
        lati = np.argmin(a)
        loni = np.argmin(b)
        refel = float(els[lati,loni]) + init_elev_m;
        if not i%10:
            print str(i) + ' of 360'
        
        for j in range(1,h):          # SHOULD NOT BE HARDCODED TO 2500
            if j == 510 and i == 64:
                tpry = 0

            lat = mylist[j,0,i]
            lon = mylist[j,1,i]
            
            a = np.abs(lats - lat)
            b = np.abs(lons - lon)
            lati = np.argmin(a)
            loni = np.argmin(b)
            
            elevations[j,i] = els[lati,loni]
            tp = els[lati,loni]

            temp2 = res_m*(j)

            if temp2 == 0:
                temp2 = 2**-1024;
            
            ea[j] = np.arctan((tp-refel)/temp2)

        
        m,indexes[i] = np.max(ea),np.argmax(ea)
        elAngles[i] = m*180/np.pi
        
    maxcoords = np.zeros([d,2])
        
    for l in np.arange(d):
        maxcoords[l,0] = mylist[int(indexes[l]),0,l]
        maxcoords[l,1] = mylist[int(indexes[l]),1,l]
   
    filename = filenamekml + '.kml'

    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    abs_file_path = os.path.join(script_dir, filename)
    path_only = os.path.dirname(abs_file_path)# abspath(abs_file_path)
    if not os.path.exists(path_only):
        os.makedirs(path_only)     

    fileID = open(abs_file_path,'w')


    
    fileID.write(
                   '<?xml version="1.0" encoding="utf-8"?>\n' 
                   '    <kml xmlns="http://www.opengis.net/kml/2.2">\n' 
                   '        <Document>\n' 
                   '            <name>%s</name>\n' 
                   % filename)

    for j in np.arange(d):
        fileID.write(
                   '            <Placemark>\n' 
                   '                <Snippet maxLines="0"> </Snippet>\n' 
                   '                <description> </description>\n' 
                   '                <name>Point %d</name>\n' 
                   '                <Point>\n' 
                   '                    <coordinates> '
                   % j)
        fileID.write(
                   '%0.15f,%0.16f,0 '
                   % (maxcoords[j,1],maxcoords[j,0]))
        fileID.write(
                   '</coordinates>\n' 
                   '                </Point>\n' 
                   '            </Placemark>\n')
    fileID.write(
                   '        </Document>\n' 
                   '    </kml>')
    fileID.close()
     
     
    azimuths = np.arange(d)     
    elevations = elAngles
    d = {'azimuth':azimuths,'elevation':elevations}
    ar = pd.DataFrame.from_dict(d)
    
    filename = filenamecsv + '.csv'
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    abs_file_path = os.path.join(script_dir, filename)
    path_only = os.path.dirname(abs_file_path)# abspath(abs_file_path)
    if not os.path.exists(path_only):
        os.makedirs(path_only)
    ar.to_csv(abs_file_path)
    
    retval = collections.namedtuple('retval',['az','el'])
    pback = retval(azimuths,elevations)
    return(pback)



