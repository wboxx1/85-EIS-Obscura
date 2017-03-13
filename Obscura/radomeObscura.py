import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os

def run(**kwargs):
    '''
    #    function [azimuths,elevations] = radomeObscura(d,h,r,azC,fn,th,res)
    #%RADOMEOBSCURA      Create azimuth and elevations for radome obscura
    #%   [AZIMUTHS,ELEVATIONS] = RADOMEOBSCURA(d,h,r,azC,fn) Creates azimuth and
    #%   elevation arrays for a radome obscura defined by d,h,r and azC. d is
    #%   the distance from the transmitter to the radome center.  h is the
    #%   height difference between center of beam and radome center (value is
    #%   negative if radome center is above center of beam).  r is the radius of
    #%   the obscura radome.  fn is the filename used for the .csv file that is
    #%   created.  All units must be the same.
    #%
    #%   [AZIMUTHS ELEVATIONS] = RADOMEOBSCURA(d,h,r,azC,fn,th) You can specify
    #%   a threshold value instead of the default one of 5.  This value is
    #%   used as a error correction because we are using discrete values.
    #%
    #%   [AZIMUTHS ELEVATIONS] = RADOMEOBSCURA(d,h,r,azC,fn,th,res)  You can
    #%   specify a resolution value instead of the default one of 360.  This
    #%   value determines the number of points in the meshgrid as res^2.
    #%
    #%   Example: Find the obscura azimuth and elevations vectors for a 20 meter
    #%               radius radome that is 60 meters away from the radiating 
    #%               source and has a center 6 meters above it.  The azimuth
    #%               correction is 300 degrees and filename will be radometest.
    #%
    #%       [azimuths,elevations] = radomeobscura(60,-6,20,300,'radometest');
    '''

    if 'dist' in kwargs:
        x0 = kwargs['dist']

    if 'height' in kwargs:
        z0 = kwargs['height']

    if 'radius' in kwargs:
        R = kwargs['radius']

    if 'azcorr' in kwargs:
        azCorrection = kwargs['azcorr']

    if 'fname' in kwargs:
        fn = kwargs['fname']

    if 'threshold' in kwargs:
        thresh = kwargs['threshold']
    else:
        thresh = 5
        
    if 'res' in kwargs:
        N = kwargs['res']
    else:
        N = 360

    y0 = 0
    filename = fn  + '.csv'

    k = 0
    firstRun = True

    phi = np.linspace(0.,np.pi/2.,N)
    theta = np.linspace(-np.pi/2.,np.pi/2.,N)

    el = []
    az = []
    point = []

    for i in range(N):
        for j in range(N):
            x1 = R*np.sin(phi[i])*np.cos(theta[j])
            y1 = R*np.sin(phi[i])*np.sin(theta[j])
            z1 = R*np.cos(phi[i])

            #evaluate gradient at point
            dx = 2*x1
            dy = 2*y1
            dz = 2*z1

            #evaluate dot product
            dp = x1*dx + y1*dy + z1*dz

            #evaluate transmit point
            tp = x0*dx + y0*dy + z0*dz

            if abs(dp-tp) < thresh:
                #create array of points
                if firstRun:
                    point = [x1,y1,z1]
                    firstRun = False
                else:
                    point = np.vstack((point,[x1,y1,z1]))

                #establish piont 1 and point 2
                P1 = {'x':x0,'y':y0,'z':z0}
                P2 = {'x':x1,'y':y1,'z':z1}

                #find diff of x and z
                dx = P1['x'] - P2['x']
                dz = P1['z'] - P2['z']
                dy = P1['y'] - P2['y']

                #for elevation, take projections on xz plane
                #calculate arctangent of z/x
                el.append(-np.arctan2(dz,np.sqrt(dx*dx + dy*dy))*180/np.pi)

                #for azimuth - take projections on xy plane
                #calculate arctangent of y/x
                #include azimuth correction
                az.append(-np.arctan2(dy,np.sqrt(dx*dx + dz*dz))*180/np.pi + azCorrection)

                if az[k] > 359:
                    az[k] = az[k] - 360
                k = k + 1

    phi = np.linspace(0,np.pi,N)

    phi,theta = np.meshgrid(phi,theta)
    X = R*np.sin(phi)*np.cos(theta)
    Y = R*np.sin(phi)*np.sin(theta)
    Z = R*np.cos(phi)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X,Y,Z)
    
    
    S = np.array(point).shape
    h = S[0]

    for i in range(h):
        ax.plot([point[i,0]],[point[i,1]],[point[i,2]],marker='x',mew=3,ms=6)

    ax.plot([x0],[y0],[z0],marker='x',mew=3,ms=6)

    #if x0 > 2*R:
    #    upp = np.ceil(x0+10)
    #    if np.mod(upp,2) > 0:
    #        upp = upp + 1
    #    ax.axis([0, upp, -upp/2, upp/2, -upp/2, upp/2])

    plt.show()

    #find values at integer azimuths by spline interpolation
    low = int(np.floor(np.min(az)))
    high = int(np.ceil(np.max(az)))
    azimuths = range(low,high+1)    
    I = np.argsort(az)
    azi = np.squeeze(np.array(az)[I])
    eli = np.squeeze(np.array(el)[I])
    elevations = np.interp(azimuths,azi,eli)

    csvdata =  pd.DataFrame([az,el],index=['az','el']).T

    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    abs_file_path = os.path.join(script_dir, filename)
    path_only = os.path.dirname(abs_file_path)# abspath(abs_file_path)
    if not os.path.exists(path_only):
        os.makedirs(path_only)

    csvdata.to_csv(abs_file_path,index=False)
    
    return azimuths,elevations