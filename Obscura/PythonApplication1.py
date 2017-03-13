import azimuthLines
import readhgt
import sys
import pandas as pd
import radomeObscura
import readsrtm


def radome(X,Z,rad,beam,azcorr,name):
    #MUE params and buffers
    meters = 1
    feet = 1
    degrees = 1
    beamwidth = (beam*meters)*(3.28125*feet/meters)
    humanbuffer = 10*feet
    totalbuffer = 0.5*beamwidth + humanbuffer

    #User variables
    x0 = X*feet             #distance from transmitter to radome center
    y0 = 0                    #will always be zero for simplicity
    z0 = Z*feet        #Height difference between center of beam and radome center
    R = rad + totalbuffer*feet  #Radius of obscura radome
    azCorrection = azcorr*degrees #azimuth to obscura radome
    filename = name
    
    az,el = radomeObscura.run(dist=x0,height=z0,radius=R,azcorr=azCorrection,fname=filename)

    return az,el


def main():
    lat = 21.564736;
    lon = -158.239802;
    #coord = azimuthLines.run(25,0.01,lat,lon,r'TDY\kirt-terra-azlines')
    #cdf = pd.Panel(data=coord)
    #cdf.to_pickle(path='.\coord.pkl')
    cdf = pd.read_pickle(path='.\coord.pkl')

    #X = readhgt.run(lat=cdf[0,0,0,],lon=cdf[0,1,0],options=['srtm3'])
    #X.to_pickle('.\X.pkl')
    #X = pd.read_pickle(path='.\X.pkl')
    
    azimuths,elevations = readsrtm.run(cdf,0.01,0.01,r'TDY\testkml',r'TDY\testcsv')

    azA,elA = radome(1160,30.5-107.5,47,13.56,223.75,r'TDY\radomeA')


    print 'done'
    kill = raw_input('Hit anything')

if __name__ == "__main__":
    sys.exit(int(main() or 0))

#%%

#
#