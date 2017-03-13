import numpy as np
import os

#Created on Fri Feb 03 14:16:52 2017

#@author: Boxx

#@profile
def run(d_km,res_km,lat,lon,filename):
    '''
    AZIMUTHLINES   Creates 360 azimuth lines radiating from a single point
    [COORD] = AZIMUTHLINES(d_km,res_km,lat1,lon1,fn)  This function will
    create 360 lines, each one representing a different azimuth, radiating
    from a source point.  It will also create a kml file with these lines.
    d_km is the distance of each radiated line in kilometers (also the 
    radius of circle). res_km is the resolution or distance between points 
    in kilometers.  lat1 is the latitude of the source point in decimal 
    degrees.  lat2 is the longitude of the source point in decimal degrees.  
    fn is the filename for the generated kml file.
    
    Example:  Generate azimuth lines from source point with latitude 
                21.564736 degrees and longitude -158.239802 degrees.  Use a
                resolution of 10 meters and a radiating distance of 25
                kilometers.  Use a filename of testazlines.
    
    lat = 21.564736;
    lon = -158.239802;
    [coord] = azimuthLines(25,0.01,lat,lon,'testazlines');
    
    Auther: William Boxx, Electronics Engineer, 85 EIS USAF
     
    Formula used to find next lat long point
    lat2 = asin(sin(lat1)*cos(d/R) + cos(lat1)*sin(d/R)*cos(theta))
    lon2 = lon1 + atan2(sin(theta)*sin(d/R)*cos(lat1), cos(d/R)-sin(lat1)*sin(lat2))
    
    asin          = arc sin()   
    d             = distance (in any unit)   
    R             = Radius of the earth (in the same unit as above)  
    and hence d/r = is the angular distance (in radians)  
    atan2(a,b)    = arc tan(b/a)  
    is the bearing (in radians, clockwise from north); 
    '''
    print 'Start azimuthLines.py'
    fname = filename + '.kml'
    km = 1
    rad = 1
    #d = 25*km;
    R =	6378.1*km
    lat1 =	lat*np.pi/180.0*rad    #21.564736*pi/180*rad;
    lon1 =	lon*np.pi/180.0*rad    #-158.239802*pi/180*rad;
    theta = np.arange(360)*np.pi/180
    #res_km = 0.01*km;
    numres = np.floor(d_km/res_km)
    coord = np.zeros([int(numres),2,360])
    index = 0
    
    print 'Start loop'
    #k = np.arange(res_km,d_km+res_km,res_km)
    for x in xrange(int(numres)):
        print x
        k = res_km + x*res_km
        lat2 = np.arcsin(np.sin(lat1)*np.cos(k/R) + np.cos(lat1)*np.sin(k/R)*np.cos(theta));
        lon2 = lon1 + np.arctan2(np.sin(theta)*np.sin(k/R)*np.cos(lat1),np.cos(k/R)-np.sin(lat1)*np.sin(lat2))
      
        #coord = np.matlib.repmat(np.array(lat2*180/np.pi , lon2*180/np.pi),1,360)
        coord[index,:,:] = [lat2*180/np.pi,lon2*180/np.pi]
        index = index + 1
    
    
    #
    # - Create kml
    #
    
    #filename = 'azLines'

    print 'Finished loop.'
    print 'Writing to file: ' + filename + '.kml'

    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    abs_file_path = os.path.join(script_dir, fname)
    path_only = os.path.dirname(abs_file_path)# abspath(abs_file_path)
    if not os.path.exists(path_only):
        os.makedirs(path_only)

    fileID = open(abs_file_path,'w')
    fileID.write(  
                    '<?xml version="1.0" encoding="utf-8"?>\n'
                    '    <kml xmlns="http://www.opengis.net/kml/2.2">\n'
                    '        <Document>\n'
                    '            <name>%s</name>\n'
                    % fname)
    for j in xrange(0,360):
        fileID.write(
                '            <Placemark>\n'
                '                <Snippet maxLines="0"> </Snippet>\n'
                '                <description> </description>\n'
                '                <name>Line %d</name>\n'
                '                <LineString>\n'
                '                    <coordinates>' % j)
        for k in xrange(0,int(numres)):
            fileID.write(
                "%0.15f,%0.16f,0 " % (coord[k,1,j],coord[k,0,j]))

        fileID.write(
                "</coordinates>\n"
                "                </LineString>\n" 
                "            </Placemark>\n")

    fileID.write(
                "        </Document>\n" 
                "    </kml>")
    fileID.close
    

    return coord

