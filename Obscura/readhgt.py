
"""
Created on Mon Feb 06 10:51:33 2017

@author: Boxx
"""
import numpy as np
import pandas as pd
import Tkinter
import tkFileDialog
import os
from io import BytesIO
from urllib import urlopen
from zipfile import ZipFile

def run(**kwargs):
    '''
    function varargout = readhgt(varargin)
    %READHGT Import/download NASA SRTM data files (.HGT).
    %	READHGT(area=) where AREA is a 4-element vector [LAT1,lAT2,LON1,LON2]
    %	downloads the SRTM data and plots a map corresponding to the geographic
    %	area defined by latitude and longitude limits (in decimal degrees). If 
    %	the needed SRTM .hgt files are not found in the current directory (or 
    %	in the path), they are downloaded from the USGS data server (needs an 
    %	Internet connection and a companion file "readhgt_srtm_index.txt"). For
    %	better plot results, it is recommended to install DEM personal function
    %	available at author's Matlab page. 
    %
    %	READHGT(lat=,lon=) reads or downloads the SRTM tiles corresponding to LAT
    %	and LON (in decimal degrees) coordinates (lower-left corner).
    %
    %	LAT and/or LON can be vectors: in that case, tiles corresponding to all
    %	possible combinations of LAT and LON values will be downloaded, and
    %	optional output structure X will have as much elements as tiles.
    %
    %	READHGT(lat=,lon=,outdir=) specifies output directory OUTDIR to write
    %	downloaded files.
    %
    %	READHGT(lat=,lon=,outdir=,url=) specifies the URL address to find HGT 
    %	files (default is USGS).
    %
    %	READHGT(filename=) reads HGT data file FILENAME, must be in the form
    %	"[N|S]yy[E|W]xxx.hgt[.zip]", as downloaded from SRTM data servers.
    %
    %
    %	X=READHGT(...) returns a structure X containing: 
    %		lat: coordinate vector of latitudes (decimal degree)
    %		lon: coordinate vector of longitudes (decimal degree)
    %		  z: matrix of elevations (meters, INT16 class)
    %		hgt: downloaded filename(s)
    %
    %
    % X=READHGT(...,crop=) if value is left blank, crops the resulting map 
    % around existing land (reduces any sea or novalue areas at the borders).
    % If crop value is [LAT1,LAT2,LON1,LON2], crops the map using latitude/lon
    % limits.  perfer to use new syntax READHGT(AREA=).
    %
    %
    %	--- Additionnal options ---
    %
    % 'plot'
    %    Also plots the tile(s).
    %
    %	'tiles'
    %	   Imports and plots individual tiles instead of merging them (default 
    %	   behavior if adjoining values of LAT and LON).
    %
    %	'interp'
    %	   Linearly interpolates missing data.
    %
    %	'crop'
    %	   crops the resulting map around existing land (reduces any sea or 
    %	   novalue areas at the borders).
    %
    %	'crop',[lat1,lat2,lon1,lon2]
    %	   Former syntax that crops the map using latitude/longitude limits. 
    %	   Prefer the new syntax READHGT(AREA).
    %
    %	'srtm1'
    %	   Downloads SRTM1 tiles which are 9 times bigger than default SRTM3 
    %	   ! EXPERIMENTAL ! since the used URL from USGS seems temporary.
    %
    %	'srtm3'
    %	   Forces SRTM3 download (by default, SRTM1 tile is downloaded only for
    %	   USA, if exists).
    %
    %
    %	--- Examples ---
    %
    %	- to plot a map of the Paris region, France (single tile):
    %		readhgt(48,2)
    %
    %	- to plot a map of Flores volcanic island, Indonesia (5 tiles):
    %		readhgt(-9,119:123)
    %
    %	- to plot a map of the Ubinas volcano, Peru (SRTM1 cropped tile):
    %	   readhgt([-16.4,-16.2,-71.5,-71.3],'srtm1','interp')
    %
    %	- to download SRTM1 data of Cascade Range (27 individual tiles):
    %		X=readhgt(40:48,-123:-121,'tiles');
    %
    %
    %	--- Informations ---
    %
    %	- each file corresponds to a tile of 1x1 degree of a square grid
    %	  1201x1201 of elevation values (SRTM3 = 3 arc-seconds), and for USA  
    %	  territory or when using the 'srtm1' option, at higher resolution 
    %	  3601x3601 grid (SRTM1 = 1 arc-second). Note that SRTM1 and SRTM3 
    %	  files have the same syntax names; only the size differs.
    %
    %	- elevations are of class INT16: sea level values are 0, unknown values
    %	  equal -32768 (there is no NaN for INT class), use 'interp' option to
    %	  fill the gaps.
    %
    %	- note that borders are included in each tile, so to concatenate tiles
    %	  you must remove one row/column in the corresponding direction (this
    %	  is made automatically by READHGT when merging tiles).
    %
    %	- downloaded file is written in the current directory or optional  
    %	  OUTDIR directory, and it remains there. Take care that mixed SRTM1
    %	  and SRTM3 files may lead to fail to merge.
    %
    %	- NASA Shuttle Radar Topography Mission [February 11 to 22, 2000] 
    %	  produced a near-global covering on Earth land, but still limited to 
    %	  latitudes from 60S to 60N. Offshore tiles will be output as flat 0
    %	  value grid.
    %
    %	- if you look for other global topographic data, take a look to ASTER
    %	  GDEM, worldwide 1 arc-second resolution (from 83S to 83N): 
    %	  http://gdex.cr.usgs.gov/gdex/ (free registration required)
    %
    %
    %	Author: Fran?ois Beauducel <beauducel@ipgp.fr>
    %		Institut de Physique du Globe de Paris
    %
    %	References:
    %		http://dds.cr.usgs.gov/srtm/version2_1
    %
    %	Acknowledgments: Yves Gaudemer, Jinkui Zhu, Greg
    %
    %	Created: 2012-04-22 in Paris, France
    %	Updated: 2016-05-06
    
    %	Copyright (c) 2016, Fran?ois Beauducel, covered by BSD License.
    %	All rights reserved.
    %
    %	Redistribution and use in source and binary forms, with or without 
    %	modification, are permitted provided that the following conditions are 
    %	met:
    %
    %	   * Redistributions of source code must retain the above copyright 
    %	     notice, this list of conditions and the following disclaimer.
    %	   * Redistributions in binary form must reproduce the above copyright 
    %	     notice, this list of conditions and the following disclaimer in 
    %	     the documentation and/or other materials provided with the distribution
    %	                           
    %	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
    %	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
    %	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
    %	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
    %	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
    %	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
    %	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
    %	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
    %	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
    %	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
    %	POSSIBILITY OF SUCH DAMAGE.
    '''

    fidx = 'readhgt_srtm_index.txt'
    # ATTENTION: this file must exist in the Matlab path to use default SRTM3 tiles
    # since USGS delivers data continent-by-continent with nominative directories,
    # this index file is needed to know the full path name of each tile.
    sz1 = 3601   # SRTM1 tile size
    sz3 = 1201   # SRTM3 tile size
    novalue =  -32768
    n = 1

    if 'options' in kwargs:
        options = kwargs['options']
        makeplot = True if 'plot' in options else False
        merge = True if 'merge' in options else False
        tiles = True if 'tiles' in options else False
        inter = True if 'interp' in options else False
        
    if 'srtm1' in options:
        # EXPERIMENTAL: SRTM1 full resolution tiles available here (2016):
        url = 'http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11'
        srtm1 = True
        srtm3 = False
        sz = 3601
    else:
        # official USGS SRTM3 tiles (and SRTM1 for USA):
        url = 'http://dds.cr.usgs.gov/srtm/version2_1'
        srtm3 = True
        srtm1 = False
        sz = 1201
        
    
    

    if 'crop' in kwargs:
        if kwargs['crop'] == '':
            cropflag = 1
        else:
            crop = kwargs['crop']
            lat = np.floor(np.arange(min(crop[0],crop[1]),max(crop[0],crop[1])))
            lon = np.floor(np.arange(min(crop[2],crop[3]),max(crop[2],crop[3])))
            cropflag = 2
    else:
        cropflag = 0        
              

    # syntax READHGT without argument: opens the GUI to select a file
    if not kwargs:
        # Make a top-level instance and hide since it is ugly and big.
        root = Tkinter.Tk()
        root.withdraw()
        
        # Make it almost invisible - no decorations, 0 size, top left corner.
        root.overrideredirect(True)
        root.geometry('0x0+0+0')
        
        # Show window again and lift it to top so it can get focus,
        # otherwise dialogs will end up behind the terminal.
        root.deiconify()
        root.lift()
        root.focus_force()
        
        options = {}    
        options['filetypes'] = [('HGT files','.hgt'),('HGT zip','.hgt.zip')]
        f = tkFileDialog.askopenfilenames(parent=root,**options) # Or some other dialog
        
        # Get rid of the top-level instance once to make it actually invisible.
        root.destroy()
        


    # syntax READHGT(FILENAME, ...)
    if 'filename' in kwargs:
        f = kwargs['filename']
        try:
            test = open(f,'r')
            test.close()
        except IOError:
            print('Cannot open file: ' + f)
    
            
    if 'lat' in kwargs and 'lon' in kwargs:
        lat = kwargs['lat']
        lon = kwargs['lon']
        
    if 'url' in kwargs:
        url = kwargs['url']
        
    try:
        vl = len(lat)**2
    except TypeError:
        vl = 1
        
#        if ~isnumeric(lon) || ~isnumeric(lat) || any(abs(lat) > 60) || any(lon < -180) || any(lon > 179) || isempty(lat) || isempty(lon)
#            error('LAT and LON must be numeric and in valid SRTM interval (abs(LAT)<60).');
#        end
#        if ~tiles && (any(diff(lat) ~= 1) || any(diff(lon) ~= 1))
#            fprintf('READHGT: Warning! LAT and LON vectors do not define adjoining tiles. Cannot merge and force TILES option.');
#            tiles = 1;
#        end
            
    if 'outdir' in kwargs:
        outdir = kwargs['outdir']
        if not os.path.isdir(outdir):
            print('Out directory is invalid. Using default')
            outdir = '.'
    else:
        outdir = '.'
                

    lon = int(np.floor(lon))
    lat = int(np.floor(lat))
    
    lat,lon = np.meshgrid(lat,lon)
    f = []
    
    for i in range(vl):
        for j in range(vl):
            if lat[i][j] < 0:
                slat = 'S%02d' % -lat[i][j]
            else:
                slat = 'N%02d' % lat[i][j]
                
            if lon[i][j] < 0:
                slon = 'W%03d' % -lon[i][j]
            else:
                slon = 'E%03d' % lon[i][j]

            f.append('%s/%s%s.hgt' % (outdir, slat, slon) )   
            
            try:
                ofile = open(f[i+j])
                ofile.close()

            except IOError:
                if srtm1:
                    ff = '/%s%s.SRTMGL1.hgt.zip' % (slat,slon)
                else:
                    fsrtm = fidx
                    try:
                        fid = open(fsrtm,'rU')
                        k = []
                        for lines in fid:
                            if slat+slon in lines:
                                k.append(lines)
                        if not k:
                            print('READHGT: Warning! Cannot find %s%s tile in SRTM database. Consider it offshore...\n' % (slat,slon))                           
                        else:
                            if srtm3:
                                ff = k[1]
                            else:
                                ff = k[len(k)]
                    except IOError:
                        print('Cannot find "%s" index file to parse SRTM database. Please download HGT file manually.' % fsrtm)            				
                    print('Download %s%s ... ' % (url,ff))
                    
                zipurl = url + ff
                zipresp = urlopen(zipurl)
                zfile = ZipFile(BytesIO(zipresp.read()))
                f[i+j] = zfile.namelist()
                zfile.extractall(outdir)
#                with urlopen(zipurl) as zipresp:
#                    with ZipFile(BytesIO(zipresp.read())) as zfile:
#                        f[i+j] = zfile.namelist()
#                        zfile.extractall(outdir)
#                zipresp = urlopen(zipurl)
#                tempzip = open('tzip.zip','wb')
#                tempzip = write(zipresp.read())
#                tempzip.close()
#                zf = ZipFile('tzip.zip')
#                zf.extractall(path=out)
#                zf.close()
                
                print('done.\n')
    
    
    #pre-allocates X structure (for each file/tile)
    X = pd.DataFrame(index=np.arange(vl),columns=['hgt','lat','lon','z'])
    
    tiles = False if vl == 1 else True
    
    for i in np.arange(vl):
        #unzip HGT file if needed
        if '.zip' in f[i]:
            with ZipFile(f[i]) as zfile:
                    X.loc[i]['hgt'] = zfile.namelist()
                    zfile.extractall(outdir)
        else:
            X.loc[i]['hgt'] = f[i]
                        
        if f[i] == '':
            #offshore empty tile...
            X.loc[i]['z'] = []
        else:
                       
            #loads data from HGT file
            fl = open(X.loc[i]['hgt'],'rb')
            data = np.fromfile(fl,np.dtype('>i2'))
            numel = int(np.sqrt(len(data)))
                        
            if numel == sz1:
                if srtm3:                    
                    z = data.reshape(sz1,sz1)
                    X.loc[i]['z'] = z[::3,::3]  #select every 3rd row
                    
                    X.loc[i]['lon'] = np.linspace(lon[i],lon[i]+1,sz3)
                    X.loc[i]['lat'] = np.linspace(lat[i],lat[i]+1,sz3)
                elif srtm1:
                    X.loc[i]['z'] = data.reshape(sz1,sz1)
            
                    X.loc[i]['lon'] = np.linspace(lon[i],lon[i]+1,sz1)
                    X.loc[i]['lat'] = np.linspace(lat[i],lat[i]+1,sz1)
            elif numel == sz3:
                X.loc[i]['z'] = data.reshape(sz3,sz3)
            
                X.loc[i]['lon'] = np.linspace(lon[i],lon[i]+1,sz3)
                X.loc[i]['lat'] = np.linspace(lat[i],lat[i]+1,sz3)
                
    
    return X
    

#if ~tiles
#	% NOTE: cannot merge mixted SRTM1 / SRTM3 or discontiguous tiles
#	Y.lat = linspace(min(lat(:)),max(lat(:))+1,size(lat,1)*(sz(1)-1)+1)';
#	Y.lon = linspace(min(lon(:)),max(lon(:))+1,size(lon,2)*(sz(2)-1)+1);
#	Y.z = zeros(length(Y.lat),length(Y.lon),'int16');
#	for n = 1:numel(X)
#		if ~isempty(X(n).z)
#			Y.z((sz(1)-1)*(X(n).lat(1)-Y.lat(1)) + (1:sz(1)),(sz(2)-1)*(X(n).lon(1)-Y.lon(1)) + (1:sz(2))) = X(n).z;
#		end
#	end
#
#	if cropflag
#		if cropflag == 1 || isempty(crop)
#			klat = firstlast(any(Y.z ~= 0 & Y.z ~= novalue,2));
#			klon = firstlast(any(Y.z ~= 0 & Y.z ~= novalue,1));
#		else
#			crop = [minmax(crop(1:2)),normlon(minmax(crop(3:4)))];
#			klat = find(Y.lat >= crop(1) & Y.lat <= crop(2));
#			klon = find(Y.lon >= crop(3) & Y.lon <= crop(4));
#		end			
#		Y.lat = Y.lat(klat);
#		Y.lon = Y.lon(klon);
#		Y.z = Y.z(klat,klon);
#	end
#	
#	if inter
#		Y.z = fillgap(Y.lon,Y.lat,Y.z,novalue);
#	end
#end
#
#if nargout == 0 || makeplot
#	if ~tiles
#		fplot(Y.lon,Y.lat,Y.z,decim,url,novalue)
#	else
#		for n = 1:numel(X)
#			fplot(X(n).lon,X(n).lat,X(n).z,decim,url,novalue)
#		end
#	end
#end
#
#if nargout == 3 % for backward compatibility...
#	varargout{1} = X(1).lon;
#	varargout{2} = X(1).lat;
#	varargout{3} = X(1).z;
#elseif nargout > 0
#	if tiles
#		varargout{1} = X;
#	else
#		varargout{1} = Y;
#	end
#	if nargout == 2
#		varargout{2} = f{1}; % for backward compatibility...
#	end
#end
#
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#function fplot(x,y,z,decim,url,novalue)
#%FPLOT plot the data using DEM function if exists, or IMAGESC
#
#demoptions = {'latlon','legend','lake','nodecim'};
#
#figure
#if decim
#	n = decim;
#else
#	n = ceil(sqrt(numel(z))/1201);
#end
#if n > 1
#	x = x(1:n:end);
#	y = y(1:n:end);
#	z = z(1:n:end,1:n:end);
#	fprintf('READHGT: In the figure data has been decimated by a factor of %d...\n',n);
#end
#
#if exist('dem','file')
#	dem(x,y,z,demoptions{:})
#else
#	warning('For better results you might install the function dem.m from http://www.ipgp.fr/~beaudu/matlab.html#DEM')
#	z(z==novalue) = 0;
#	imagesc(x,y,z);
#	if exist('landcolor','file')
#		colormap(landcolor(256).^1.3)
#	else
#		colormap(jet)
#	end
#	% aspect ratio (lat/lon) is adjusted with mean latitude
#	xyr = cos(mean(y)*pi/180);
#	set(gca,'DataAspectRatio',[1,xyr,1])
#
#	orient tall
#	axis xy, axis tight
#end
#
#title(sprintf('Data SRTM/NASA from %s',url),'FontSize',10,'Interpreter','none')
#
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#function y = firstlast(x)
#
#k = find(x);
#y = k(1):k(end);
#
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#function y = minmax(x)
#
#y = [min(x(:)),max(x(:))];
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#function z = fillgap(x,y,z,novalue)
#% GRIDDATA is not efficient for large arrays, but has great advantage to be
#% included in Matlab core functions! To optimize interpolation, we
#% reduce the number of relevant data by building a mask of surrounding
#% pixels of novalue areas... playing with linear index!
#
#sz = size(z);
#k = find(z == novalue);
#k(k == 1 | k == numel(z)) = []; % removes first and last index (if exist)
#if ~isempty(k)
#	[xx,yy] = meshgrid(x,y);
#	mask = zeros(sz,'int8');
#	k2 = ind90(sz,k); % k2 is linear index in the row order
#	% sets to 1 every previous and next index, both in column and row order
#	mask([k-1;k+1;ind90(fliplr(sz),[k2-1;k2+1])]) = 1; 
#	mask(k) = 0; % removes the novalue index
#	kb = find(mask); % keeps only border values
#	z(k) = int16(griddata(xx(kb),yy(kb),double(z(kb)),xx(k),yy(k)));
#end
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#function k2 = ind90(sz,k)
#
#[i,j] = ind2sub(sz,k);
#k2 = sub2ind(fliplr(sz),j,i); % switched i and j: k2 is linear index in row order
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#function y = normlon(x)
#% normalize longitude between -180 and 180
#y = mod(x+180,360) - 180;
