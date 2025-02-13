import os
import scipy.io as sio
from astropy.io import fits
import numpy as np
from datetime import datetime
# from matplotlib import pyplot as plt
# from pylab import *
import time
from astropy import constants as const


# https://docs.astropy.org/en/stable/io/fits/usage/table.html
# https://docs.astropy.org/en/stable/io/fits/


'''
Prototype conversion from ALFALFA IDL data structures into a FITS file.
'''

def create_grid(dir_suffix, file_suffix, path_prefix, grid_file_name):

    start = time.time()

    if not os.path.exists(dir_suffix):
        os.makedirs(dir_suffix)

    #dir_suffix  = '1220+09'
    #file_suffix = '1220+09a'
    #path_prefix = '/home/adjuntas2/galaxy/grids/{!s}/'.format(dir_suffix)
    #grid_file_name = path_prefix + 'gridbf_{!s}.sav'.format(file_suffix)

    dictnames = {}
    gridfile = sio.readsav(grid_file_name, idict=dictnames,python_dict=True)

    grid = gridfile['grid'][0]

    name           = grid[0].decode("utf-8")
    ramin          = grid[1]
    decmin         = grid[2]
    epoch          = grid[3]
    deltara        = grid[4]
    deltadec       = grid[5]
    nx             = grid[6]
    ny             = grid[7]
    map_projection = grid[8].decode("utf-8")
    czmin          = grid[9]
    nz             = grid[10]
    velarr         = grid[11]
    wf_type        = grid[12].decode("utf-8")
    wf_fwhm        = grid[13]
    han            = grid[14]
    medsubtract    = grid[15]
    baseline       = grid[16]
    calib_facs     = grid[17]
    grms           = grid[18]
    date           = grid[19].decode("utf-8")
    who            = grid[20].decode("utf-8")
    pos            = grid[21]
    drift_list     = grid[22]
    grid_makeup    = grid[23]
    d              = grid[24]
    w              = grid[25]
    cont           = grid[26]
    cw             = grid[27]

    # ion()

    #---------------------------------------------------------
    # Reshape data cube - 1st polarization
    cube1 =  d[:,:,0,0]
    image = d[:,:,0,1]
    cube1 = np.stack((cube1, image))
    for i in range(2, nz):
        image = d[:,:,0,i]
        cube1 = np.concatenate((cube1, [image]), axis=0)
    

    # Reshape data cube - 2nd polarization
    cube2 =  d[:,:,1,0]
    image = d[:,:,1,1]
    cube2 = np.stack((cube2, image))
    for i in range(2, nz):
        image = d[:,:,1,i]
        cube2 = np.concatenate((cube2, [image]), axis=0)

    cube3 = np.stack((cube1,cube2))

    # cube = np.reshape(cube, (nz, nx, ny))

    #----------------------------------------------------------
    # Reshape data cube weights - 1st polarization
    weights1 =  w[:,:,0,0]
    image = w[:,:,0,1]
    weights1 = np.stack((weights1, image))
    for i in range(2, nz):
        image = w[:,:,0,i]
        weights1 = np.concatenate((weights1, [image]), axis=0)
    

    # Reshape data cube weights - 2nd polarization
    weights2 =  w[:,:,1,0]
    image = w[:,:,1,1]
    weights2 = np.stack((weights2, image))
    for i in range(2, nz):
        image = w[:,:,1,i]
        weights2 = np.concatenate((weights2, [image]), axis=0)

    weights3 = np.stack((weights1,weights2))

    #----------------------------------------------------------
    # Reshape continuum map 
    cont1 = np.stack((cont[:,:,0], cont[:,:,1]))

    #----------------------------------------------------------
    # Reshape continuum weights
    cw1 = np.stack((cw[:,:,0], cw[:,:,1]))

    # Create PRIMARY FITS header
    utcnow = datetime.utcnow().strftime("%B %d %Y %H:%M:%S")
    restfreq = 1420.405751e6   # Rest frequency of neutral hydrogen
    hdr = fits.Header()
    hdr['NAXIS']    = (len(d.shape), 'Number of dimensions')
    hdr['NAXIS1']   = (nx, 'Length of x axis')
    hdr['NAXIS2']   = (ny, 'Length of y axis')
    hdr['NAXIS3']   = (nz, 'Length of velocity axis')
    hdr['NAXIS4']   = (d.shape[2], 'Length of polarization axis')
    hdr['BMAJ']     = (3.8/60.0, 'Major axis ALFALFA beam degrees')
    hdr['BMIN']     = (3.3/60.0, 'Minor axis ALFALFA beam degrees')

    hdr['OBSERVER'] = ('ALFALFA Survey Team', 'Observers')
    hdr['OBJECT']   = (name, 'ALFALFA data cube and grid center')
    hdr['TELESCOP'] = ('Arecibo', 'Telescope')
    hdr['INSTRUME'] = ('ALFA', 'Receiver instrument')
    hdr['BSCALE']   = 1.0
    hdr['BZERO']    = 0.0
    hdr['BUNIT']    = ('mJY/BEAM', 'Units of flux')
    hdr['EPOCH']    = (epoch, 'Epoch of RA DEC')
    hdr['VELREF']   = (2, 'HELIO')  
    hdr['RESTFREQ'] = (restfreq, 'Rest frequency')

    hdr['CTYPE1']   = ('RA---CAR', 'X-axis type')
    hdr['CRVAL1']   = (ramin*15.0, 'Reference pixel value RA Degrees')
    hdr['CDELT1']   = ((deltara / np.cos(decmin*np.pi/180.0))*15.0/3600.0, 'RA Degrees/pixel')
    hdr['CRPIX1']   = (0.0000, 'Reference pixel')
    hdr['CROTA1']   = 0.0

    hdr['CTYPE2']   = ('DEC--CAR', 'Y-axis type')
    hdr['CRVAL2']   = (decmin, 'Reference pixel value Dec Degrees')
    hdr['CDELT2']   = (deltadec/60.0, 'Dec Degrees/pixel')
    hdr['CRPIX2']   = (0.0000, 'Reference pixel')
    hdr['CROTA2']   = 0.0

    hdr['CTYPE3']   = ('FREQ-HEL', 'Z-axis type')

    freqmin = restfreq/((velarr[0]*1000.0 / const.c.value)+1)/1.e6
    freqmin2 = restfreq/((velarr[1]*1000.0 / const.c.value)+1)/1.e6
    deltafreq = freqmin2-freqmin  # Delta MHz

    hdr['CRVAL3']   = (freqmin, 'Reference pixel value frequency (MHz)')
    hdr['CDELT3']   = (deltafreq, 'Channel width (MHz)')
    hdr['CRPIX3']   = (0.0000, 'Reference pixel')
    hdr['CROTA3']   = 0.0

    hdr['CTYPE4']   = ('STOKES', 'Pol-axis type')
    hdr['CRVAL4']   = (1.0, '')
    hdr['CDELT4']   = (1.0, '')
    hdr['CRPIX4']   = (1.0, '')
    hdr['CROTA4']   = 0.0

    hdr['EQUINOX']  = (epoch, 'Equinox of coordinates')
    hdr['ORIGIN']   = ('Cornell University', 'File creation location')
    hdr['SURVEY']   = ('ALFALFA', 'Survey Name')
    hdr['COMMENT']  = ('ALFALFA Dual Polarization Data Cube')
    # hdr['COMMENT']  = ('Intensity units: mJy km/s beam^-1')
    hdr['COMMENT']  = ('FITS header created: {!s}'.format(utcnow))
    hdr['COMMENT']  = ('IDL data cube created: {!s}'.format(date))
    hdr['COMMENT']  = ('Website: http://egg.astro.cornell.edu/alfalfa/')

    hdr['HISTORY'] = ('Data cube name: {!s}'.format(name))
    hdr['HISTORY'] = ('Map projection: {!s}'.format(map_projection))
    hdr['HISTORY'] = ('wf_type: {!s}'.format(wf_type))
    hdr['HISTORY'] = ('wf_fwhm: {!s}'.format(wf_fwhm))
    hanning = 'No'
    if han == 1: hanning = 'Yes'
    hdr['HISTORY'] = ('Hanning? {!s}'.format(hanning))
    mediansubtract = 'No'
    if medsubtract == 1: mediansubtract = 'Yes'
    hdr['HISTORY'] = ('Median subtraction? {!s}'.format(mediansubtract))
    hdr['HISTORY'] = ('GRMS: {!s}'.format(grms))

    hdr['HISTORY'] = 'File list of drift scans'
    hdr['HISTORY'] = '  used to create this data cube:'

    for rec in drift_list:
        hdr['HISTORY'] = 'File: {!s} '.format(rec[0].decode("utf-8"))

    #--------------------------------------------------------------------------
    # CREATE CONTINUUM IMAGE FITS HEADER
    restfreq = 1420.405751e6   # Rest frequency of neutral hydrogen
    chdr = fits.Header()
    chdr['NAXIS']    = (len(cont.shape), 'Number of dimensions')
    chdr['NAXIS1']   = (nx, 'Length of x axis')
    chdr['NAXIS2']   = (ny, 'Length of y axis')
    chdr['NAXIS3']   = (cont.shape[2], 'Length of polarization axis')

    chdr['OBSERVER'] = ('ALFALFA Survey Team', 'Observers')
    chdr['OBJECT']   = (name, 'ALFALFA grid center')
    chdr['TELESCOP'] = ('Arecibo', 'Telescope')
    chdr['INSTRUME'] = ('ALFA', 'Receiver instrument')
    chdr['EPOCH']    = (epoch, 'Epoch of RA DEC')
    chdr['RESTFREQ'] = (restfreq, 'Rest frequency')

    chdr['CTYPE1']   = ('RA---CAR', 'X-axis type')
    chdr['CRVAL1']   = (ramin*15.0, 'Reference pixel value RA Degrees')
    chdr['CDELT1']   = ((deltara / np.cos(decmin*np.pi/180.0))*15.0/3600.0, 'RA Degrees/pixel')
    chdr['CRPIX1']   = (0.0000, 'Reference pixel')
    chdr['CROTA1']   = 0.0

    chdr['CTYPE2']   = ('DEC--CAR', 'Y-axis type')
    chdr['CRVAL2']   = (decmin, 'Reference pixel value Dec Degrees')
    chdr['CDELT2']   = (deltadec/60.0, 'Dec Degrees/pixel')
    chdr['CRPIX2']   = (0.0000, 'Reference pixel')
    chdr['CROTA2']   = 0.0

    chdr['CTYPE3']   = ('STOKES', 'Pol-axis type')
    chdr['CRVAL3']   = (1.0, '')
    chdr['CDELT3']   = (1.0, '')
    chdr['CRPIX3']   = (1.0, '')
    chdr['CROTA3']   = 0.0

    chdr['EQUINOX']  = (epoch, 'Equinox of coordinates')
    chdr['ORIGIN']   = ('Cornell University', 'File creation location')
    chdr['SURVEY']   = ('ALFALFA', 'Survey Name')
    chdr['COMMENT']  = ('ALFALFA Dual Polarization Continuum Map')
    chdr['COMMENT']  = ('Intensity units: mJy beam^-1')
    chdr['COMMENT']  = ('FITS header created: {!s}'.format(utcnow))
    chdr['COMMENT']  = ('IDL continuum map created: {!s}'.format(date))
    chdr['COMMENT']  = ('Website: http://egg.astro.cornell.edu/alfalfa/')

    chdr['HISTORY'] = ('Data cube name: {!s}'.format(name))
    chdr['HISTORY'] = ('Map projection: {!s}'.format(map_projection))
    chdr['HISTORY'] = ('wf_type: {!s}'.format(wf_type))
    chdr['HISTORY'] = ('wf_fwhm: {!s}'.format(wf_fwhm))
    hanning = 'No'
    if han == 1: hanning = 'Yes'
    chdr['HISTORY'] = ('Hanning? {!s}'.format(hanning))
    mediansubtract = 'No'
    if medsubtract == 1: mediansubtract = 'Yes'
    chdr['HISTORY'] = ('Median subtraction? {!s}'.format(mediansubtract))
    chdr['HISTORY'] = ('GRMS: {!s}'.format(grms))


    # --------------------------------------------------------------
    # Write out a set of cubes

    # Primary Spectral data cube
    spectralfitsfile = '{!s}/{!s}_spectral.fits'.format(dir_suffix,file_suffix)
    print('    Writing spectral cube: {!s}'.format(spectralfitsfile))
    os.system('rm -rf ' + spectralfitsfile)
    spectral_hdu = fits.PrimaryHDU(cube3, header = hdr)
    spectral_hdul = fits.HDUList([spectral_hdu])
    spectral_hdul.writeto(spectralfitsfile)
    # No.    Name      Ver    Type      Cards   Dimensions   Format
    #   0  PRIMARY       1 PrimaryHDU     102   (144, 144, 1024, 2)   float32 

    # Spectral weights HDU
    spweightsfitsfile = '{!s}/{!s}_spectralweights.fits'.format(dir_suffix,file_suffix)
    print('    Writing spectral weights: {!s}'.format(spweightsfitsfile))
    os.system('rm -rf ' + spweightsfitsfile)
    weights_hdu = fits.PrimaryHDU(weights3, header=hdr)
    weights_hdul = fits.HDUList([weights_hdu])
    weights_hdul.writeto(spweightsfitsfile)
    # No.    Name      Ver    Type      Cards   Dimensions   Format
    #   0  PRIMARY       1 PrimaryHDU     102   (144, 144, 1024, 2)   float32 

    # Continuum Map
    contfitsfile = '{!s}/{!s}_continuum.fits'.format(dir_suffix,file_suffix)
    print('    Writing continuum map: {!s}'.format(contfitsfile))
    os.system('rm -rf ' + contfitsfile)
    continuum_hdu = fits.PrimaryHDU(cont1, header=chdr)
    continuum_hdul = fits.HDUList([continuum_hdu])
    continuum_hdul.writeto(contfitsfile)
    # No.    Name      Ver    Type      Cards   Dimensions   Format
    #   0  PRIMARY       1 PrimaryHDU      42   (144, 144, 2)   float32

    # Continuum Weights
    cwfitsfile = '{!s}/{!s}_continuumweights.fits'.format(dir_suffix,file_suffix)
    print('    Writing continuum weights: {!s}'.format(cwfitsfile))
    os.system('rm -rf ' + cwfitsfile)
    continuum_weights_hdu = fits.PrimaryHDU(cw1, header=chdr)
    continuum_weights_hdul = fits.HDUList([continuum_weights_hdu])
    continuum_weights_hdul.writeto(cwfitsfile)
    # No.    Name      Ver    Type      Cards   Dimensions   Format
    #   0  PRIMARY       1 PrimaryHDU      42   (144, 144, 2)   float32

    end = time.time()
    print("Seconds to complete: {!s} seconds".format(end - start))

#------------------------------------------------------------------------------
# dir_suffix  = '1220+09'
# file_suffix = '1220+09a'
# path_prefix = '/home/adjuntas2/galaxy/grids/{!s}/'.format(dir_suffix)
# grid_file_name = path_prefix + 'gridbf_{!s}.sav'.format(file_suffix)

# Example grid name directories
#gridnames = ['0140+03','0740+15','0940+27','1116+33','1300+01',
#              '1436+07','1628+01', '2316+17','0004+03','0804+11',
#              '0852+13']


filename = 'gridlist00a'

with open(filename, 'r') as f:
    allgrids = [line.rstrip() for line in f]

#gridnames = allgrids[::35]
gridnames = allgrids
gridnames = [gridname.strip('/') for gridname in gridnames]


#gridnames = ['1220+09']
#gridnames = ['0732+01']
#gridnames = ['1004+01']
#gridnames = ['1204+27']
grid_labels = ['a', 'b', 'c', 'd']

count = 0

for dir_suffix in gridnames:
    print(dir_suffix+'... ('+str(count)+')')
    for grid_label in grid_labels:
        file_suffix = dir_suffix + grid_label
        path_prefix = '/home/adjuntas2/galaxy/grids/{!s}/'.format(dir_suffix)
        grid_file_name = path_prefix + 'gridbf_{!s}.sav'.format(file_suffix)

        try:        
            create_grid(dir_suffix, file_suffix, path_prefix, grid_file_name)
        except Exception as e:
            print(e)
    count = count + 1



# Open and manipulate the images
# from astropy.io.fits import getheader
# primary_hdr = getheader(spectralfitsfile, 0)
# weights_hdr = getheader(spweightsfitsfile, 0)
# continuum_hdr = getheader(contfitsfile, 0)
# contweights_hdr = getheader(cwfitsfile, 0)

# from astropy.io.fits import getdata
# datacube, primary_hdr = getdata(spectralfitsfile, 0, header=True)
# weights, weights_hdr = getdata(spweightsfitsfile, 0, header=True)
# cont, continuum_hdr = getdata(contfitsfile, 0, header=True)
# contweights, contweights_hdr = getdata(cwfitsfile, 0, header=True)







