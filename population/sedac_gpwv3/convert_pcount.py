# -*- coding: utf-8 -*- 
'''
Translating SEDAC GPW v3 Population Count Grid ascii files into netCDF format
(http://sedac.ciesin.columbia.edu/data/set/gpw-v3-population-count, http://dx.doi.org/10.7927/H4639MPP)

Code is written and tested under Python 2.7 

@license: This work is licensed under a Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
@author: Mate Rozsai
'''
import sys
import datetime
import csv
import gzip
import numpy
import netCDF4
    
def convertSEDACpcountAscii2nc(asciiGzFile, ncFile):
    def log(message):
        sys.stdout.write(message)

    def createAllMaskedArray(dimA, dimB, defValue):
        res = numpy.ma.array(numpy.zeros([dimA, dimB])) 
        res[:,:] = defValue
        return res

    def saveSEDACncFile(lats, lons, variable, ncFile, NODATA_value):
        #write result to netCDF file
        nc = netCDF4.Dataset(ncFile, 'w', format='NETCDF4')
        
        nc.createDimension('lat', lats.size)
        nc.createDimension('lon', lons.size)
        
        rvLat = nc.createVariable('lat','f8',('lat',))
        rvLat.setncattr('standard_name', 'latitude')
        rvLat.setncattr('long_name', 'latitude')
        rvLat.setncattr('axis', 'Y')
        rvLat.units = 'degrees_north'
        rvLat[:] = lats
        
        rvLon = nc.createVariable('lon','f8',('lon',))
        rvLon.setncattr('standard_name', 'longitude')
        rvLon.setncattr('long_name', 'longitude')
        rvLon.setncattr('axis', 'X')
        rvLon.units = 'degrees_east'
        rvLon[:] = lons
        
        rvVar = nc.createVariable('pcount','f4',('lat','lon',), fill_value=1e+20)
        rvVar.setncattr('long_name', 'Population counts in 2000 adjusted to match UN totals (SEDAC GPWv3)')
        rvVar.units = 'persons'
        rvVar.setncattr('standard_name', 'population')
        rvVar[:,:] = variable[:,:]
        nc.close()
        
    def loadAsciiFile(asciiGzFile, dtype):
        # ****** loading the ascii file   
        log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   reading sedac ascii file (can take a while).')
        # checking header rows
        with gzip.open(asciiGzFile, 'rb') as fcsv:
            reader = csv.reader(fcsv, delimiter=' ', quoting=csv.QUOTE_NONE)
            row = next(reader)
            if 'ncols'!=row[0]: raise Exception('IO error', 'Invalid 1st line: does not start with "ncols".')
            else: ncols = int(row[-1])
            row = next(reader)
            if 'nrows'!=row[0]: raise Exception('IO error', 'Invalid 2nd line: does not start with "nrows".')
            else: nrows = int(row[-1])
            row = next(reader)
            if 'xllcorner'!=row[0]: raise Exception('IO error', 'Invalid 3rd line: does not start with "xllcorner".')
            else: xllcorner = float(row[-1])
            if -180!=xllcorner: raise Exception('IO error', 'Invalid 3rd line: "xllcorner" must be -180 (correct the script otherwise).')
            row = next(reader)
            if 'yllcorner'!=row[0]: raise Exception('IO error', 'Invalid 4th line: does not start with "yllcorner".')
            else: yllcorner = float(row[-1])
            if 0<yllcorner: raise Exception('IO error', 'Invalid 4th line: "yllcorner" must be negative (correct the script otherwise).')
            row = next(reader)
            if 'cellsize'!=row[0]: raise Exception('IO error', 'Invalid 5th line: does not start with "cellsize".')
            row = next(reader)
            if 'NODATA_value'!=row[0]: raise Exception('IO error', 'Invalid 6th line: does not start with "NODATA_value".')
            else: NODATA_value = row[-1]
        # reading data
        log('.')
        var_asc = numpy.ma.masked_values(numpy.loadtxt(asciiGzFile, dtype=dtype, skiprows=6), int(NODATA_value))
        log('.\n')
        return ncols, nrows, yllcorner, NODATA_value, var_asc 

    def createNcFile(ncols, nrows, yllcorner, NODATA_value, var_asc, defValue, ncFile):
        # ****** generating nc file
        log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+"]   creating 2.5' nc file...")
        # calculating the lat and lon coordinates
        lons = numpy.zeros([ncols])
        lats = numpy.zeros([int(ncols/2)])
        gridSize = 360.0 / ncols             
        for i in range(lons.size): lons[i] = -180.0 + gridSize/2 + i*gridSize
        for i in range(lats.size): lats[i] = 90.0 - gridSize/2 - i*gridSize
        # reordering ascii variable according to lats and lons
        var_nc = createAllMaskedArray(lats.size, lons.size, defValue)
        startLat = int(round((90.0-(yllcorner+gridSize*nrows))/gridSize-1))
        stopLat = int(round((90.0-yllcorner)/gridSize-1))
        for idxLat in range(startLat,stopLat): 
            var_nc[idxLat,:]=var_asc[idxLat-startLat,:]
            if 0 == idxLat % int((stopLat-startLat)/20): log('.')
        saveSEDACncFile(lats, lons, var_nc, ncFile + '_25.nc', int(NODATA_value))
        log('..\n')
        
    def createHalfDegreeNcFile(ncols, nrows, yllcorner, NODATA_value, glbnds_asc, ncFile):
        # create another resolution - half degree grid - from the original data
        log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   creating 1/2° nc file...')
        # calculating the lat and lon coordinates
        gridSize = 0.5
        lons = numpy.zeros([int(360/gridSize)])
        lats = numpy.zeros([int(180/gridSize)])
        for i in range(lons.size): lons[i] = -180 + gridSize/2 + i*gridSize
        for i in range(lats.size): lats[i] = 90 - gridSize/2 - i*gridSize
        # creating the result array with lower resolution
        glbnds_half_nc = createAllMaskedArray(lats.size, lons.size, 0.0)
        startLat = int(round((90-(yllcorner+nrows*(360.0/ncols)))/gridSize-1))
        stopLat = int(round((90-yllcorner)/gridSize-1))
        dv = (int)(ncols/lons.size)                  # reshaping factor
        for idxLat in range(startLat,stopLat):
            idxOriLat = (idxLat-startLat)*dv
            for idxLon in range(lons.size):
                idxOriLon = idxLon*dv
                # the new value is the sum of the dv*dv square from the original data 
                glbnds_half_nc[idxLat,idxLon]=numpy.sum(glbnds_asc[idxOriLat:idxOriLat+dv,idxOriLon:idxOriLon+dv])
            if 0 == idxLat % int((stopLat-startLat)/20): log('.')
        saveSEDACncFile(lats, lons, glbnds_half_nc, ncFile+'_half.nc', int(NODATA_value))
        log('.\n')

    # Converting the SEDAC ascii file to netCDF file
    ncols, nrows, yllcorner, NODATA_value, var_asc = loadAsciiFile(asciiGzFile, 'float')
    createNcFile(ncols, nrows, yllcorner, NODATA_value, var_asc, 0.0, ncFile)
    createHalfDegreeNcFile(ncols, nrows, yllcorner, NODATA_value, var_asc, ncFile)
    log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   pcount translation ready\n')

if __name__ == "__main__":
    convertSEDACpcountAscii2nc('gl_gpwv3_pcount_00_ascii_25/glp00ag.asc.gz', 'results/pcount')
