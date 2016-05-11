# -*- coding: utf-8 -*- 
'''
Sample code of  SEDAC GPW v3 dataset translation aggregations

Code is written and tested under Python 2.7 

@license: This work is licensed under a Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
@author: Mate Rozsai
'''
import csv
import dbf
import netCDF4
import numpy
import sys
import datetime

def log(message):
    sys.stdout.write(message)

def loadRegions(tsvFile):
    # Load regions into a list of tuples, where each tuple has two elements: first the region name, second a list of countries that belong to the region
    log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   loading regions...')
    regions=list()
    with open(tsvFile, 'rb') as fcsv:
        reader = csv.reader(fcsv, delimiter='\t', quoting=csv.QUOTE_NONE)
        for row in reader:
            regions.append((row[0], row[1:]))
    log('............\n')
    return regions

def createCountryMapping(dbfFile, regions):
    # Load SEDAC country list and make a numpy array whose index refers to the country's sedaccode, and the value is the index of the region that it belongs to
    def getRegionIndex(iso3v10, regions):
        # returns -1 if code cannot be found in regions
        i, res = 0, -1
        while i<len(regions) and res<0:
            if iso3v10 in regions[i][1]: res=i
            i+=1
        return res
    log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   creating country mapping...')
    # loading the table
    table = dbf.Table(dbfFile)
    table.open()
    maxid=0
    countries = list()
    for row in table:
        countries.append((row.value, row.iso3v10))
        maxid=max(maxid,row.value)
    table.close()
    # assigning region codes
    cntMapping = numpy.empty([maxid+1], dtype='int') 
    cntMapping[:] = numpy.nan
    for country in countries:
        cntMapping[country[0]] = getRegionIndex(country[1], regions)
    log('...\n')
    return cntMapping

def getTotalPopulation(glbndsFile, pcountFile, regions, cntMapping):
    def loadNc(ncFile, var):
        nc = netCDF4.Dataset(ncFile, 'r')
        ncVar = nc.variables[var][:]
        lons = nc.variables['lon'][:]
        lats = nc.variables['lat'][:]
        return ncVar, lats, lons
        nc.close()
    
    log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   loading nc files...')
    glbnds, latsG, lonsG = loadNc(glbndsFile, 'glbnds')
    pcount, latsP, lonsP = loadNc(pcountFile, 'pcount')
    # check if we use the same grid resolution and ordering
    if not(numpy.array_equal(latsG, latsP)): raise Exception('Different lat coordinates!')
    if not(numpy.array_equal(lonsG, lonsP)): raise Exception('Different lon coordinates!')
    # calculate the total population
    totPop = numpy.zeros(len(regions))
    for idxLat in range(latsG.shape[0]):
        for idxLon in range(lonsG.shape[0]):
            if not glbnds.mask[idxLat,idxLon]:
                regionIdx = cntMapping[glbnds[idxLat,idxLon]]
                if 0<=regionIdx: totPop[regionIdx] += pcount[idxLat, idxLon]
        if 0 == idxLat % int(latsG.shape[0]/10): log('.')
    log('.\n')
    return totPop

if __name__ == "__main__":
    regions = loadRegions('regions.tsv')
    cntMapping = createCountryMapping('../gl_gpwv3_ntlbndid_ascii_25/bndsg.dbf', regions)
    totPop = getTotalPopulation('../results/glbnds_25.nc', '../results/pcount_25.nc', regions, cntMapping)
    log('results:\n')
    for i in range(totPop.shape[0]):
        log('{0}: {1}\n'.format(regions[i][0], totPop[i]))