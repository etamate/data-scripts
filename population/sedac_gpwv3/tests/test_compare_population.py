# -*- coding: utf-8 -*- 
'''
Testing the SEDAC GPW v3 dataset translations

Code is written and tested under Python 2.7 

@license: This work is licensed under a Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
@author: Mate Rozsai
'''
import csv
import dbf
import netCDF4
import numpy
import sys
import collections
import datetime

Country = collections.namedtuple('Country', 'name iso3v10 unsdcode sedaccode')
UNData = collections.namedtuple('UNData', 'name unsdcode population')

def log(message):
    sys.stdout.write(message)

def loadCountries(dbfFile):
    log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   loading country dbf...')
    table = dbf.Table(dbfFile)
    table.open()
    maxid=0
    countries = dict()
    for row in table:
        countries[row.iso3v10]=Country(name=row.countryeng.strip().replace('\t', ' '), iso3v10=row.iso3v10, unsdcode=row.unsdcode, sedaccode=row.value)
        maxid=max(maxid,row.value)
    table.close()
    log('.............\n')
    return countries, maxid

def loadUNPop(tsvFile):
    log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   loading UN population data...')
    countries=dict()
    with open(tsvFile, 'rb') as fcsv:
        reader = csv.reader(fcsv, delimiter='\t', quoting=csv.QUOTE_NONE)
        rownum=0
        for row in reader:
            rownum+=1
            #check header rows 
            if 1==rownum:
                if 'Country code'!=row[1] or '2000'!=row[2]: raise Exception('IO error', 'Invalid header in UN population file!.')
            else:
                countries[int(row[1])] = UNData(name=row[0], unsdcode=int(row[1]), population=float(row[2])*1000)
    log('......\n')
    return countries
    
def getTotalPopulation(glbndsFile, pcountFile, ftype, maxid):
    def loadNc(ncFile, var):
        nc = netCDF4.Dataset(ncFile, 'r')
        ncVar = nc.variables[var][:]
        lons = nc.variables['lon'][:]
        lats = nc.variables['lat'][:]
        return ncVar, lats, lons
        nc.close()
    
    log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   loading '+ftype+' nc files...')
    glbnds, latsG, lonsG = loadNc(glbndsFile, 'glbnds')
    pcount, latsP, lonsP = loadNc(pcountFile, 'pcount')
    # check if we use the same grid resolution and ordering
    if not(numpy.array_equal(latsG, latsP)): raise Exception('Different lat coordinates!')
    if not(numpy.array_equal(lonsG, lonsP)): raise Exception('Different lon coordinates!')
    # calculate the total population
    totPop = numpy.zeros([maxid+1])
    for idxLat in range(latsG.shape[0]):
        for idxLon in range(lonsG.shape[0]):
            if not glbnds.mask[idxLat,idxLon]:
                totPop[glbnds[idxLat,idxLon]] += pcount[idxLat, idxLon]
        if 0 == idxLat % int(latsG.shape[0]/10): log('.')
    log('.\n')
    return totPop

if __name__ == "__main__":
    # load country data, UN population data and calculate total populations from the two resolutions 
    countries, maxid = loadCountries('../gl_gpwv3_ntlbndid_ascii_25/bndsg.dbf')
    unPopData = loadUNPop('../un_population/un_pop_2000.tsv')
    totPopHalf = getTotalPopulation('../results/glbnds_half.nc', '../results/pcount_half.nc', '1/2°', maxid)
    totPop25 = getTotalPopulation('../results/glbnds_25.nc', '../results/pcount_25.nc', "2.5'", maxid)
    # combine the country list with the calculated total populations and the UN dataset
    log('['+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+']   creating comparison table...')
    resTable=dict()
    for country in countries.itervalues():
        unData=None
        if country.unsdcode in unPopData: unData = unPopData[country.unsdcode]
        resTable[country.iso3v10] = [country, unData, totPopHalf[country.sedaccode], totPop25[country.sedaccode]]
    # writing a comparison table to a tsv file (by population in decreasing order)
    resultFile = 'results/population_comparison.tsv'
    with open(resultFile, 'w') as f:
        f.write("iso3v10\tname\tunsdcode\tUN name\tUN population\tSEDAC 1/2° population\tSEDAC 2.5' population\n")
        for code in sorted(resTable.keys(), key=lambda x: resTable[x][-1], reverse=True):
            country=resTable[code]
            f.write('{0}\t{1}\t{2}\t'.format(country[0].iso3v10, country[0].name, country[0].unsdcode))
            if None==country[1]: f.write('<missing>\t\t')
            else: f.write('{0}\t{1}\t'.format(country[1].name, country[1].population))
            f.write('{0}\t{1}\n'.format(country[2], country[3]))
    log('.......\n                           check results in ' + resultFile)

