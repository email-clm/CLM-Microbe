#!/usr/bin/python

import os, csv, time, math
from Scientific.IO import NetCDF
from Scientific import MPI
import numpy

site_code  = 'US-MOz'  #AmeriFlux/FLUXNET identifier


#Function for cacluating saturation specific humidity
def qsat(t,pres):
    a = [6.107799961, 4.436518521e-01, 1.428945805e-02, 2.650648471e-04, \
             3.031240396e-06, 2.034080948e-08, 6.136820929e-11]
    b = [6.109177956, 5.034698970e-01, 1.886013408e-02, 4.176223716e-04, \
             5.824720280e-06, 4.838803174e-08, 1.838826904e-10]

    if (t > 150):
        t = t-273.15  #convert K to C
    if (t > 0):
        esat = 100.*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*(a[4]+t*(a[5]+t*a[6]))))))
    else:
        esat = 100.*(b[0]+t*(b[1]+t*(b[2]+t*(b[3]+t*(b[4]+t*(b[5]+t*b[6]))))))
    myqsat = 0.622 * esat / (pres - 0.378*esat)
    return myqsat


#-----------------get relevant tower information -----------------------
tower_data = open('./Towerinfo.txt','r')
for line in tower_data:
    line_data = line.split()
    if (line_data[0] == site_code):
        lst = int(line_data[6])
        site_lat = float(line_data[7])
        site_lon = float(line_data[8])
        print 'Site code is '+site_code
        print 'Site latitude is  '+str(site_lat)
        print 'Site longitude is '+str(site_lon)
        print 'Local time offset is UTC'+str(lst)


#-----------------Get data from CDIAC (Level 2, with gaps, netcdf)  ---------------
inputdata = '/home/zdr/models/ccsm_inputdata/atm/datm7/CLM1PT_data'
redownload = False

getdata = True

if (os.access(inputdata, os.W_OK)):
    AmeriFluxdir = inputdata+'/AmeriFlux/'+site_code
else:
    print 'No write access to input data directory.  Putting data in ./temp'
    os.system('rm -rf ./temp')
    inputdata = './temp'
    AmeriFluxdir = './temp'

if (os.path.exists(AmeriFluxdir)):
    print 'AmeriFlux data already downloaded.'
    if (redownload):
        print '--redownload set.  Downloading data again.'
        os.system('rm -rf '+AmeriFluxdir)
    else:
        getdata = False
else:
    'Creating directory '+AmeriFluxdir
os.system('mkdir -p '+AmeriFluxdir)

if (getdata):
    print 'Getting data from CDIAC'
    workdir = os.path.abspath('.')
    os.chdir(AmeriFluxdir)
    if ('UMB' in site_code):
        getdir = 'http://cdiac.ornl.gov/ftp/ameriflux/data/Level2/Sites_ByID/'+site_code+ \
            '/with_gaps/hourly/'
    else:
        getdir = 'http://cdiac.ornl.gov/ftp/ameriflux/data/Level2/Sites_ByID/'+site_code+ \
            '/with_gaps/'
    os.system('wget -q --mirror --continue -nd '+getdir)
    os.chdir(workdir)

startyear = 0
endyear   = 0

#figure out which years and which version number we are using
for y in range(1990,2015):
    for v in range(0,20):
        vst = str(1000+v)
        fname = AmeriFluxdir+'/AMF_'+site_code[0:2]+site_code[3:6]+'_'+str(y)+'_L2_WG_V'+vst[1:4]+'.nc'
        if (os.path.isfile(fname)):
            if (startyear == 0):
                startyear = y
            endyear = y
            thisversion = v

print "loading "+site_code+' for years '+str(startyear)+' to '+str(endyear)
print "AMF file version "+str(thisversion)            
        
#---------------load the data -----------------------------------------

vars_in = ['DOY', 'HRMIN', 'TA', 'Rg', 'PAR', 'RH', 'PREC', 'WS', 'Rgl', 'PRESS', 'CO2height']
vars_out = ['TBOT','RH','WIND','FSDS','FLDS','PSRF','PRECTmms','ZBOT']
#vars_ncep = ['air', 'dlwrf', 'dswrf', 'prate', 'pres', 'shum', 'uwnd', 'vwnd']
vars_ncep = ['air', 'shum', 'uwnd', 'dswrf', 'dlwrf', 'pres', 'prate', 'vwnd']
#levs_ncep = ['2m', 'sfc', 'sfc', 'sfc', 'sfc', '2m', '10m', '10m']
levs_ncep = ['2m', '2m', '10m', 'sfc', 'sfc', 'sfc', 'sfc', '10m']
vars_out_i = [0, 3, 5, 2, 6, 7, 4, 8]  #indices mapping output to input vars
vars_out_units = ['K','kg/kg','m/s','W/m2','W/m2','Pa','kg/m2/s','m']
long_names = []
long_names.append('temperature at the lowest atm level (TBOT)')
long_names.append('humidity at the lowest atm level (RH)')
long_names.append('wind at the lowest atm level (WIND)')
long_names.append('incident solar (FSDS)')
long_names.append('incident longwave (FLDS)')
long_names.append('pressure at the lowest atm level (PSRF)')
long_names.append('precipitation (PRECTmms)')
long_names.append('observational height')

diurmean     = numpy.zeros((9,366,24),dtype=numpy.float)
diurmean_ct  = numpy.zeros((9,366,24),dtype=numpy.int)
alldata      = numpy.zeros((9,1000000), dtype=numpy.float)-9999
alldata_ncep = numpy.zeros((8,1000000), dtype=numpy.float)-0000

for y in range(startyear,endyear+1):
    vst = str(1000+thisversion)
    fname =  AmeriFluxdir+'/AMF_'+site_code[0:2]+site_code[3:6]+'_'+str(y)+'_L2_WG_V'+vst[1:4]+'.nc'
    data_in = NetCDF.NetCDFFile(fname,'r')
    vnum=0

    npd=24
    for v in vars_in:
        vst = str(1000+thisversion)
        myvar = data_in.variables[v]
        myvar_vals = myvar.getValue()
        if (v == 'DOY'):
            mydoy   = myvar_vals
        elif (v == 'HRMIN'): 
            myhrmin = myvar_vals
            if ((str(myhrmin[1])[-2:] == '30') or (str(myhrmin[0])[-2:] == '30')):
                npd=48
            if (y == startyear):
                print 'Number of timesteps per day: '+str(npd)
            starti = (y-startyear)*365*npd+int(mydoy[0]-1)*npd+int(myhrmin[0])/100 + \
                (str(myhrmin[0])[-2:] == '30')
        ngood=0
        ntot=0
        if ((y % 4) == 0):
            mvar_vals = myvar_vals[:-1*npd]
        if (vnum >= 2):
            alldata[vnum-2][starti:starti+len(myvar_vals)] = myvar_vals
 
        for thisval in myvar_vals:
            if (thisval > -900 and vnum >=2):
                ngood = ngood+1
                #print mydoy[ntot]-1, vnum-2
                diurmean[vnum-2][int(mydoy[ntot])-1][int(myhrmin[ntot])/100] = \
                    diurmean[vnum-2][int(mydoy[ntot])-1][int(myhrmin[ntot])/100] + float(thisval)
                diurmean_ct[vnum-2][int(mydoy[ntot])-1][int(myhrmin[ntot])/100] = \
                    diurmean_ct[vnum-2][int(mydoy[ntot])-1][int(myhrmin[ntot])/100] + 1
            ntot = ntot+1
        vnum = vnum+1      

for v in range(0,9):
    for d in range(0,365):
        for h in range(0,24):
            if (diurmean_ct[v][d][h] > 0):
                diurmean[v][d][h] = diurmean[v][d][h]/diurmean_ct[v][d][h]
            else:
                diurmean[v][d][h] = -9999


npoints = (endyear-startyear+1)*365*npd

ctbad = numpy.zeros(9, numpy.int)
for i in range(0,npoints):
    for v in range(0,9):
        myday = (i % (365*npd))/npd
        myhour = ((i % (365*npd)) % npd) / (npd/24)
        if (alldata[v][i] < -900):
            alldata[v][i] = diurmean[v][myday][myhour]
            if (diurmean[v][myday][myhour] < -900):
                ctbad[v]=ctbad[v]+1

#-------------get reanalysis data (NCEP for now) ----------------
reanal_dir = inputdata+'/reanalysis/'+site_code
os.system('mkdir -p '+reanal_dir)
#startyear=1948
for v in range(0,8):
    reanal_file = site_code+'_'+str(startyear)+'-'+str(endyear)+'_'+vars_ncep[v]+'.nc'
    file_list = ''
    for y in range(startyear,endyear+1): #startyear,
        if (y == startyear):
            fname_first = '/home/zdr/models/reanalysis/NCEP/'+vars_ncep[v]+'.'+levs_ncep[v]+'/' \
                +vars_ncep[v]+'.'+levs_ncep[v]+'.gauss.'+str(y)+'.nc'
        file_list = file_list+ ' /home/zdr/models/reanalysis/NCEP/'+vars_ncep[v]+'.'+levs_ncep[v]+'/' \
            +vars_ncep[v]+'.'+levs_ncep[v]+'.gauss.'+str(y)+'.nc'
    if (v == 0):
        data_in = NetCDF.NetCDFFile(fname_first,'r')
        myvar = data_in.variables['lat']
        ncep_lat = myvar.getValue()
        myvar = data_in.variables['lon']
        ncep_lon = myvar.getValue()
        data_in.close()
        mindist = 99999.
        i=0
        for mylon in ncep_lon:
            j=0
            for mylat in ncep_lat:
                if (mylon > 180):
                    mylon = mylon -360.
                thisdist = (mylon - site_lon)**2 + math.cos(mylat*math.pi/180.) \
                    *(mylat - site_lat)**2
                if (thisdist < mindist):
                    xind = i
                    yind = j
                    mindist = thisdist
                j=j+1
            i=i+1
        print ncep_lon[xind], ncep_lat[yind], site_lon, site_lat
    if (not os.path.isfile(reanal_dir+'/'+reanal_file)):
        print ('Creating reanalysis file for '+vars_ncep[v])
        os.system('/home/zdr/nco-4.0.0/gcc/bin/ncrcat -d lat,'+str(yind)+','+str(yind)+' -d lon,' + \
                      str(xind)+','+str(xind)+' -v '+vars_ncep[v]+' '+file_list+' -o '+ \
                      reanal_dir+'/'+reanal_file)
    else:
        print ('Reanalysis file already exists for '+vars_ncep[v])

    data_in = NetCDF.NetCDFFile(reanal_dir+'/'+reanal_file ,'r')
    myvar = data_in.variables[vars_ncep[v]]
    myvar_vals_ncep = myvar.getValue()
    data_in.close()
    starti=0
    starti_ly=0
    for y in range(startyear,endyear+1):
        for h in range(0,1460):
            alldata_ncep[v,starti+h] = myvar_vals_ncep[starti_ly+h,0,0]
        starti = starti+1460
        starti_ly = starti_ly+1460
        if (y % 4 == 0):
            starti_ly = starti_ly+4

for v in range(0,9):
  print 'Number of bad points for '+vars_in[v+2]+': '+str(ctbad[v])

# ------------create CLM-style output variables -----------------

toffset = lst*(npd/24)-1  # conversion from LST to UTC, CLM units
ndaysm = [31,28,31,30,31,30,31,31,30,31,30,31]  #number of days per month
#output meteorogical data to netcdf file(s)

#unit conversions
vnum=0
for v in vars_out:
    if (vnum == 0):     #convert temperature from C to K
        alldata[vars_out_i[vnum],0:npoints] = \
        alldata[vars_out_i[vnum],0:npoints] + 273.15
    if (vnum == 1):
        #NCEP specific humidity to relative humidity
        for i in range(0,npoints):
            alldata_ncep[1,i] = (alldata_ncep[1,i]/ \
                qsat(alldata_ncep[0,i], alldata_ncep[5,i]))*100.
            alldata_ncep[1,i] = max( min(alldata_ncep[1,i], 100.), 0.)
    elif (vnum == 3 and vars_out_i[3] == 2):  #convert PAR to W/m2
        alldata[vars_out_i[vnum],0:npoints] = \
        alldata[vars_out_i[vnum],0:npoints] / 0.48 / 4.6 
    elif (vnum == 5):   #convert pressure from kpa to pa
        alldata[vars_out_i[vnum],0:npoints] = \
        alldata[vars_out_i[vnum],0:npoints] * 1000.0
    elif (vnum == 6):   #convert precip from mm to kg/m2/s
        for i in range(0,npoints):
            if (alldata[vars_out_i[vnum],i] > -9000):
                    alldata[vars_out_i[vnum],i] = \
                    alldata[vars_out_i[vnum],i] / (3600.0 * 24.0/npd)
    vnum = vnum+1

for ftype in range(0,2):
    if (ftype == 0):     #standard monthly files
        nfiles = (endyear-startyear+1)*12
    else:                #CPL_BYPASS style files 
        nfiles = 1
    starti = 0
    outdirname = inputdata+'/1x1pt_'+site_code
    os.system('mkdir -p '+outdirname)
    for n in range(0,nfiles):
        if (ftype == 0):
            yst = str(startyear+n/12)
            mst = str((n % 12)+101)
            out_nc = NetCDF.NetCDFFile(outdirname+'/'+yst+'-'+mst[1:]+'.nc','w')
            out_nc.createDimension('lon', 1)
            out_nc.createDimension('lat', 1)
            npoints_out =  ndaysm[n % 12]*npd
            timedim = 'time'
        else:
            yst = str(startyear)
            mst = str(101)
            out_nc = NetCDF.NetCDFFile(outdirname+'/all_hourly.nc','w')
            out_nc.createDimension('gridcell',1)
            npoints_out = npoints
            timedim = 'DTIME'
        out_nc.createDimension('scalar', 1)
        out_nc.createDimension(timedim, npoints_out)

        #create variables and attributes
        time = out_nc.createVariable(timedim,'d',(timedim,))
        time.units = 'days since '+yst+'-'+mst[1:]+'-01 00:00:00'
        time.long_name = 'Time axis'
        time.calendar = 'noleap'
        time[:] = numpy.arange(npoints_out)*(1.0/npd)
        if (ftype == 0):
            longxy = out_nc.createVariable('LONGXY','d',('lon',))
            latixy = out_nc.createVariable('LATIXY','d',('lat',))
        else:
            longxy = out_nc.createVariable('LONGXY','d',('gridcell',))
            latixy = out_nc.createVariable('LATIXY','d',('gridcell',))
        longxy.long_name = 'longitude'
        longxy.units = 'degrees E'
        longxy[:] = site_lon
        latixy.long_name = 'latitude'
        latixy.units = 'degrees N'
        latixy[:] = site_lat
        edgew = out_nc.createVariable('EDGEW','d',('scalar',))
        edgew.long_name = 'western edge in atmospheric data'
        edgew.units = 'degrees E'
        edgew[:] = site_lon-0.05
        edgee = out_nc.createVariable('EDGEE','d',('scalar',))
        edgee.long_name = 'eastern edge in atmospheric data'
        edgee.units = 'degrees E'
        edgee[:] = site_lon+0.05
        edgen = out_nc.createVariable('EDGEN','d',('scalar',))
        edgen.long_name = 'northern edge in atmospheric data'
        edgen.units = 'degrees N'
        edgen [:]= site_lat+0.05
        edges = out_nc.createVariable('EDGES','d',('scalar',))
        edges.long_name = 'southern edge in atmospheric data'
        edges.units = 'degrees N'
        edges[:] = site_lat-0.05
        start_year = out_nc.createVariable('start_year','i',('scalar',))
        start_year[:] = startyear
        end_year = out_nc.createVariable('end_year','i',('scalar',))
        end_year[:] = endyear

        vnum = 0
        for v in vars_out:
            if (ftype == 1):
                thisvar = out_nc.createVariable(v,'s',('gridcell',timedim,))
                #figure out scale factors and add_offsets for compressed data
                if (v == 'PRECTmms'):
                    if (numpy.max(alldata[vars_out_i[vnum], 0:npoints]) > -9000):
                        scale_factor = (numpy.max(alldata[vars_out_i[vnum], \
                                      0:npoints]))*1.1 / 2.0**15
                    else:
                        scale_factor = (numpy.max(alldata_ncep[vnum, \
                                      0:npoints/(npd/4)]))*1.1 / 2.0**15
                    add_offset = scale_factor*2.0**14
                else:
                    if (numpy.max(alldata[vars_out_i[vnum], 0:npoints]) > -9000):
                        add_offset = numpy.mean(alldata[vars_out_i[vnum], \
                                                        0:npoints])
                        scale_factor = (numpy.max(alldata[vars_out_i[vnum], \
                          0:npoints]) - numpy.min(alldata[vars_out_i[vnum], \
                          0:npoints])) / 2.0**15
                    else:
                        add_offset = numpy.mean(alldata_ncep[vnum,0:npoints/(npd/4)])
                        scale_factor = (numpy.max(alldata_ncep[vnum,0:npoints/(npd/4)]) \
                                   - numpy.min(alldata_ncep[vnum,0:npoints/(npd/4)])) / 2.0**15
                scale_factor = max(1e-9, scale_factor)
                thisvar.add_offset = add_offset
                thisvar.scale_factor = scale_factor
            else:
                thisvar = out_nc.createVariable(v,'f',(timedim,'lat','lon',))

            thisvar.long_name = long_names[vnum]
            thisvar.units = vars_out_units[vnum]
            thisvar.mode = 'time-dependent'

            for np in range(starti,starti+npoints_out):
                if (alldata[vars_out_i[vnum], max(min(np+toffset,npoints-1),0)] < -900):
                    #replace still-missing data with NCEP, use linear interpolation
                    wt1 = 1 - float(np % (npd/4))/(npd/4)
                    wt2 = float(np % (npd/4))/(npd/4)
                    myval = alldata_ncep[vnum, int(np/(npd/4))]*wt1 + \
                         alldata_ncep[vnum, min(int(np/(npd/4))+1, npoints/(npd/4)-1)]*wt2
                else:
                    myval = alldata[vars_out_i[vnum], max(min(np+toffset,npoints-1),0)]
                if (ftype == 1):      #Convert to short int format
                    thisvar[0,np-starti] = int((myval - add_offset) / scale_factor)
                else:
                    thisvar[np-starti,0,0] = myval
            vnum = vnum+1
        out_nc.close()
        starti = starti+ndaysm[n % 12]*npd
  
