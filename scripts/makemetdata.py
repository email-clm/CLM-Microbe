#!/usr/bin/env python

import os, sys
from optparse import OptionParser
from scipy.io import netcdf
#from scipy.interpolate import interp1d
import pp
import numpy

def getvar(fname, prod, year, thisx, thisy, nx, ny, nhours, timeres, offset, var):

    f = scipy.io.netcdf.netcdf_file(fname,'r')
   
    #interpoloate values to hourly resolution
    xin  = numpy.cumsum(numpy.ones(nhours/timeres))*timeres + offset  #variables for interpolation
    xout = numpy.cumsum(numpy.ones(nhours))
    print fname
    myvar     = f.variables[var]
    thisvar = myvar[0:(nhours/timeres),thisy,thisx]*myvar.scale_factor+myvar.add_offset
    f_lin = numpy.interp(xout, xin, thisvar)
    f.close()
    return f_lin


parser = OptionParser()
parser.add_option("--site", dest="site", default='MYTEST', \
                  help = '6-character site code to run')
parser.add_option("--latitude", dest="lat", default='25', \
                  help = 'Site Latitude (deg N)')
parser.add_option("--longitude", dest="lon", default='-90', \
                  help = 'Site Longitude (deg E)')
parser.add_option("--height", dest="zbot", default='10', \
                  help = 'Observation height (meters)')
parser.add_option("--startyear", dest="startyear", default='1990', \
                  help = 'Start year')
parser.add_option("--endyear", dest="endyear", default='2010', \
                  help = 'End year')
parser.add_option("--netcdf_path", dest='nc_path', default='', \
                  help = 'netcdf path')
parser.add_option("--reanalysis_path", dest='rean_path', default='/home/zdr/data/Reanalysis', \
                  help = 'Path to reanalysis datasets')
parser.add_option("--ccsm_input", dest='ccsm_input', default='../../ccsm_inputdata', \
                  help = 'Input data directory')
parser.add_option("--metdir", dest='metdir', default='', \
                  help = 'Subdirectory for met data')
parser.add_option("--product", dest='prod', default='NCEP', \
                  help = 'Reanalysis product')
parser.add_option("--numxpts", dest="numx", default=1, \
                  help = "number of points in x")
parser.add_option("--numypts", dest="numy", default=1, \
                  help = "number of points in y")
parser.add_option("--resx", dest="resx", default=0.1, \
                  help = "resolution (deg) in x")
parser.add_option("--resy", dest="resy", default=0.1, \
                  help = "resolution (deg) in y")


(options, args) = parser.parse_args()

#set up parallel
ppservers = ("*",)
job_server = pp.Server(ppservers=ppservers)
print ("Starting pp with", job_server.get_ncpus(), "workers")


print('Creating meteorlogy data for '+options.site+' using '+options.prod)
ccsm_input = os.path.abspath(options.ccsm_input)
print(ccsm_input)
if (options.nc_path == ''):
    print('Using netcdf path from environment:')
    os.system('echo $NETCDF_PATH')
    options.nc_path = os.getenv('NETCDF_PATH')

#os.chdir('../..')
os.system('pgf90 PTCLM_files/writevar.f90 -I$NETCDF_PATH/include -L$NETCDF_PATH/lib -lnetcdf -lnetcdff')

#parameters for saturation vapor pressure equation (for converting qair to RH)
dim = (2,7)
a = numpy.zeros(dim)
a[:,0] = [6.10779961, 6.109177956]
a[:,1] = [4.436518521e-1, 5.034698970e-1]
a[:,2] = [1.428945805e-2, 1.886013408e-2]
a[:,3] = [2.650648471e-4, 4.176223716e-4]
a[:,4] = [3.031240396e-6, 5.824720280e-6]
a[:,5] = [2.034080948e-8, 4.838803174e-8]
a[:,6] = [6.136920929e-11, 1.838826904e-10]


if (options.prod == 'NARR'):
    #NARR-specific information
    nx = 349
    ny = 277
    dirs_in = ['air.2m', 'apcp', 'dlwrf', 'dswrf', 'pres.sfc', 'shum.2m', 'uwnd.10m', 'vwnd.10m']
    vars_in = ['air', 'apcp', 'dlwrf', 'dswrf', 'pres', 'shum', 'uwnd', 'vwnd']
    timeres = 3   
if (options.prod == 'NCEP'):
    dirs_in = ['air.2m.gauss', 'prate.sfc.gauss', 'dlwrf.sfc.gauss', 'dswrf.sfc.gauss', 'pres.sfc.gauss', \
                   'shum.2m.gauss', 'uwnd.10m.gauss', 'vwnd.10m.gauss']
    vars_in = ['air', 'prate', 'dlwrf', 'dswrf', 'pres', 'shum', 'uwnd', 'vwnd']
    timeres = 6
    nx = 192
    ny =  94
if (options.prod == 'Princeton'):
    dirs_in = ['tas', 'prcp', 'dlwrf', 'dswrf', 'pres', 'shum', 'wind']
    vars_in = ['tas', 'prcp', 'dlwrf', 'dswrf', 'pres', 'shum', 'wind']
    timeres = 3
    nx = 360
    ny = 180

offsets = [0, 0, int(-timeres/2), int(-timeres/2), 0, 0, 0, 0]

#CLM variables to output
vars_out = ['TBOT','PRECTmms','FLDS','FSDS','PSRF','RH','WIND']

ndaysm=[31,28,31,30,31,30,31,31,30,31,30,31]

nyears = (int(options.endyear) - int(options.startyear)+ 1 )
nhours = 8760
dim=(8,nhours*nyears)
thisvar_out = numpy.zeros(dim)
dim=(7,12)
monthly_stats = numpy.zeros(dim)
dim=7
annual_stats = numpy.zeros(dim)
startyear = int(options.startyear)
endyear = int(options.endyear)

ynum = 0
jobs=[]


#-------------------------get closest point---------------------------------------------------

if (options.prod == 'NARR' or options.prod == 'NCEP'):
    fname = options.rean_path+'/'+options.prod+'/'+dirs_in[0]+'/'+dirs_in[0]+'.'+str(startyear)+'.nc'
if (options.prod == 'Princeton'):
    fname = options.rean_path+'/'+options.prod+'/'+dirs_in[0]+'/'+dirs_in[0]+'_3hourly_'+str(startyear)+'-'+ \
        str(startyear)+'.nc'

f = netcdf.netcdf_file(fname,'r')
if (options.prod == 'NARR' or options.prod == 'NCEP'):
    lat = f.variables['lat']
    lon = f.variables['lon']
if (options.prod == 'Princeton'):
    lat = f.variables['latitude']
    lon = f.variables['longitude']
dist=999999
if (options.prod == 'NCEP' or options.prod == 'Princeton'):
    if (float(options.lon) < 0):
        options.lon = str(float(options.lon)+360)
for i in range(0,nx):
    for j in range(0,ny):
        if (options.prod == 'NARR'):
            thisdist = ((110*abs(float(options.lat) - lat[j][i]))**2 + \
                            (110*abs(float(options.lon) - lon[j][i])* \
                                 numpy.cos(lat[j][i]*numpy.pi/180.))**2)**(0.5)
        else:
            thisdist = ((110*abs(float(options.lat) - lat[j]))**2 + \
                            (110*abs(float(options.lon) - lon[i])* \
                                 numpy.cos(lat[j]*numpy.pi/180.))**2)**(0.5)
        if (thisdist < dist):
            dist=thisdist
            thisx = i
            thisy = j
print thisx, thisy, lat[thisy], lon[thisx]
f.close()


#----------------get data (paralllelized) ---------------------------------------------------------
for y in range(startyear,endyear+1):
    #print('Reading '+str(y))
    vind=0

    #get reanalysis data, extract requested point and convert to hourly resolution
    for v in dirs_in:
        #get filename
        if (options.prod == 'NARR' or options.prod == 'NCEP'):
            fname = options.rean_path+'/'+options.prod+'/'+v+'/'+v+'.'+str(y)+'.nc'
        if (options.prod == 'Princeton'):
            fname = options.rean_path+'/'+options.prod+'/'+v+'/'+v+'_3hourly_'+str(y)+'-'+str(y)+'.nc'
        offset = 0
        jobs.append(job_server.submit(getvar,(fname, options.prod, y, thisx, thisy, nx, ny, nhours, \
                                            timeres,offsets[vind],vars_in[vind]),(), ("scipy.io.netcdf","numpy",)))
        vind = vind+1

jnum=0
for job in jobs:
    yind = jnum/len(dirs_in)
    vind = jnum % len(dirs_in)
    thisvar_out[vind,(yind*nhours):(yind*nhours+nhours)] = job()
    jnum = jnum+1 
        
#convert shum.2m to RH
RH_temp = thisvar_out[5,:]
for i in range(0,nhours*nyears):
    ice=0
    myt = thisvar_out[0,i]-273.15
    if (myt < 0):
        ice=1
    myesat = a[ice,0] + myt*a[ice,1] + myt**2*a[ice,2] + myt**3*a[ice,3] + myt**4*a[ice,4] + \
        myt**5*a[ice,5]+myt**6*a[ice,6]
    RH_temp[i] = 100.0 * thisvar_out[5,i]*thisvar_out[4,i]/(myesat*100*(0.622+thisvar_out[5,i]))
    thisvar_out[5,i] = max(min(RH_temp[i],100.0),0.0)

    #convert U, V wind to wind speed
    if (options.prod == 'NARR' or options.prod == 'NCEP'):
        thisvar_out[6,i] = (thisvar_out[6,i]**2 + thisvar_out[7,i]**2)**0.5
    #convert precip units to mm s-1
    if (options.prod == 'NARR'):
        thisvar_out[1,i] = thisvar_out[1,i]/(3600*timeres)
        
#create CLM-style monthly files
starti = 0
for y in range(0,nyears):
    print 'Writing '+str(startyear+y)
    for m in range(0,12):
        os.system('mkdir -p '+ccsm_input+'/atm/datm7/CLM1PT_data/'+ \
                      options.metdir+'/1x1pt_'+options.site)
        if (m+1 < 10):
            fname_out = ccsm_input+'/atm/datm7/CLM1PT_data/'+ \
                options.metdir+'/'+str(options.numx)+'x'+str(options.numy)+'pt_'+ \
                options.site+'/'+ str(startyear+y)+'-0'+str(m+1)+'.nc'
            os.system('cp ./PTCLM_files/sample_metdata/1999-0'+str(m+1)+'.nc '+fname_out)
        else:
            fname_out = ccsm_input+'/atm/datm7/CLM1PT_data/'+options.metdir+ \
                '/'+str(options.numx)+'x'+str(options.numy)+'pt_'+ options.site+ \
                '/'+str(startyear+y)+'-'+str(m+1)+'.nc'
            os.system('cp ./PTCLM_files/sample_metdata/1999-'+str(m+1)+'.nc '+fname_out)
        output = open('data.in','w')                #write for fortran input file:
        output.write(str(startyear+y)+' '+str(m+1)+'\n')        #year and month
        output.write(options.lat+ ' '+options.lon+' '+options.zbot+'\n') #lat, lon, ht

        output.write(fname_out+'\n')                #filename
        output.write(str(len(vars_out))+'\n')       #number of vars
        outstr=''
        for i in range(0,len(vars_out)):
            outstr=outstr+str(vars_out[i]+' ')      #variable names
        output.write(outstr+'\n')              
        output.write(str(ndaysm[m]*24)+'\n')        #number of time points
        output.write(str(options.numx)+' '+str(options.numy)+' '+str(options.resx)+ \
                         ' '+str(options.resy)+'\n')     #regional grid info
        for i in range(0,ndaysm[m]*24):
            outstr=''                               #data values
            for j in range(0,len(vars_out)):
                outstr = outstr+str(thisvar_out[j,starti+i])+' '
            output.write(outstr+'\n')
        output.close()
        os.system('./a.out')                        #execute fortran code to modify netcdf file
        #os.system('rm data.in')
        for v in range(0,len(vars_out)):
            monthly_stats[v,m] += numpy.mean(thisvar_out[v,starti:starti+ndaysm[m]*24]/nyears)

        starti = starti + ndaysm[m]*24
    #annual statistics
    for v in range(0,len(vars_out)):
        annual_stats[v] +=  numpy.mean(thisvar_out[v,y*nhours:y*nhours+nhours]/nyears)
    

os.system('rm a.out')
monthly_stats[0,:] -= 273.15   #C
annual_stats[0] -= 273.15
for m in range(0,12):
    monthly_stats[1,m] *=3600*24*ndaysm[m]   #mm/month
annual_stats[1] *= 3600*24*365

output = open(ccsm_input+'/atm/datm7/CLM1PT_data/'+options.metdir+'/'+str(options.numx)+'x' \
                  +str(options.numy)+'pt_'+options.site+'/stats.out','w')
for m in range(0,12):
    outstr=''
    for v in range(0,len(vars_out)):
        outstr += str(monthly_stats[v,m])+' '
    output.write(outstr+'\n')
outstr=''
for v in range(0,len(vars_out)):
    outstr+= str(annual_stats[v])+' '
output.write(outstr+'\n')

output.close()

        

    
        
                
            


