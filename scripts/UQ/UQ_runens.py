#!/usr/bin/env python

import os, sys, csv, time, math, numpy, getpass
from optparse import OptionParser

#Create, run and process a CLM/ALM model ensemble member
#  given specified case and parameters (see parm_list and parm_data files)
#  Parent case must be pre-built and all namelists in run directory.
#  Post-processing calcuates normalized sum of squared errors (SSE) given
#  data constraints specified in "constraints" directory"
#  DMRicciuto 12/1/2015
#
#  Note:  This will only work for single-point CLM/ALM compiled with MPI_SERIAL

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--runroot", dest="runroot", default="", \
                  help="Directory where the run would be created")
parser.add_option("--ens_num", dest="ensnum", default=1, \
                  help="Ensemble member number")
parser.add_option("--parm_list", dest="parm_list", default="", \
                  help="File containing parameter names/pfts to modify")
parser.add_option("--parm_data", dest="parm_data", default="", \
                  help="File containing parameter values and ranges")
parser.add_option("--constraints", dest="constraints", default="", \
                  help="Directory containing constraining variables")
parser.add_option("--norun", dest="norun", default=False, action="store_true", \
                  help="Don't run model (use for testing purposes)")
parser.add_option("--machine", dest="machine", default="cades", \
                  help="My machine")

(options, args) = parser.parse_args()


#================= netcdf manipulation functions ===============================#
def getvar(fname, varname):
    usescipy = False
    try:
    	import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"r",mmap=False)
        var = nffile.variables[varname]
        varvals = var[:].copy()    #works for vector only?
        nffile.close()
    else:    
    	nffile = netcdf.NetCDFFile(fname,"r")
    	var = nffile.variables[varname]
    	varvals = var.getValue()
    	nffile.close()
    return varvals

def putvar(fname, varname, varvals):
    usescipy = False
    try:
        import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"a",mmap=False)
        var = nffile.variables[varname]
        var[:] = varvals
        nffile.close()
    else:
        nffile = netcdf.NetCDFFile(fname,"a")
        var = nffile.variables[varname]
        var.assignValue(varvals)
        nffile.close()
    ierr = 0
    return ierr

#======================================================================
UQdir = os.getcwd()

parm_names=[]
parm_indices=[]
parm_values=[]
myinput = open(options.parm_list, 'r')
lnum=0

username = getpass.getuser()
if (options.machine == 'cades' and options.runroot == ''):
    options.runroot = '/lustre/or-hydra/cades-ccsi/scratch/'+username
elif (options.runroot == ''):
    options.runroot = '../../run'


#get parameter names and PFT information
casenames = []
for s in myinput:
    if (lnum == 0):
        if ('20TR' in s or '1850' in s):
            casenames.append(s.strip())
            isfullrun = False
        else:
            casenames.append(s.strip()+'_I1850CLM45CN_ad_spinup')
            casenames.append(s.strip()+'_I1850CLM45CN')
            casenames.append(s.strip()+'_I20TRCLM45CN')
            #for now, hard-code the number of years for ad_spinup and final spinup
            isfullrun = True
    else:
        pdata = s.split()
        #print pdata
        parm_names.append(pdata[0])
        parm_indices.append(int(pdata[1]))
    lnum = lnum+1
myinput.close()

#get parameter values
myinput = open(options.parm_data, 'r')
for s in myinput:    
    parm_values.append(float(s))
myinput.close()

n_parameters = len(parm_names)
gst=str(100000+int(options.ensnum))

#create ensemble directories from original case(s)
isfirstcase = True
workdir = os.path.abspath('.')
for casename in casenames: 
    orig_dir = str(os.path.abspath(options.runroot)+'/'+casename+'/run')
    ens_dir  = os.path.abspath(options.runroot)+'/UQ/'+casename+'/g'+gst[1:]
    print 'test1'
    os.system('mkdir -p '+options.runroot+'/UQ/'+casename+'/g'+gst[1:]+'/timing/checkpoints')
    os.system('cp '+orig_dir+'/*_in* '+ens_dir)
    os.system('cp '+orig_dir+'/*nml '+ens_dir)
    os.system('cp '+orig_dir+'/*stream* '+ens_dir)
    #os.system('cp '+orig_dir+'/*.r.*.nc '+ens_dir)
    os.system('cp '+orig_dir+'/*.rc '+ens_dir)
    os.system('cp '+orig_dir+'/*para*.nc '+ens_dir)
    os.system('cp '+orig_dir+'/*initial* '+ens_dir)
    os.system('cp '+orig_dir+'/*pftdyn* '+ens_dir)
    print casename, 'test2'
    os.system('cp '+ens_dir+'/microbepar_in '+ens_dir+'/microbepar_in_orig')
    username = getpass.getuser()
    os.system('cp /home/'+username+'/models/CLM_SPRUCE/inputdata/lnd/clm2/paramdata/'+ \
               'clm_params_spruce_calveg.nc '+ens_dir+'/clm_params_'+ \
               gst[1:]+".nc")	
    if (isfullrun):     #Full spinup simulation
        mydrv_in = open(orig_dir+'/drv_in','r')
        for s in mydrv_in:
            if ('stop_n' in s and 'ad_spinup' in casename):
                nyears_ad_spinup = int(s.split()[2])
            elif ('stop_n' in s and '1850' in casename and \
                  'ad_spinup' not in casename):
                nyears_final_spinup =  int(s.split()[2])
        mydrv_in.close()

        inifile=''
        if ('1850' in casename and 'ad_spinup' not in casename):
            yst = str(10000+nyears_ad_spinup+1)
            inifile = ens_dir_last+'/'+casename_last+'.clm2.r.'+yst[1:]+'-01-01-00000.nc'
        if ('20TR' in casename):
            yst = str(10000+nyears_final_spinup+1)
            inifile = ens_dir_last+'/'+casename_last+'.clm2.r.'+yst[1:]+'-01-01-00000.nc'
    else:                #Trasient case only
        inifile = ens_dir+'/'+casename+'.clm2.r.1974-01-01-00000.nc'
    casename_last = casename
    ens_dir_last = ens_dir


#inifile =  ens_dir+'/SPRUCE-finalspinup-peatland-carbon-initial.nc'
#cwdc=getvar(inifile, 'cwdc_vr')
#ierr = putvar(inifile, 'cwdc_vr', cwdc*0.0)
#cwdn=getvar(inifile, 'cwdn_vr')
#ierr = putvar(inifile, 'cwdn_vr', cwdn*0.0)
#lit1 = getvar(inifile, 'litr1c_vr')
#ierr = putvar(inifile, 'litr1c_vr', lit1*0.0)
#lit1 = getvar(inifile, 'litr1n_vr')
#ierr = putvar(inifile, 'litr1n_vr', lit1*0.0)
#lit2 = getvar(inifile, 'litr2c_vr')
#ierr = putvar(inifile, 'litr2c_vr', lit2*0.0)
#lit2 = getvar(inifile, 'litr2n_vr')
#ierr = putvar(inifile, 'litr2n_vr', lit2*0.0)
#lit3 = getvar(inifile, 'litr3c_vr')
#ierr = putvar(inifile, 'litr3c_vr', lit3*0.0)
#lit3 = getvar(inifile, 'litr3n_vr')
#ierr = putvar(inifile, 'litr3n_vr', lit3*0.0)

    #loop through all filenames, change directories in namelists, change parameter values
    for f in os.listdir(ens_dir):
        if (os.path.isfile(ens_dir+'/'+f) and (f[-2:] == 'in' or f[-3:] == 'nml' or 'streams' in f)):
            myinput=open(ens_dir+'/'+f)
            myoutput=open(ens_dir+'/'+f+'.tmp','w')
            for s in myinput:
                if ('finidat = ' in s):
                    myoutput.write(" finidat = '"+inifile+"'\n")
                elif ('fpftcon' in s):
                    est = str(100000+int(options.ensnum))
                    myoutput.write(" fpftcon = './clm_params_"+est[1:]+".nc'\n")
                    #Hard-coded parameter file
                    pftfile = ens_dir+'/clm_params_'+est[1:]+'.nc'
                    microbefile = ens_dir+'/microbepar_in'
                    pnum = 0
                    for p in parm_names:
                        if (p[0:2] == 'm_' or p[0:2] == 'k_'):   #Assume this is a microbe_par parameter
                            moutput = open(microbefile,'w')
                            minput = open(microbefile+'_orig','r')
                            for s2 in minput:
                                if (p.lower() in s2.lower()):
                                    moutput.write(p+'   '+str(parm_values[pnum]) \
                                                  +'\n')
                                else:
                                    moutput.write(s2)
                            minput.close()
                            moutput.close()
                            os.system('cp '+microbefile+' '+microbefile+'_orig')
                        else:
                            if (pnum == 0):
                                stem_leaf = getvar(pftfile, 'stem_leaf')
                                stem_leaf[2:5]=-1
                                ierr = putvar(pftfile, 'stem_leaf', stem_leaf)
                            param = getvar(pftfile, p)
                            if (parm_indices[pnum] > 0):
                                param[parm_indices[pnum]-1] = parm_values[pnum]
                            elif (parm_indices[pnum] == 0):
                                param = parm_values[pnum]
                            else:
                                param[:] = parm_values[pnum]
                        ierr = putvar(pftfile, p, param)
                        pnum = pnum+1
                #elif ('logfile =' in s):
                #    myoutput.write(s.replace('`date +%y%m%d-%H%M%S`',timestr))
                else:
                    myoutput.write(s.replace(orig_dir,ens_dir))
            myoutput.close()
            myinput.close()
            os.system(' mv '+ens_dir+'/'+f+'.tmp '+ens_dir+'/'+f)

    os.chdir(ens_dir)
    if (isfirstcase):
        exedir = os.path.abspath(orig_dir+'/../bld/')
    if (options.norun == False):
        os.system(exedir+'/cesm.exe > ccsm_log.txt')
    isfirstcase=False

#---------  code to post-process ensebmle member and cacluate total normalized SSE ----------
sse=0
myoutput = open('myoutput_sse.txt','w')
myind = 0
for p in parm_names:
        myoutput.write(str(parm_names[myind])+' '+str(parm_indices[myind])+' '+str(parm_values[myind])+'\n')
	myind = myind+1

for filename in os.listdir(UQdir+'/'+options.constraints):
  if (not os.path.isdir(filename)):
    myinput = open(UQdir+'/'+options.constraints+'/'+filename,'r')
    myvarname = filename.split('.')[0]  #model variable is filename
    #code to deal with special variables and/or aggregation
    #-------------
    lnum = 0
    year = 0
    for s in myinput:
        if (lnum == 0):
            header = s.split()
        else:
            hnum = 0
            PFT=-1      #default:  don't use PFT-specific info 
                        #  if specified, use h1 file (PFT-specific)
            doy=-1      #default:  annual average
            month=-1    #default:  don't use monthly data
            depth=-1
            unc = -999
            for h in header:
                if (h.lower() == 'year'):
                    year_last = year
                    year = int(s.split()[hnum])
                if (h.lower() == 'doy'):
                    doy = int(s.split()[hnum])
                if (h.lower() == 'month'):
                    month = int(s.split()[hnum])
                if (h.lower() == 'pft'):
                    PFT = int(s.split()[hnum])
                if (h.lower() == 'value'):
                    value = float(s.split()[hnum])
                if (h.lower() == 'depth'):
                    depth = float(s.split()[hnum])
                if ('unc' in h.lower()):
                    unc   = float(s.split()[hnum])
                hnum = hnum+1
            #get the relevant variable/dataset
            #Assumes annual file with daily output
            if (year != year_last and year <= 2013):
                if (PFT == -1):
                    myfile = casename+'.clm2.h0.'+str(year)+'-01-01-00000.nc'
                else:
                    myfile = casename+'.clm2.h1.'+str(year)+'-01-01-00000.nc'
                #post processing of model output with nco to match constraining variables 
                if (myvarname == 'STEMC'):
                    os.system('ncap -s "STEMC=DEADSTEMC+LIVESTEMC" '+myfile+' '+myfile+'.tmp')
                    os.system('mv '+myfile+'.tmp '+myfile)
                if (myvarname == 'WTHT'):
                    os.system('ncap -s "WTHT=(ZWT*-1+0.3)*1000 " '+myfile+' '+myfile+'.tmp')
                    os.system('mv '+myfile+'.tmp '+myfile)         
                if (myvarname == 'AGBIOMASS'):
                    os.system('ncap -s "AGBIOMASS=DEADSTEMC+LIVESTEMC+LEAFC" '+myfile+' '+myfile+'.tmp')
                    os.system('mv '+myfile+'.tmp '+myfile)

                myvals = getvar(myfile, myvarname)
            if (doy > 0 and value > -900):
                if (myvarname == 'WTHT'):
                    unc = 30.0   #no uncertainty given for water table height.
                if (PFT > 0):
                    #PFT-specific constraints (daily)
                    if (myvarname == 'AGBIOMASS' and PFT == 3):
                        #Both tree types - Larch and spruce, hummock+hollow
                        model_val = (myvals[doy,PFT-1]*0.25+myvals[doy,PFT]*0.25)*0.75 \
                                   +(myvals[doy,PFT+16]*0.25+myvals[doy,PFT+17]*0.25)*0.25
                    else:
                        #use hummock+hollow
                        model_val = myvals[doy,PFT-1]*0.25*0.75 + myvals[doy,PFT+16]*0.25*0.25
                    if (unc < 0):
    	              unc = value*0.25 #default uncertainty set to 25%
                    sse = sse + ((model_val-value) /unc)**2
                    myoutput.write(str(myvarname)+' '+str(year)+' '+str(doy)+' '+str(PFT)+' '+ \
                                       str(model_val)+' '+str(value)+' '+str(unc)+' '+str(sse)+'\n')
                elif (depth > 0):
                    #depth-specific constraint (relative to hollow)
                    layers = [0,1.8,4.5,9.1,16.6,28.9,49.3,82.9,138.3,229.6,343.3]
                    for l in range(0,10):
                        if (depth >= layers[l] and depth < layers[l+1]):
                            thislayer = l
                            model_val = myvals[doy,thislayer,1]   #Holow 
                            sse = sse + ((model_val-value) / unc )**2        
                            myoutput.write(str(myvarname)+' '+str(year)+' '+str(doy)+' '+str(depth)+' '+ \
                                               str(model_val)+' '+str(value)+' '+str(unc)+' '+str(sse)+'\n')
                else:
                    #Column-level constraint (daily)
                    #Water table, column-level (no PFT info), use hummock only
                    model_val = myvals[doy,0]
                    sse = sse + ((model_val-value) / unc )**2        
                    myoutput.write(str(myvarname)+' '+str(year)+' '+str(doy)+' '+str(PFT)+' '+ \
                                       str(model_val)+' '+str(value)+' '+str(unc)+' '+str(sse)+'\n')
            elif (value > -900):
                #Annual constraints.  Only one is AGNPP, moss only.  Need to convert units
                model_val = sum(myvals[0:365,PFT-1]*0.25*0.75)*24*3600 + \
                            sum(myvals[0:365,PFT+16]*0.25*0.25)*24*3600
                sse = sse + ((model_val-value) / unc )**2      
                myoutput.write(myvarname+' '+str(year)+' '+str(doy)+' '+str(PFT)+' '+ \
                                   str(model_val)+' '+str(value)+' '+str(unc)+' '+str(sse)+'\n')

        lnum = lnum+1
myoutput.close()

myoutput = open(workdir+'/ssedata/mysse_'+gst[1:]+'.txt','w')
myoutput.write(str(sse))
myoutput.close()



