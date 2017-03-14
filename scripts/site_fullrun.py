#!/usr/bin/env python

import os, sys, csv
from optparse import OptionParser

parser = OptionParser();

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup", default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--machine", dest="machine", default = 'generic_linux_pgi', \
                  help = "machine to use (default = generic_linux_pgi)")
parser.add_option("--compiler", dest="compiler", default = 'gnu', \
                  help = "compiler to use (pgi*, gnu)")
parser.add_option("--mpilib", dest="mpilib", default="mpi-serial", \
                      help = "mpi library (openmpi*, mpich, ibm, mpi-serial)")
parser.add_option("--csmdir", dest="csmdir", default='..', \
                  help = "base CESM directory (default = ../)")
parser.add_option("--runroot", dest="runroot", default="../run", \
                  help="Directory where the run would be created")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='../inputdata', \
                  help = "input data directory for CESM (required)")
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--nyears_final_spinup", dest="nyears_final_spinup", default='1000', \
                  help="base no. of years for final spinup")
parser.add_option("--clean_build", action="store_true", default=False, \
                  help="Perform a clean build")
parser.add_option("--hist_mfilt", dest="hist_mfilt", default="365", \
                  help = 'number of output timesteps per file (transient only)')
parser.add_option("--hist_nhtfrq", dest="hist_nhtfrq", default="-24", \
                  help = 'output file timestep (transient only)')
parser.add_option("--parm_file", dest="parm_file", default="", \
                  help = 'parameter file to use')
parser.add_option("--nofire", dest="nofire", default=False, \
                  action="store_true", help='Turn off fire algorithms')
parser.add_option("--regional", action="store_true", \
                   dest="regional", default=False, \
                   help="Flag for regional run (2x2 or greater)")
parser.add_option("--np", dest="np", default=1, \
                  help = 'number of processors')
parser.add_option("--tstep", dest="tstep", default=0.5, \
                  help = 'CLM timestep (hours)')
parser.add_option("--co2_file", dest="co2_file", default="fco2_datm_1765-2007_c100614.nc", \
                  help = 'CLM timestep (hours)')
parser.add_option("--nyears_ad_spinup", dest="ny_ad", default=600, \
                  help = 'number of years to run ad_spinup')
parser.add_option("--metdir", dest="metdir", default="none", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--C13", dest="C13", default=False, \
                      help = 'Switch to turn on C13', action="store_true")
parser.add_option("--C14", dest="C14", default=False, \
                  help = 'Use C14 as C13 (no decay)', action="store_true")
parser.add_option("--ninst", dest="ninst", default=1, \
                      help = 'number of land model instances')
parser.add_option("--npoolmod", action="store_true", dest="npoolmod", default=False, \
                    help="To turn on nitrogen pool modifications")
parser.add_option("--cpoolmod", action="store_true", dest="cpoolmod", default=False, \
                    help="To turn on carbon pool modifications")
parser.add_option("--q10wbmod", action="store_true", dest="q10wbmod", default=False, \
                    help="To turn on Woodrow-Berry Q10 curve")
parser.add_option("--tfmod", action="store_true", dest="tfmod", default=False, \
                    help="To set temperature threshold (0 degC) for plant wilting factor")
parser.add_option("--harvmod", action="store_true", dest='harvmod', default=False, \
                    help="turn on harvest modification:  All harvest at first timestep")
parser.add_option("--humhol", action="store_true", dest="humhol", \
                      default=False, help = "SPRUCE Hummock/Hollow modification")
parser.add_option("--vertsoilc", dest="vsoilc", default=False, \
                  help = 'To turn on CN with multiple soil layers, excluding CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--centbgc", dest="centbgc", default=False, \
                  help = 'To turn on CN with multiple soil layers, CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--CH4", dest="CH4", default=False, \
                  help = 'To turn on CN with CLM4me', action="store_true")
parser.add_option("--MICROBE", dest="MICROBE", default=False, \
                  help = 'To turn on MICROBE with CN', action="store_true")
parser.add_option("--makemetdata", action="store_true", dest="makemet", default=False, \
                    help="generate site meteorology")
parser.add_option("--refcase", dest="refcase" , default='none', \
                  help = 'Use already compiled CLM case')
parser.add_option("--xpts", dest="xpts", default=1, \
                      help = 'for regional runs: xpts')
parser.add_option("--ypts", dest="ypts", default=1, \
                      help = 'for regional runs: ypts')
parser.add_option("--cruncep", dest="cruncep", default=False, \
                  action="store_true", help = 'Use CRU-NCEP meteorology')
parser.add_option("--surfdata_grid", dest="surfdata_grid", default=False, \
                  help = 'Use gridded surface data instead of site data', action="store_true")
parser.add_option("--no_siteparms",dest = "no_siteparms", default=False, \
                  action="store_true", help = 'Use default PFT parameters')
parser.add_option("--cpl_bypass", dest = "cpl_bypass", default=False, \
                  help = "Bypass coupler", action="store_true")
parser.add_option("--spinup_vars", dest = "spinup_vars", default=False, help = "limit output variables for spinup", action="store_true")

(options, args) = parser.parse_args()

ccsm_input = options.ccsm_input
mycaseid   = options.mycaseid
srcmods    = options.srcmods_loc
csmdir     = os.path.abspath(options.csmdir)

#get start and year of input meteorology from site data file
fname = './PTCLM_files/PTCLM_sitedata/'+ \
    options.sitegroup+'_sitedata.txt'
AFdatareader = csv.reader(open(fname, "rb"))


for row in AFdatareader:
    if row[0] == options.site or (options.site == 'all' and row[0] !='site_code' \
                                      and row[0] != ''):
        site      = row[0]
        if (options.cruncep):
                startyear = 1901
                endyear = 1920
        else:
            startyear = int(row[6])
            endyear   = int(row[7])

        site_endyear = int(row[7])
        ncycle   = endyear-startyear+1   #number of years in met cycle
        translen = endyear-1850+1        #length of transient run

        #use site parameter file if it exists
        if (not options.no_siteparms):
            if (os.path.exists('./PTCLM_files/parms_'+site)):
                print ('Using parameter file PTCLM_Files/parms_'+site)
                options.parm_file = './PTCLM_files/parms_'+site
            else:
                options.parm_file = ''

        for i in range(0,ncycle+1):  #figure out length of final spinup run
            fsplen = int(options.nyears_final_spinup)+i
            if ((fsplen+translen) % ncycle == 0):
                break

        #get align_year for transient run
        year_align = (startyear-1850) % ncycle
        #print year_align, fsplen
        basecmd = 'python runCLM.py --site '+site+' --ccsm_input '+ \
            os.path.abspath(ccsm_input)+' --rmold --no_submit --sitegroup ' + \
            options.sitegroup+' --xpts '+str(options.xpts)+' --ypts '+str(options.ypts)
        basecmd = basecmd+' --machine '+options.machine
        basecmd = basecmd+' --runroot '+options.runroot
        if (srcmods != ''):
            srcmods    = os.path.abspath(srcmods)
            basecmd = basecmd+' --srcmods_loc '+srcmods
        if (mycaseid != ''):
            basecmd = basecmd+' --caseidprefix '+mycaseid
        #if (options.parm_file != ''):
        #    basecmd = basecmd+' --parm_file '+options.parm_file
        if (options.clean_build):
            basecmd = basecmd+' --clean_build '
        if (options.regional):
            basecmd = basecmd+' --xpts 2 --ypts 1 '
        if (options.metdir !='none'):
            basecmd = basecmd+' --metdir '+options.metdir
        if (options.C13):
            basecmd = basecmd+' --C13 '
        if (options.C14):
            basecmd = basecmd+' --C14 '
        if (options.ninst > 1):
            basecmd = basecmd+' --ninst '+str(options.ninst)
        if (options.npoolmod):
            basecmd = basecmd+' --npoolmod '
        if (options.cpoolmod):
            basecmd = basecmd+' --cpoolmod '
        if (options.q10wbmod):
            basecmd = basecmd+' --q10wbmod '
        if (options.tfmod):
            basecmd = basecmd+' --tfmod '
        if (options.refcase != 'none'):
            basecmd = basecmd+' --refcase '+options.refcase
        if (options.nofire):
            basecmd = basecmd+' --nofire'
        if (options.humhol):
            basecmd = basecmd+' --humhol'
        if (options.harvmod):
            basecmd = basecmd+' --harvmod'
        if (options.vsoilc):
            basecmd = basecmd+' --vertsoilc'
        if (options.centbgc):
            basecmd = basecmd+' --centbgc'
        if (options.CH4):
            basecmd = basecmd+' --CH4'
        if (options.MICROBE):
	    basecmd = basecmd+' --MICROBE'
        if (options.cruncep):
            basecmd = basecmd+' --cruncep'
        if (options.surfdata_grid):
            basecmd = basecmd+' --surfdata_grid'
        if (options.cpl_bypass):
            basecmd = basecmd+' --cpl_bypass'
        basecmd = basecmd + ' --np '+str(options.np)
        basecmd = basecmd + ' --tstep '+str(options.tstep)
        basecmd = basecmd + ' --co2_file '+options.co2_file
        basecmd = basecmd + ' --compiler '+options.compiler
        basecmd = basecmd + ' --mpilib '+options.mpilib

#----------------------- build commands for runCLM.py -----------------------------

        #AD spinup
        cmd_adsp = basecmd+' --ad_spinup --nyears_ad_spinup '+ \
            str(options.ny_ad)+' --hist_mfilt 1 --hist_nhtfrq -8760'
        if (options.makemet):
            cmd_adsp = cmd_adsp+' --makemetdat'
        if (options.spinup_vars):
            cmd_adsp = cmd_adsp+' --spinup_vars'
        ad_case = site+'_I1850CLM45CN_ad_spinup'
        if (mycaseid != ''):
            ad_case = mycaseid+'_'+ad_case
        #final spinup
        if mycaseid !='':
            basecase=mycaseid+'_'+site+'_I1850CLM45CN'
        else:
            basecase=site+'_I1850CLM45CN'
        cmd_fnsp = basecmd+' --finidat_case '+basecase+'_ad_spinup '+ \
            '--finidat_year '+str(int(options.ny_ad)+1)+' --run_units nyears --run_n '+ \
            str(fsplen)+' --hist_mfilt 1 --hist_nhtfrq -8760 --exeroot_case '+ad_case
        if (options.spinup_vars):
		cmd_fnsp = cmd_fnsp+' --spinup_vars'
        #transient
        cmd_trns = basecmd+' --finidat_case '+basecase+ \
            ' --finidat_year '+str(fsplen+1)+' --run_units nyears' \
            +' --run_n '+str(translen)+' --align_year '+ \
            str(year_align+1850)+' --hist_nhtfrq '+ \
            options.hist_nhtfrq+' --hist_mfilt '+options.hist_mfilt + \
            ' --compset I20TRCLM45CN --exeroot_case '+ad_case
        if (options.spinup_vars):
               cmd_trns = cmd_trns + ' --spinup_vars'
        #transient phase 2 (CRU-NCEP only)
        if (options.cruncep):
            basecase=basecase.replace('1850','20TR')+'_phase1'
            thistranslen = site_endyear - 1921 + 1
            cmd_trns2 = basecmd+' --trans2 --finidat_case '+basecase+ \
                ' --finidat_year 1921 --run_units nyears --branch ' \
                +' --run_n '+str(thistranslen)+' --align_year 1921'+ \
                ' --hist_nhtfrq '+options.hist_nhtfrq+' --hist_mfilt '+ \
                options.hist_mfilt+' --compset I20TRCLM45CN'

#---------------------------------------------------------------------------------

        #build cases
        print('\nSetting up ad_spinup case\n')
        os.system(cmd_adsp)
        print('\nSetting up final spinup case\n')
        os.system(cmd_fnsp)
        print('\nSetting up transient case\n')
        os.system(cmd_trns)
        if (options.cruncep):
             print('\nSetting up transient case phase 2\n')
             os.system(cmd_trns2)
                
        #Create the PBS script
        output = open('./site_fullrun.pbs','w')
        input = open(ad_case+'/'+ad_case+'.run')
        for s in input:
            if ("#!" in s or '#PBS' in s):
                output.write(s)
        input.close()
        output.write("\n")
        output.write("cd "+csmdir+'/scripts/'+ad_case+'\n')
        output.write("./"+ad_case+'.run\n')
        output.write("cd ../"+basecase+'\n')
        output.write("./"+basecase+'.run\n')
        #initialize (add the code here)
        output.write("cd ../"+basecase.replace('1850','20TR')+'\n')
        output.write("./"+basecase.replace('1850','20TR')+'.run\n')
        output.close()

        #submit
        os.system('qsub site_fullrun.pbs')

