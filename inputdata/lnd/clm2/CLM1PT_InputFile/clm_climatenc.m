% import ascii data of CLM climate and convert to .nc format in NACP

txtindir = 'monthly_txt';
ncoutdir = 'monthly_nc';

ncinfo_clm_jan = ncinfo('example-jan.nc')
ncinfo_clm_feb = ncinfo('example-feb.nc')
ncinfo_clm_apr = ncinfo('example-apr.nc')

fileheader = 'y';
fileender ='';

yrst = 2012;
yrend= 2012;
no_leapyear = 1;

% one station data: NOAA-ESRL-GMD, Barrow Baseline Observatory (AK)
LONGXY = 203.3886;
LATIXY = 71.3230;
ELEVXY = 11.0;
EDGEE  = 203.2886;
EDGEW  = 203.4886;
EDGES  = 71.2230;
EDGEN  = 71.4230;

clmdata.yyyy = [];
clmdata.mm   = [];
clmdata.time = [];
clmdata.ZBOT = [];
clmdata.WIND = [];
clmdata.PSRF = [];
clmdata.TBOT = [];
clmdata.RH   = [];
clmdata.PRECTmms = [];
clmdata.FSDS = [];
clmdata.FLDS = [];

for iyr=yrst:yrend;
    for im = 1:12
        %
        infile = fullfile(txtindir, [fileheader,num2str(iyr),'m',num2str(im), fileender]); 
        [time,ZBOT,WIND,PSRF,TBOT,RH,PRECTmms,FSDS,FLDS] = importclmtxtfile(infile);
        RH(RH>100.)=100.;
        WIND(WIND<0.01)=0.01;
        PSRF=PSRF*100.;
        
        clmdata.yyyy = [clmdata.yyyy; ones(size(time))*iyr];
        clmdata.mm   = [clmdata.mm; ones(size(time))*im];
        clmdata.time = [clmdata.time; time];      % day-time of the month
        clmdata.ZBOT = [clmdata.ZBOT; ZBOT];
        clmdata.WIND = [clmdata.WIND; WIND];
        clmdata.PSRF = [clmdata.PSRF; PSRF];
        clmdata.TBOT = [clmdata.TBOT; TBOT+273.15];   % CLM: tair in Kevin
        clmdata.RH   = [clmdata.RH; RH];
        clmdata.PRECTmms = [clmdata.PRECTmms; PRECTmms];
        clmdata.FSDS = [clmdata.FSDS; FSDS];
        clmdata.FLDS = [clmdata.FLDS; FLDS];
        
        %
        ofile = fullfile(ncoutdir, [num2str(iyr),'-',num2str(im, '%02d'), '.nc']);
        if (im == 2)
            ncsch = ncinfo_clm_feb;
            txt= ['days since ',num2str(iyr), '-',num2str(im,'%02d'),'-01 00:00:00'];
            ncsch.Variables(1,1).Attributes(1,2).Value=txt;

            ncwriteschema(ofile, ncsch);
            
            
            if (no_leapyear==1 && length(time)>1344)
               time(1345:end)=[];
               ZBOT(1345:end)=[];
               WIND(1345:end)=[];
               PSRF(1345:end)=[];
               TBOT(1345:end)=[];
               RH(1345:end)=[];
               PRECTmms(1345:end)=[];
               FSDS(1345:end)=[];
               FLDS(1345:end)=[];
               
            end
        elseif (im==1 | im==3 | im==5 | im==7 | im==8 | im==10 | im==12) 
            ncsch = ncinfo_clm_jan;
            txt= ['days since ',num2str(iyr), '-',num2str(im,'%02d'),'-01 00:00:00'];
            ncsch.Variables(1,1).Attributes(1,2).Value=txt;

            ncwriteschema(ofile, ncsch);
        elseif (im==4 | im==6 | im==9 | im==11) 
            ncsch = ncinfo_clm_apr;
            txt= ['days since ',num2str(iyr), '-',num2str(im,'%02d'),'-01 00:00:00'];
            ncsch.Variables(1,1).Attributes(1,2).Value=txt;

            ncwriteschema(ofile, ncsch);
        else
            display('month index is wrong!');
        end
        
        ncwrite(ofile, 'time', time);
        ncwrite(ofile, 'LONGXY', LONGXY);
        ncwrite(ofile, 'LATIXY', LATIXY);
        ncwrite(ofile, 'EDGEW', EDGEW);
        ncwrite(ofile, 'EDGEE', EDGEE);
        ncwrite(ofile, 'EDGES', EDGES);
        ncwrite(ofile, 'EDGEN', EDGEN);
        
        ncwrite(ofile, 'ZBOT', reshape(ZBOT, 1, 1, []));
        ncwrite(ofile, 'TBOT', reshape(TBOT, 1, 1, []));
        ncwrite(ofile, 'RH', reshape(RH, 1, 1, []));
        ncwrite(ofile, 'WIND', reshape(WIND, 1, 1, []));
        ncwrite(ofile, 'FSDS', reshape(FSDS, 1, 1, []));
        ncwrite(ofile, 'FLDS', reshape(FLDS, 1, 1, []));
        ncwrite(ofile, 'PSRF', reshape(PSRF, 1, 1, []));
        ncwrite(ofile, 'PRECTmms', reshape(PRECTmms, 1, 1, []));
        
    end

end