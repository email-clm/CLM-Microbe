# Xiaofeng Xu created this file in May 2020, and updated in Sep 2020, and finalized with 14C on Nov 11 2020
library("RNetCDF")
library("akima")
rm(list=ls())
setwd("/Users/xxuadmin/BUSINESS/PUBLICATIONS/WorkingOn_Xu_Isotope/SPRUCE_2020Sep")
# 13C is based on the PD standard Pee Dee Belemnite (PDB). 13C:12C = 0.0112372.
# 13C/12C(sample) / 13C/12C(standard - 1)
# 13C = (delta / 1000 + 1) * 0.0112372 * 12C

# for the spruce site, the 14C = 12C * ((1 + delta14C / 1000) * (1.176e-12/((1+delta13C)*(1+delta13C)/(1-0.025)/(1-0.025)))
# the 13Cdealta is -25 per mil for standard sample; this is for SPRUCE site based on the data from Erik Hobbie at UNH

data <- read.csv("Peat_Characteristics_T0_20180425_plot.csv", header=TRUE)
O13C <- data[1:16,3:4]
O14C <- data[1:16,5]
#depth 25, 15, 5, -5, -25, -35, -45, -55, -65, -85, -113, -163, -188, -225, -255
#CLM depth 0.71, 2.8, 6.2, 11.9, 21.2, 36.6, 61.97, 103.8, 172.8, 286.5, 473.9, 782.97, 1292.5, 2132.6, 3517.8 

depth <- c(0.0071, 0.0028, 0.0623, 0.1189, 0.2122, 0.366, 0.6197, 1.038, 1.72, 2.86)

#Hommock
#25 for layer 1, 2, 3
#15 for layer 4
# 5 for layer 5, 
#-5 for layer 6,
#-25, -35, -45 for layer 7
#-55, -65, and -85 for layer 8
# -113, -163, and -188 for layer 9
# -255 for layer 10 2.86

#Hollow
#-5 for layer 1, 2, 3
# -15 for layer 4, 
#-25 for layer 5,
#-35 for layer 6
#-45 -55, and -65 for layer 7 0.6197
#-85 and -133 for layer 8 1.038
# -163 and -188 for layer 9 1.72
# -225 and -275 for layer 10 2.86
SO13C = array(0, dim=c(15,2))
#L13C = O13C[3,1]
SO13C[1,1] = O13C[1,1]
SO13C[2,1] = O13C[1,1]
SO13C[3,1] = O13C[1,1]
SO13C[4,1] = O13C[2,1]
SO13C[5,1] = O13C[3,1]
SO13C[6,1] = O13C[4,1]
SO13C[7,1] = mean(O13C[5:7,1])
SO13C[8,1] = mean(O13C[8:10,1])
SO13C[9,1] = mean(O13C[11:13,1])
SO13C[10,1] = mean(O13C[14:15,1])

SO13C[1,2] = mean(O13C[3:4,1])
SO13C[2,2] = mean(O13C[3:4,1])
SO13C[3,2] = mean(O13C[3:4,1])
SO13C[4,2] = mean(O13C[4:5,1])
SO13C[5,2] = mean(O13C[6:7,1])
SO13C[6,2] = mean(O13C[8:9,1])
SO13C[7,2] = mean(O13C[8:9,1])
SO13C[8,2] = mean(O13C[10:11,1])
SO13C[9,2] = mean(O13C[12:13,1])
SO13C[10,2] = mean(O13C[14:16,1])

#plot(SO13C[1:10,1],-depth,type="l",xlab="soil 13CH4 (mmol/L)", ylab="soil depth (m)",cex=1,las=1, cex.lab = 1.5, cex.axis = 1.5,lwd=2)

SO14C = array(0, dim=c(15,2))
SO14C[1,1] = O14C[1]
SO14C[2,1] = O14C[1]
SO14C[3,1] = O14C[1]
SO14C[4,1] = O14C[2]
SO14C[5,1] = O14C[3]
SO14C[6,1] = O14C[4]
SO14C[7,1] = mean(O14C[5:7])
SO14C[8,1] = mean(O14C[8:10])
SO14C[9,1] = mean(O14C[11:13])
SO14C[10,1] = mean(O14C[14:15])

SO14C[1,2] = O14C[4]
SO14C[2,2] = O14C[4]
SO14C[3,2] = O14C[4]
SO14C[4,2] = O14C[5]
SO14C[5,2] = O14C[6]
SO14C[6,2] = mean(O14C[7:9])
SO14C[7,2] = mean(O14C[7:9])
SO14C[8,2] = mean(O14C[10:11])
SO14C[9,2] = mean(O14C[12:13])
SO14C[10,2] = mean(O14C[14:16])

#13C_ave	Depth_ave
#-29.058	25
#-28.34	15
#-28.1075	5
#-28.77631579	-5
#-27.79684211	-15
#-27.006	-25
#-26.72888889	-35
#-26.45055556	-45
#-26.60555556	-55
#-26.36611111	-65
#-26.05333333	-85
#-25.93	-113
#-26.39176471	-163
#-26.27	-188
#-26.63857143	-225
#-25.75	-275

#depth	14C_mean	14Cto12C
#25	40.65714286	1.23406E-12
#15	54.37692308	1.24849E-12
#5	112.4714286	1.31665E-12
#-5	68.17647059	1.26596E-12
#-15	140.0176471	1.34839E-12
#-25	114.2	        1.31571E-12
#-35	-99.17647059	1.06314E-12
#-45	-180.2058824	9.66953E-13
#-55	-199.9941176	9.43913E-13
#-65	-286.2058824	8.41779E-13
#-85	-358.8117647	7.55669E-13
#-113	-392.4058824	7.15896E-13
#-163	-486.13125	6.06039E-13
#-188	-350.8	        7.65452E-13
#-225	-571.1	        5.06086E-13
#-275	-703.75	        3.48927E-13

computC13 <- function(C12, deltaC)
{
        (deltaC / 1000.0 + 1.0) * 0.0112372 * C12
}

computC14 <- function(C12, deltaC13, deltaC14)
{
        C12 * ((1.0 + deltaC14 / 1000) * (1.176e-12/((1+deltaC13)*(1+deltaC13)/(1-0.000025)/(1-0.000025))))
}
        
#14C = 12C * ((1 + delta14C / 1000) * (1.176e-12/((1+delta13C)*(1+delta13C)/(1-0.025)/(1-0.025)))
             
#notice: please copy a file with the name of fileoutput in the following code for variable to be incorporated
fileinput <- open.nc("SPRUCE-finalspinup-peatland-carbon-initial.nc")  # 2012 year climate data
#filerestart <- open.nc("SPR_I20TRCLM45CN_adspinup.clm2.r.1201-01-01-00000.nc")
fileoutput <- open.nc("SPRUCE-finalspinup-peatland-ISOcarbon-initial.nc",write=TRUE)

#depth <- c(0.0071, 0.0028, 0.0623, 0.1189, 0.2122, 0.366, 0.6197, 1.038, 1.72, 2.86)
#thick <- c(0.0175128179, 0.0275789693, 0.0454700332, 0.0749674110, 0.1236003651, 0.2037825510, 0.3359806264, 0.5539384054, 0.9132900316, 1.5057607014)
#realcarbondensity = c(22385, 22385, 22385, 30294, 43978, 95584, 92983, 83658, 83799, 90616)

cwdc_vr <- var.get.nc(fileinput, "cwdc_vr")
litr1c_vr <- var.get.nc(fileinput, "litr1c_vr")
litr2c_vr <- var.get.nc(fileinput, "litr2c_vr")
litr3c_vr <- var.get.nc(fileinput, "litr3c_vr")
soil1c_vr <- var.get.nc(fileinput, "soil1c_vr")
soil2c_vr <- var.get.nc(fileinput, "soil2c_vr")
soil3c_vr <- var.get.nc(fileinput, "soil3c_vr")
soil4c_vr <- var.get.nc(fileinput, "soil4c_vr")

cwdn_vr <- var.get.nc(fileinput, "cwdn_vr")
litr1n_vr <- var.get.nc(fileinput, "litr1c_vr")
litr2n_vr <- var.get.nc(fileinput, "litr2n_vr")
litr3n_vr <- var.get.nc(fileinput, "litr3n_vr")
soil1n_vr <- var.get.nc(fileinput, "soil1n_vr")
soil2n_vr <- var.get.nc(fileinput, "soil2n_vr")
soil3n_vr <- var.get.nc(fileinput, "soil3n_vr")
soil4n_vr <- var.get.nc(fileinput, "soil4n_vr")

#rc13_canair <- var.get.nc(filerestart, "rc13_canair")
#rc13_psnsun <- var.get.nc(filerestart, "rc13_psnsun")
#rc13_psnsha <- var.get.nc(filerestart, "rc13_psnsha")
#leafc_13 <- var.get.nc(filerestart, "leafc_13")
#leafc_storage_13 <- var.get.nc(filerestart, "leafc_storage_13")
#leafc_xfer_13 <- var.get.nc(filerestart, "leafc_xfer_13")
#frootc_13 <- var.get.nc(filerestart, "frootc_13")
#frootc_storage_13 <- var.get.nc(filerestart, "frootc_storage_13")
#frootc_xfer_13 <- var.get.nc(filerestart, "frootc_xfer_13")
#livestemc_13 <- var.get.nc(filerestart, "livestemc_13")
#livestemc_storage_13 <- var.get.nc(filerestart, "livestemc_storage_13")
#livestemc_xfer_13 <- var.get.nc(filerestart, "livestemc_xfer_13")
#deadstemc_13 <- var.get.nc(filerestart, "deadstemc_13")
#deadstemc_storage_13 <- var.get.nc(filerestart, "deadstemc_storage_13")
#deadstemc_xfer_13 <- var.get.nc(filerestart, "deadstemc_xfer_13")
#livecrootc_13 <- var.get.nc(filerestart, "livecrootc_13")
#deadstemc_13 <- var.get.nc(filerestart, "deadstemc_13")

cwdc_13_vr = cwdc_vr
litr1c_13_vr = litr1c_vr
litr2c_13_vr = litr2c_vr
litr3c_13_vr = litr3c_vr
soil1c_13_vr = soil1c_vr
soil2c_13_vr = soil2c_vr
soil3c_13_vr = soil3c_vr
soil4c_13_vr = soil4c_vr

# C14 
cwdc_14_vr = cwdc_vr
litr1c_14_vr = litr1c_vr
litr2c_14_vr = litr2c_vr
litr3c_14_vr = litr3c_vr
soil1c_14_vr = soil1c_vr
soil2c_14_vr = soil2c_vr
soil3c_14_vr = soil3c_vr
soil4c_14_vr = soil4c_vr

bacteriac_vr = litr1c_vr
fungic_vr = litr1c_vr
domc_vr = litr1c_vr
seedc = litr1c_vr
col_ctrunc_vr = litr1c_vr

totlitc = litr1c_vr
totcolc = litr1c_vr
prod10c = litr1c_vr
prod100c = litr1c_vr

bacteriac_13_vr = bacteriac_vr
fungic_13_vr = fungic_vr
domc_13_vr = domc_vr
seedc_13 = seedc
col_ctrunc_13_vr = col_ctrunc_vr

totlitc_13 = totlitc
totcolc_13 = totcolc
prod10c_13 = prod10c
prod100c_13 = prod100c

#C14
bacteriac_14_vr = bacteriac_vr
fungic_14_vr = fungic_vr
domc_14_vr = domc_vr
seedc_14 = seedc
col_ctrunc_14_vr = col_ctrunc_vr

totlitc_14 = totlitc
totcolc_14 = totcolc
prod10c_14 = prod10c
prod100c_14 = prod100c

for(i in 1:10)
{
cwdc_13_vr[i,1] = computC13(cwdc_vr[i,1], SO13C[i])
litr1c_13_vr[i,1] = computC13(litr1c_vr[i,1], SO13C[i])
litr2c_13_vr[i,1] = computC13(litr2c_vr[i,1], SO13C[i])
litr3c_13_vr[i,1] = computC13(litr3c_vr[i,1], SO13C[i])
soil1c_13_vr[i,1] = computC13(soil1c_vr[i,1], SO13C[i])
soil2c_13_vr[i,1] = computC13(soil2c_vr[i,1], SO13C[i])
soil3c_13_vr[i,1] = computC13(soil3c_vr[i,1], SO13C[i])
soil4c_13_vr[i,1] = computC13(soil4c_vr[i,1], SO13C[i])

bacteriac_13_vr[i,1] = computC13(bacteriac_vr[i,1], SO13C[i])
fungic_13_vr[i,1] = computC13(fungic_vr[i,1], SO13C[i])
domc_13_vr[i,1] = computC13(domc_vr[i,1], SO13C[i])
seedc_13[i,1] = computC13(seedc[i,1], SO13C[i])
col_ctrunc_13_vr[i,1] = computC13(col_ctrunc_vr[i,1], SO13C[i])
totlitc_13[i,1] = computC13(totlitc[i,1], SO13C[i])
totcolc_13[i,1] = computC13(totcolc[i,1], SO13C[i])
prod10c_13[i,1] = computC13(prod10c[i,1], SO13C[i])
prod100c_13[i,1] = computC13(prod100c[i,1], SO13C[i])

cwdc_13_vr[i,2] = computC13(cwdc_vr[i,2], SO13C[i])
litr1c_13_vr[i,2] = computC13(litr1c_vr[i,2], SO13C[i])
litr2c_13_vr[i,2] = computC13(litr2c_vr[i,2], SO13C[i])
litr3c_13_vr[i,2] = computC13(litr3c_vr[i,2], SO13C[i])
soil1c_13_vr[i,2] = computC13(soil1c_vr[i,2], SO13C[i])
soil2c_13_vr[i,2] = computC13(soil2c_vr[i,2], SO13C[i])
soil3c_13_vr[i,2] = computC13(soil3c_vr[i,2], SO13C[i])
soil4c_13_vr[i,2] = computC13(soil4c_vr[i,2], SO13C[i])

bacteriac_13_vr[i,2] = computC13(bacteriac_vr[i,2], SO13C[i])
fungic_13_vr[i,2] = computC13(fungic_vr[i,2], SO13C[i])
domc_13_vr[i,2] = computC13(domc_vr[i,2], SO13C[i])
seedc_13[i,2] = computC13(seedc[i,2], SO13C[i])
col_ctrunc_13_vr[i,2] = computC13(col_ctrunc_vr[i,2], SO13C[i])
totlitc_13[i,2] = computC13(totlitc[i,2], SO13C[i])
totcolc_13[i,2] = computC13(totcolc[i,2], SO13C[i])
prod10c_13[i,2] = computC13(prod10c[i,2], SO13C[i])
prod100c_13[i,2] = computC13(prod100c[i,2], SO13C[i])
}

for(i in 11:15)
{
cwdc_13_vr[i,1] = 0.0
litr1c_13_vr[i,1] = 0.0
litr2c_13_vr[i,1] = 0.0
litr3c_13_vr[i,1] = 0.0
soil1c_13_vr[i,1] = 0.0
soil2c_13_vr[i,1] = 0.0
soil3c_13_vr[i,1] = 0.0
soil4c_13_vr[i,1] = 0.0
        
bacteriac_13_vr[i,1] = 0.0
fungic_13_vr[i,1] = 0.0
domc_13_vr[i,1] = 0.0
seedc_13[i,1] = 0.0
col_ctrunc_13_vr[i,1] = 0.0

totlitc_13[i,1] = 0.0
totcolc_13[i,1] = 0.0
prod10c_13[i,1] = 0.0
prod100c_13[i,1] = 0.0
        
cwdc_13_vr[i,2] = 0.0
litr1c_13_vr[i,2] = 0.0
litr2c_13_vr[i,2] = 0.0
litr3c_13_vr[i,2] = 0.0
soil1c_13_vr[i,2] = 0.0
soil2c_13_vr[i,2] = 0.0
soil3c_13_vr[i,2] = 0.0
soil4c_13_vr[i,2] = 0.0
        
bacteriac_13_vr[i,2] = 0.0
fungic_13_vr[i,2] = 0.0
domc_13_vr[i,2] = 0.0
seedc_13[i,2] = 0.0
col_ctrunc_13_vr[i,2] = 0.0
        
totlitc_13[i,2] = 0.0
totcolc_13[i,2] = 0.0
prod10c_13[i,2] = 0.0
prod100c_13[i,2] = 0.0
}

for(i in 1:10)
{
        cwdc_14_vr[i,1] = computC14(cwdc_vr[i,1], cwdc_13_vr[i,1], SO14C[i])
        litr1c_14_vr[i,1] = computC14(litr1c_vr[i,1], litr1c_13_vr[i,1], SO14C[i])
        litr2c_14_vr[i,1] = computC14(litr2c_vr[i,1], litr2c_13_vr[i,1], SO14C[i])
        litr3c_14_vr[i,1] = computC14(litr3c_vr[i,1], litr3c_13_vr[i,1], SO14C[i])
        soil1c_14_vr[i,1] = computC14(soil1c_vr[i,1], soil1c_13_vr[i,1], SO14C[i])
        soil2c_14_vr[i,1] = computC14(soil2c_vr[i,1], soil2c_13_vr[i,1], SO14C[i])
        soil3c_14_vr[i,1] = computC14(soil3c_vr[i,1], soil3c_13_vr[i,1], SO14C[i])
        soil4c_14_vr[i,1] = computC14(soil4c_vr[i,1], soil4c_13_vr[i,1], SO14C[i])
        
        bacteriac_14_vr[i,1] = computC14(bacteriac_vr[i,1], bacteriac_13_vr[i,1], SO14C[i])
        fungic_14_vr[i,1] = computC14(fungic_vr[i,1], fungic_13_vr[i,1], SO14C[i])
        domc_14_vr[i,1] = computC14(domc_vr[i,1], domc_13_vr[i,1], SO14C[i])
        seedc_14[i,1] = computC14(seedc[i,1], seedc_13[i,1], SO14C[i])
        col_ctrunc_14_vr[i,1] = computC14(col_ctrunc_vr[i,1], col_ctrunc_13_vr[i,1], SO14C[i])
        totlitc_14[i,1] = computC14(totlitc[i,1], totlitc_13[i,1], SO14C[i])
        totcolc_14[i,1] = computC14(totcolc[i,1], totcolc_13[i,1], SO14C[i])
        prod10c_14[i,1] = computC14(prod10c[i,1], prod10c_13[i,1], SO14C[i])
        prod100c_14[i,1] = computC14(prod100c[i,1], prod100c_13[i,1], SO14C[i])
        
        cwdc_14_vr[i,2] = computC14(cwdc_vr[i,2], cwdc_13_vr[i,2], SO14C[i])
        litr1c_14_vr[i,2] = computC14(litr1c_vr[i,2], litr1c_13_vr[i,2], SO14C[i])
        litr2c_14_vr[i,2] = computC14(litr2c_vr[i,2], litr2c_13_vr[i,2], SO14C[i])
        litr3c_14_vr[i,2] = computC14(litr3c_vr[i,2], litr3c_13_vr[i,2], SO14C[i])
        soil1c_14_vr[i,2] = computC14(soil1c_vr[i,2], soil1c_13_vr[i,2], SO14C[i])
        soil2c_14_vr[i,2] = computC14(soil2c_vr[i,2], soil2c_13_vr[i,2], SO14C[i])
        soil3c_14_vr[i,2] = computC14(soil3c_vr[i,2], soil3c_13_vr[i,2], SO14C[i])
        soil4c_14_vr[i,2] = computC14(soil4c_vr[i,2], soil4c_13_vr[i,2], SO14C[i])
        
        bacteriac_14_vr[i,2] = computC14(bacteriac_vr[i,2], bacteriac_13_vr[i,2], SO14C[i])
        fungic_14_vr[i,2] = computC14(fungic_vr[i,2], fungic_13_vr[i,2], SO14C[i])
        domc_14_vr[i,2] = computC14(domc_vr[i,2], domc_13_vr[i,2], SO14C[i])
        seedc_14[i,2] = computC14(seedc[i,2], seedc_13[i,2], SO14C[i])
        col_ctrunc_14_vr[i,2] = computC14(col_ctrunc_vr[i,2], col_ctrunc_13_vr[i,2], SO14C[i])
        totlitc_14[i,2] = computC14(totlitc[i,2], totlitc_13[i,2], SO14C[i])
        totcolc_14[i,2] = computC14(totcolc[i,2], totcolc_13[i,2], SO14C[i])
        prod10c_14[i,2] = computC14(prod10c[i,2], prod10c_13[i,2], SO14C[i])
        prod100c_14[i,2] = computC14(prod100c[i,2], prod100c_13[i,2], SO14C[i])
}

for(i in 11:15)
{
        cwdc_14_vr[i,1] = 0.0
        litr1c_14_vr[i,1] = 0.0
        litr2c_14_vr[i,1] = 0.0
        litr3c_14_vr[i,1] = 0.0
        soil1c_14_vr[i,1] = 0.0
        soil2c_14_vr[i,1] = 0.0
        soil3c_14_vr[i,1] = 0.0
        soil4c_14_vr[i,1] = 0.0
        
        bacteriac_14_vr[i,1] = 0.0
        fungic_14_vr[i,1] = 0.0
        domc_14_vr[i,1] = 0.0
        seedc_14[i,1] = 0.0
        col_ctrunc_14_vr[i,1] = 0.0
        
        totlitc_14[i,1] = 0.0
        totcolc_14[i,1] = 0.0
        prod10c_14[i,1] = 0.0
        prod100c_14[i,1] = 0.0
        
        cwdc_14_vr[i,2] = 0.0
        litr1c_14_vr[i,2] = 0.0
        litr2c_14_vr[i,2] = 0.0
        litr3c_14_vr[i,2] = 0.0
        soil1c_14_vr[i,2] = 0.0
        soil2c_14_vr[i,2] = 0.0
        soil3c_14_vr[i,2] = 0.0
        soil4c_14_vr[i,2] = 0.0
        
        bacteriac_14_vr[i,2] = 0.0
        fungic_14_vr[i,2] = 0.0
        domc_14_vr[i,2] = 0.0
        seedc_14[i,2] = 0.0
        col_ctrunc_14_vr[i,2] = 0.0
        
        totlitc_14[i,2] = 0.0
        totcolc_14[i,2] = 0.0
        prod10c_14[i,2] = 0.0
        prod100c_14[i,2] = 0.0
}


# 
#rc13_canair = c(0.005582864, 0.002791432)
#rc13_psnsha = c(0.005544063, 0.002771556)
#rc13_psnsun = c(0.0055442, 0.002771556)
#

#levgrnd = 15
#column = 2
#var.def.nc(fileoutput,"levgrnd","NC_FLOAT","levgrnd")
#var.def.nc(fileoutput,"column","NC_FLOAT","column")

var.def.nc(fileoutput,"cwdc_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"litr1c_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"litr2c_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"litr3c_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"soil1c_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"soil2c_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"soil3c_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"soil4c_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"bacteriac_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"fungic_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"domc_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"seedc_13","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"col_ctrunc_13_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"totlitc_13","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"totcolc_13","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"prod10c_13","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"prod100c_13","NC_FLOAT",c("levgrnd","column"))

var.def.nc(fileoutput,"cwdc_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"litr1c_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"litr2c_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"litr3c_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"soil1c_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"soil2c_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"soil3c_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"soil4c_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"bacteriac_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"fungic_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"domc_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"seedc_14","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"col_ctrunc_14_vr","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"totlitc_14","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"totcolc_14","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"prod10c_14","NC_FLOAT",c("levgrnd","column"))
var.def.nc(fileoutput,"prod100c_14","NC_FLOAT",c("levgrnd","column"))

var.put.nc(fileoutput, "cwdc_vr",cwdc_vr)
var.put.nc(fileoutput, "litr1c_vr",litr1c_vr)
var.put.nc(fileoutput, "litr2c_vr",litr2c_vr)
var.put.nc(fileoutput, "litr3c_vr",litr3c_vr)
var.put.nc(fileoutput, "soil1c_vr",soil1c_vr)
var.put.nc(fileoutput, "soil2c_vr",soil2c_vr)
var.put.nc(fileoutput, "soil3c_vr",soil3c_vr)
var.put.nc(fileoutput, "soil4c_vr",soil4c_vr)

var.put.nc(fileoutput, "cwdn_vr",cwdn_vr)
var.put.nc(fileoutput, "litr1n_vr",litr1n_vr)
var.put.nc(fileoutput, "litr2n_vr",litr2n_vr)
var.put.nc(fileoutput, "litr3n_vr",litr3n_vr)
var.put.nc(fileoutput, "soil1n_vr",soil1n_vr)
var.put.nc(fileoutput, "soil2n_vr",soil2n_vr)
var.put.nc(fileoutput, "soil3n_vr",soil3n_vr)
var.put.nc(fileoutput, "soil4n_vr",soil4n_vr)

att.put.nc(fileoutput, "cwdc_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "cwdc_13_vr",cwdc_13_vr)

att.put.nc(fileoutput, "litr1c_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "litr1c_13_vr",litr1c_13_vr)

att.put.nc(fileoutput, "litr2c_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "litr2c_13_vr",litr2c_13_vr)

att.put.nc(fileoutput, "litr3c_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "litr3c_13_vr",litr3c_13_vr)

att.put.nc(fileoutput, "soil1c_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "soil1c_13_vr",soil1c_13_vr)

att.put.nc(fileoutput, "soil2c_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "soil2c_13_vr",soil2c_13_vr)

att.put.nc(fileoutput, "soil3c_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "soil3c_13_vr",soil3c_13_vr)

att.put.nc(fileoutput, "soil4c_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "soil4c_13_vr",soil4c_13_vr)

att.put.nc(fileoutput, "bacteriac_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "bacteriac_13_vr",bacteriac_13_vr)

att.put.nc(fileoutput, "fungic_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "fungic_13_vr",fungic_13_vr)

att.put.nc(fileoutput, "domc_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "domc_13_vr",domc_13_vr)

att.put.nc(fileoutput, "seedc_13","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "seedc_13",seedc_13)

att.put.nc(fileoutput, "col_ctrunc_13_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "col_ctrunc_13_vr",col_ctrunc_13_vr)

att.put.nc(fileoutput, "totlitc_13","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "totlitc_13",totlitc_13)

att.put.nc(fileoutput, "totcolc_13","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "totcolc_13",totcolc_13)

att.put.nc(fileoutput, "prod10c_13","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "prod10c_13",prod10c_13)

att.put.nc(fileoutput, "prod100c_13","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "prod100c_13",prod100c_13)

# 14C
att.put.nc(fileoutput, "cwdc_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "cwdc_14_vr",cwdc_14_vr)

att.put.nc(fileoutput, "litr1c_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "litr1c_14_vr",litr1c_14_vr)

att.put.nc(fileoutput, "litr2c_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "litr2c_14_vr",litr2c_14_vr)

att.put.nc(fileoutput, "litr3c_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "litr3c_14_vr",litr3c_14_vr)

att.put.nc(fileoutput, "soil1c_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "soil1c_14_vr",soil1c_14_vr)

att.put.nc(fileoutput, "soil2c_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "soil2c_14_vr",soil2c_14_vr)

att.put.nc(fileoutput, "soil3c_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "soil3c_14_vr",soil3c_14_vr)

att.put.nc(fileoutput, "soil4c_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "soil4c_14_vr",soil4c_14_vr)

att.put.nc(fileoutput, "bacteriac_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "bacteriac_14_vr",bacteriac_14_vr)

att.put.nc(fileoutput, "fungic_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "fungic_14_vr",fungic_14_vr)

att.put.nc(fileoutput, "domc_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "domc_14_vr",domc_14_vr)

att.put.nc(fileoutput, "seedc_14","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "seedc_14",seedc_14)

att.put.nc(fileoutput, "col_ctrunc_14_vr","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "col_ctrunc_14_vr",col_ctrunc_14_vr)

att.put.nc(fileoutput, "totlitc_14","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "totlitc_14",totlitc_14)

att.put.nc(fileoutput, "totcolc_14","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "totcolc_14",totcolc_14)

att.put.nc(fileoutput, "prod10c_14","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "prod10c_14",prod10c_14)

att.put.nc(fileoutput, "prod100c_14","missing_value", "NC_FLOAT", -9999)
var.put.nc(fileoutput, "prod100c_14",prod100c_14)

# new file
#var.def.nc(fileoutput,"rc13_canair","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "rc13_canair","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "rc13_canair",rc13_canair)

#var.def.nc(fileoutput,"rc13_psnsun","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "rc13_psnsun","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "rc13_psnsun",rc13_psnsun)

#var.def.nc(fileoutput,"rc13_psnsha","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "rc13_psnsha","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "rc13_psnsha",rc13_psnsha)

#var.def.nc(fileoutput,"leafc_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "leafc_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "leafc_13",leafc_13)

#var.def.nc(fileoutput,"leafc_storage_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "leafc_storage_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "leafc_storage_13",leafc_storage_13)

#var.def.nc(fileoutput,"leafc_xfer_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "leafc_xfer_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "leafc_xfer_13",leafc_xfer_13)

#var.def.nc(fileoutput,"frootc_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "frootc_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "frootc_13",frootc_13)

#var.def.nc(fileoutput,"frootc_storage_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "frootc_storage_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "frootc_storage_13",frootc_storage_13)

#var.def.nc(fileoutput,"frootc_xfer_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "frootc_xfer_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "frootc_xfer_13",frootc_xfer_13)

#var.def.nc(fileoutput,"livestemc_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "livestemc_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "livestemc_13",livestemc_13)

#var.def.nc(fileoutput,"livestemc_storage_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "livestemc_storage_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "livestemc_storage_13",livestemc_storage_13)

#var.def.nc(fileoutput,"livestemc_xfer_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "livestemc_xfer_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "livestemc_xfer_13",livestemc_xfer_13)

#var.def.nc(fileoutput,"deadstemc_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "deadstemc_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "deadstemc_13",deadstemc_13)

#var.def.nc(fileoutput,"deadstemc_storage_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "deadstemc_storage_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "deadstemc_storage_13",deadstemc_storage_13)

#var.def.nc(fileoutput,"deadstemc_xfer_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "deadstemc_xfer_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "deadstemc_xfer_13",deadstemc_xfer_13)

#var.def.nc(fileoutput,"livecrootc_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "livecrootc_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "livecrootc_13",livecrootc_13)

#var.def.nc(fileoutput,"deadstemc_13","NC_FLOAT",c("pft"))
#att.put.nc(fileoutput, "deadstemc_13","missing_value", "NC_FLOAT", -9999)
#var.put.nc(fileoutput, "deadstemc_13",deadstemc_13)

att.put.nc(fileoutput,"NC_GLOBAL","title","NC_CHAR","isotopic 13 carbon pools and 14 carbon pools")

close.nc(fileinput)
close.nc(fileoutput)

