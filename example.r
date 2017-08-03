# R utilities created for data analysis of "Evaluating Canopy 
# Spectral Reflectance Vegetation Indices to Estimate Nitrogen 
# Use Traits in Hard Winter Wheat"
# Copyright (C) 2017  Brian Joyce
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

# This is an example script of how the utilities were used to calculate
# reflectance indices
#

setwd('C:/Users/brian/Dropbox/Manuscripts/Github/')
source('DirtyMySQL.r')
source('reflectance_equation.r')
source('derivative_reflectance_equation.r')
source('remove_bumple.r')

# initialize DMSR so it can be used by the functions
dmsr.init(
  mysql_path = 'C:\\Program Files\\MySQL\\MySQL Server 5.7\\bin\\mysql.exe',
  hostname = 'localhost',
  username = 'bdjoyce',
  passwd = 'unsecure',
  database = 'agr'
);


# basic reflectance indices use the reflectance data and therefore the reflectance_equation.r function

#Maccioni	(R780-R710)/(R780-R680)
Maccioni = reflectance_equation(" ( median(r1)-median(r2) ) / ( median(r1)-median(r3) ) ", 780, 710, 680)

#mSR2	(R750/R705)-1/SQRT((R750-R705)+1)
mSR2 = reflectance_equation("median(r1)/median(r2) - 1/sqrt(median(r1)+median(r2)) + 1", 750, 705, 705)

#EVI	2.5*(R800-R670)/(R800-(6*R670)-(7.5*R475)+1))
EVI = reflectance_equation("2.5 * ( (median(r1) - median(r2)) / (median(r1) - (6*median(r2)) - (7.5*median(r3)) + 1) )", 800, 670, 475)



# derivative reflectance indices use the derivative_reflectance_equation.r function

#Boochs2	D720
Boochs2 = derivative_reflectance_equation("dr1", 720, 720, 720, 720)

#Vogelmann3	D715/D705
Vogelmann3 = derivative_reflectance_equation("dr1/dr2", 715, 705, 705, 705)



# the rededge position indices are a bit more complicated and so don't fit the same structure

# REIP	Wavelength for maximum value of the 1st derivative in red-edge region
#luckily this is something the function calculates automatically and returns for use as reip
REIP = derivative_reflectance_equation("reip", 0, 0, 0, 0)

# REIPLE	Linear Extrapolation from Cho and skidmore
# this fits in the derivative_reflectance_equation but the equations are gnarly
#
#y1 = 679.65
#y2 = 694.30
#y3 = 732.46
#y4 = 760.41
#m1 = dr2-dr1 / w2-w1
#c1 = dr1 - (m1*w1)
#m2 = dr4-dr3 / w4-w3
#c2 = dr3 - (m2*w3)
#reiple = -1*(c2-c1)/(m2-m1) 
#
REIPLE = derivative_reflectance_equation("-1*((dr3 - ((dr4-dr3 / w4-w3)*w3))-(dr1 - ((dr2-dr1 / w2-w1)*w1)))/((dr4-dr3 / w4-w3)-(dr2-dr1 / w2-w1))", 679.65, 694.30, 732.46, 760.41)

# mREIP  inverted Gaussian rededge position
# This one was too complex to fit in the reflectance_equation function so a modified fucntion was 
# used and saved it as an external sheet <ig_rededge_reflectance.r>
source("ig_rededge_reflectance.r")


