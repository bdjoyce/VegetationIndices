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

library("signal")


# derivative_reflectance_equation
# ================================
# Function to generalize execution of reflectance equations relating to the derivative of the
# spectrum to make large batch analysis easier allows a reflectance equation to be listed in 
# expression and the function queries the database to find the data needed to evaluate the equation.
# The derivative, noise removal functions, and filtering all happen automatically. Up to 4 
# wavelengths can be used in the calculation.
#
# NOTE:  function assumes dmsr library is available to run the queries
# NOTE:  function assumes remove_bumple function is sourced and available to run the filtering
# NOTE:  function assumes database structure of at least ...
# +--------------------------------------------------------------------+
# | date | plot | wavelength | reflectance | serial_number_1 | exclude |
# +--------------------------------------------------------------------+
# NOTE:  the following analysis performed on the data
# 1:  Any plot/day measurement made with the same serial_number_1 multiple times is averaged
# 2:  Any plot/day measurement made with different serial_number_1 will be averaged after the median is calculated
#
# Usage:
# expression 	 = (string) 	reflectance equation to evaluate
#
#								y1, y2, y3, y4 correspond to wavelength and reflectance variables that can be used in the expression
#								w1, w2, w3, y4 correspond to the exact wavelengths available in the dataset for the y1 y2 y3 values
#								dr1, dr2, dr3, and dr4 are the derivative reflectance values that correspond to the wavelengths (y1 y2 y3 y4)
#									the r's used are already filtered so single values are returned
#								reip corresponds to the wavelength for maximum value of the 1st derivative in red-edge region
#
#							Example expressions
#         						Vogelmann3: "dr1/dr2"
#
# y1, y2, y3, y4 = (float) 		4 wavelengths to pull data for use in expression
# table			 = (string)  	table to use when executing queries
#
# Returns:
# result		 = (data frame) a data frame containing each plot/day data along with the wavelengths 
#								pulled, reflectances for those wavelengths, and the results of the 
#								calculated reflectance equation
#
# Example:
# Maccioni = reflectance_equation("(median(r1)-median(r2))/(median(r1)-median(r3))", 780, 710, 680, "2014_CSRV")
#
derivative_reflectance_equation <- function(expression, y1, y2, y3, y4, table="2014_CSRV") {
  
  #run some setup queries to get information used in looping
  serials = dmsr.query(paste('SELECT DISTINCT serial_number_1 FROM', table ,' WHERE exclude IS NULL', sep=" "))$serial_number_1
  ys = dmsr.query(paste('SELECT distinct(wavelength) FROM ', table, ' WHERE exclude IS NULL', sep=" "))$wavelength
  dates = dmsr.query(paste('SELECT distinct(date) FROM ', table, ' WHERE exclude IS NULL ORDER by date ASC', sep=" "))$date
  plots = dmsr.query(paste('SELECT distinct(plot) FROM ',table,' WHERE exclude IS NULL ORDER BY plot ASC', sep=" "))$plot
  
  #initialize data structures so we can push data onto them later in the fuction
  dr1c = c()
  dr2c = c()
  dr3c = c()
  dr4c = c()
  reipc = c()
  datec = c()
  plotc = c()
  eqc = c()
  idx = 1
    
  for (serial in serials) {
        
    #get all the wavelengths all the data
	#more data is more filtering and more processor cycles, the 500 and 900 limits here speed up the calcuation but can also be removed is a wider range is needed
    query = paste('SELECT date, plot, wavelength, AVG(reflectance) AS avg_ref FRoM ', table, ' WHERE serial_number_1=\'', serial, '\' AND wavelength>500 AND wavelength<900 AND exclude IS NULL GROUP BY date, plot, wavelength order by date, plot, wavelength ASC;', sep='');
    res = dmsr.query(query)
    
    #iterate over each day/plot
    for (day in dates) {
      res_day = subset(res, date==day)
      
      #date sub
      for (displot in plots) {
        res_day_plot = subset(res_day, plot==displot)
        
        #skip empty datas
        if (length(res_day_plot$wavelength) == 0) {
          next
        }
        
		#store the reflectance data in 'r_rem' variable anticipating that eventually this will store the "remove_bumple" data
		#r:	reflectance
		#rem: remove bumple function run
        r_rem = res_day_plot$avg_ref
        
		#upper and lower limits for the "remove_bumple" function since this was the region the noise was found in
		#so we only want to remove noise for this region
        u = 800
        l = 600
        
        il = which(res_day_plot$wavelength<u & res_day_plot$wavelength>l)
        il_w = res_day_plot$wavelength[il]
        
        removed = remove_bumple(res_day_plot$wavelength[il], res_day_plot$avg_ref[il], iterations=5, interpolate=FALSE)
        
        while (is.na(removed[length(removed)])) {
          u = u + 5
          il = which(res_day_plot$wavelength<u & res_day_plot$wavelength>l)
          il_w = res_day_plot$wavelength[il]
          
          removed = remove_bumple(res_day_plot$wavelength[il], res_day_plot$avg_ref[il], iterations=5, interpolate=FALSE)
        }
        
		#now, for the bumple region, have reflectnace with bumple removed
        r_rem[il] = remove_bumple(res_day_plot$wavelength[il], res_day_plot$avg_ref[il], iterations=5, interpolate=TRUE)
        
		#take derivative
		#d: derivative
        d_r_rem = dx(res_day_plot$wavelength, r_rem)
		
		#filter with 3rd order sgolay
		#f: filter
		#3: 3rd order
        d_r_rem_f3 = sgolayfilt(d_r_rem, p=3, n=41)
		
		#filter with 0th order sgolay (moving average)
		#0: 0th order
        d_r_rem_f30 = sgolayfilt(d_r_rem_f3, p=0, n=41)
        
        #find the REIP max
        reip_max_i = which.max(d_r_rem_f30[il])
        reip = il_w[reip_max_i]
        
        
        #find the indicies that correspond to the desired wavelengths
        idx1 = max(which(res_day_plot$wavelength<=y1))
        idx2 = max(which(res_day_plot$wavelength<=y2))
        idx3 = max(which(res_day_plot$wavelength<=y3))
        idx4 = max(which(res_day_plot$wavelength<=y4))
        
		#find the data points
        dr1c[idx] = d_r_rem_f30[idx1]
        dr2c[idx] = d_r_rem_f30[idx2]
        dr3c[idx] = d_r_rem_f30[idx3]
        dr4c[idx] = d_r_rem_f30[idx4]

        #store everything into nice variables for expression to use
        dr1 = res_day_plot$avg_ref[idx1]
        dr2 = res_day_plot$avg_ref[idx2]
        dr3 = res_day_plot$avg_ref[idx3]
		dr4 = res_day_plot$avg_ref[idx4]
        w1 = res_day_plot$wavelength[idx1]
        w2 = res_day_plot$wavelength[idx2]
        w3 = res_day_plot$wavelength[idx3]
		w4 = res_day_plot$wavelength[idx4]
		
		#now finally evaluate the mathametically expression passed in
		#at the begining of the function
        eqc[idx] = eval(parse(text=expression))
        
		#max rededge
        reipc[idx] = reip
		
		#save drs
		dr1c[idx] = dr1
		dr2c[idx] = dr2
		dr3c[idx] = dr3
		dr4c[idx] = dr4
		
        #information for analysis
        datec[idx] = day
        plotc[idx] = displot
        
        #index
        idx = idx+1
        
      }
    }
  }
  
  #return the dataframe
  foo = data.frame(date=datec, plot=plotc, reip=reipc, eq=eqc, dr1=dr1c, dr2=dr2c, dr3=dr3c, dr4=dr4c)
  return(foo)
  
}


# dx
# ================================
# Takes a discrete time 2 sided derivative of x with respect to t
# 
# Usage:
# t			= (list) data to be used in the calucations
# x			= (list) data to be used in the calucations
#
# Returns:
# dx			= (list) list of first derivative of x with respect to t
dx <- function(t, x) {
  
  #
  dx = vector(length =length(t))
  
  for (i in 1:length(t)) {
    
    #can't take 2 sided derivitive on the endpoints
    if (i==1) {
      #right sided derative
      dx[i] = (x[i+1]-x[i])/(t[i+1]-t[i])
    } else if (i==length(x)) {
      #left sided dsa;rafosi
      dx[i] = (x[i]-x[i-1])/(t[i]-t[i-1])
    } else {
      #balance deravitive 
      dx[i] = (x[i+1]-x[i-1]) / (t[i+1]-t[i-1])
      
    }
  }
  
  return(dx)  
}
