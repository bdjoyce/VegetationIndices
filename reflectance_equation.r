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



# reflectance_equation
# ================================
# Function to generalize execution of reflectance equations to make large batch analysis easier
# allows a reflectance equation to be listed in expression and the function queries the database
# to find the data needed to evaluate the equation.  This allows up to 4 wavelengths to be used
# in the calculation.
#
# NOTE:  function assumes dmsr library is available to run the queries
# NOTE:  function assumes database structure of at least ...
# +--------------------------------------------------------------------+
# | date | plot | wavelength | reflectance | serial_number_1 | exclude |
# +--------------------------------------------------------------------+
# NOTE:  the following analysis performed on the data
# 1:  Any plot/day measurement made with the same serial_number_1 multiple times is averaged
# 2:  Any plot/day measurement made with different serial_number_1 will be averaged after the median is calcuated
#
# Usage:
# expression 	= (string) 	reflectance equation to evaluate
#
#								y1, y2, y3 correspond to wavelength and reflectance variables that can be used in the expression
#								w1, w2, w3 correspond to the exact wavelengths available in the dataset for the y1 y2 y3 values
#								r1, r2, r3 are the reflectance values that correspond to the wavelengths (y1 y2 y3)
#									the r's used are a list of the 20 closest reflectances to the wavelength
#									so using a math function like median or mean around r1 r2 r3 is necessary
#
#							Example expressions
#         						NDVI14: "(median(r1)-median(r2))/(median(r1)+median(r2))"
#        			 			Boochs2: "mean(dx(w1, r1))"
#
# y1, y2, y3	= (float) 	3 wavelengths to pull data for use in expression
# table			= (string)  table to use when executing queries
#
# Returns:
# result		= (data frame) a data frame containing each plot/day data along with with the wavelengths 
#								pulled, reflectances for those wavelengths, and the results of the 
#								calculated reflectance equation
#
# Example:
# Vogelmann3 = rededge_reflectance_equation("dr1/dr2", 715, 705, 705, 705, "2014_CSRV")
#
reflectance_equation <- function(expression, y1, y2, y3, table="2014_CSRV") {
  
  #run some setup queries to get information used in looping
  serials = dmsr.query(paste('SELECT DISTINCT serial_number_1 FROM', table ,' WHERE exclude IS NULL', sep=" "))$serial_number_1
  ys = dmsr.query(paste('SELECT distinct(wavelength) FROM ', table, 'WHERE exclude IS NULL', sep=" "))$wavelength
  dates = dmsr.query(paste('SELECT distinct(date) FROM ', table, ' WHERE exclude IS NULL ORDER by date ASC', sep=" "))$date
  plots = dmsr.query(paste('SELECT distinct(plot) FROM ',table,' WHERE exclude IS NULL ORDER BY plot ASC', sep=" "))$plot
  
  #initialize data structures so we can push data onto them later in the fuction
  r1c = c()
  r2c = c()
  r3c = c()
  datec = c()
  plotc = c()
  eqc = c()
  idx = 1
  
  for (serial in serials) {
    
    #find the ys we want for this wavelength = y1
    wavelength = y1
    for(i in 1:length(ys)){
      if (ys[i] >= wavelength) {
        break
      }
    }
    if (ys[i] == wavelength) {
      upper1 = ys[i+25]
    } else {
      upper1 = ys[i+25]
    }
    lower1 = ys[i-25]
    
    #find the ys we want for this wavelength = y2
    wavelength = y2
    for(i in 1:length(ys)){
      if (ys[i] >= wavelength) {
        break
      }
    }
    if (ys[i] == wavelength) {
      upper2 = ys[i+25]
    } else {
      upper2 = ys[i+25]
    }
    lower2 = ys[i-25]
    
    #find the ys we want for this wavelength = y3
    wavelength = y3
    for(i in 1:length(ys)){
      if (ys[i] >= wavelength) {
        break
      }
    }
    if (ys[i] == wavelength) {
      upper3 = ys[i+25]
    } else {
      upper3 = ys[i+25]
    }
    lower3 = ys[i-25]
    
    #Run queries to get the wavelengths and reflectances
    query1 = paste('SELECT date, plot, wavelength, AVG(reflectance) AS avg_ref FROM ', table, ' where wavelength >=', lower1, ' AND wavelength <=', upper1, ' AND serial_number_1=\'', serial, '\' AND exclude IS NULL GROUP BY date, plot, wavelength order by date, plot, wavelength ASC;', sep='');
    res1 = dmsr.query(query1)
    query2 = paste('SELECT date, plot, wavelength, AVG(reflectance) AS avg_ref FROM ', table, ' where wavelength >=', lower2, ' AND wavelength <=', upper2, ' AND serial_number_1=\'', serial, '\' AND exclude IS NULL GROUP BY date, plot, wavelength order by date, plot, wavelength ASC;', sep='');
    res2 = dmsr.query(query2)
    query3 = paste('SELECT date, plot, wavelength, AVG(reflectance) AS avg_ref FROM ', table, ' where wavelength >=', lower3, ' AND wavelength <=', upper3, ' AND serial_number_1=\'', serial, '\' AND exclude IS NULL GROUP BY date, plot, wavelength order by date, plot, wavelength ASC;', sep='');
    res3 = dmsr.query(query3)
    w = unique(res1$wavelength)
    
    #Now parse by date and plot 
    for (day in dates) {
      res_day1 = subset(res1, date==day)
      res_day2 = subset(res2, date==day)
      res_day3 = subset(res3, date==day)
      
      for (displot in plots) {
        #Make subsets
        res_day_plot1 = subset(res_day1, plot==displot)
        res_day_plot2 = subset(res_day2, plot==displot)
        res_day_plot3 = subset(res_day3, plot==displot)
        
		    #if there is no data no need to keep going
        if (length(res_day_plot1$wavelength) == 0 && length(res_day_plot2$wavelength) == 0 && length(res_day_plot3$wavelength) == 0) {
          next
        }
        
		    #Find the indexes that correspond to wavelength = y1
        #take the 20 closest wavelengths to y1 (10 above and 10 below)
        #unless there is an exact match in the database then take 19 closest (9 below, 9 ablove and matched)
        wavelength = y1
        ys1 = res_day_plot1$wavelength
        for(i in 1:length(ys1)){
          if (ys1[i] >= wavelength) {
            break
          }
        }
        if (ys1[i] == wavelength) {
          ui1 = i+10
        } else {
          ui1 = i+9
        }
        li1 = i-10
        
        #Find the indexes that correspond to wavelength = y2
        wavelength = y2
        ys2 = res_day_plot2$wavelength
        for(i in 1:length(ys2)){
          if (ys2[i] >= wavelength) {
            break
          }
        }
        if (ys2[i] == wavelength) {
          ui2 = i+10
        } else {
          ui2 = i+9
        }
        li2 = i-10
        
        #Find the indexes that correspond to wavelength = y3
        wavelength = y3
        ys3 = res_day_plot3$wavelength
        for(i in 1:length(ys3)){
          if (ys3[i] >= wavelength) {
            break
          }
        }
        if (ys3[i] == wavelength) {
          ui3 = i+10
        } else {
          ui3 = i+9
        }
        li3 = i-10
        
        
        #store everything into nice variables for expression to use
        r1 = res_day_plot1$avg_ref[li1:ui1]
        r2 = res_day_plot2$avg_ref[li2:ui2]
        r3 = res_day_plot3$avg_ref[li3:ui3]
        w1 = res_day_plot1$wavelength[li1:ui1]
        w2 = res_day_plot2$wavelength[li2:ui2]
        w3 = res_day_plot3$wavelength[li3:ui3]
        r1c[idx] = median(res_day_plot1$avg_ref)
        r2c[idx] = median(res_day_plot2$avg_ref)
        r3c[idx] = median(res_day_plot3$avg_ref)
		
		#now finally evaluate the mathametically expression passed in
		#at the begining of the function
        eqc[idx] = eval(parse(text=expression))
        
        #information for analysis
        datec[idx] = day
        plotc[idx] = displot
        
        #index
        idx = idx+1
      }
    }
    
  }
  
  r1 = r1c
  r2 = r2c
  r3 = r3c
  
  #return the dataframe
  result = data.frame(date=datec, plot=plotc, r1=r1c, r2=r2c, r3=r3c, eq=eqc)
  return(result)
  
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
