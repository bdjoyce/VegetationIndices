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

# ig_rededge_reflectance_equation
# ================================
# Script to calculate the inverted gaussian rededge position reflectance index
# this was a little too complex to fit inside an expression for use in the reflectance fucntions
# 
# NOTE:  script assumes dmsr library is available to run the queries
# NOTE:  script assumes database structure of at least ...
# +--------------------------------------------------------------------+
# | date | plot | wavelength | reflectance | serial_number_1 | exclude |
# +--------------------------------------------------------------------+
# NOTE:  the following analysis performed on the data
# 1:  Any plot/day measurement made with the same serial_number_1 multiple times is averaged
# 2:  Any plot/day measurement made with different serial_number_1 will be averaged after the median is calculated
#
# for basic structure of reflectance scripts see reflectance_equation.r or rededge_reflectance_equation.r
#

  library("optimx")
  library("signal")

  #returns a data frame
  table="2014_CSRV"
  
  serials = dmsr.query(paste('SELECT DISTINCT serial_number_1 FROM', table ,' WHERE exclude IS NULL', sep=" "))$serial_number_1
  ys = dmsr.query(paste('SELECT distinct(wavelength) FROM ', table, ' WHERE exclude IS NULL', sep=" "))$wavelength
  dates = dmsr.query(paste('SELECT distinct(date) FROM ', table, ' WHERE exclude IS NULL ORDER by date ASC', sep=" "))$date
  plots = dmsr.query(paste('SELECT distinct(plot) FROM ',table,' WHERE exclude IS NULL ORDER BY plot ASC', sep=" "))$plot
  
  iG_REc = c()
  R0c = c()
  Rsc = c()
  W0c = c()
  Sc = c()
  datec = c()
  plotc = c()
  eqc = c()
  idx = 1
    
  for (serial in serials) {
        
    #get all the wavelengths in the red edge region
	#we can filter later the regions we need to speed up analysis
    query = paste('SELECT date, plot, wavelength, AVG(reflectance) AS avg_ref FROM ', table, ' WHERE serial_number_1=\'', serial, '\' AND wavelength>500 AND wavelength<900 AND exclude IS NULL GROUP BY date, plot, wavelength order by date, plot, wavelength ASC;', sep='');
    res = dmsr.query(query)
    
    #iterate over each day/plot
    for (day in dates) {
      res_day = subset(res, date==day)
      
      #date subroutine
      for (displot in plots) {
        res_day_plot = subset(res_day, plot==displot)
        
        #skip emptys
        if (length(res_day_plot$wavelength) == 0) {
          next
        }
        
        rem = res_day_plot$avg_ref
        
		#upper and lower boung for red edge smoothing and analysis
        u = 810
        l = 660
        
        il = which(res_day_plot$wavelength<u & res_day_plot$wavelength>l)
        il_w = res_day_plot$wavelength[il]
        
		#if the upper frequency range seen is removed by the bumple smoothing
		#increase the frequency range because we need an upper point to use for
		#interpolation
        removed = remove_bumple(res_day_plot$wavelength[il], res_day_plot$avg_ref[il], iterations=5, interpolate=FALSE)
 
		while (is.na(removed[length(removed)])) {
          u = u + 5
          il = which(res_day_plot$wavelength<u & res_day_plot$wavelength>l)
          il_w = res_day_plot$wavelength[il]
 
          removed = remove_bumple(res_day_plot$wavelength[il], res_day_plot$avg_ref[il], iterations=5, interpolate=FALSE)
          
        }
		
		#now that the frequency range is OK remove_bumple with interpolation
        rem[il] = remove_bumple(res_day_plot$wavelength[il], res_day_plot$avg_ref[il], iterations=5, interpolate=TRUE)
        
        
        # use these to do the optimization
        r = rem[il]
        w = res_day_plot$wavelength[il]
        
        
        #find the initial values
        Rsr = max(r)
        R0r = min(r)
        W0r = w[which.min(r)]
        Sg = 30
        
        #R0 is non-interated
        R0 = mean(r[which(w>(W0r-10) & w<=W0r)])
        
        #plotting
        iG_start = Rsr - (Rsr-R0) * exp( -1*((W0r-w)**2) / (2*Sg**2)  )
        
        #test function
        test.f <- function(v,R0,w,r){
          Rs = v[1] 
          W0 = v[2]
          S = v[3]
          R0 = R0
          
          #create the iG vector
          iG = Rs - (Rs-R0) * exp( -1*((W0-w)**2) / (2*S**2)  )
          
          #now check how close to the actual data it is
          #via RMSE for all of the data points
          rmse = sum(sqrt((iG-r)**2))
          return(rmse)
        }
        
        #initial condition vector
        initial = c(Rsr, W0r, Sg)
        opti = optimx(initial, test.f, R0=R0, w=w, r=r)
                
        Rs = opti$p1[1]
        W0 = opti$p2[1]
        S = opti$p3[1]
        
        iG_RE = W0+S
        
        iG_REc[idx] = iG_RE
        Rsc[idx] = Rs
        W0c[idx] = W0
        Sc[idx] = S
        R0c[idx] = R0
        
        #info
        datec[idx] = day
        plotc[idx] = displot
        
        #index
        idx = idx+1
        
      }
    }
  }
  
  mREIP = data.frame(date=datec, plot=plotc, iG_RE=iG_REc, Rs=Rsc, W0=W0c, S=Sc, R0=R0c)
  



