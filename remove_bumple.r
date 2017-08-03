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


# Functions use sgolayfilt from signal package in R to do analysis
library("signal")




# remove_bumple
# ================================
# Function designed to remove peak in reflectance spectra due to instrument error characterized by
# large shift up or down in reflectance across a narrow wavelength band that eventually settles
# after approximately 15nm
#
#
# Usage:
# w 			= (list) wavelength data to be used in the calculations
# r 			= (list) reflectance data to be used in the calculations (must match size of w[])
# width 		= (float) distance in units of w[] (ns for this paper) used to detect whether a bumple existed
# iterations 	= (int)   allows the function to run recursively in case multiple bumples exist in 1 dataset
#					 	  1 or 0 means no recursive functionality
# interpolate 	= (bool)  when true a high order polynomial is fit to the remaining, non-bumple data and 
#               		  used to fill in the points which were removed by this 
#
# Returns:
# r				= (list) reflectance data with the bumple removed (if found) and interpolated data (if interpolate=TRUE)
#
#
# Example: 
# 			u = 810
#			l = 660   
#	        il = which(results_frame$wavelength<u & results_frame$wavelength>l)
#			w_removed = remove_bumple(results_frame$wavelength[il], results_frame$avg_ref[il], iterations=5, interpolate=FALSE)
#
#
remove_bumple <- function(w, r, width=5, iterations=1, interpolate=FALSE) {  
  
  #find the first derivative
  dr = dx(w, r)
  
  #get the min and the  max
  maxi = which.max(dr)
  mini = which.min(dr)
  
  #if the min and the max are too close together then eliminate those points
  if (abs(w[maxi]-w[mini]) < width) {
    #found a bumple
    
	#eliminate the space between the min and max which identify the bumple
	dist = abs(maxi-mini)
    dr[(maxi-dist):(maxi+dist)] = NA
    dr[(mini-dist):(mini+dist)] = NA
    r[(maxi-dist):(maxi+dist)] = NA
    r[(mini-dist):(mini+dist)] = NA
	
	#eliminate data a bit outside of the min and max as well as these are also affected by the bumple
    r[mini:(mini+30)] = NA
    r[maxi:(maxi+30)] = NA
    r[mini:(mini-10)] = NA
    r[maxi:(maxi-10)] = NA
  }
  
  #do recursion if needed
  if (iterations > 1) {
    r = remove_bumple(w, r, width, iterations-1)
  }
  
  # do interpolations if needed
  if (interpolate) {
    fit = lm(formula = r ~ base_r_poly(w, 17))
    rp = predict(fit, data.frame(x=r))
    
    r[which(is.na(r))] = rp[which(is.na(r))]
  }
  
  return(r)
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
  
  #hold this for later
  dx = vector(length =length(t))
  
  for (i in 1:length(t)) {
    
    #can't take 2 sided derivative on the endpoints
    if (i==1) {
      #right sided derivative
      dx[i] = NA
    } else if (i==length(x)) {
      #left sided derivative
      dx[i] = NA
    } else {
      #balance derivative 
      dx[i] = (x[i+1]-x[i-1])/(t[i+1]-t[i-1])
    }
  }
  
  return(dx)
}


# Utilizes sgolayfilt function but doesn't spit an error if there is internal NAs
#
# https://cran.r-project.org/web/packages/signal/signal.pdf
# 
# All function arguments are passed directly to sgolayfilt, see sgolayfilt documentation for details
sgolayfilt_na <- function(x, p = 3, n = p + 3 - p%%2, m = 0, ts = 1) {
  last = 1
  ret = c()
  for (i in c(which(is.na(x)), length(x)+1)) {
    if (i-last >= n) {
      ret[last:(i-1)] = sgolayfilt(x[last:(i-1)], p=p, n=n, m=m, ts=ts)
    } else {
      ret[last:(i-1)] = NA
    }
    
    last = i+1
  }
  return(ret)
}





# Below is the base R poly function
# When sourcing the filter package or some other package in the analysis the base R poly function was overridden
# so this function was retyped under the base_r_poly function we can still call it
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/poly.html


#  =================================================================
#  File src/library/stats/R/contr.poly.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2015 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

base_r_poly = function (x, ..., degree = 1, coefs = NULL, raw = FALSE) 
{
  dots <- list(...)
  if (nd <- length(dots)) {
    if (nd == 1 && length(dots[[1L]]) == 1L) 
      degree <- dots[[1L]]
    else return(polym(x, ..., degree = degree, raw = raw))
  }
  if (is.matrix(x)) {
    m <- unclass(as.data.frame(cbind(x, ...)))
    return(do.call("polym", c(m, degree = degree, raw = raw)))
  }
  if (degree < 1) 
    stop("'degree' must be at least 1")
  if (anyNA(x)) 
    stop("missing values are not allowed in 'poly'")
  n <- degree + 1
  if (raw) {
    Z <- outer(x, 1L:degree, "^")
    colnames(Z) <- 1L:degree
    attr(Z, "degree") <- 1L:degree
    class(Z) <- c("poly", "matrix")
    return(Z)
  }
  if (is.null(coefs)) {
    if (degree >= length(unique(x))) 
      stop("'degree' must be less than number of unique points")
    xbar <- mean(x)
    x <- x - xbar
    X <- outer(x, seq_len(n) - 1, "^")
    QR <- qr(X)
    if (QR$rank < degree) 
      stop("'degree' must be less than number of unique points")
    z <- QR$qr
    z <- z * (row(z) == col(z))
    raw <- qr.qy(QR, z)
    norm2 <- colSums(raw^2)
    alpha <- (colSums(x * raw^2)/norm2 + xbar)[1L:degree]
    Z <- raw/rep(sqrt(norm2), each = length(x))
    colnames(Z) <- 1L:n - 1L
    Z <- Z[, -1, drop = FALSE]
    attr(Z, "degree") <- 1L:degree
    attr(Z, "coefs") <- list(alpha = alpha, norm2 = c(1, 
                                                      norm2))
    class(Z) <- c("poly", "matrix")
  }
  else {
    alpha <- coefs$alpha
    norm2 <- coefs$norm2
    Z <- matrix(, length(x), n)
    Z[, 1] <- 1
    Z[, 2] <- x - alpha[1L]
    if (degree > 1) 
      for (i in 2:degree) Z[, i + 1] <- (x - alpha[i]) * 
      Z[, i] - (norm2[i + 1]/norm2[i]) * Z[, i - 1]
    Z <- Z/rep(sqrt(norm2[-1L]), each = length(x))
    colnames(Z) <- 0:degree
    Z <- Z[, -1, drop = FALSE]
    attr(Z, "degree") <- 1L:degree
    attr(Z, "coefs") <- list(alpha = alpha, norm2 = norm2)
    class(Z) <- c("poly", "matrix")
  }
  Z
}