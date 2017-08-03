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


# Dirty MySQL
# A better option would be to use RMySQL.  But since I can't get RMySQL 
# running, this is my workaround
#
# This tool works by invoking the mysql executable using the 
# R system2 commands.  Data is piped into files and then parsed
# into a R using read.delim.  THIS METHOD IS INHERENTLY INSECURE
# because the database password will be called from the command line
# and partial results will be stored with whatever permissions
# are given to R tempfiles.  This is truly a quick-and-dirty way
# to get MySQL data into R hence the name Dirty MySQL.
# 


# Setup environment to hold dmsr
dmsr <- new.env()
dmsr$initalized = FALSE


# dmsr.init
# ================================
# initializes the dmsr library with the information needed to execute queries
# provide the connection information
#
# Usage:
# mysql_path 	= (string) path to mysql on the machine
#  							eg.  C:\\Program Files\\MySQL\\MySQL Server 5.6\\bin\\mysql.exe
# hostname 		= (string) hostname of MySQL server
# username 		= (string) username to use when authenticating to the server
# passwd 		= (string) password to use when authenticating to the server
# database 		= (string) MySQL database to select from
# na.strings 	= (list) argument passed to read.delim as na.strings when parsing output
#  					eg.  c('NULL')  would make null retuns NA in R
#
#
# Example: 
#  dmsr.init(
#  		mysql_path 	= 'C:\\Program Files\\MySQL\\MySQL Server 5.6\\bin\\mysql.exe',
#  		hostname 	= 'localhost',
#  		username 	= 'root',
#  		passwd 		= 'unsecure',
#  		database 	= 'agr');
#
dmsr.init <- function(mysql_path, hostname, username, passwd, database, na.strings) {
  dmsr$initalized = TRUE
  dmsr$host = hostname
  dmsr$path = mysql_path
  dmsr$user = username
  dmsr$passwd = passwd
  dmsr$db = database
  if (missing(na.strings)) {
    dmsr$na.strings = c()
  } else {
    dmsr$na.strings = na.strings
  }
}



# dmsr.query
# ================================
# Dirty mysql query function, passes the query to the DB
# and returns the results as an R data frame
#
# Usage:
# query 			= (string) mysql query text
#
# Example:
# data = dmsr.query('SELECT wavelength, reflectance FROM 2014_CSRV WHERE date=20130605 AND plot=1010 ORDER by wavelength ASC')
#
dmsr.query <- function(query) {
  
  if (!dmsr$initalized) {
    stop('ERROR: dmsr must be initalized with dmsr.init() prior to running query')
  }
  
  
  args = c(
    paste('-h', dmsr$host, sep=''),
    paste('-u', dmsr$user, sep=''),
    paste('-p', dmsr$passwd, sep=''),
    paste('-D', dmsr$db, sep=''),
    '-e',
    paste('\"', query, '\"', sep="")
  )
  
  #swap filenames
  stdtemp = tempfile(pattern = "dirtyMySQLstd", fileext = ".txt")
  errtemp = tempfile(pattern = "dirtyMySQLerr", fileext = ".txt")
  
  
  tryCatch({
    system2(dmsr$path, args=args, stdout=stdtemp, stderr=errtemp)
  }, warning = function(w) {
    #warning-handler-code
    warning(w)
  }, error = function(e) {
    #error-handler-code
    dmsr$parse.stderr(errtemp)
    file.remove(stdtemp)
    file.remove(errtemp)
    stop(e)
  }, finally = {
    #cleanup-code
  })
  
  #get the results
  dmsr$parse.stderr(errtemp)
  results = dmsr$parse.stdout(stdtemp)
  
  #leave everything like we found it
  if (file.exists(stdtemp)) {
    file.remove(stdtemp)
  }
  if (file.exists(errtemp)) {
    file.remove(errtemp)
  }
  
  #be free
  return(results)
}

# internal function used by the dmsr library
dmsr$parse.stderr <- function(filename) {
  #no file = no errors
  if (!file.exists(filename)) {
    return(FALSE)
  }
  
  #file = something to think about perhaps
  lines = readLines(filename)
  for (i in 1:length(lines)) {
    write(lines[i], stderr())
  }
  return(TRUE)
}

# internal function used by the dmsr library
dmsr$parse.stdout <- function(filename) {
  #no file = no output
  if (!file.exists(filename)) {
    return(c())
  }
  results <- read.delim(filename, na.strings=dmsr$na.strings)
  return(results)
}
