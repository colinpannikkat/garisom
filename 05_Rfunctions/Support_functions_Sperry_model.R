# Useful functions for supporting parameterization of the Sperry model
# Martin Venturas, 03/December/2020
# Modification by Antoine (17 feb 2021) in function Get_weather_NLDAS: pressure calculated from elevation expressed in Pa instead of kPa

# DESCRIPTION: This script contains several functions that facilitate parameterizing the
# Sperry et al. 2017 model. Details of each function are provided in the anotations.

# LIST OF FUNCTIONS:
# Calculate_Kmax
# USDA_soil_texture
# Root_beta
# Solar_noon_correction
# Calculate_D
# Calculate_leap_years
# thermaltime
# SSE_GS
# Calculate_growing_season
# Calculate_calendar_days_since
# Convert_DOY_to_Month_Day
# Convert_UTC
# Convert_to_DOY
# Get_weather_NLDAS
# Get_grid_NLDAS
# Divide_timesteps_weather_repetition
# Divide_timesteps_weather_interpolation

#----------------------------------------------------------------------------------------------------------
#### Calculate_Kmax ####
# Estimates the Kmax of a tree for running the Sperry model.
# This script calculates the Kmax of the tree when you know the tree leaf
# area to basal area (LABA), leaf specific conductance (LSC), and the
# percent resistance in the leaf at saturation (Rleaf).
# INPUTS
# LABA (m^2 m^-2), usually obtained from allometry
# LSC (mmol s^-1 m^-2 MPa^-1), where the m^2 is of leaf area
# Rleaf (%), assuming the resistance is at 0 MPa pressure
# OUTPUT
# Kmax (kg h^-1 MPa-1 m^2), where the m^2 is of basal area
Calculate_Kmax<-function(LABA,LSC,Rleaf){
  Kleaf=LABA*LSC; # Calculate the total canopy leaf conductance (mmol s^-1 m^-2 MPa^-1) in relation to BA
  # Unit conversions of Kleaf
  # 1 mmol H2O = 18.01528*10^-6 kg H2O
  # 1 s = 1/3600 h
  Kleaf=Kleaf*(18.01528*10^-6)*3600; # Kleaf (kg h^-1 m^-2 MPa^-1)
  # Now we can solve for tree Kmax knowing that K=1/Resistance
  # As the resistances are in series there is the following equation that has
  # to be satisfied
  # Rtree=Rroot+Rstem+Rleaf, thus: 1/Kmax=1/Kroot+1/Kstem+1/Kleaf
  Rtree=(1/Kleaf)*100/Rleaf; # For calcutating total tree resistance
  Kmax=1/Rtree; # Kmax of the tree (kg h^-1 MPa-1 m^2)
  return(Kmax)
}

#----------------------------------------------------------------------------------------------------------
#### USDA_soil_texture ####
# Script written to obtain the USDA soil texture classification from the percent
# clay, sand and silt.
# INPUTS
# Percent clay
# Pencent sand
# Percent silt
# OUTPUT
# Soil type (string)
USDA_soil_texture<-function(clay,sand,silt){
  if ((silt + 1.5*clay) < 15){
    texture="sand"
  } else if ((silt + 1.5*clay >= 15) && (silt + 2*clay < 30)){
    texture="loamy sand"
  } else if ((clay >= 7 && clay < 20) && (sand > 52) && ((silt + 2*clay) >= 30) || (clay < 7 && silt < 50 && (silt+2*clay)>=30)){
    texture="sandy loam"
  } else if ((clay >= 7 && clay < 27) && (silt >= 28 && silt < 50) && (sand <= 52)){
    texture="loam"
  } else if ((silt >= 50 && (clay >= 12 && clay < 27)) || ((silt >= 50 && silt < 80) && clay < 12)){
    texture="silt loam"
  } else if (silt >= 80 && clay < 12) {
    texture="silt"
  } else if ((clay >= 20 && clay < 35) && (silt < 28) && (sand > 45)) {
    texture="sandy clay loam"
  } else if ((clay >= 27 && clay < 40) && (sand > 20 && sand <= 45)){
    texture="clay loam"
  } else if ((clay >= 27 && clay < 40) && (sand  <= 20)){
    texture="silty clay loam"
  } else if (clay >= 35 && sand > 45){
    texture="sandy clay"
  } else if (clay >= 40 && silt >= 40){
    texture="silty clay"
  } else if (clay >= 40 && sand <= 45 && silt < 40){
    texture="clay"
  } else {
    texture="Soil Texture Error"
  }
  return(texture)
}

#----------------------------------------------------------------------------------------------------------
#### Root_beta ####
# Calculate rooting depth parameters
# INPUT
# Rdmax, maximum rooting depth (m). It can be an array with multiple values
# OUTPUT
# Beta coefficient. This determines the root biomass distribution in
# relation to depth.
Root_beta<-function(Rdmax){
  # Root biomass distribution is allocated based on the equation provided in
  # Love et al (2019):
  # M = 1 - Beta^d, where M is the fraction of biomass above depth d
  # expressed in cm. We find the Beta that provides an M of 0.995 for the
  # maximum rooting depth.
  d=Rdmax*100; # Converts bedrock rooting depth from m to cm
  Beta=(1-0.995)^(1/d); # Calculates Beta
  return(Beta)
}
# REFERECE
#  Love DM, Venturas MD, Sperry JS, Brooks PD, Pettit JL, Wang Y, Anderegg WLR, Tai X, Mackay DS (2019) Dependence of aspen stands on a subsurface water subsidy: implications for climate change impacts. Water Resources Research 55: 1833-1848

#----------------------------------------------------------------------------------------------------------
#### Solar_noon_correction
# Calculate the UTC time zone and solar noon correction
# INPUT
# Longitude (+ East; - West)
# UTC.weather, time zone of the weather file in UTC time zone (optional). E.g., -7
# OUTPUT
# Solar_noon_correction, for inputing in the Sperry model. It is mainly used to calculate
# diffuse versus direct radiation.
# UTC, UTC only based on longitude to coordinate sun hour if UTC.weather is not provided, esle UTC=UTC.weather
Solar_noon_correction<-function(longitude,UTC.weather){
  UTC<-trunc(longitude/15) #rounds to the nearest integer in the direction of 0
  if (is.na(UTC.weather)) {
    Solar_noon_correction<-(longitude/15)-UTC; # What fraction of an hour it is away from noon
  } else {
    Solar_noon_correction<-((longitude/15)-UTC)+(UTC-UTC.weather); # What fraction of an hour it is away from noon
    UTC<-UTC.weather
  }
  Output<-c(Solar_noon_correction,UTC)
  return(Output)
}

#----------------------------------------------------------------------------------------------------------
#### Calculate_D #####
# Function for calculating the vapor pressure defficit of the atmosphere
# from specific humidity (shum), air temperature (tas) and atmospheric
# pressure (pres)
# INPUTS
# Air temperature, tas (C)
# Specific humidity, shum (kg kg-1, ratio of mass of water vator to the
# total mass of moist air)
# Pressure, pres (Pa)
# OUTPUT
# Vapor pressure deficit, D (kPa)
# REFERENCE
# Campbell GA, Norman JM (1998) An Introduction to Environmental Biophysics. Second Edition. Springer-Verlag.
Calculate_D<-function(tas,shum,pres){
  q<-shum #specific humidity kg kg-1
  T<-tas # temperature (C)
  P<-pres/1000 # pressure (kPa)
  # Saturation water vapor pressure from temperature (es, in kPa) using the Tetens formula
  # from Campbell & Norman (1998) equation 3.8
  a<-0.611 #kPa
  b<-17.502
  c<-240.97 #Celsius
  es<-a*exp((b*T)/(T+c)) # saturation vapor pressure in kPa
  
  # Saturation water vapor mixing ratio ws
  ws<-0.622*es/P
  
  # Relative humidigy (hr) is the ratio of water vapor mixing ratio (which is
  # approximately equal to q) divided by ws. 
  hr<-q/ws
  
  # D = es(T)-ea = es(Ta)(1-hr), Campbell & Norman (1998), eq. 3.12
  D=es*(1-hr)
  D[D[]<0]<-0 # Any D value that is smaller than 0 due to rounding I change it to 0
  return(D)
}

#----------------------------------------------------------------------------------------
##### Calculate_leap_years ####
# Determine if a year is a leap or not a leap year using the division method
# INPUT
# year, array with n years
# OUTPUT
# out, dataframe with two columns: out$year, out$days.in.year (365 normal, 366 leap)
Calculate_leap_years<-function(year){
  nyr<-length(year)
  days.in.year<-rep(0,nyr)
  for (i in 1:length(year)){
    c1<-(year[i]/4)-floor(year[i]/4);
    if (c1!=0){ # not a leap year
      days.in.year[i]<-365
    } else {
      c2<-(year[i]/100)-floor(year[i]/100)
      if (c2!=0){ # this is a leap year)
        days.in.year[i]<-366
      } else {
        c3<-(year[i]/400)-floor(year[i]/400)
        if (c3==0){ # this is a leap year
          days.in.year[i]<-366
        } else{
          days.in.year[i]<-365
        }
      }
    }
  }
  out<-cbind(year,days.in.year)
  return(out)
}

#----------------------------------------------------------------------------------------
#### themaltime ####
# Fucntion for calculating thermal degree days for each day in the year
# Using the same algorith as for the Sperry et al. weather curator used 
# for Dave et al.(2019) Water Research Resources paper.
# INPUTS
# year, array with the year (yyyy)
# day, array with the day of year (1-365 or 1-366)
# temperature, array with the corresponding air temperature in Celsius (mean daily, hourly or any other format)
# It requires the full year data.
# OUTPUT
# sol, a dataframe with the first column being the year, second day of year
# and third the thermaltime above 5 C
thermaltime<-function(year,day,temperature){
  X<-as.data.frame(cbind(year,day,temperature))
  colnames(X)<-c({'year'},{'day'},{'temperature'})
  minyr<-min(year)
  maxyr<-max(year)
  tb<-5; # Critical temperature threshold in degrees Celsius
  ct<-1
  # Run throught years
  for (y in minyr:maxyr){
    ttime<-0; # set themal time to zero
    # Check if the year is normal or gap
    max_doy<-Calculate_leap_years(y)
    # For storing the outputs
    year1<-rep(y,max_doy[2])
    day1<-c(1:max_doy[2])
    tt1<-rep(0,max_doy[2])
    X1<-as.data.frame(cbind(year1,day1,tt1))
    for (i in 1:max_doy[2]) {# Run through all days of the year
      dif <- mean(X[X$year==y & X$day==i,3])-tb
      if (dif<0){
        dif<-0
      }
      ttime<-ttime+dif
      X1[i,3]<-ttime
    }
    if (ct==1){
      out<-X1
    } else {
      out<-rbind(out,X1)
    }
    ct<-ct+1
  }
  colnames(out)<-c({'year'},{'day'},{'ttime'})
  return(out)
}
#---------------------------------------------------------------------------------------------
#### SSE_GS ####
# Function for calculating the start day and end day of growing a
# season. It is based on the the sum of squared errors
# INPUTS
# X, an array with six values with the initial guesses
#       X[1]=start day
#       X[2]=ttime of the start day
#       X[3]=end day
#       X[4]=ttime of the end day
#       X[5]=slope of the first line
#       X[6]=slope of the third line
# day, array with the day of year
# ttime, array with the ttime for each day
# OUTPUT
# sseGS, squared mean standard error for the X guess
SSE_GS<-function(X,day,ttime){
  # X is the vector containing starting points and slopes (6 parameters)
  # X=[X0,Y0,X1,Y1,K0,K2]; X0=Start day; X1=End day;
  X0<-X[1]
  Y0<-X[2]
  X1<-X[3]
  Y1<-X[4]
  K0<-X[5]
  K2<-X[6]
  sseGS<-0 # Initialize the value of sum squared means
  dstart<-32 # day of year to start the calculations. This only considers from February 1st onwards
  for (i in dstart:length(day)){
    if (day[i]<= X0){# For the line prior to start of growing season start
      sseGS<-sseGS+(ttime[i]-(Y0-(X0-day[i])*K0))^2;
    } else if (day[i]> X0 & day[i]<= X1){# For main slope
      sseGS<-sseGS+(ttime[i]-(Y0+((Y1-Y0)/(X1-X0))*(day[i]-X0)))^2;
    } else {# Third line  day[i]> X1
      sseGS<-sseGS+(ttime[i]-(Y1+(day[i]-X1)*K2))^2;
    }
  }
  return(sseGS)
}
#----------------------------------------------------------------------------------------
#### Calculate_growing_season ####
# Function for calculating the start and end date of the growing season for
# each year in the weather files based on thermal time. Equal to the approach used
# in Sperry et al. 2019 PNAS. We fit a three lines to the sections of the cumulative
# thermal degree days above 5 C to identify the inflexion points as those where the
# three lines that intersect provide the best fit based on mean squared errors.
# INPUTS
# Year, array with year (yyyy)
# Day, array with day of year (1-365 or 1-366) 
# Tair, array with air temperature (C)
# OUTPUT
# GS, dataframe with three columns (GS$Year, GS$start_day, GS$end_day)
Calculate_growing_season<-function(Year,Day,Tair){
  start_yr<-min(Year)
  end_yr<-max(Year)
  nyr=end_yr-start_yr+1
  # Calculate the thermal time for all years at once
  ttime<-thermaltime(Year,Day,Tair)
  # Create dataframe to save the year, start and end day
  ttime$year<-as.integer(ttime$year)
  g1<-rep(0,nyr) 
  g2<-rep(0,nyr)
  g3<-rep(0,nyr)
  GS<-as.data.frame(cbind(g1,g2,g3))
  colnames(GS)<-c({'Year'},{'Start_day'},{'End_day'})
  yr<-start_yr
  for (i in 1:nyr){
    GS$Year[i]<-yr
    xday<-ttime[ttime$year==yr,2]
    yttime<-ttime[ttime$year==yr,3]
    # Initial guess for calculating GS days: X=[X0,Y0,X1,Y1,K0,K2], where
    # X0=Start day; Y0=Start day ttime; X1=End day; Y1=End day ttime;
    # K0=Slope of the first line; K2=Slope of the third line
    # As the sseGSval function is not linear, depending on the initial
    # guesses the fminsearch can run into problems. A good way to constrain
    # it initially is to get guesses within the range
    X<-c(120,yttime[120],270,yttime[270],yttime[120]/120,(yttime[365]-yttime[270])/(365-270))
    Fit<-nlminb( X , objective=SSE_GS, day = xday, ttime = yttime)
    GS$Start_day[i]<-round(Fit$par[1])
    GS$End_day[i]<-round(Fit$par[3])
    if (GS$End_day[i]>max(xday)){ # Just in case the fitting reaches a weird solution
      GS$End_day[i]<-max(xday)
    }
    yr<-yr+1
  }
  return(GS)
}
#----------------------------------------------------------------------------------------
#### Convert_calendar_days_since ####
# Obtains the Year, Day (DOY) and Hour values from days since year-01-01
# INPUTS
# time, array with the timestamp as days since yyyy-01-01
# calendar, string with the calendar used
# reference_year, year since days are counted
# OUTPUT
# out, matrix with three columns (Year, Day, Hour)
Convert_calendar_days_since<-function(time,calendar,reference_year){
  dst<-reference_year
  nt<-length(time)
  Year<-time # Generates the output array of the correct size
  Day<-time # Generates the output array of the correct size
  Hour<-time # Generates the output array of the correct size
  if (calendar=="noleap" | calendar=="365_day"){ # For no leap calendars
    for (i in 1:nt){
      Year[i]<-trunc(dst+time[i]/365)
      Day[i]<-(((dst+time[i]/365)-Year[i])*365)+1 # Day of year (1-365) with fractions for h
      Hour[i]<-(Day[i]-trunc(Day[i]))*24 # Hour with fractions if step is smaller than one hour
      Day[i]<-trunc(Day[i]) # Day of year (1-365) as an integer
    }
  } else if (calendar=="gregorian" | calendar=="standard" | calendar=="proleptic_gregorian"){ # For calendars with leap years
    yrs<-seq(from=dst,to=dst+400,by=1)
    DaysYear<-Calculate_leap_years(yrs) # I calculate it for 400 years from days since to make sure it covers my years
    # Calculate the bounds for days since per year
    day.limit<-0
    for (i in 2:length(DaysYear[,1])){ 
      day.limit[i]<-day.limit[i-1]+DaysYear[i-1,2]
    }
    for (i in 1:nt){
      for (j in 1:length(DaysYear[,1])){
        if (time[i]>=day.limit[j] & time[i]<day.limit[j+1]) {
          Year[i]<-DaysYear[j,1]
          Day[i]<-time[i]-day.limit[j]+1 # Day with fraction hour
          Hour[i]<-(Day[i]-trunc(Day[i]))*24 # Hour with fractions if step is smaller than one hour
          Day[i]<-trunc(Day[i]) # Day of year (1-366) as an integer
        }       
      }
      
    }
  }
  out<-cbind(Year,Day,Hour)
  return(out)
}
#----------------------------------------------------------------------------------------
##### Convert_DOY_to_Month_Day ######
# Function for converting day of year (DOY) to Month and Day
# INPUTS
# Year, array with the year
# DOY, array with the day of year
# leap, string with "yes" if the calendar has leap years or "no" if it does not
# OUTPUT
# Out, dataframe with four columns (Year, DOY, Month, Day)
Convert_DOY_to_Month_Day<-function(Year,DOY,leap){
  n<-length(DOY)
  Month<-rep(0,n)
  Day<-rep(0,n)
  if (leap=="no") { # calendars with no leap year
    MoDa<-c(31,28,31,30,31,30,31,31,30,31,30,31)
    for (i in 1:n){
      if (DOY[i]<=31){
        Month[i]<-1
        Day[i]<-DOY[i]
      } else {
        for (j in 1:11){
          if (DOY[i]>sum(MoDa[1:j]) & DOY[i]<=sum(MoDa[1:(j+1)])){
            Month[i]<-j+1
            Day[i]<-DOY[i]-sum(MoDa[1:j])
          }
        }
      }
    }
  } else { # calendars with leap years
    yrs<-unique(Year)
    DaysYear<-Calculate_leap_years(yrs)
    for (i in 1:n){
      if (DaysYear[DaysYear[,1]==Year[i],2]==365){
        MoDa<-c(31,28,31,30,31,30,31,31,30,31,30,31)
      } else {
        MoDa<-c(31,29,31,30,31,30,31,31,30,31,30,31)
      }
      if (DOY[i]<=31){
        Month[i]<-1
        Day[i]<-DOY[i]
      } else {
        for (j in 1:11){
          if (DOY[i]>sum(MoDa[1:j]) & DOY[i]<=sum(MoDa[1:(j+1)])){
            Month[i]<-j+1
            Day[i]<-DOY[i]-sum(MoDa[1:j])
          }
        }    
      }
    }
  }
  out<-as.data.frame(cbind(Year,DOY,Month,Day))
  return(out)
}

#----------------------------------------------------------------------------------------
#### Convert_UTC ####
# Function that converts UTC+0 time zone (Zulu time zone) to the UTC zone we want as + or - nhours, e.g., UTC-8
# This function does not take into consideration day time savings, i.e., it allways substracts or adds the same
# number of hours. It does deal with leap / non leap years.
# INPUTS
# time, dataframe with three columns time$year, time$DOY, time$hour
# UTC, the time zone to which you want to convert the time as + or - number of hours, e.g. Mountain Standard Time
# (MST) is -7, Pacific Standard Time (PST) is -8
# OUTPUT
# out, dataframe with the same three columns but with the timezone corrected.
# Adjust to UTC time zone
Convert_UTC<-function(time,UTC){
  # First add or substract the number of hours
  time$hour<-time$hour+UTC
  # Change year, DOY and hour as needed so hour (0-23)
  for (i in 1:length(time$hour)){
    if (time$hour[i]<0){ # Happens when UTC is negative
      time$hour[i]<-time$hour[i]+24
      time$DOY[i]<-time$DOY[i]-1 # The day before
      if (time$DOY[i]==0){ # Change year
        time$year[i]<-time$year[i]-1
        leap<-Calculate_leap_years(time$year[i])
        if (leap[1,2]==365){
          time$DOY[i]<-365
        } else {
          time$DOY[i]<-366
        }
      }
    } else if (time$hour[i]>23) {# Happens when UTC is positive
      time$hour[i]<-time$hour[i]-24
      time$DOY[i]<-time$DOY[i]+1 # The next day
      if (time$DOY[i]==366){ # Decide if we change year
        leap<-Calculate_leap_years(time$year[i])
        if (leap[1,2]==365){ # non leap year
          time$DOY[i]<-1 # Change to first day in the next year
          time$year[i]<-time$year[i]+1 # Advance year
        }
      } else if (time$DOY[i]==367){
        time$DOY[i]<-1 # Change to first day in the next year
        time$year[i]<-time$year[i]+1 # Advance year
      }
    }
  }
  return(time)
}

#----------------------------------------------------------------------------------------
#### Convert_to_DOY #####
# Convert year, month, day to day of year (DOY)
# INPUTS
# year, array with year
# month, array with month
# day, array with day
# OUTPUT
# out, dataframe with four columns out$year, out$month, out$day, out$DOY 
Convert_to_DOY<-function(year,month,day){
  minyr<-min(year)
  maxyr<-max(year)
  yr<-c(minyr:maxyr)
  leap<-Calculate_leap_years(yr)
  MoDa<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  MoDaleap<-c(31,29,31,30,31,30,31,31,30,31,30,31)
  DOY<-NA
  for (i in 1:length(year)){
    if (leap[leap[,1]==year[i],2]==365){ # Non leap years
      if (month[i]==1){ # January
        DOY[i]<-day[i]
      } else { # Other months
        DOY[i]<-sum(MoDa[1:(month[i]-1)])+day[i]
      }
    } else { # leap year
      if (month[i]==1){ # January
        DOY[i]<-day[i]
      } else { # Other months
        DOY[i]<-sum(MoDaleap[1:(month[i]-1)])+day[i]
      }
    }
  }
  out<-as.data.frame(cbind(year,month,day,DOY))
  return(out)
}

#----------------------------------------------------------------------------------------
##### Get_weather_NLDAS #####
# This function downloads NLDAS weather hourly weather files for a location and
# formats them for using them as an input in the Sperry model.
# INPUTS
# site, string array with the site name or identifiyer. It is used for saving files with that name
# lat, latitude of the site as degrees North (South is -)
# lon, longitude of the site as degrees East (West is -)
# elev, elevation above sea level (m)
# startTime, string with "yyyy-mm-dd-Thh" format (e.g. "1980-01-01T08" for strarting at 8 am) in Zulu time (UTC+0)
# endTime, string with "yyyy-mm-dd-Thh" format (e.g. "2019-12-31T07" for ending at 7 am) in Zulu time (UTC+0)
# timeZone, UTC for which we calculate local time (e.g. MST would be -7). No daytime savings.
# OUTPUT
# out, list with three elements out$metadata, out$site, out$weather
#      out$metadata, contains the metadata information of the files downloaded
#      out$siteinfo, contains a datafame with columns: site, lat, lon, elev, lat_grid, lon_grid, elev_grid
#      out$weather, contains a dataframe with columns: Year, Day, Hour, Solar_Wm2, Rain_mm, Wind_ms-1, Tair_C, Tsoil_C, D_kPa
# dowloaded weather files are saved in the work directory
Get_weather_NLDAS<-function(site,lat,lon,elev,startTime,endTime,timeZone){
  out<-list() # List where the output data will be included
  # Air Temperature
  varname <- "NLDAS:NLDAS_FORA0125_H.002:TMP2m"
  URL <- sprintf("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=%s&location=GEOM:POINT(%f,%s%f)&startDate=%s&endDate=%s&type=asc2", varname, lon, "%20", lat, startTime, endTime) # %20 is the URL version of a space
  #URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:TMP2m&location=GEOM:POINT(-108.18458,%2037.47872)&startDate=2014-01-01T07&endDate=2020-01-01T06&type=asc2"
  fn<-paste(site,"Tair.txt",sep='_')
  download.file(URL, destfile = fn, method="curl")
  a <- read.table(fn, sep="\t",  stringsAsFactors = FALSE)
  # First get the grid info
  # latitude grid
  eval<-0
  ct<-0
  while (ct<length(a[,1])){
    ct<-ct+1
    eval<-grep("grid_x",a[ct,1],fixed=TRUE) # to identify the row where this is written in the csv file
    if(length(eval)!=0){break} # ct will be the row with the txt
  }
  lat_grid<-strsplit(as.character(a[ct,1])," ",fixed=TRUE) # The time array is in one txt line
  lat_grid<-lat_grid[[1]][1]
  lat_grid<-strsplit(lat_grid,"=")
  lat_grid<-as.numeric(lat_grid[[1]][2])
  # longitude grid
  eval<-0
  ct<-0
  while (ct<length(a[,1])){
    ct<-ct+1
    eval<-grep("grid_y",a[ct,1],fixed=TRUE) # to identify the row where this is written in the csv file
    if(length(eval)!=0){break} # ct will be the row with the txt
  }
  lon_grid<-strsplit(as.character(a[ct,1])," ",fixed=TRUE) # The time array is in one txt line
  lon_grid<-lon_grid[[1]][1]
  lon_grid<-strsplit(lon_grid,"=")
  lon_grid<-as.numeric(lon_grid[[1]][2])
  # elevation grid
  eval<-0
  ct<-0
  while (ct<length(a[,1])){
    ct<-ct+1
    eval<-grep("elevation",a[ct,1],fixed=TRUE) # to identify the row where this is written in the csv file
    if(length(eval)!=0){break} # ct will be the row with the txt
  }
  elev_grid<-strsplit(as.character(a[ct,1]),"=",fixed=TRUE) # The time array is in one txt line
  elev_grid<-as.numeric(elev_grid[[1]][2])
  # Join
  siteinfo<-data.frame(site,lat,lon,elev,lat_grid,lon_grid,elev_grid,stringsAsFactors = FALSE)
  
  metadata<-a[1:35,1]
  a<-a[37:(nrow(a)-1),1]
  b<-strsplit(a," ")
  ymd<-NA
  year<-NA
  month<-NA
  day<-NA
  hour<-NA
  data_Tair<-NA
  for (i in 1:length(b)) {
    ymd[i] <- b[[i]][7]
    y1<-strsplit(ymd[i],"-")
    year[i]<-as.integer(y1[[1]][1])
    month[i]<-as.integer(y1[[1]][2])
    day[i]<-as.integer(y1[[1]][3])
    hour[i]<-as.integer(substr(b[[i]][8],1,2))
    data_Tair[i]<-as.double(b[[i]][13])
  }
  time<-Convert_to_DOY(year,month,day)
  time<-cbind(time[,c(1,4)],hour) # Stay with year, DOY and add h
  time<-Convert_UTC(time,timeZone)
  data_Tair<-data_Tair-273.15 # Convert from Kelvin to Celsius
  
  # Tsoil
  varname <- "NLDAS:NLDAS_NOAH0125_H.002:TSOIL0-10cm"
  URL <- sprintf("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=%s&location=GEOM:POINT(%f,%s%f)&startDate=%s&endDate=%s&type=asc2", varname, lon, "%20", lat, startTime, endTime) # %20 is the URL version of a space
  #URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:DSWRFsfc&location=GEOM:POINT(-108.18458,%2037.47872)&startDate=2014-01-01T07&endDate=2020-01-01T06&type=asc2"
  fn<-paste(site,"Tsoil.txt",sep='_')
  download.file(URL, destfile = fn, method="curl")
  a <- read.table(fn, sep="\t",  stringsAsFactors = FALSE)
  metadata<-rbind(metadata,a[1:35,1])
  a<-a[37:(nrow(a)-1),1]
  b<-strsplit(a," ")
  data_Tsoil<-NA
  for (i in 1:length(b)) {
    data_Tsoil[i]<-as.double(b[[i]][13])
  }
  data_Tsoil<-data_Tsoil-273.15 # Convert from Kelvin to Celsius
  
  # Solar
  varname <- "NLDAS:NLDAS_FORA0125_H.002:DSWRFsfc"
  URL <- sprintf("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=%s&location=GEOM:POINT(%f,%s%f)&startDate=%s&endDate=%s&type=asc2", varname, lon, "%20", lat, startTime, endTime) # %20 is the URL version of a space
  #URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:DSWRFsfc&location=GEOM:POINT(-108.18458,%2037.47872)&startDate=2014-01-01T07&endDate=2020-01-01T06&type=asc2"
  fn<-paste(site,"Solar.txt",sep='_')
  download.file(URL, destfile = fn, method="curl")
  a <- read.table(fn, sep="\t",  stringsAsFactors = FALSE)
  metadata<-rbind(metadata,a[1:35,1])
  a<-a[37:(nrow(a)-1),1]
  b<-strsplit(a," ")
  data_solar<-NA
  for (i in 1:length(b)) {
    data_solar[i]<-as.double(b[[i]][13])
  }
  
  # Precipitation
  varname <- "NLDAS:NLDAS_FORA0125_H.002:APCPsfc"
  URL <- sprintf("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=%s&location=GEOM:POINT(%f,%s%f)&startDate=%s&endDate=%s&type=asc2", varname, lon, "%20", lat, startTime, endTime) # %20 is the URL version of a space
  #URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:APCPsfc&location=GEOM:POINT(-108.18458,%2037.47872)&startDate=2014-01-01T07&endDate=2020-01-01T06&type=asc2"
  fn<-paste(site,"Precipitation.txt",sep='_')
  download.file(URL, destfile = fn, method="curl")
  a <- read.table(fn, sep="\t",  stringsAsFactors = FALSE)
  metadata<-rbind(metadata,a[1:35,1])
  a<-a[37:(nrow(a)-1),1]
  b<-strsplit(a," ")
  data_precip_kgm2 <- NA
  for (i in 1:length(b)) {
    data_precip_kgm2[i]<-as.double(b[[i]][13])
  }
  
  # Wind
  varname <- "NLDAS:NLDAS_FORA0125_H.002:UGRD10m"
  URL <- sprintf("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=%s&location=GEOM:POINT(%f,%s%f)&startDate=%s&endDate=%s&type=asc2", varname, lon, "%20", lat, startTime, endTime) # %20 is the URL version of a space
  #URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:UGRD10m&location=GEOM:POINT(-108.18458,%2037.47872)&startDate=2014-01-01T07&endDate=2020-01-01T06&type=asc2"
  fn<-paste(site,"Wind1.txt",sep='_')
  download.file(URL, destfile = fn, method="curl")
  a <- read.table(fn, sep="\t",  stringsAsFactors = FALSE)
  metadata<-rbind(metadata,a[1:35,1])
  a<-a[37:(nrow(a)-1),1]
  b<-strsplit(a," ")
  data_wind_x <- 0
  for (i in 1:length(b)) {
    data_wind_x[i]<-as.double(b[[i]][13])
  }
  
  varname <- "NLDAS:NLDAS_FORA0125_H.002:VGRD10m"
  URL <- sprintf("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=%s&location=GEOM:POINT(%f,%s%f)&startDate=%s&endDate=%s&type=asc2", varname, lon, "%20", lat, startTime, endTime) # %20 is the URL version of a space
  #URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:VGRD10m&location=GEOM:POINT(-108.18458,%2037.47872)&startDate=2014-01-01T07&endDate=2020-01-01T06&type=asc2"
  fn<-paste(site,"Wind2.txt",sep='_')
  download.file(URL, destfile = fn, method="curl")
  a <- read.table(fn, sep="\t",  stringsAsFactors = FALSE)
  metadata<-rbind(metadata,a[1:35,1])
  a<-a[37:(nrow(a)-1),1]
  b<-strsplit(a," ")
  data_wind_y <- 0
  for (i in 1:length(b)) {
    data_wind_y[i]<-as.double(b[[i]][13])
  }
  data_wind_x[is.na(data_wind_x)] <- 0
  data_wind_y[is.na(data_wind_y)] <- 0
  data_wind_diagonal <- (data_wind_x^2 + data_wind_y^2)^0.5
  
  # Specific humidity
  varname <- "NLDAS:NLDAS_FORA0125_H.002:SPFH2m"
  URL <- sprintf("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=%s&location=GEOM:POINT(%f,%s%f)&startDate=%s&endDate=%s&type=asc2", varname, lon, "%20", lat, startTime, endTime) # %20 is the URL version of a space
  #URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:SPFH2m&location=GEOM:POINT(-108.18458,%2037.47872)&startDate=2014-01-01T07&endDate=2020-01-01T06&type=asc2"
  fn<-paste(site,"Specific_humidity.txt",sep='_')
  download.file(URL, destfile = fn, method="curl")
  a <- read.table(fn, sep="\t",  stringsAsFactors = FALSE)
  metadata<-rbind(metadata,a[1:35,1])
  a<-a[37:(nrow(a)-1),1]
  b<-strsplit(a," ")
  data_shum <- NA
  for (i in 1:length(b)) {
    data_shum[i]<-as.double(b[[i]][13])
  }
  
  # # Pressure
  # varname <- "NLDAS:NLDAS_FORA0125_H.002:PRESsfc"
  # URL <- sprintf("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=%s&location=GEOM:POINT(%f,%s%f)&startDate=%s&endDate=%s&type=asc2", varname, lon, "%20", lat, startTime, endTime) # %20 is the URL version of a space
  # #URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:PRES2m&location=GEOM:POINT(-108.18458,%2037.47872)&startDate=2014-01-01T07&endDate=2020-01-01T06&type=asc2"
  # fn<-paste(site,"Pressure.txt",sep='_')
  # download.file(URL, destfile = fn, method="curl")
  # a <- read.table(fn, sep="\t",  stringsAsFactors = FALSE)  
  # metadata<-rbind(metadata,a[1:35,1])
  # b <- a[37:nrow(a),1]
  # c <- strsplit(b, "     ")
  # date <- 0
  # data_pres <- 0
  # for (x in 1:length(c)) {
  #   dateVal <- c[[x]]
  #   date[x] <- trimws(dateVal[2], which="left")
  #   data_pres[x] <- as.numeric(dateVal[3])
  # }
  
  pressure <- matrix(101325*exp(-elev/8200), nrow=length(data_shum))
  
  data_VPD <- Calculate_D(data_Tair,data_shum,pressure)
  weather <- cbind(time,data_solar,data_precip_kgm2,data_wind_diagonal,data_Tair,data_Tsoil,data_VPD)
  colnames(weather)<-c("Year", "Day", "Hour", "Solar_Wm2", "Rain_mm", "Wind_ms-1", "Tair_C", "Tsoil_C", "D_kPa")
  
  out$siteinfo<-siteinfo
  out$metadata<-metadata
  out$weather<-weather
  
  return(out)
}

#----------------------------------------------------------------------------------------
#### Get_grid_NLDAS ####
# Function that checks in what is the NLDAS grid for a site (latitude, longitude) 
# INPUTS
# lat, latitude of the site (degrees North)
# lon, longitude of the site (degrees East, negative for West)
# OUTPUTS
# gridID, dataframe with columns: gridID ("longrid_latgrid"), lat_grid, lon_grid, elev_grid, lat, lon
Get_grid_NLDAS<-function(lat,lon){
  out<-list() # List where the output data will be included
  # Download a file with two hours data to read the metadata
  startTime<-"1980-01-01T01"
  endTime<-"1980-01-01T02"
  # Air Temperature
  varname <- "NLDAS:NLDAS_FORA0125_H.002:TMP2m"
  URL <- sprintf("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=%s&location=GEOM:POINT(%f,%s%f)&startDate=%s&endDate=%s&type=asc2", varname, lon, "%20", lat, startTime, endTime) # %20 is the URL version of a space
  #URL = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi?variable=NLDAS:NLDAS_FORA0125_H.002:TMP2m&location=GEOM:POINT(-108.18458,%2037.47872)&startDate=2014-01-01T07&endDate=2020-01-01T06&type=asc2"
  fn<-"temp.txt"
  download.file(URL, destfile = fn, method="curl")
  a <- read.table(fn, sep="\t",  stringsAsFactors = FALSE)
  # First get the grid info
  # latitude grid
  eval<-0
  ct<-0
  while (ct<length(a[,1])){
    ct<-ct+1
    eval<-grep("grid_x",a[ct,1],fixed=TRUE) # to identify the row where this is written in the csv file
    if(length(eval)!=0){break} # ct will be the row with the txt
  }
  lat_grid<-strsplit(as.character(a[ct,1])," ",fixed=TRUE) # The time array is in one txt line
  lat_grid<-lat_grid[[1]][1]
  lat_grid<-strsplit(lat_grid,"=")
  lat_grid<-as.numeric(lat_grid[[1]][2])
  # longitude grid
  eval<-0
  ct<-0
  while (ct<length(a[,1])){
    ct<-ct+1
    eval<-grep("grid_y",a[ct,1],fixed=TRUE) # to identify the row where this is written in the csv file
    if(length(eval)!=0){break} # ct will be the row with the txt
  }
  lon_grid<-strsplit(as.character(a[ct,1])," ",fixed=TRUE) # The time array is in one txt line
  lon_grid<-lon_grid[[1]][1]
  lon_grid<-strsplit(lon_grid,"=")
  lon_grid<-as.numeric(lon_grid[[1]][2])
  # elevation grid
  eval<-0
  ct<-0
  while (ct<length(a[,1])){
    ct<-ct+1
    eval<-grep("elevation",a[ct,1],fixed=TRUE) # to identify the row where this is written in the csv file
    if(length(eval)!=0){break} # ct will be the row with the txt
  }
  elev_grid<-strsplit(as.character(a[ct,1]),"=",fixed=TRUE) # The time array is in one txt line
  elev_grid<-as.numeric(elev_grid[[1]][2])
  # Join
  gridID<-paste(as.character(lon_grid),as.character(lat_grid),sep='_')
  out<-data.frame(gridID,lat_grid,lon_grid,elev_grid,lat,lon,stringsAsFactors = FALSE)
  return(out)
}

#----------------------------------------------------------------------------------------
#### Divide_timesteps_weather_repetition ####
# Function that downscales the weather timesteps to smaller intervals by repeting the same value
# for each timespet (except for rain (mm) which is kept in one single timestep)
# INPUT
# weather, matrix or dataframe with the structure that the Sperry model reads.
#          Columns: 1.Year, 2.Day, 3.Hour, 4.Solar_Wm2, 5.Rain_mm, 6.Wind_ms-1, 7.Tair_C, 8.Tsoil_C, 9.D_kPa
# n_div, number (integer) determining into how many timesteps we want to divide the original timestep
# n_hours, original timestep length in hours
# OUTPUT
# weather, dataframe with the same structure as the original weather file but with smaller timesteps
Divide_timestep_weather_repetition<-function(weather,n_div,n_hours){
  T<-as.matrix(weather)
  # We just repeat the hourly weather variable for smaller timesteps
  n1<-length(T[,1]) # Original length of the dataset
  n2<-n_div # How many intervals we divide each original timestep. Has to be integer
  n3<-n1*n2 # Final length of the datafile
  T2<-matrix(NA,nrow=n3,ncol=length(T[1,]))
  ct<-1
  T2[,5]<-0 # Set rain to zero
  for (i in 1:n1){
    T2[ct,5]<-T[i,5] # Rain only in one of the steps
    for(j in 0:(n2-1)){
      T2[ct,1]<-T[i,1]
      T2[ct,2]<-T[i,2]
      T2[ct,3]<-T[i,3]+j*n_hours/n2 # Add time difference
      T2[ct,4]<-T[i,4]
      T2[ct,6]<-T[i,6]
      T2[ct,7]<-T[i,7]
      T2[ct,8]<-T[i,8]
      T2[ct,9]<-T[i,9]
      ct<-ct+1
    }
  }
  T2<-as.data.frame(T2)
  colnames(T2)<-colnames(weather)
  return(T2)
}

#----------------------------------------------------------------------------------------
#### Divide_timesteps_weather_interpolation ####
# Function that downscales the weather timesteps to smaller intervals by interpolation between timesteps,
# except for rain (mm) which is kept in one single timestep or split for the timesteps
# INPUT
# weather, matrix or dataframe with the structure that the Sperry model reads.
#          Columns: 1.Year, 2.Day, 3.Hour, 4.Solar_Wm2, 5.Rain_mm, 6.Wind_ms-1, 7.Tair_C, 8.Tsoil_C, 9.D_kPa + optional: 10.Pressure
# n_div, number (integer) determining into how many timesteps we want to divide the original timestep
# n_hours, original timestep length in hours
# divide_rain, string. Set to "yes" if the total rain of one timestep is divided into fraction for the smaller timesteps.
#              Set to "no" for just keeping the total rain in one single timestep
# OUTPUT
# weather, dataframe with the same structure as the original weather file but with smaller timesteps
Divide_timestep_weather_interpolation<-function(weather,n_div,n_hours,divide_rain){
  M1<-as.matrix(weather)
  n<-length(M1[,1])
  M2<-matrix(0,nrow=n*n_div,ncol=length(M1[1,])) # Matrix where smaller timesteps are stored
  for (i in 1:n){
    j<-(i*n_div)-(n_div-1)
    if (i<n) {
      M2[j,]<-M1[i,] # The same values for the same hour. Does not split rain
      if (divide_rain=="yes"){
        M2[j,5]<-M1[i,5]/n_div
      }
      for (s in 1:(n_div-1)){ # for the remaining subdivision hour steps
        # Year
        M2[j+s,1]<-M1[i,1] # Initially equal. Later we check if it has to be adjusted
        # Day
        M2[j+s,2]<-M1[i,2] # Initially equal. Later we check if it has to be adjusted
        # Hour
        M2[j+s,3]<-M1[i,3]+s*n_hours/n_div # Add s*n_hours/n_div to the timestep
        # Now interpolate for solar (column 4)
        M2[j+s,4]<-M1[i,4]+s*((M1[i+1,4]-M1[i,4])/n_div) # solar
        # Rain
        if (divide_rain=="yes"){
          M2[j+s,5]<-M1[i,5]/n_div
        } else {
          M2[j+s,5]<-0 # All the rain is stored in the first timestep
        }
        # Now interpolate the wind (column 6)
        M2[j+s,6]<-M1[i,6]+s*((M1[i+1,6]-M1[i,6])/n_div) # wind 
        # Now interpolate the Tair (column 7)
        M2[j+s,7]<-M1[i,7]+s*((M1[i+1,7]-M1[i,7])/n_div); # Tair
        # Assing Tsoil equal to Tair
        M2[j+s,8]<-M1[i,8]+s*((M1[i+1,8]-M1[i,8])/n_div); # Tsoil
        # Interpolate D (column 9)
        M2[j+s,9]<-M1[i,9]+s*((M1[i+1,9]-M1[i,9])/n_div) # D
        # Interpolate Pressure (column 10), optional
        if (length(M1[1,])==10){
                  M2[j+s,10]<-M1[i,10]+s*((M1[i+1,10]-M1[i,10])/n_div) # Pressure
        }
      } 
    } else { # For the last values of the matrix
      M2[j,]<-M1[i,] # The same values
      for (s in 1:(n_div-1)){
        M2[j+s,c(1:2,4,6:9)]<-M1[i,c(1:2,4,6:9)] # Same value
        M2[j+s,3]<-M1[i,3]+s*n_hours/n_div # Add time to the timestep
        if (divide_rain=="yes"){
          M2[j+s,5]<-M1[i,5]/n_div
        } else {
          M2[j+s,5]<-0 # All the rain is stored in the first timestep
        }
        if (length(M1[1,])==10){
          M2[j+s,10]<-M1[i,10] # Pressure
        }
      }
    }
  }
  colnames(M2)<-colnames(M1) # For renaming the columns
  # Depending on how the hours are reported (e.g., if original 3h steps go as 1.5, 4.5,...,22.5) you can have instances where we change days
  # Check if that is a possibility
  if (max(M2[,3])>=24){
    print("Adjusted days and years during interpolation of timestep")
    M2[M2[,3]>=24,2]<-M2[M2[,3]>=24,2]+1 # Change day to following day
    M2[M2[,3]>=24,3]<-M2[M2[,3]>=24,3]-24 # Change hour to hour the following day
    # For last hours of a leap year
    M2[M2[,2]>366,1]<-M2[M2[,2]>366,1]+1 # Pass to next year as this only happens at the end of leap years
    M2[M2[,2]>366,2]<-1 # Set the day to the first day of the year
    # For non-leap years
    leap<-Calculate_leap_years(unique(M2[,1]))
    for (i in 1:length(M2[,1])){
      if(M2[i,2]==366 & leap[leap[,1]==M2[i,1],2]==365){# non-leap year
        M2[i,2]<-1 # Pass to first day of next year
        M2[i,1]<-M2[i,1]+1 # Pas to next year
      }
    }
  }
  M2<-as.data.frame(M2)
  return(M2)
}
