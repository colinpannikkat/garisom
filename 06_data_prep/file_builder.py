"""
Author: Colin Pannikkat
Date: 02/12/25
Description: This program is a generator for configuration, parameter, dataset, 
             and growing season files used as input to the stomatal optimzation 
             model garisomv3.
"""

# For data
import pandas as pd
from pandas import DataFrame

# For program arguments
import argparse
import os

# For calculations, fitting growing season and visualization
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt

# For fetching NLDAS data and timezone conversion
import requests
import io
from datetime import datetime, timedelta
from timezonefinder import TimezoneFinder
import pytz
from tqdm import tqdm

# Types
from typing import Callable


def _calc_thermal_time(climate_data: DataFrame) -> DataFrame:
    '''
    Calculate the thermal time for every day in the climate data.

    Inputs:
    - 'climate_data': DataFrame containing the climate data with (atleast) 
                      columns 'Year', 'Day', and 'Tair_C'.
    '''
    x = climate_data[['Year', 'Day', 'Tair_C']]
    x = x.dropna()
    minyr = int(min(x['Year']))
    maxyr = int(max(x['Year']))
    tb = 5 # critical temperature threshold in degrees celsius

    # days_in_years = _calc_leap_years(x['Year'])
    ttime_dt = DataFrame(columns=('Year', 'Day', 'ThermalTime'))
    first = True

    for y in range(minyr, maxyr + 1):
        ttime = 0
        # max_day = days_in_years.loc['Year' == y]
        for _, group in x.loc[x['Year'] == y].groupby("Day"):
            day = int(group['Day'].values[0])
            dif = group['Tair_C'].mean() - tb
            if (dif < 0):
                dif = 0
            ttime += dif
            if (first):
                ttime_dt = DataFrame([[y, day, ttime]], columns=ttime_dt.columns)
                first = False
            else:
                ttime_dt = pd.concat([DataFrame([[y, day, ttime]], columns=ttime_dt.columns), ttime_dt], ignore_index=True)

    return ttime_dt.sort_values(by=['Year', 'Day']).reset_index(drop=True)

def _sigmoid(x, L, x0, k, b):
    '''
    Sigmoid function for fitting to thermal-time data.

    Inputs:
    - `x`: x-input
    - `L`: The curve's maximum value
    - `x0`: The x-value of the sigmoid's midpoint
    - `k`: The steepness of the curve
    - `b`: The curve's minimum value
    '''
    return L / (1 + np.exp(-k * (x - x0))) + b

def _fit_sigmoid(x, y):
    '''
    Function to fit data to a sigmoid function.

    Inputs:
    - `x`: Day of the year
    - `y`: Thermal-time/degree-days

    Outputs:
    - `param`: Fitted sigmoid curve parameters
    '''
    p0 = [max(y), np.median(x), 1, min(y)]
    param, _ = curve_fit(_sigmoid, x, y, p0, maxfev=10000)
    return param

def _fit_three_segments(x, y):
    '''
    Fits three lines with lowest MSE to sigmoidal plot of degree-days vs.
    day of the year.

    Inputs:
    - `x`: day of the year
    - `y`: fitted thermal-time/degree-days

    Outputs:
    - `best_fit`: Best-fit line parameters
    '''
    min_mse = float('inf')
    best_fit = None
    
    # Iterates through all possible start and end growing season days
    # Assumes growing season won't start or end halfway through the year.
    for i in range(1, len(x) // 2):
        for j in range(len(x) // 2, len(x)):
            if i >= j:
                continue
            
            x1, y1 = x[:i], y[:i]
            x2, y2 = x[i:j], y[i:j]
            x3, y3 = x[j:], y[j:]
            
            reg1 = LinearRegression().fit(x1.reshape(-1, 1), y1)
            reg2 = LinearRegression().fit(x2.reshape(-1, 1), y2)
            reg3 = LinearRegression().fit(x3.reshape(-1, 1), y3)
            
            y_pred1 = reg1.predict(x1.reshape(-1, 1))
            y_pred2 = reg2.predict(x2.reshape(-1, 1))
            y_pred3 = reg3.predict(x3.reshape(-1, 1))
            
            mse = (np.mean((y1 - y_pred1) ** 2) + 
                   np.mean((y2 - y_pred2) ** 2) + 
                   np.mean((y3 - y_pred3) ** 2))
            
            if mse < min_mse:
                min_mse = mse
                best_fit = (reg1, reg2, reg3)
    
    return best_fit

def _find_intersection(reg1, reg2):
    '''
    Given two lines, finds the intersection between them.

    Inputs:
    - `reg1`: Line parameters for the first line
    - `reg2`: Line parameters for the second line
    '''
    m1, b1 = reg1.coef_[0], reg1.intercept_
    m2, b2 = reg2.coef_[0], reg2.intercept_
    x_intersect = (b2 - b1) / (m1 - m2)
    return x_intersect

def _fit_lines(xday, yt_time, show_plt=False):
    '''
    Fits the thermal-time data to a sigmoid curve, and then finds the optimal
    segmentation lines to determine the start and end of the growing season
    data.
    '''
    # Only fit from day 32 onwards
    x = np.array(xday[31:])
    y = np.array(yt_time[31:])
    
    param = _fit_sigmoid(x, y)
    y_fit = _sigmoid(x, *param)
        
    best_fit = _fit_three_segments(x, y_fit)
    reg1, reg2, reg3 = best_fit
    
    gs_start = round(_find_intersection(reg1, reg2))
    gs_end = round(_find_intersection(reg2, reg3))
    
    if show_plt:
        plt.plot(xday, yt_time, label='Cumulative Degree Days', color='black')
        plt.plot(x, y_fit, label='Fitted Sigmoid', linestyle='dashed', color='gray')
        plt.plot(xday[:gs_start], reg1.predict(np.array(xday[:gs_start]).reshape(-1, 1)), label='Segment 1', color='blue')
        plt.plot(xday[gs_start:gs_end], reg2.predict(np.array(xday[gs_start:gs_end]).reshape(-1, 1)), label='Segment 2', color='green')
        plt.plot(xday[gs_end:], reg3.predict(np.array(xday[gs_end:]).reshape(-1, 1)), label='Segment 3', color='orange')
        plt.axvline(gs_start, color='blue', linestyle='dashed', label='GS Start')
        plt.axvline(gs_end, color='red', linestyle='dashed', label='GS End')
        plt.legend()
        plt.xlabel('Day of Year')
        plt.ylabel('Cumulative Degree Days')
        plt.show()
    
    return gs_start, gs_end

def _calculate_growing_season(climate_data: DataFrame) -> DataFrame:
    """
    Function for calculating the start and end date of the growing season for
    each year in the weather files based on thermal time. Equal to the approach 
    used in Sperry et al. 2019 PNAS. We fit three lines to the sections of the 
    cumulative thermal degree days above 5°C to identify the inflexion points as 
    those where the three lines that intersect provide the best fit based on mean 
    squared errors.
    
    Input:
    - `climate_data`: DataFrame containing the climate data with (atleast) 
                      columns 'Year', 'Day', and 'Tair_C'.

    Output:
    - `DataFrame`: DataFrame with three columns ('Year', 'start_day', 'end_day') 
                   indicating the start and end day of the growing season for 
                   each year.
    """
    start = int(climate_data['Year'].min())
    end = int(climate_data['Year'].max())
    num_years = end - start + 1

    thermal_time = _calc_thermal_time(climate_data)
    thermal_time['Year'] = thermal_time['Year'].astype(int)
    gs = DataFrame(columns=['Year', 'Start_day', 'End_day'])

    yr = start
    for i in range(1, num_years+1):
        xday = thermal_time.loc[thermal_time['Year'] == yr]['Day'].to_list()
        yt_time = thermal_time.loc[thermal_time['Year'] == yr]['ThermalTime'].to_list()

        start_day, end_day = _fit_lines(xday, yt_time, show_plt=False)

        if (end_day > max(xday)): # if end_day reaches the max_day
            end_day = max(xday)
        if (i == 1):
            gs = DataFrame([[yr, start_day, end_day]], columns=gs.columns)
        else:
            gs = pd.concat([DataFrame([[yr, start_day, end_day]], columns=gs.columns), gs], ignore_index=True)

        yr += 1
    
    return gs.sort_values(by='Year').reset_index(drop=True)


def build_growing_season_data(noaa_data: DataFrame, 
                              climate_data: tuple[DataFrame, str], 
                              out_data_path: str, 
                              initial_year: int, 
                              final_year: int) -> None:
    '''
    This function takes in NOAA atmospheric CO2 trend data, climate_data,
    and builds the growing season file "seasonlimits.csv" and the limited
    climate data file "dataset.csv", outputting to `out_data_path/site/`.

    - `noaa_data`: DataFrame with NOAA annual atmospheric Carbon Dioxide 
                   concentration records
    - `climate_data`: Tuple with DataFrame of climate data (Year, Day, Tair) and
                      site.
    - `out_data_path`: Output data path
    - `initial_year`: Initial year to begin pulling CO2 data from
    - `final_year`: Final year to pull CO2 data
    '''
    noaa_data = noaa_data.loc[(noaa_data['Year'] >= initial_year) & (noaa_data['Year'] <= final_year)]
    for site_data, site in climate_data:
        site_data = site_data.loc[(site_data['Year'] >= initial_year) & (site_data['Year'] <= final_year)]
        gs_data = _calculate_growing_season(site_data)
        gs_data['ca_ppm'] = noaa_data.set_index('Year').loc[gs_data['Year'], 'ca_ppm'].values
        gs_data.to_csv(f"{out_data_path}/{site}/seasonlimits.csv", index=False)
        print(f"{site} growing season data saved to {out_data_path}/{site}/seasonlimits.csv.")
        
        # Add empty last row to site_data
        site_data = pd.concat([site_data, pd.DataFrame([[""] * len(site_data.columns)], columns=site_data.columns)], ignore_index=True)
        site_data.to_csv(f"{out_data_path}/{site}/dataset.csv", index=False, columns=['Year', 'Day', 'Hour', 'Solar_Wm2', 'Rain_mm', 'Wind_ms.1', 'Tair_C', 'Tsoil_C', 'D_kPa'])
        print(f"{site} hourly data saved to {out_data_path}/{site}/dataset.csv")

def _convert_time_to_utc(start_date, end_date, latitude, longitude):
    '''
    Converts the start_date and end_date to UTC time depending on local time
    determined by latitude and longitude of site.

    Inputs:
    - `start_date`: Start date in the format 'YYYY-MM-DDT%H'
    - `end_date`: End date in the format 'YYYY-MM-DDT%H'
    - `latitude`: Latitude of the location
    - `longitude`: Longitude of the location

    Outputs:
    - `start_date_utc`: Start date in UTC time in the format 'YYYY-MM-DDT%H'
    - `end_date_utc`: End date in UTC time in the format 'YYYY-MM-DDT%H'
    '''
    tf = TimezoneFinder()
    timezone_str = tf.timezone_at(lat=latitude, lng=longitude)
    if timezone_str is None:
        raise ValueError("Could not determine the timezone for the given coordinates.")
    
    local_tz = pytz.timezone(timezone_str)
    start_date_local = local_tz.localize(datetime.strptime(start_date, "%Y-%m-%dT%H"))
    end_date_local = local_tz.localize(datetime.strptime(end_date, "%Y-%m-%dT%H"))
    
    start_date_utc = start_date_local.astimezone(pytz.utc)
    end_date_utc = end_date_local.astimezone(pytz.utc)
    
    return start_date_utc.strftime("%Y-%m-%dT%H"), end_date_utc.strftime("%Y-%m-%dT%H")

''''
The following two functions are used/modified from:
https://github.com/nasa/gesdisc-tutorials/blob/main/notebooks/How_to_Access_the_Hydrology_Data_Rods_API_Using_Python.ipynb
'''

def _get_time_series(start_date, end_date, latitude, longitude, variable):
    """
    Calls the data rods service to get time series data.

    Inputs:
    - `start_date`: Start date in the format 'YYYY-MM-DDT00'
    - `end_date`: End date in the format 'YYYY-MM-DDT00'
    - `latitude`: Latitude of the location
    - `longitude`: Longitude of the location
    - `variable`: Variable to retrieve from the API

    Outputs:
    - `r.text`: Response text containing the time series data
    """
    base_url = "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi"
    query_parameters = {
        "variable": variable,
        "startDate": start_date,
        "endDate": end_date,
        "location": f"GEOM:POINT({longitude}, {latitude})",
        "type": "asc2"
    }

    full_url = base_url+"?"+ \
         "&".join(["{}={}".format(key,query_parameters[key]) for key in query_parameters])

    session = requests.Session()  # Reuse session for efficiency
    retries = 5
    for attempt in range(retries):
        try:
            r = session.get(full_url, timeout=10)  # Set timeout to prevent hanging
            if r.status_code == 200:
                return r.text
        except requests.RequestException as e:
            print(f"Attempt {attempt+1} failed: {e}")

    try:
        raise Exception(f"Error code {r.status_code} from url {full_url} : {r.text}")
    except:
        raise Exception(f"Could not connect to {full_url}. Please try again later.")

def _parse_time_series(ts_str, conv_func: Callable = None):
    """
    Parses the response from data rods into a series.

    Inputs:
    - `ts_str`: String from NLDAS API with data
    - `conv_func`: A conversion function if needed

    Outputs:
    - `parameters`: All parameters involved with data (header)
    - `df`: DataFrame with data
    """
    lines = ts_str.split("\n")
    parameters = {}
    for line in lines[2:11]:
        key,value = line.split("=")
        parameters[key] = value
    
    df = pd.read_table(io.StringIO(ts_str),
                       sep="\t",
                       names=["time","data"],
                       header=10,
                       parse_dates=["time"])
    
    if conv_func:
        df['data'] = df['data'].apply(conv_func)

    return parameters, df

def _calc_wind_from_nldas(start_date, end_date, latitude, longitude):
    '''
    Calculates wind speed from NLDAS zonal and meridional wind.
    '''
    x_wind = _parse_time_series(
                _get_time_series(
                        start_date=start_date, 
                        end_date=end_date,
                        latitude=latitude,
                        longitude=longitude,
                        variable="NLDAS2:NLDAS_FORA0125_H_v2.0:Wind_E"
                )
            )[1]['data'].fillna(0)
    
    y_wind = _parse_time_series(
                _get_time_series(
                        start_date=start_date, 
                        end_date=end_date,
                        latitude=latitude,
                        longitude=longitude,
                        variable="NLDAS2:NLDAS_FORA0125_H_v2.0:Wind_N"
                )
            )[1]['data'].fillna(0)
    
    return np.sqrt(np.pow(x_wind, 2) + np.pow(y_wind, 2))

def _calc_vpd_from_nldas(start_date, end_date, latitude, longitude):
    '''
    Calculates VPD from NLDAS air temp, specific humidity, and air pressure.
    
    Modified from `Support_functions_Sperry_model.R`
    '''

    # Get temperature of air
    tair = _parse_time_series(
        _get_time_series(
                start_date=start_date, 
                end_date=end_date,
                latitude=latitude,
                longitude=longitude,
                variable="NLDAS2:NLDAS_FORA0125_H_v2.0:Tair"
        ), lambda x: x - 273.15 # convert K to C
    )[1]['data']

    # Get specific humidity
    shum = _parse_time_series(
        _get_time_series(
                start_date=start_date, 
                end_date=end_date,
                latitude=latitude,
                longitude=longitude,
                variable="NLDAS2:NLDAS_FORA0125_H_v2.0:Qair"
        )
    )[1]['data']

    # Get pressure data
    # maybe swap out with pressure <- matrix(101325*exp(-elev/8200), nrow=length(data_shum))?
    pair = _parse_time_series(
        _get_time_series(
                start_date=start_date, 
                end_date=end_date,
                latitude=latitude,
                longitude=longitude,
                variable="NLDAS2:NLDAS_FORA0125_H_v2.0:PSurf"
        ), lambda x : x / 1000 # convert from Pa to kPa
    )[1]['data']

    # Saturation water vapor pressure from temperature (es, in kPa) using the Tetens formula
    # from Campbell & Norman (1998) equation 3.8
    a = 0.661 # kPa
    b = 17.502
    c = 240.97 # C
    es = a * np.exp((b * tair)/(tair + c))

    # Saturation water vapor mixing ratio ws
    ws = 0.622 * es / pair

    # Relative humidity is the ratio of water vapor mixing ratio (which is
    # approximately equal to q) divided by ws. 
    hr = shum / ws

    # D = es(T)-ea = es(Ta)(1-hr), Campbell & Norman (1998), eq. 3.12
    vpd = es * (1 - hr)

    vpd = vpd.apply(lambda x: max(x, 0))

    return vpd

def readin_nldas_data(climate_data: list[tuple[DataFrame, str]], 
                      missing_vars: list[str], lat: float, lon: float):
    '''
    This function retrieves weather data for a site from the NLDAS Hydrology Data 
    Rods Time Series API and fills in any missing data.

    Information on data can be found here: https://disc.gsfc.nasa.gov/information/tools?title=Hydrology%20Data%20Rods
    
    Inputs:
    - `lat`: Latitude 
    - `lon`: Longitude
    - `elev`: Elevation of site

    Outputs:
    - `climate_data`: Tuple with DataFrame and name of site with NLDAS climate 
                      data imputed for missing columns:
        * Solar_Wm2, Rain_mm, Wind_ms.1, Tair_C, Tsoil_C, D_kPa.
    '''
    
    variables = {
        "Solar_Wm2" : ("NLDAS2:NLDAS_FORA0125_H_v2.0:SWdown", None),
        "Rain_mm" : ("NLDAS2:NLDAS_FORA0125_H_v2.0:Rainf", None),
        "Wind_ms.1" : ("", None), # see calc_wind_from_nldas
        "Tair_C" : ("NLDAS2:NLDAS_FORA0125_H_v2.0:Tair", lambda x: x - 273.15), # convert K to C
        "Tsoil_C" : ("NLDAS2:NLDAS_NOAH0125_H_v2.0:SoilT_0_10cm", lambda x: x - 273.15), # convert K to C
        "D_kPa" : ("", None) # see calc_vpd_from_nldas
    }
    
    for i, site_missing_vars in enumerate(tqdm(missing_vars, desc="Sites")):

        c_data, _ = climate_data[i]

        latitude = lat.iloc[i]
        longitude = lon.iloc[i]

        # This assumes data is sorted oldest to newest, which it should be if
        # loaded via readin_climate_data()
        start_date = datetime.strptime(f"{int(c_data.iloc[0]['Year'])}-{int(c_data.iloc[0]['Day'])}-{int(c_data.iloc[0]['Hour'])}", "%Y-%j-%H").strftime("%Y-%m-%dT%H")
        if c_data.iloc[-1]['Day'] == 366: # leap year, for some reason datetime strptime can't handle 366 day
            end_date = (datetime(int(c_data.iloc[-1]['Year']), 1, 1) + timedelta(days=int(c_data.iloc[-1]['Day']) - 1, hours=int(c_data.iloc[-1]['Hour']))).strftime("%Y-%m-%dT%H")
        else:
            end_date = datetime.strptime(f"{int(c_data.iloc[-1]['Year'])}-{int(c_data.iloc[-1]['Day'])}-{int(c_data.iloc[-1]['Hour'])}", "%Y-%j-%H").strftime("%Y-%m-%dT%H")

        start_date, end_date = _convert_time_to_utc(start_date, end_date, latitude, longitude)

        for var in tqdm(site_missing_vars, leave=False, desc="Missing vars"):

            if var == "D_kPa":
                # VPD calculations from air temp, specific humidity, and pressure
                vpd = _calc_vpd_from_nldas(start_date, end_date, latitude, longitude)
                c_data[var] = vpd
                continue

            if var == "Wind_ms.1":
                # Need to calculate diagonal wind from both directions
                wind = _calc_wind_from_nldas(start_date, end_date, latitude, longitude)
                c_data[var] = wind
                continue
            
            # Check that there isn't a weird column
            if var not in variables:
                print(f"Provided var {var} is not one of the variables that can be imputed.")
                continue
            else:
                variable, conv_func = variables[var]

            # Get the time series data for that variable and parse.
            _, df_var = _parse_time_series(
                _get_time_series(
                        start_date=start_date, 
                        end_date=end_date,
                        latitude=latitude,
                        longitude=longitude,
                        variable=variable
                ),
                conv_func
            )

            c_data[var] = df_var['data']

    return climate_data

def readin_noaa_data(noaa_url: str) -> DataFrame:
    '''
    This function reads in NOAA atmospheric CO2 trend data.

    - `noaa_url`: URL to NOAA annual atmospheric Carbon Dioxide concentration 
                  records, check here: https://gml.noaa.gov/ccgg/trends/gl_data.html
    '''
    try:
        noaa_data = pd.read_csv(noaa_url, skiprows=range(37), header=0, names=['Year', 'ca_ppm', 'sd_year'])
    except Exception as e:
        raise(Exception(f"ERROR: Unable to fetch NOAA CO2 data: {e}"))
            
    return noaa_data

def readin_climate_data(climate_file: list[str], sites: list[str]) -> list[tuple[DataFrame, str]]:
    '''
    Returns climate data for every specified row in param_data, along with the
    site.

    Climate_file must be a csv with at least 'Year', 'Day', and 'Hour'. Ensure
    that the time is local to the site as it will be converted to UTC time if
    data is needed to be fetched from NLDAS.
    '''
    dt = []
    seen_sites = set() # ensure no duplicate files for same site as each site has same weather
    for url, site in zip(climate_file, sites):
        if site not in seen_sites:
            try:
                df = pd.read_csv(url, on_bad_lines="error", index_col=False)

                # Sort dataframe and ensure year, day, and hour are ints
                df = df.sort_values(by=['Year', 'Day', 'Hour'], ascending=True)
                df['Year'] = df['Year'].astype(int)
                df['Day'] = df['Day'].astype(int)
                df['Hour'] = df['Hour'].astype(int)

                dt.append((df, site))
                seen_sites.add(site)
            except:
                raise(Exception("File unable to be opened. Check file format and columns."))
    return dt

def check_climate_data(climate_data: list[tuple[DataFrame, str]]) -> list[list]:
    '''
    Checks for missing climate_data columns and inputs NLDAS data if missing.

    Input:
    - `climate_data`: Tuple with DataFrame that has bare minimum: year, date, 
                      and hour, and the site corresponding to that climate_data.

    Output:
    - `missing_columns`: List with Year, Day, Hour, Solar_Wm2, Rain_mm, 
    Wind_ms.1, Tair_C, Tsoil_C, D_kPa if missing columns.
    '''
    columns = {'Year', 'Day', 'Hour', 'Solar_Wm2', 'Rain_mm', 'Wind_ms.1', 
                'Tair_C', 'Tsoil_C', 'D_kPa'}
    missing_columns = []
    for data, _ in climate_data:
        mis_col = []
        for column in columns:
            if column not in data.columns:
                mis_col.append(column)
        missing_columns.append(mis_col)

    return missing_columns

def readin_data(filename: str) -> DataFrame:
    '''
    Given a filename to csv file, reads in the data file and stores in 
    a dataframe. Used for building configuration and parameter files.
    
    The parameter data file can contain multiple rows for different species.
    Same sites will be grouped into the same file. Different sites will be
    outputted in unique files

    Parameter file should contain the following columns below. If data is
    unavailable for a certain column, supply it anyways and it will either be
    filled with the default or an error will be thrown. Recommend using
    data_template.csv to start:

    - `sp`: Species
    - `region`: Geographical region of simulations
    - `site`: Site/simulation ID
    - `latitude`: Latitude in degree fraction north
    - `longitude`: Longitude in degree fraction west
    - `elevation`: Site elevation in meters above sea level
    - `slopeI`: Slope inclination; degrees from horizontal
    - `slopeA`: Slope aspect; counterclockwise degrees from south
    - `gWaterP`: Ground water pressure
    - `gWaterDist`: Distance to ground water source in meters from the bottom
                    of the root system.
    - `atmTrans`: Atmospheric transmittance from weather data
    - `emiss`: Long wave emissivity
    - `co2AmbPPM`: Atmospheric/experiment CO2 ppm (can be empty if multiple 
                   years are supplied later when building season_limits)
    - `layers`: Number of soil layers.
    - `fieldCapFrac`: Fraction that field capacity is of saturation 
                      (minus residual).
    - `fieldCapPercInit`: Percent field capacity for starting the season.
    - `rockFrac`: Fraction of soil volume as rocks (0-1).
    - `soilAbsSol`: Absorptivity of soil surface for solar.
    - `rhizoPer`: Average percent of whole plant resistance in rhizosphere 
                  (maximum soil limitation)
    - `texture`: USDA soil texture category (equal for all layers but could 
                 also be determined per layer)
    - `baPerGa`: Basal area per ground area (m2 ha-1)
    - `leafAreaIndex`: Canopy leaf area index (m2 m-2)
    - `soilXHeight`: Height above soil surface for understory wind and gh 
                     in meters
    - `height`: Average tree height in meters
    - `leafPerBasal`: Initial leaf area per basal area (m2 m-2)
    - `leafWidth`: Leaf width in meters
    - `leafAngleParam`: Leaf angle parameter; CN 15.4
    - `aspect`: Max radius of root system per max depth
    - `rootBeta`: Root biomass distribution is allocated based on the equation 
                  reported in Love et al (2019): M = 1 - Beta^d, where M is 
                  the fraction of biomass above depth d expressed in cm. We 
                  find the Beta that provides an M of 0.995 for the maximum 
                  rooting depth.
    - `leafPercRes`: Saturated % of tree resistance in leaves
    - `pInc`: Pressure increment for curve generation, (MPa) - higher is 
              faster, but less accurate (setting too high can cause 
              Newton-Rhapson root pressure solving failure)
    - `lsc`: Leaf specific conductance in mmol m-2 w-1 MPa-1
    - `cr`: Root element Weibull parameter c
    - `br`: Root element Weibull parameter b
    - `cs`: Stem element Weibull parameter c
    - `bs`: Stem element Weibull parameter b
    - `cl`: Leaf element Weibull parameter c
    - `bl`: Leaf element Weibull parameter b
    - `sapwoodT`: Change in sapwood per change in diameter at breast height
    - `conduitDiam`: Vessel or tracheid diameter in um
    - `qMax`: Quantum yield of electron transport; moles e per mols photons
    - `vmax25`: Maximum carboxylation rate (vmax) at 25C (umol m-2 s-1)
    - `jmax25`: Maximum electron transport rate (jmax) at 25C (umol m-2 s-1), 
                can be assumed to be Vmax25 * 1.67
    - `kc25`: Michaelis-Menten constant for CO2 in mole fraction at 25C. 
              Bernacchi T response
    - `ko25`: Michaelis-Menten constant for O2 in mole fraction at 25C. 
              Bernacchi T response
    - `comp25`: Photorespiratory compensation point in mole fraction at 25C. 
                Bernacchi T response
    - `thetaC`: Shape factor for A-ci colimitation
    - `havmax`: Temp-dependency parameters from Leunig 2002 (J mol-1)
    - `hdvmax`: Temp-dependency parameters from Leunig 2002 (J mol-1)
    - `svvmax`: Temp-dependency parameters from Leunig 2002 (J mol-1)
    - `lightCurv`: Temp-dependency parameters from Leunig 2002
    - `lightComp`: Light compensation point in ppfd
    - `hajmax`: Temp-dependency parameters from Leunig 2002 (J mol-1)
    - `hdjmax`: Temp-dependency parameters from Leunig 2002 (J mol-1)
    - `svjmax`: Temp-dependency parameters from Leunig 2002 (J mol-1 K-1)
    - `gWaterEnable`: Turns on/off groundwater flow
    - `soilRedEnable`: Turns on/off soil redistribution routine
    - `soilEvapEnable`: Turns on/off soil evaporation routine
    - `rainEnable`: Turns on/off rain inputs
    - `useGSDataStress`: Turns on/off growing season data for multiple year
                         modeling
    - `useGSDataOpt`: Turns on/off growing season for BAGA optimization
    - `refilling`: Turns on/off xylem refilling.
    - `predawnsMode`: Turns on/off if model should consider pre-dawn water
                      potential values.
    - `cavitFatigue`: Turns on/off xylem stress hystersis.
    - `stemOnly`: Turns on/off xylem stress hystersis only in stem.
    - `multipleSP`: Turns on/off if model config has 1 species per site or 
                    multiple
    - `speciesN`: Number of species to run the model
    - `climateData`: Path to file with climate forcing variables `dataset.csv`
    '''
    try:
        df = pd.read_csv(filename, on_bad_lines="error", index_col=False)
    except:
        print("File unable to be opened. Check file format and columns.")
    return df

# Not used here, calculated in model if LSC/LSC_pref is provided
def _calculate_kmax(la_ba: pd.Series, 
                    lsc: pd.Series, 
                    r_leaf: pd.Series) -> pd.Series:
    '''
    Calculates the estimated maxmimum conductance of the whole tree given
    leaf area to basal area ratio, leaf specific conductance, and the resistance
    of the leaf.
    '''
    k_leaf = la_ba * lsc; # Calculate the total canopy leaf conductance (mmol s^-1 m^-2 MPa^-1) in relation to BA
    # Unit conversions of Kleaf
    # 1 mmol H2O = 18.01528*10^-6 kg H2O
    # 1 s = 1/3600 h
    k_leaf = k_leaf * (18.01528 * 1e-6) * 3600; # Kleaf (kg h^-1 m^-2 MPa^-1)
    # Now we can solve for tree Kmax knowing that K=1/Resistance
    # As the resistances are in series there is the following equation that has
    # to be satisfied
    # Rtree = Rroot+Rstem+Rleaf, thus: 1/Kmax = 1/Kroot+1/Kstem+1/Kleaf
    r_tree = (1/k_leaf) * 100/r_leaf; # For calcutating total tree resistance
    k_max = 1/r_tree; # Kmax of the tree (kg h^-1 MPa-1 m^2)
    return k_max

def _solar_noon_correction(longitude: pd.Series, utc_weather: int = 8):
    '''
    Calculate the UTC time zone and solar noon correction.
    '''
    utc = longitude/15 # round to nearest integer in direction of 0
    
    # Calculate fraction of an hour away from noon
    if (utc.isna().any()):
        corrected = longitude/15 - utc 
    else:
        corrected = (longitude/15 - utc) + (utc - utc_weather)

    return corrected

def build_config_param_files(data: DataFrame, out_data_path: str) -> None:
    '''
    Builds configuration and parameter files for separate sites and different
    species based on a given data DF.

    Input:
    - `data`: Data file that has all of the information needed for building
              configuration and parameter files.
    - `out_data_path`: Path for outputting

    Output:
    - `None`: Configuration and parameter files are outputted and built in
              `out_data_path/site/configuration.csv` and 
              `out_data_path/site/parameters.csv`
    '''

    # define data defaults
    # Default values for various parameters
    co2 = 410               # Atmospheric CO2 concentration in ppm
    leaf_angle = 1          # Leaf angle parameter
    leaf_width = 0.04096    # Leaf width in meters
    soil_Xheight = 1        # Height above soil surface for understory wind and gh in meters
    soil_layers = 5         # Number of soil layers
    field_cap_frac = 0.5    # Fraction that field capacity is of saturation (minus residual)
    field_cap_init = 100    # Percent field capacity for starting the season
    soil_abs_sol = 0.94     # Absorptivity of soil surface for solar
    p_inc = 0.00075         # Pressure increment for curve generation in MPa
    rhizo_per = 50          # Average percent of whole plant resistance in rhizosphere
    leaf_per = 25           # Saturated % of tree resistance in leaves
    kmax = 0                # if kmax is not provided, set to 0 and model calculates
                            # kmax based on LSC parameters (lsc/lsc_pref)
    lsc = 0                 
    lsc_pref = 0
    root_aspect = 1         # Max radius of root system per max depth
    root_beta = 0.95
    emiss = 0.97            # Long wave emissivity
    atm_trans = 0.75        # Atmospheric transmittance from weather data
    light_curv = 0.9        # 0.9 is curvature of light per Medlyn 2002
    q_max = 0.3             # "The value of α (q_max) was fixed at 0.3 mol electrons 
                            # mol−1 photon, based on an average C3 photosynthetic 
                            # quantum yield of 0.093 and a leaf absorptance of 0.8 
                            # (Long, Postl & Bolharnordenkampf 1993)." -- Medyln 2002
    mole_frac = 0.21        # Mole fraction (not used in the provided code)
    ground_water_p = 0      # Ground water pressure
    ground_water_d = 1      # Distance to ground water source in meters from the bottom of the root system
    ground_water = "n"      # Turns on/off groundwater flow
    soil_evap = "n"         # Turns on/off soil evaporation routine
    refilling = "n"         # Turns on/off xylem refilling
    rainfall = "y"          # Turns on/off rain inputs
    gs = "n"                # Turns on/off growing season data for multiple year modeling
    predawn = "n"           # Turns on/off if model should consider pre-dawn water potential values

    site_param_data = data.groupby(['site'])
    for site, group in site_param_data:
        df_config = DataFrame()
        df_param = DataFrame()
        df_config['i_gWaterEnable'] = group['gWaterEnable'].fillna(ground_water)
        df_config['i_soilRedEnable'] = group['soilRedEnable']
        df_config['i_soilEvapEnable'] = group['soilEvapEnable'].fillna(soil_evap)
        df_config['i_rainEnable'] = group['rainEnable'].fillna(rainfall)
        df_config['i_useGSDataStress'] = group['useGSDataStress'].fillna(gs)
        df_config['i_useGSDataOpt'] = group['useGSDataOpt']
        df_config['i_refilling'] = group['refilling'].fillna(refilling)
        df_config['i_predawnsMode'] = group['predawnsMode'].fillna(predawn)
        df_config['i_cavitFatigue'] = group['cavitFatigue']
        df_config['i_stemOnly'] = group['stemOnly']
        df_config['i_multipleSP'] = group['multipleSP']
        df_config['i_speciesN'] = group['speciesN']

        # Get current working directory and provide relative file path to
        # dataset.csv and seasonlimits.csv
        cwd = os.getcwd()
        if out_data_path != ".":
            if out_data_path[0] == '/':
                cwd = cwd + out_data_path
            else:
                cwd = cwd + "/" + out_data_path
                
        df_config['i_ClimateData'] = f"{cwd}/{site[0]}/dataset.csv"
        df_config['i_GSData'] = f"{cwd}/{site[0]}/seasonlimits.csv"
        df_param['i_sp'] = group['sp']
        df_param['i_region'] = group['region']
        df_param['i_site'] = group['site']
        df_param['i_latitude'] = group['latitude']
        df_param['i_longitude'] = group['longitude']
        df_param['i_elevation'] = group['elevation']
        df_param['i_slopeI'] = group['slopeI']
        df_param['i_slopeA'] = group['slopeA']
        df_param['i_gWaterP'] = group['gWaterP'].fillna(ground_water_p)
        df_param['i_gWaterDist'] = group['gWaterDist'].fillna(ground_water_d)
        df_param['i_atmTrans'] = group['atmTrans'].fillna(atm_trans)
        df_param['i_solarNoon'] = 0 # _solar_noon_correction(df_param['i_longitude']), your mileage may vary by using this
        df_param['i_emiss'] = group['emiss'].fillna(emiss)
        df_param['i_co2AmbPPM'] = group['co2AmbPPM'].fillna(co2)
        df_param['i_layers'] = group['layers'].fillna(soil_layers)
        df_param['i_fieldCapFrac'] = group['fieldCapFrac'].fillna(field_cap_frac)
        df_param['i_fieldCapPercInit'] = group['fieldCapPercInit'].fillna(field_cap_init)
        df_param['i_rockFrac'] = group['rockFrac']
        df_param['i_soilAbsSol'] = group['soilAbsSol'].fillna(soil_abs_sol)
        df_param['i_rhizoPer'] = group['rhizoPer'].fillna(rhizo_per)
        df_param['i_texture'] = group['texture']
        df_param['i_baperga'] = group['baPerGa']
        df_param['i_leafAreaIndex'] = group['leafAreaIndex']
        df_param['i_soilXHeight'] = group['soilXHeight'].fillna(soil_Xheight)
        df_param['i_height'] = group['height']
        df_param['i_leafPerBasal'] = group['leafPerBasal']
        df_param['i_leafWidth'] = group['leafWidth'].fillna(leaf_width)
        df_param['i_leafAngleParam'] = group['leafAngleParam'].fillna(leaf_angle)
        df_param['i_aspect'] = group['aspect'].fillna(root_aspect)
        df_param['i_rootBeta'] = group['rootBeta'].fillna(root_beta)
        df_param['i_leafPercRes'] = group['leafPercRes'].fillna(leaf_per)
        df_param['i_kmaxTree'] = group['kmax_tree'].fillna(kmax)
        df_param['i_pinc'] = group['p_inc'].fillna(p_inc)
        df_param['i_LSC'] = group['lsc'].fillna(lsc)
        df_param['i_LSCpref'] = group['lsc_pref'].fillna(lsc_pref)
        df_param['i_cr'] = group['cr']
        df_param['i_br'] = group['br']
        df_param['i_cs'] = group['cs']
        df_param['i_bs'] = group['bs']
        df_param['i_cl'] = group['cl']
        df_param['i_bl'] = group['bl']
        df_param['i_sapwoodT'] = group['sapwoodT']
        df_param['i_conduitDiam'] = group['conduitDiam']
        df_param['i_qMax'] = group['qMax'].fillna(q_max)
        df_param['i_vmax25'] = group['vmax25']
        df_param['i_jmax25'] = group['jmax25']
        df_param['i_kc25'] = group['kc25']
        df_param['i_ko25'] = group['ko25']
        df_param['i_comp25'] = group['comp25']
        df_param['i_thetaC'] = group['thetaC']
        df_param['i_havmax'] = group['havmax']
        df_param['i_hdvmax'] = group['hdvmax']
        df_param['i_svvmax'] = group['svvmax']
        df_param['i_lightCurv'] = group['lightCurv'].fillna(light_curv)
        df_param['i_lightComp'] = group['lightComp']
        df_param['i_hajmax'] = group['hajmax']
        df_param['i_hdjmax'] = group['hdjmax']
        df_param['i_svjmax'] = group['svjmax']

        # Create path, default data_path = "." (current directory)
        if not os.path.exists(f"{out_data_path}/{site[0]}"):
            os.makedirs(f"{out_data_path}/{site[0]}")
        df_config.to_csv(f"{out_data_path}/{site[0]}/configuration.csv", index=False)
        print(f"{site[0]} configuration file saved to {out_data_path}/{site[0]}/configuration.csv")
        df_param.to_csv(f"{out_data_path}/{site[0]}/parameters.csv", index=False)
        print(f"{site[0]} parameter file saved to {out_data_path}/{site[0]}/parameters.csv")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", help="File path to data file. See data_template.csv. `file_builder.py [path/to/file/file_builder.csv]`")
    parser.add_argument("-i", "--initial-year", help="Initial year to pull climate data from.", required=True, type=int)
    parser.add_argument("-f", "--final-year", help="Final year to pull climate data from.", required=True, type=int)
    parser.add_argument("-p", "--path", help="Optional data path to write parameter files in.", default=".")
    args = parser.parse_args()
    
    # Read in data and build config and param files.
    print("Reading in data and building the configuration and parameter files...")
    param_data = readin_data(args.data_file)
    build_config_param_files(param_data, args.path)

    # Get NOAA data for CO2 levels
    print("Fetching NOAA data for CO2 levels...")
    noaa_data = readin_noaa_data("https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_gl.csv")

    # Readin climate data specified in configuration
    print("Reading in climate data specified in configuration...")
    climate_data = readin_climate_data(param_data['climateData'], param_data['site'])

    # Check if there is climate data within the provided range
    for c_data in climate_data:
        assert(args.initial_year >= c_data[0]['Year'].min() and args.final_year <= c_data[0]['Year'].max())

    # Check climate data for missing columns and if there are the missing values
    # will be imputed from NLDAS site data.
    print("Checking climate data for missing columns...")
    missing_columns = check_climate_data(climate_data)
    print("Imputing missing values from NLDAS site data...")
    climate_data = readin_nldas_data(climate_data, missing_columns, param_data['latitude'], param_data['longitude'])

    # Calculate the output the growing season range, as well as shortened dataset
    # within initial_year and final_year.
    print("Calculating and outputting the growing season range and shortened dataset...")
    build_growing_season_data(noaa_data, climate_data, args.path, args.initial_year, args.final_year)

if __name__ == "__main__":
    main()