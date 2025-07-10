# Carbon Gain vs Hydraulic Risk Stomatal Optimization Model V3.0

***Ga****in* ***RI****sk* ***S****tomatal* ***O****ptimization* ***M****odel* (**GARISOM**)

__Coding language:__ C++

__Authors:__ Colin Pannikkat, German Vargas G. & William R.L. Anderegg

__Contact:__ vargasgg (at) oregonstate (dot) edu

------------

## Introduction:

The model uses a stomatal carbon gain vs. hydraulic risk optimization (Sperry et al. 2017 -- see references below) combined with a soil water budget to run a continuous growing season simulation, producing stand outputs such as net carbon assimilation (Anet), internal [CO2] (Ci), transpiration (E), total evapotranspiration (ET), and element conductances (k) on an hourly and summary basis (see output details below).

------------

## Basic usage instructions:

First, fork the repository to your own account. If you don't have a Github account, just download the repository. Provide plant and site traits through the parameters file, then supply hourly weather data to drive model. The model will expect these files to be located in the current working directory. See the included examples and the details below for formatting and units. These example files contain all inputs necessary to test-run the model immediately after building.

__a)__ To run, build and execute the model program (no command line arguments -- runs from files in the working directory). Building requires C++17, GNU example:

```{}
g++ -std=c++17
```

The -O2 ~~and -ffast-math~~ optimizations are recommended with GNU compilers:

```{}
g++ -std=c++17 -O2
```

Using ffast-math does not guarantee reproducibility as this model uses a lot of floating point operations and is incredibly sensitive to numerical instabilities. See [here](https://stackoverflow.com/questions/7420665/what-does-gccs-ffast-math-actually-do) for more information. There is currently a header guard in `01Utilities.h` that prevents the use of it, but if you would like to compile with it, then feel free to comment that out.

__b)__ To build and run from the terminal in a Linux/Mac OS (*nix) system:

__b.1)__ Clone the repository in a desired location in your system.

```
git clone https://github.com/gevargu/garisom.git
```

__b.2)__ Navigate to the repository by usin the command cd:

```
cd garisom/02_program_code
```

__b.3)__ Run the following command will compile the code and build an executable:

```
make
```

__b.4)__ Run this program from the same folder with this command:

```
./run [param_data_file] [config_data_file] [config_setting] [output_path]
```

The `param_data_file`, `config_data_file`, and `config_setting` are optional command line arguments. If not provided, the default file paths will be used which are defined in `01Utilities.h`.

```c++
#define CONFIG_FILE_PATH "../03_test_data/configuration.csv" // path to configuration file
#define PARAMETER_FILE_PATH "../03_test_data/parameters.csv" // path to parameter file
#define CONFIG_SETTING 0                                     // which configuration to use (this is 0-indexed)
```
You may also run `./run --help` for more explicit usage instructions.

__b.5)__ Press <kbd>Command</kbd> + <kbd>C</kbd> if you want to stop the model before it completes (<kbd>Ctrl</kbd> + <kbd>C</kbd> in Windows/Linux)

The following files (included in this repository) should be located in the working directory (normally the same directory as the executable) before running:
	
  - __parameters.csv__ (site, atmospheric, soils, stand, plant, hydraulics, and carbon assimilation parameters).
  - __configuration.csv__ (model controls)
  - __dataset.csv__ (hourly weather drivers).
  - __dataheader.csv__ (a header row for the hourly data output).
  - __seasonlimits.csv__ (growing season limits and yearly atmospheric CO2, only required if using "sequential year mode" described below)

You can also provide a relative directory path to **parameters.csv** and **configuration.csv** for input to the model. If the files were built using `file_builder.py`, the seasonlimits and dataset files will automatically be inputted with absolute file paths.

Upon completion, two output files are produced:
	
  - The __timesteps_output.csv__ output contains all of the hourly model outputs corresponding to the weather inputs.
  - The __sum_output.csv__ output includes various total values for each year (for example, net growing season productivity Anet, net transpiration E, etc.)

------------

## Generating Model Input Files

Rather than generating the input files by hand, users can use the `file_builder.py` script. This generates the parameter, configuration, climate data, and growing season files.

**a)** User's may want to isolate their local Python installation, of which they can use [**venv**](https://docs.python.org/3/tutorial/venv.html) or equivalent to create a virtual environment to install the requirements into. After doing so (or not), users can then run the following commands to install all necessary Python packages.

**Activate the environment.**
```
source env/bin/activate
```
**Then install the required packages given the requirements file.**
```
pip install -r requirements.txt
```

**b)** Users should then prepare an initial data file using `data_template.csv`, some values may be optional and will be replaced by defaults if left empty.

```python
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
```

The only other file that is required is a `dataset.csv` file with at least **Year**, **Day**, and **Hour**. Put this in the `climateData` column.

**c)** Then, run the file builder script.

```
python file_builder.py [-h] -i INITIAL_YEAR -f FINAL_YEAR [-p PATH] data_file
```

Inputs:
- `INITIAL_YEAR`: Initial year to include in output dataset.csv and calculate seasonlimits for.
- `FINAL_YEAR`: Final year to include in output, this is inclusive.
- `PATH`: Optional path to output all the files to, by default this is set to `"."` (current working directory).
- `data_file`: The file path to the full data file (**data_template.csv**).

Since the seasonlimits and climate data are homogenous for sites, output files will be split via the specified site. Therein, configuration and parameter files are outputted and built in `PATH/site/configuration.csv` and `PATH/site/parameters.csv`.

If a column is missing from `datasets.csv`, it will fetch the NLDAS equivalent based on plot location.

------------

## Model Input Files

### Model Parameters

Configure plant traits and other parameters in __parameters.csv__ (expected input units are indicated).

| Group    		| Parameter       	    	| Description										|
| ---------------------	| ----------------------------- | -------------------------------------------------------------------------------------	|
| Site			| __i_sp__			| Species or PFT represented in parameter data	|
| Site			| __i_region__			| Geographical region of the simulations.	|
| Site			| __i_site__			| Site/simulation ID. This ID is used for naming the result files that are exported	|
| Site			| __i_latitude__		| Latitude in degree fraction north.	|
| Site			| __i_longitude__		| Longitude in degree fraction west.	|
| Site			| __i_elevation__		| Site elevation in m above sea level.	|
| Site			| __i_slopeI__			| Slope inclination; degrees from horizontal. (not used in model)|
| Site			| __i_slopeA__			| Slope aspect; counterclockwise degrees from south. (not used in model)|
| Site			| __i_gWaterP__			| Ground water pressure (MPa)	|
| Site			| __i_gWaterDist__		| Distance to ground water source in m from the bottom of the rootsystem.	|
| Atmosphere		| __i_atmTrans__		| Atmospheric transmittance from weather data (set to 0.65 as default if no data available).	|
| Atmosphere		| __i_solarNoon__		| Solar noon correction from weather data in hours. (calculated in `file_builder.py`)	|
| Atmosphere		| __i_emiss__			| Long wave emissivity.	|
| Atmosphere		| __i_co2AmbPPM__		| Atmospheric/experiment CO2 ppm, it will update if working with multiple years.	|
| Soil			| __i_layers__			| Number of soil layers (select 1-5). Can technically do greater than 5, but would need to add that pressure column to the data header. |
| Soil			| __i_fieldCapFrac__		| Fraction that field capacity is of saturation (minus residual).	|
| Soil			| __i_fieldCapPercInit__	| Percent field capacity for starting the season.	|
| Soil			| __i_rockFrac__		| Fraction of soil volume as rocks (0-1).	|
| Soil			| __i_soilAbsSol__		| Absorptivity of soil surface for solar.	|
| Soil			| __i_rhizoPer__		| Average percent of whole plant resistance in rhizosphere (maximum soil limitation)	|
| Soil			| __i_texture__			| USDA soil texture category (equal for all layers but could also be determined per layer)	|
| Stand			| __i_baperga__			| Basal area per ground area (m2 ha-1)	|
| Stand			| __i_leafAreaIndex__		| Canopy lai (m2 m-2)	|
| Stand			| __i_soilXHeight__		| Height above soil surface for understory wind and gh in m	|
| Stand			| __i_height__			| Average tree height in m	|
| Stand			| __i_treeToPhotoLAI__		| Deprecated since V3.0 |
| Stand			| __i_leafPerBasal__		| Initial leaf area per basal area per individual tree; (m2 m-2)	|
| Tree			| __i_leafWidth__		| Leaf width in m	|
| Tree			| __i_leafAngleParam__		| Leaf angle parameter; CN 15.4	|
| Tree			| __i_aspect__			| Max radius of root system per max depth	|
| Tree			| __i_rootBeta__		| Root biomass distribution is allocated based on the equation reported in Love et al (2019): M = 1 - Beta^d, where M is the fraction of biomass above depth d expressed in cm. We find the Beta that provides an M of 0.995 for the maximum rooting depth.	|
| Hydraulics		| __i_leafPercRes__		| Saturated % of tree resistance in leaves	|
| Hydraulics		| __i_kmaxTree__		| Kmax of tree in kg hr-1 m-2 MPa-1 per basal area (can include this, or LSC/LSCpref)|
| Hydraulics		| __i_pinc__			| Pressure increment for curve generation, (MPa) - higher is faster, but less accurate (setting too high can cause Newton-Rhapson root pressure solving failure)	|
| Hydraulics		| __i_LSC__			| Leaf specific conductance in mmol m-2 s-1 MPa-1 (per leaf area), used in calculated kmaxTree |
| Hydraulics		| __i_LSCpref__			| Water potential for LSC |
| Hydraulics		| __i_cr__			| Root element Weibull parameter c	|
| Hydraulics		| __i_br__			| Root element Weibull parameter b	|
| Hydraulics		| __i_cs__			| Stem element Weibull parameter c	|
| Hydraulics		| __i_bs__			| Stem element Weibull parameter b	|
| Hydraulics		| __i_cl__			| Leaf element Weibull parameter c	|
| Hydraulics		| __i_bl__			| Leaf element Weibull parameter b	|
| Hydraulics		| __i_sapwoodT__		| Change in sapwood per change in diameter at breast height	|
| Hydraulics		| __i_conduitDiam__		| Vessel or tracheid diameter in um	|
| Photosynthesis	| __i_qMax__			| Quantum yield of electron transport; moles e per mols photons	|
| Photosynthesis	| __i_vmax25__			| Maximum carboxylation rate (vmax) at 25C (umol m-2 s-1)	|
| Photosynthesis	| __i_jmax25__	 		| Maximum electron transport rate (jmax) at 25C (umol m-2 s-1), can be assumed to be Vmax25 * 1.67	|
| Photosynthesis	| __i_kc25__			| Michaelis-Menten constant for CO2 in mole fraction at 25C. Bernacchi T response	|
| Photosynthesis	| __i_ko25__			| Michaelis-Menten constant for O2 in mole fraction at 25C. Bernacchi T response	|
| Photosynthesis	| __i_comp25__			| Photorespiratory compensation point in mole fraction at 25C. Bernacchi T response	|
| Photosynthesis	| __i_thetaC__			| Shape factor for A-ci colimitation	|
| Photosynthesis	| __i_havmax__			| Temp-dependency parameters from Leunig 2002 (J mol-1)	|
| Photosynthesis	| __i_hdvmax__			| Temp-dependency parameters from Leunig 2002 (J mol-1)	|
| Photosynthesis	| __i_svvmax__			| Temp-dependency parameters from Leunig 2002 (J mol-1)	|
| Photosynthesis	| __i_lightCurv__		| Temp-dependency parameters from Leunig 2002	|
| Photosynthesis	| __i_lightComp__		| Light compensation point in ppfd	|
| Photosynthesis	| __i_hajmax__			| Temp-dependency parameters from Leunig 2002 (J mol-1)	|
| Photosynthesis	| __i_hdjmax__			| Temp-dependency parameters from Leunig 2002 (J mol-1)	|
| Photosynthesis	| __i_svjmax__			| Temp-dependency parameters from Leunig 2002 (J mol-1 K-1)	|

### Model Configuration

All model configuration parameters in **configuration.csv**.

| Group         | Model Control         | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| ------------- | --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Soil          | __igWaterEnable__     | Turns on/off groundwater flow. Values: n (off); y (on). It provides an unlimited source of water at a set potential and distance below the root layers. This water will flow up into the soil layers, and potentially allow layers to fill above field capacity (from the bottom layer up). When disabled (default), the only sources of water input will be the initial fraction of field capacity and observed rainfall (and any water over field capacity will become "drainage").                                                                      |
| Soil          | __i_soilRedEnable__   | Turns on/off soil redistribution routine. Values: n (off); y (on). It allows water to flow between soil layers.                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| Soil          | __i_soilEvapEnable__  | Turns on/off soil evaporation routine. Values: n (off); y (on). It enables simulation of water evaporation from the surface soil layer.                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Climate       | __i_rainEnable__      | Turns on/off rain inputs. Values: n (off); y (on). It allows for precipitation events. Weather data rainfall will be ignored if disabled.                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Climate       | __i_useGSDataStress__ | Turns on/off growing season data for multiple year modeling. Vakyes: n (off); y (on). If enabled, multiple years will be run "sequentially" with on and off seasons defined in seasonlimits_2.0.0.csv and a continuous water budget. When disabled (default), all weather timesteps provided are treated as part of the growing season and the user is expected to truncate individual years to their start/end days. Water budget is reset between years when disabled, treating years as totally independent.                                            |
| Climate       | __i_useGSDataOpt__    | Turns on/off growing season data for multiple year modeling during BAGA optimization. Values: n (off); y (on).                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Hydraulics    | __i_refilling__       | Turns on/off xylem refilling within a growing season. Values: n (off); y (on). It allows trees to restore lost conductance, however the refilling model is not sufficient to simulate authentic xylem refilling behavior and has __not been thoroughly tested in the current version of the code__.                                                                                                                                                                                                                                                        |
| Hydraulics    | __i_predawnsMode__    | Turns on/off if model should consider measured pre-dawn water potential values. Values: n (off); y (on). If set to 'y', disables soil simulation and runs from hourly inputs of canopy predawn water potential. These are read from the "rain" column (rain is not used with the soil sim disabled) in MPa. See "dataset - predawns example.csv". This mode is especially useful when comparing to other models which run from canopy predawn measurements, and generally runs significantly faster as it does not need to solve for root layer pressures. |
| Hydraulics    | __i_cavitFatigue__    | Turns on/off xylem stress hysteresis to carry effects from previous growing season. Values: n (off); y (on). It allows for a weighted estimation of xylem vulnerability to embolism                                                                                                                                                                                                                                                                                                                                                                        |
| Hydraulics    | __i_stemOnly__        | Turns on/off xylem stress hysteresis only in stem xylem. Values: n (off); y (on). When disabled it allows for a weighted estimation of xylem vulnerability to embolism for both stem and roots                                                                                                                                                                                                                                                                                                                                                             |
| Community     | __i_multipleSP__      | Turns on/off whether our model configuration has 1 species per site (monodominant) or multiple species per site (diverse). Values: n (off); y (on).                                                                                                                                                                                                                                                                                                                                                                                                        |
| Community     | __i_speciesN__        | Number of species/PFT to run the model. The species number indicated should correspond to the row in the parameter file.                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Forcing files | __i_ClimateData__     | Path to file with climate forcing variables __dataset.csv__                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| Forcing files | __i_GSData__          | Path to file with growing season data __seasonlimits.csv__                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |

### Weather Data (`dataset.csv`)

See example data for formatting. Weather drivers should be in hourly timesteps and can include multiple years of data. 

Inputs:
- `Year`: Year.
- `Day`: Julian Day (1-366).
- `Hour`: Hour (0-23).
- `Solar_Wm2`: Obs. Solar (W m-2).
- `Rain_mm`: Rain (mm).
- `Wind_ms.1`: Wind (m s-1).
- `Tair_C`: Tair (C).
- `Tsoil_C`: Tsoil (C).
- `D_kPa`: D (kPa).

If one of these input parameters is not available, `file_builder.py` is able to use the coordinates specified in the parameter file to provide NLDAS data for the specific plot (see above section).

### Growing Season (`seasonlimits.csv`)

| Column    | Description                                |
| --------- | ------------------------------------------ |
| Year      | Year for growing season and CO2 levels     |
| Start_day | Initial day of growing season              |
| End_day   | Last day of growing season                 |
| ca_ppm    | Air CO2 levels near plot according to NOAA |

To estimate the growing season, following [Sperry et al. 2019 PNAS](https://www.pnas.org/doi/abs/10.1073/pnas.1913072116), the script fits three lines to the sections of the cumulative thermal degree days above 5°C to identify the inflexion points as those where the three lines that intersect provide the best fit based on mean squared errors.

------------

## Outputs:

- Autosave is always enabled regardless of the setting, as this version of the model has no alternative output method. Output files will be generated in the working directory when the run completes.

### Hourly Outputs

| Column                | Unit                     | Description                                                                                                          |
| --------------------- | ------------------------ | -------------------------------------------------------------------------------------------------------------------- |
| P0                    | -MPa                     | Soil pressure for layer 0 (soil top layer)                                                                           |
| P1                    | -MPa                     | Soil pressure for layer 1                                                                                            |
| P2                    | -MPa                     | Soil pressure for layer 2                                                                                            |
| P3                    | -MPa                     | Soil pressure for layer 3                                                                                            |
| P4                    | -MPa                     | Soil pressure for layer 4                                                                                            |
| P5                    | -MPa                     | Soil pressure for layer 5                                                                                            |
| P-PD                  | -MPa                     | Predawns sunlit canopy pressure                                                                                      |
| P-MD                  | -MPa                     | Midday sunlit canopy pressure                                                                                        |
| E-MD                  | mmol s-1 m-2 (LA)        | Midday sunlit plant transpiration per leaf area                                                                      |
| GW                    | mmol s-1 m-2 (LA)        | Sunlit canopy stomatal conductance per leaf area                                                                     |
| leaf-air-vpd          | kPa                      | Sunlit leaf to air vapor pressure deficit. Difference in water vapor pressure between a leaf and the surrounding air |
| leaftemp              | C                        | Sunlit leaf temperature in sunlit canopy layer                                                                       |
| Anet-la               | $\mu$-mol s-1 m-2 (LA)   | Sunlit net carbon/photosynthetic assimilation per leaf area                                                          |
| ci                    | Pa                       | Sunlit canopy partial pressure of CO2                                                                                |
| PPFD                  | $\mu$-mol s-1m-2         | Sunlit photon flux density in sunlit canopy layer                                                                    |
| S-P-MD                | -MPa                     | Midday sunlit canopy pressure                                                                                        |
| S-E-MD                | mmol s-1 m-2 (LA)        | Midday shaded plant transpiration per leaf area                                                                      |
| S-GW                  | mmol s-1 m-2 (LA)        | Shaded canopy stomatal conductance per leaf area                                                                     |
| S-leaf-air-vpd        | kPa                      | Shaded leaf to air vapor pressure deficit.                                                                           |
| S-Anet-la             | $\mu$-mol s-1 m-2 (LA)   | Shaded net carbon/photosynthetic assimilation per leaf area                                                          |
| S-ci                  | Pa                       | Shaded canopy partial pressure of CO2                                                                                |
| S-PPFD                | $\mu$-mol s-1m-2         | Shaded photon flux density in sunlit canopy layer                                                                    |
| S-E-Tree              | mmol s-1 m-2 (LA)        | Weighted mean of transpiration in tree per leaf area                                                                 |
| Anet-tree             | $\mu$-mol s-1 m-2 (LA)   | Total carbon assimilation in tree per leaf area                                                                      |
| Pcrit                 | MPa                      | Critical pressure for tree                                                                                           |
| Ecrit                 | $\mu$-mol h-1 m-2 (LA)   | Transpiration at Pcrit per leaf area                                                                                 |
| P-leaf                | MPa                      | Composite leaf pressure                                                                                              |
| P-stem                | MPa                      | Composite stem pressure                                                                                              |
| P-root                | MPa                      | Composite root pressure                                                                                              |
| K-stem                | kg hr-1 m-2 (BA)         | Composite stem hydraulic conductance per basal area                                                                  |
| K-leaf                | kg hr-1 m-2 (BA)         | Composite leaf hydraulic conductance per basal area                                                                  |
| K-plant               | kg hr-1 m-2 (BA)         | Composite plant hydraulic conductance per basal area                                                                 |
| K-xylem               | kg hr-1 m-2 (BA)         | Total xylem hydraulic conductance per basal area                                                                     |
| K-root-1              | kg hr-1 m-2 (BA)         | Hydraulic conductance for root layer 1 per basal area                                                                |
| K-root-2              | kg hr-1 m-2 (BA)         | Hydraulic conductance for root layer 2 per basal area                                                                |
| K-root-3              | kg hr-1 m-2 (BA)         | Hydraulic conductance for root layer 3 per basal area                                                                |
| K-root-4              | kg hr-1 m-2 (BA)         | Hydraulic conductance for root layer 4 per basal area                                                                |
| K-root-5              | kg hr-1 m-2 (BA)         | Hydraulic conductance for root layer 5 per basal area                                                                |
| K-root-all            | kg hr-1 m-2 (BA)         | Total hydraulic conductance of roots per basal area at midday                                                        |
| E-root-1              | mmol s-1 m-2 (LA)        | Total root uptake for layer 1 per leaf area                                                                          |
| E-root-2              | mmol s-1 m-2 (LA)        | Total root uptake for layer 2 per leaf area                                                                          |
| E-root-3              | mmol s-1 m-2 (LA)        | Total root uptake for layer 3 per leaf area                                                                          |
| E-root-4              | mmol s-1 m-2 (LA)        | Total root uptake for layer 4 per leaf area                                                                          |
| E-root-5              | mmol s-1 m-2 (LA)        | Total root uptake for layer 5 per leaf area                                                                          |
| water-content         | mm                       | Root zone water content in mm (= m^3/m^2 (GA))                                                                       |
| water-content-delta   | mm timestep-1            | Change in water-content over previous time-step                                                                      |
| end-rain              | mm timestep-1            | Rain input per previous time-step                                                                                    |
| end-ground-water      | mm timestep-1            | Groundwater input per previous time-step                                                                             |
| end-E                 | mm timestep-1            | Transpiration per time-step                                                                                          |
| end-drainage          | mm timestep-1            | Total drainage per time-step                                                                                         |
| end-soil-evap         | mm timestep-1            | Evaporative water loss per time-step                                                                                 |
| end-ET                | mm timestep-1            | Total evaporation (E + soil-evap) per time-step                                                                      |
| end-Anet-la           | mmol timestep-1 m-2 (LA) | Net carbon assimilation per time-step per leaf area                                                                  |
| end-total-water-input | mm timestep-1            | Total water input (rain + ground-water) per time-step                                                                |
| end-PLC-plant         | %                        | Percent loss of conductivity for plant in previous time-step                                                         |
| end-PLC-xylem         | %                        | Percent loss of conductivity for xylem in previous time-step                                                         |
| end-runoff            | mm timestep-1            | Excess root zone water per time-step                                                                                 |

### Summary Outputs
- Total Anet (mmol yr-1 m-2(leaf area))
- Total E (mm = mm3/mm2(ground area))
- Minimum whole plant conductance during the growing season (kghr-1m-2)
- Percent Loss Conductance (PLC, percent, relative to a reference conductance at field capacity)
- Mean Ci/Ca (+ A weighted Ci/Ca) 
- Water summary (start/end content, total growing season input (mm))

------------

## Year Modes

### Independent year mode (default)

- Set "Use GS Data" to "n" under "Program Options"
- Use growing season trimmed data (see the example: "dataset.csv").
- Ensure that the growing season limits are defined in "seasonlimits.csv"

The default setting is to reset the tree hydraulics and reset the soil water content to the specified percent of field capacity every year. The years are completely independent, only run in a single dataset for convenience. In this mode, weather data should be trimmed to __only the growing season days__ as in the included dataset.csv (so that the last day of one growing season is followed immediately by the first day of the next). Year values are used for output, and each must be unique and optimally sequential (to determine when new years begin), but the values are otherwise unused by the model and thus do not need to be meaningful values. When running in this mode the growing season limits (seasonlimits.csv) will not be used; All days in the dataset will be considered to be in the growing season.

### Sequential year mode:

In this mode, plant hydraulics will reset between seasons and plant transpiration/productivity will be disabled during the off-season, but soil water budget will continue to be computed. Soil may or may not be refilled to field capacity depending on the availability of off-season precipitation.
- Set "Use GS Data" to "y" under "Program Options".
- Use full-year data (see the example: "dataset - full year example.csv").
- Ensure that the growing season limits are defined in "seasonlimits.csv".
  - __Note:__ The year values in "seasonlimits.csv" are for reference only; The first row of start/end days will be used for the first year of data, etc.
- Soil surface evaporation will also be disabled during the off-season. This is not particularly realistic, but the functionality was intended to answer the question: Is there at minimum enough recorded rainfall to refill the soil? A more robust off-season water simulation would require additional data (snow pack) and simulation of soil behavior under snow and is not provided here.

-------------

## References:

Describing the version 1.0 of the C++ code:
- Venturas MD, JS Sperry, DM Love, EH Frehner, MG Allred, Y Wang, and WRL Anderegg. (2017). A Stomatal Control Model Based on Optimization of Carbon Gain versus Hydraulic Risk Predicts Aspen Sapling Responses to Drought. New Phytologist 220: 836–50.

Describing the gain/risk algorithm used in the model:
- Sperry JS, MD Venturas, WRL Anderegg, M Mencucinni, DS Mackay, Y Wang, and DM Love. (2017). Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment 40: 816-830

Describing the original hydraulic model the gain-risk optimization was based on:
- Sperry JS, and DM Love (2015) Tansley Review: What plant hydraulics can tell us about plant responses to climate-change droughts. New Phytologist 207: 14-17.
- Sperry JS, Y Wang, BT Wolfe, DS Mackay, WRL Anderegg, NG McDowell, and WT Pockman. (2016). Pragmatic hydraulic theory predicts stomatal responses to climatic water deficits. New Phytologist 212: 577-589
-------------
