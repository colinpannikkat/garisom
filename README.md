# Carbon Gain vs Hydraulic Risk Stomatal Optimization Model V 2.0

__Coding language:__ C++

__Authors:__ German Vargas G. & William R.L. Anderegg

__Contact:__ german.vargas@utah.edu

------------

## Introduction:

The model uses a stomatal carbon gain vs. hydraulic risk optimization (Sperry et al. 2017 -- see references below) combined with a soil water budget to run a continuous growing season simulation, producing stand outputs such as net carbon assimilation (Anet), internal [CO2] (Ci), transpiration (E), total evapotranspiration (ET), and element conductances (k) on an hourly and summary basis (see output details below).

------------

## Basic usage instructions:

First, fork the repository to your own account. If you don't have a github account, just download the repository. Provide plant and site traits through the parameters file, then supply hourly weather data to drive model. The model will expect these files to be located in the current working directory. See the included examples and the details below for formatting and units. These example files contain all inputs necessary to test-run the model immediately after building.

__a)__ To run, build and execute the model program (no command line arguments -- runs from files in the working directory). Building requires c++11, GNU example:

```{}
g++ -std=c++11
```

The -O3 and -ffast-math optimizations are recommended with GNU compilers:

```{}
g++ -std=c++11 -O3 -ffast-math
```

__b)__ To build and run from the terminal in a Mac OS system:

__b.1)__ Clone the repository in a desired location in your system.

```
git clone https://github.com/gevargu/garisom.git
```

__b.2)__ Navigate to the repository by usin the command cd:

```
cd garisom
```

__b.3)__ Run the following command will compile the code and build an executable:

```
make
```

__b.4)__ Run this program from the same folder with this command:

```
./garisom
```

__b.5)__ Press <kbd>Command</kbd> + <kbd>C</kbd> if you want to stop the model before it completes (<kbd>Ctrl</kbd> + <kbd>C</kbd> in Windows/Linux)

This version has also been tested with Visual Studio 2017's compiler with similar optimizations (floating point mode fast, maximum optimization preferring speed). The following files (included in this repository) should be located in the working directory (normally the same directory as the executable) before running:
	
  - __parameters_2.0.0.csv__ (site, atmospheric, soils, stand, plant, hydraulics, and carbon assimilation parameters).
  - __configuration_2.0.0.csv__ (model controls)
  - __dataset.csv__ (hourly weather drivers).
  - __dataheader.csv__ (a header row for the hourly data output).
  - __sumheader.csv__ (a header row for the summary data output).
  - __seasonlimits_2.0.0.csv__ (growing season limits and yearly atmospheric CO2, only required if using "sequential year mode" described below)

Upon completion, two output files are produced:
	
  - The __timesteps.csv__ output contains all of the hourly model outputs corresponding to the weather inputs.
  - The __summary.csv__ output includes various total values for each year (for example, net growing season productivity Anet, net transpiration E, etc.)

------------

## Model Input Files

__Model Parameters:__

Configure plant traits and other parameters in __parameters_2.0.0.csv__ (expected input units are indicated). Coupled with the parameter build function in R to produce the __paramaters_2.0.0.csv__ file.

| Group    		| Parameter       	    	| Description										|
| ---------------------	| ----------------------------- | -------------------------------------------------------------------------------------	|
| Site			| __i_sp__			| Species or PFT represented in parameter data	|
| Site			| __i_region__			| Site/simulation ID. This ID is used for naming the result files that are exported	|
| Site			| __i_latitude__		| Latitude in degree fraction north.	|
| Site			| __i_longitude__		| Longitude in degree fraction west.	|
| Site			| __i_elevation__		| Site elevation in m above sea level.	|
| Site			| __i_slopeI__			| Slope inclination; degrees from horizontal.	|
| Site			| __i_slopeA__			| Slope aspect; counterclockwise degrees from south.	|
| Site			| __i_gWaterP__			| Ground water pressure.	|
| Site			| __i_gWaterDist__		| Distance to ground water source in m from the bottom of the rootsystem.	|
| Atmosphere		| __i_atmTrans__		| Atmospheric transmittance from weather data (set to 0.65 as default if no data available).	|
| Atmosphere		| __i_solarNoon__		| Solar noon correction from weather data in hours.	|
| Atmosphere		| __i_emiss__			| Long wave emissivity.	|
| Atmosphere		| __i_co2AmbPPM__		| Atmospheric/experiment CO2 ppm, it will update if working with multiple years.	|
| Soil			| __i_layers__			| Number of soil layers (select 1-5).	|
| Soil			| __i_fieldCapFrac__		| Fraction that field capacity is of saturation (minus residual).	|
| Soil			| __i_fieldCapPercInit__	| Percent field capacity for starting the season.	|
| Soil			| __i_rockFrac__		| Fraction of soil volume as rocks (0-1).	|
| Soil			| __i_soilAbsSol__		| Absorptivity of soil surface for solar.	|
| Soil			| __i_rhizoPer__		| Average percent of whole plant resistance in rhizosphere (maximum soil limitation)	|
| Soil			| __i_texture__			| USDA soil texture category (equal for all layers but could also be determined per layer)	|
| Stand			| __i_baperga__			| Basal area per ground area m2 ha-1	|
| Stand			| __i_leafAreaIndex__		| Canopy lai (m2 m-2)	|
| Stand			| __i_soilXHeight__		| Height above soil surface for understory wind and gh in m	|
| Stand			| __i_height__			| Average tree height in m	|
| Stand			| __i_treeToPhotoLAI__		|	|	
| Stand			| __i_leafPerBasal__		| Initial leaf area per basal area per individual tree; m2 m-2	|
| Tree			| __i_leafWidth__		| Leaf width in m	|
| Tree			| __i_leafAngleParam__		| Leaf angle parameter; CN 15.4	|
| Tree			| __i_aspect__			| Max radius of root system per max depth	|
| Tree			| __i_rootDepth__		| Maximum rooting depth in m	|
| Hydraulics		| __i_leafPercRes__		| Saturated % of tree resistance in leaves	|
| Hydraulics		| __i_kmaxTree__		| Kmax of tree in kg hr-1 m-2 MPa-1 per basal area	|
| Hydraulics		| __i_pinc__			| Pressure increment for curve generation, (MPa) - higher is faster, but less accurate (setting too high can cause Newton-Rhapson root pressure solving failure)	|
| Hydraulics		| __i_rootP12__			| Root element p12	|
| Hydraulics		| __i_rootP50__			| Root element p50	|
| Hydraulics		| __i_stemP12__			| Stem p12	|
| Hydraulics		| __i_stemP50__			| Stem p50	|
| Hydraulics		| __i_leafP12__			| Leaf p12	|
| Hydraulics		| __i_leafP50__			| Leaf p50	|
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
| BAGA_optimizer	| __i_iter_gwInc__		|	|
| BAGA_optimizer	| __i_iter_gwStart__		|	|
| BAGA_optimizer	| __i_iter_gwEnd__		|	|
| BAGA_optimizer	| __i_iter_ffcInc__		| Note: If FFC start < FFC end_ will start curve gen by incrementing FFC before ground water	|
| BAGA_optimizer	| __i_iter_ffcStart__		|	|
| BAGA_optimizer	| __i_iter_ffcEnd__		|	|
| BAGA_optimizer	| __i_iter_bagaInc__		|	|
| BAGA_optimizer	| __i_iter_bagaStart__		|	|
| BAGA_optimizer	| __i_iter_bagaEnd__		|	|
| BAGA_optimizer	| __i_iter_bagaRef__		|	|
| BAGA_optimizer	| __i_iter_bagaCutoff__		| WLT K dropoff threshold (fraction of reference iteration kmin)	|

__Model Configuratiosn:__

| Group    	| Model Control       	    | Description          |
| ------------- | ------------------------- | -------------------- |
| Soil     	| __igWaterEnable__   	    | Turns on/off groundwater flow. Values: n (off); y (on). It provides an unlimited source of water at a set potential and distance below the root layers. This water will flow up into the soil layers, and potentially allow layers to fill above field capacity (from the bottom layer up). When disabled (default), the only sources of water input will be the initial fraction of field capacity and observed rainfall (and any water over field capacity will become "drainage").	|
| Soil     	| __i_soilRedEnable__  	    | Turns on/off soil redistribution routine. Values: n (off); y (on). It allows water to flow between soil layers.	|
| Soil     	| __i_soilEvapEnable__ 	    | Turns on/off soil evaporation routine. Values: n (off); y (on). It enables simulation of water evaporation from the surface soil layer.	|
| Climate  	| __i_rainEnable__  	    | Turns on/off rain inputs. Values: n (off); y (on). It allows for precipitation events. Weather data rainfall will be ignored if disabled.	|
| Climate  	| __i_useGSDataStress__     | Turns on/off growing season data for multiple year modeling. Vakyes: n (off); y (on). If enabled, multiple years will be run "sequentially" with on and off seasons defined in seasonlimits_2.0.0.csv and a continuous water budget. When disabled (default), all weather timesteps provided are treated as part of the growing season and the user is expected to truncate individual years to their start/end days. Water budget is reset between years when disabled, treating years as totally independent.	|
| Climate	| __i_useGSDataOpt__ 	    | Turns on/off growing season data for multiple year modeling during BAGA optimization. Values: n (off); y (on).	|
| Hydraulics  	| __i_refilling__ 	    | Turns on/off xylem refilling within a growing season. Values: n (off); y (on). It allows trees to restore lost conductance, however the refilling model is not sufficient to simulate authentic xylem refilling behavior and has __not been thoroughly tested in the current version of the code__.	|
| Hydraulics  	| __i_predawnsMode__ 	    | Turns on/off if model should consider measured pre-dawn water potential values. Values: n (off); y (on). If set to 'y', disables soil simulation and runs from hourly inputs of canopy predawn water potential. These are read from the "rain" column (rain is not used with the soil sim disabled) in MPa. See "dataset - predawns example.csv". This mode is especially useful when comparing to other models which run from canopy predawn measurements, and generally runs significantly faster as it does not need to solve for root layer pressures.	|
| Hydraulics  	| __i_cavitFatigue__ 	    | Turns on/off xylem stress hysteresis to carry effects from previous growing season. Values: n (off); y (on). It allows for a weighted estimation of xylem vulnerability to embolism	|
| Hydraulics  	| __i_stemOnly__ 	    | Turns on/off xylem stress hysteresis only in stem xylem. Values: n (off); y (on). When disabled it allows for a weighted estimation of xylem vulnerability to embolism for both stem and roots	|
| BAGA  	| __i_iter_gwEnable__ 	    | __Not tested in version 2.0.1__.	|
| BAGA  	| __i_iter_ffcEnable__ 	    | __Not tested in version 2.0.1__.	|
| BAGA  	| __i_iter_bagaEnable__     | Increate the BA:GA to find the basal area that puts the stand in ecohydrological equilibrium with the weather conditions. Values: n (off); y (on). __Not tested in version 2.0.1__.	|
| BAGA  	| __i_iter_useAreaTable__   | if y will pull GA:BA_ LA:BA_ LAI from AreaData table per year and per site. Values: n (off); y (on). __Not tested in version 2.0.1__.	|
| BAGA  	| __i_iter_yearsAsCount__   | The "year" values represent different data set ID's_ not actual years. In this mode the "start" year is always 0. Values: n (off); y (on). __Not tested in version 2.0.1__.	|
| BAGA  	| __i_iter_runSupplyCurve__ | Turns on/off all iteration settings. Values: n (off); y (on). __Not tested in version 2.0.1__. 	|
| Community	| __i_multipleSP__	    | Turns on/off whether our model configuration has 1 species per site (monodominant) or multiple species per site (diverse). Values: n (off); y (on). __NOT READY FOR MULTIPLE SPECIES/PFTs in version 2.0.1__.	|
| Community	| __i_speciesN__	    | Nnumber of species/PFT to run the model.	|
| Forcing files	| __i_ClimateData__	    | Path to file with climate forcing variables __dataset.csv__	|
| Forcing files	| __i_GSData__		    | Path to file with growing season data __seasonlimits.2.0.0.csv__	|

------------

## Weather Data:

See example data for formatting. Weather drivers should be in hourly timesteps and can include multiple years of data. Note that while year values are arbitrary, they should be sequential. For example, if running data for the years 1997 and 2005 these should be numbered sequentially as 1 and 2 (or 1997 and 1998, etc). Inputs:

- Year.
- Julian Day (1-366).
- Hour (0-23).
- Obs. Solar (W m-2).
- Rain (mm).
- Wind (m s-1).
- Tair (C).
- Tsoil (C, if not available substitute air temp).
- D (kPa).

------------

## Outputs:

- Autosave is always enabled regardless of the setting, as this version of the model has no alternative output method. Output files will be generated in the working directory when the run completes.

-Hourly Outputs (see dataheader.csv for full list):
	-Pressures (predawn soil layer pressures, sun and shade "mid-day" canopy pressures, MPa), 
	-Water flows (mmol m-2s-1), 
	-PS assimilation (A, umol s-1m-2 (leaf area)), 
	-Gain-risk optimized stomatal conductance to water (Gw, mmol m-2s-1),
	-Element and whole plant conductances, (k, kghr-1m-2), 
	-Water content and deltas (mm), 
	-Ci

-Summary Outputs (per year, see sumheader.csv for full list):
	-Total Anet (mmol yr-1 m-2(leaf area)),
	-Total E (mm = mm3/mm2(ground area)), 
	-Minimum whole plant conductance during the growing season (kghr-1m-2), 
	-Percent Loss Conductance (PLC, percent, relative to a reference conductance at field capacity),
	-Mean Ci/Ca (+ A weighted Ci/Ca), 
	-Water summary (start/end content, total growing season input (mm)).

------------

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
