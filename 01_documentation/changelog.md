# Changelog
Changelog as of 04/16/25

## 04/16/25 (3.0.4)
* Adding command line argument for configuration selection
* Added `file_builder.py` that can be used to generate the model input files
* Added function headers to all functions
* Fixed summary output
* Removed data and summary data header from the configuration files, scrapping summary data header entirely, and hard coded them in `01Utilities.h`
  * Summary data was redundant since 100% of its outputs were being written from gs_data, so now we just output that instead
* Added direct `()` operator for accessing data stores derived from CSVData
* Modified data_template and removed some deprecated input parameters
* Added notebook for plotting the model outputs

## 02/05/25 (3.0.3)
* Finished porting over code
  * Adding remainder of functions, including canopypressure, soilflow, soilevaporation, etc...
  * Moved leaf and carbon functions and parameters into separate files
  * Added historical curve storage and restoration
* Fixed numerical instability in updatecurves when checking if the conductance in roots were less than the layer's minimum conductance
* Added model data output (save to file)
* Updated makefile
* Added command line arguments for parameter and configuration files

## 11/09/24 (3.0.2)
* Porting over code
  * Minor changes to getsoilwetness to preserve fieldcapacity and water data
  * Ported solarcalc, predawns, newtonrhapson, and compositecurve 

## 11/02/24 (3.0.1)
* Continued porting over code
  * componentpcrits() done
  * modelTimestepIter() beginning
  * getsoilwetness() done
* Added ability for CSVData to create new columns if doesn't already exist
* Changed dataheader to be easier to use

## 10/19/24 (3.0.0)
* Began refactoring code, converted to classes with distinct header files following OOP principles
* Changed storage of parameters to hashmap for more efficient access when running the model
* Created new CSVData with data accessible via column string lookup
* Changed dataheader columns to more usable formatting
* Started porting over componentpcrits()
  * Added E(P) (water supply / transpiration flow rate) approximations to Component class

## 09/26/24 (2.0.7)
* Added a constant seed for random to allow for easier reproducibility

## 09/25/24 (2.0.6)
* Added a clean target to Makefile
* Changed config and parameter file paths to header constants
* Changed file paths within config file to relative paths
* Added a .gitignore (need to get rid of .DS_Store)
* Added more debugging flags to Makefile
* **MAJOR**: Allocated modelProgramMain class within main.


