# Changelog
Changelog as of 02/05/25

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


