# FlexiGIS_Philadelphia

Open source GIS-based plattform for the configurations of renewable power generation in cities: Philadelphia as a case study

The code remix_calculator.py is developed for Python 3. The script is published under the open source license 3-Clause BSD (Berkeley Standard Distribution) https://opensource.org/licenses/BSD-3-Clause.

Normalised load and renewable generation data are stored as an hourly timeseries in individual CSV files that are read by the provided python code. The load file provided is based on a simulated residential building in Philadelphia, as PJM's license does not allow for redistribution of the data that was used to generate the figures in the paper.

* `norm_load.csv`
* `norm_solar.csv`
* `norm_wind.csv`

The required balancing energy from storage is simulated for different scenarios of wind and solar mixtures, as well as the overall renewable energy fraction.
