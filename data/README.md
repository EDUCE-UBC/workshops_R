# Workshop data
Data sets for use in data science workshops.

## Files
### `data.csv`
Simplified data from a subset of the Saanich geochemical data from cruise 72 including season, depth in meters, oxygen (O2)  concentration in micromolar, and an indication if accompanying microbial data exists (Add_data).

These data are used in the workshops:

* [Introduction to R - 2 hour](https://github.com/EDUCE-UBC/workshops_data_science/tree/master/intro_R_2hr)
* [Reproducible research](https://github.com/EDUCE-UBC/workshops_data_science/tree/master/reproducible_research)

### `Saanich_Data.csv`
Full geochemical data from Saanich Inlet over time including longitude, latitude, cruise number, date, depth in km, cells per ml, temperature in Celsius, salinity, and density as well as mean and standard deviation of micromolar concentrations of:

* oxygen (O2)
* phosphate (PO4)
* nitrate (NO3), ammonium (NH4), nitrite (NO2), nitrogen (N2), nitrous oxide (N2O), 
* hydrogen sulfide (H2S),
* carbon dioxide (CO2), methane (CH4)

These data are used in the workshops:

* [The R tidyverse](https://github.com/EDUCE-UBC/workshops_data_science/tree/master/intro_tidyverse)
* [Intermediate R programming](https://github.com/EDUCE-UBC/workshops_data_science/tree/master/intermediate_R)


## The data
These data were collected as part of an on-going oceanographic time series program in Saanich Inlet, a seasonally anoxic fjord on the East coast of Vancouver Island, British Columbia (Figure 1). We use data from various geochemical measurements at many depths in Saanich Inlet over time (approximately monthly from 2006 to 2014) as well as microbial data from selected time frames and depths of interest.

![](https://github.com/EDUCE-UBC/workshops_data_science/blob/master/intro_tidyverse/images/Saanich.png){width=4in}

**Figure 1.** Map of Saanich Inlet indicating conventional sample collection stations (S1-S9). Data used in this workshop is sourced from S3.

Saanich Inlet is a steep sided fjord characterized by a shallow glacial sill located at the mouth of the inlet that restricts circulation in basin waters below 100 m (Figure 2).

![](https://github.com/EDUCE-UBC/workshops_data_science/blob/master/intro_tidyverse/images/Inlet_structure.png)

**Figure 2.** Structure of Saanich Inlet. The glacial sill restricts water circulation into and out of the lower depth of the inlet basin.

During spring and summer months, elevated primary production (like photosynthesis) in surface waters combined with restricted circulation results in progressive water column stratification and complete oxygen starvation (anoxia) in deep basin waters. In late summer, pulses of oxygenated nutrient-rich ocean waters upwelling from the Haro Straight cascade over the sill, displacing oxygen starved bottom waters upward. The intensity of these renewal events varies from year to year with implications for microbial ecology and biogeochemical cycles (Figure 3). 

![](https://github.com/EDUCE-UBC/workshops_data_science/blob/master/intro_tidyverse/images/oxygen_timeseries.png)

**Figure 3.** Contour plot of water column oxygen concentrations over multiple years in the time series. Warmer colors indicate high oxygen concentrations while cooler colors are low. Note the recurring pattern of oxygen decline below 100 m depth intervals followed by seasonal renewal events in late Summer into early Fall carrying more oxygenated waters into the Inlet. 

The seasonal cycle of stratification and deep water renewal enables spatial and temporal profiling across a wide range of water column energy states and nutrients, thus making Saanich Inlet a model ecosystem for studying microbial community responses to ocean deoxygenation. Ocean deoxygenation is a widespread phenomenon currently increasing due to climate change. 

For a brief introduction to the data, see Hallam SJ *et al*. 2017. Monitoring microbial responses to ocean deoxygenation in a model oxygen minimum zone. Sci Data 4: 170158 [doi:10.1038/sdata.2017.158](https://www.nature.com/articles/sdata2017158) More detailed information on the geochmeical data can be found in Torres-Beltrán M *et al*. 2017. A compendium of geochemical information from the Saanich Inlet water column. Sci Data 4: 170159. [doi:10.1038/sdata.2017.159](https://www.nature.com/articles/sdata2017159) More detailed information on the mutli-omic microbial data can be found in Hawley AK *et al*. 2017. A compendium of multi-omic sequence information from the Saanich Inlet water column. Sci Data 4: 170160. [doi:10.1038/sdata.2017.160](https://www.nature.com/articles/sdata2017160)
