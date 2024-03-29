[![DOI](https://zenodo.org/badge/173319564.svg)](https://zenodo.org/badge/latestdoi/173319564)

### A repository for:

Bates, J. M., Fidino, M., Nowak-Boyd, L., Strausberger, B. M., Schmidt, K. A., and Whelan, C. J. (in review). Climate change affects nesting phenology of midwestern birds: comparison of modern field records with historical records obtained from museum collections.

<div align="center"><img width="150" height="auto" src="./images/american_robin.png" alt="A drawing of a robin that Mason made." /></div>

<div align="center"> <h3>Scripts</h3> </div>

---

**This repository has 3 R Scripts used for the analysis. They will be in your working directory. These include:**

**calc_climate_residuals.R:** This script reads in the atmospheric CO2 data, fits a linear model to it (with year as the independent variable), and calculates the residuals from the model.
These residuals are then used in our primary analysis.

**analysis_script.R:** This script fits our robust to outlier model to bird nesting records.

**plotting.R:** This script can be used to generate the figures in the manuscript.

To conduct the analysis the scripts should be ran in the order above.

---

<div align="center"><img width="150" height="auto" src="./images/shrike.png" alt="A drawing of a shrike that Mason made." /></div>

<div align="center"> <h3>Models</h3> </div>

---

**This repository has 2 JAGS models that we used for our analysis. They should be placed within the jags_models sub-folder of the working directory. These include:**

**climate_resid_model.R:** This is the model that is called by `bates_2017_calc_climate_residuals.R`. 

**robust_t_model.R:** This is the model that is called by `bates_2017_analysis_script.R`. 

---

<div align="center"><img width="150" height="auto" src="./images/blue_jay.png" alt="A drawing of a blue jay that Mason made." /></div>

<div align="center"> <h3>Data</h3> </div>

---

**There are three data files within the data sub-folder which are used in this analysis. They include:**

**co2.csv:** This is the global atmospheric CO2 levels per year. 

| Column header | Data type | Description |
|---|---|---|
| `yr`| Integer | The year the global atmospheric CO2 level is associated to. Ranges from 1744 to 2015. |
| `co2` | Numeric | The global atmospheric CO2 level on a given year. |

Between 1744 and 1953, global CO2 levels were compiled from ice cores collected at Siple Station, West Antarctica (Neftel et al. 1994). 
For 1958 to 2015, direct observations of atmospheric CO2 levels were collected from the Mauna Loa Observatory (Keeling et al. 2008).

Keeling, RF, Piper, SC, Bollenbacher, AF, Walker, JS. 2008 Atmospheric CO2 Records from Sites in the Scripps Institution of Oceanography [SIO] Air Sampling Network [1985-2007]. Oak Ridge, TN (USA): Carbon Dioxide Information Analysis Center (CDIAC), Oak Ridge National Laboratory (ORNL).

Neftel, A, Friedli, H, Moor, E, Lötscher, H, Oeschger, H, Siegenthaler, U, Stauffer, B. 1994
Historical CO2 record from the Siple Station ice core. In Trends: A Compendium of Data on Global Change. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, Tenn., U.S.A., U.S. Department of Energy.


<br>
<br>

**migratory_status.csv:** These data relate a species to a specific migratory group as well as its American Ornithological Union (AOU) 4-letter alpha code.

| Column header | Data type | Description |
|---|---|---|
| `cmn` | Character | The common name of a given bird species. Species names are lowercase. |
| `migstat` | Character | The migratory status of a species. `long` indicates long-distance migrants (species that spend the non-breeding season primarily in the subtropics/tropics south of the United States border). `short` indicates short-distance migrants (species that spend the non-breeding season in southern temperate regions of the southern U.S.), and `resident` indicates permanant residents (species that maintain most of their populations in the study region throughout the year). |
 |`species`| Character | The AOU 4-letter alpha code of a species. |
 
 Migatory status for all species was compiled from https://www.allaboutbirds.org/
 
<br>
<br>
 
 **bird_lay_dates.csv:** This is the lay date information for midwestern birds that were used in this analysis.
 
 | Column header | Data type | Description |
|---|---|---|
|`species`| Character | The AOU 4-letter alpha code of a species. |
|`jdate`| Integer | The julian date of the first lay date of a nest (Number of days from January 1 on a given year) |
|`year` | Integer | The year the nest was found. |
| `period` | Categorical | `low` indicates historic records housed at the Field museum of egg collections. `high` indicates current records of nest phenology colleted through comprehensive field work by Chris Whelan and Bill Strausburger.|
 


<div align="center"><img width="150" height="auto" src="./images/rock_dove.png" alt="A drawing of a pigeon that Mason made." /></div>


<div align="center"> <h3>images</h3> </div>

---

This folder houses the bird drawings, which were done by Mason Fidino. They are just used in this README file to add a little bit of fun to it.
