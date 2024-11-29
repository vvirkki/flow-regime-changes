## flow-regime-changes

### Code repository for

### "Archetypal flow regime change classes and their associations with anthropogenic drivers of global streamflow alterations"

Vili Virkki, Reetik Kumar Sahu, Mikhail Smilovic,
Josias Láng-Ritter, Miina Porkka, Matti Kummu

**Published in _Environmental Research Communications_**  
[https://doi.org/10.1088/2515-7620/ad9439](https://doi.org/10.1088/2515-7620/ad9439)

**Supporting data is available in a Zenodo repository**  
[https://doi.org/10.5281/zenodo.11102422](https://doi.org/10.5281/zenodo.11102422)

**Corresponding authors of the article**  
Vili Virkki (vili.virkki@uef.fi)  
Matti Kummu (matti.kummu@aalto.fi)

**Repository author**  
Vili Virkki (vili.virkki@uef.fi)

**Please cite the _Environmental Research Communications_ publication if
using code from this repository or data from the Zenodo repository in
another publication.**

### Folders and files in repository

#### Data/

Raw data are not provided with this repository, except for a land-sea mask used
for distinguishing land cells from sea cells when using ERA5-Land data.

Raw data can be acquired from the
[ERA5-Land repository](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-monthly-means?tab=overview),
the GSIM repositories
([metadata](https://doi.pangaea.de/10.1594/PANGAEA.887477), 
[indices](https://doi.pangaea.de/10.1594/PANGAEA.887470)),
[the Huang et al. (2018) repository](https://zenodo.org/records/1209296),
[the GRDC repository](https://portal.grdc.bafg.de/applications/public.html?publicuser=PublicUser#dataDownload/StationCatalogue), and the
[GeoDAR repository](https://zenodo.org/records/6163413). Dam attributes
merged with GeoDAR records are proprietary to ICOLD and thus their redistribution
is prohibited (see GeoDAR repository for details).

The folder tree of the complete `Data` folder prior to analysis is available
in the `txt` folder.

#### R/

The R scripts are run in sequence from `01_` to `08_` and `res01_` to `res02_`.
Given a complete `Data` folder, the scripts will create an `output` folder and
collect intermediate and final outputs there.

`01_` to `08_`: sample catchments, pre-process streamflow and driver data,
determine streamflow and driver trends in sampled catchments.

`res01_`: assign flow regime change (FRC) classes and use them to aggregate
driver trends, produce outputs for figures.

`res02_`: export csv files for release in the Zenodo repository.

#### txt/

`data_folder.txt`: how the `Data` folder looks like prior to running the full
analysis. Individual GSIM catchment files are not shown.

`output_folder.txt`: how the `output` folder, created during the analysis, looks
like after running the full analysis.

`sessionInfo.txt`: R and package versions used for running the full analysis.

#### Other files

`params`: set of parameters for analysis; read in the beginning of `01_`, saved
explicitly within `output` folder and referred to throughout the analysis.

`run_analysis.sh`: bash script to execute the full analysis at once; silently
overwrites existing outputs.

### Changelog and versioning

Should the repository be updated with new versions, main changes will be briefly
summarised here. Commits describing new versions are marked respectively with
Git tags.

#### v1.0.0.

Version used in producing the results shown in the
[_Environmental Research Communications_ manuscript](https://doi.org/10.1088/2515-7620/ad9439)

### License

Attribution 4.0 International (CC BY 4.0)  
https://creativecommons.org/licenses/by/4.0/

