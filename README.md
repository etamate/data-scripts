**_Code samples of netCDF file handling_** - I've collected some useful code that I use often when working with netCDF climate data aggregations.

### converting SEDAC GPW v3 files to netCDF format
Whenever you want to create a country level aggregation of some climate indexes, you will need a gridded country map to decide which country a given cell belongs to.
Besides, you often want to weight your calculated index with the corresponding population.

E.g. when I create heating-degree-days indexes I usually weight the result with the population - this can radically change the result in big countries like Canada, China or Russia.

In the [population/sedac_gpwv3](population/sedac_gpwv3) folder I use some python code to create netCDF files from the [SEDAC Gridded Population](http://sedac.ciesin.columbia.edu/data/collection/gpw-v3) dataset, which I use later to create regional aggregations and population weighted indexes.

[Read more here about the example](population/sedac_gpwv3/README.md)

