# tres-pheno-sensitivity

A version of this repository is permanently archived at Zenodo.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17804486.svg)](https://doi.org/10.5281/zenodo.17804486)

***Goal***

This repository accompanies the manuscript "Divergent climate impacts despite similar response to temperature in a widespread aerial insectivore".

The analysis includes nesting data from tree swallows contributed by groups monitoring populations across the United States and Canada, combined with citizen science monitoring data and associated climate data. The code is organized and heavily annotated to facilitate reproducibility and reuse. In order to fully execute the code, publicly available data must be downloaded from the sources indicated in the code comments.

For questions please contact Conor Taff (cct663@gmail.com). Last updated December 3, 2025.

***Publicly available data***

- NestWatch project wide download: https://nestwatch.org/explore-data/nestwatch-open-dataset-downloads/

  Download data from project website and add to repository as indicated.

- Berkeley Earth historical climate data: http://berkeleyearth.org/data/

  Download data from project website and add to repository as indicated.
  
- NALCMS land cover data: https://www.mrlc.gov/

  Download data from project website and add to repository as indicated.
  
- eBird data (for plotting only)
  
  Downloaded inline by code, but you will need to register for a free access key and add it as indicated.
  
- Daymet & MERRA-2 climate data

  Downloaded inline by code, no action needed.
  
***Data files included with repository***

- **data_by_nest.txt**: Compiled nesting data for tree swallows with one row per observed nest for all included study sites other than NestWatch data.

- **data_by_population.txt**: Metadata about the professionally monitored sites with one row per site. Used mainly for simple plotting with some notes on each population.

- **bbs_trends_80_22.txt**: Breeding Bird Survey trend data and confidence intervals for the areas included in this study.

- **tree_swallow_2024.txt**: Arrival dates for each grid cell derived from eBird raw data as described in the manuscript.

***Detailed description of columns***

Below is a description of each column for the data files included in the repository. For the additional publicly available data, please refer to the original sources for documentation.

**data_by_nest.txt**
- *pop_id*: Factor indicating what population the nest came from.
- *sub_pop*: Factor indicating sub-population within a population, if applicable.
- *site_identifier*: Factor indicating specific site within a population, if applicable.
- *nest_identifier*: Factor indicating unique nest within a site. Note that nest identifiers are often reused across sites, so must be combined with site/population info to be unique.
- *attempt_num*: Nesting attempt number in this box/year if known.
- *year*: Year of observation.
- *first_egg_date*: Date of first laid egg in day of year (January 1 = 1).
- *est_first_egg_date*: Estimated date of first laid egg in day of year (January 1 = 1) if exact date not known.
- *egg_check_interval*: Number of days between egg checks.
- *egg_check_int_numeric*: Numeric version of number of days between egg checks.
- *clutch_size*: Clutch size if known.
- *hatch_date*: Hatch date in day of year.
- *est_hatch_date*: Estimated hatch date in day of year if exact date not known.
- *hatch_check_interval*: Number of days between hatch checks.
- *hatch_check_int_numeric*: Numeric version of number of days between hatch checks.
- *simple_fate*: Fate as failed or fledged if known.
- *full_fate*: More detailed fate if known.
- *number_hatched*: Number of eggs hatched if known.
- *number_fledged*: Number of chicks fledged if known.
- *exact_latitude*: Exact latitude of this nest box.
- *exact_longitude*: Exact longitude of this nest box.
- *approx_latitude*: Approximate latitude of this nest box if exact not known (often at site/field level).
- *approx_longitude*: Approximate longitude of this nest box if exact not known (often at site/field level).
- *laying_treatment*: Details on any manipulations that could have influenced laying date.
- *morph_treatment*: Details on any manipulations that could have influenced morphological measures.
- *experiment*: Any experimental groups that this nest was a part of.
- *exclude_fitness*: Should fitness measures from this nest be excluded (i.e., because of manipulation).
- *additional_info*: Open form with any additional information or notes.
- *number_d12*: Number of nestlings alive on day 12 if known.
- *old_site_identifier*: Original site identifier before cleanup for standardization.
- *old_nest_identifier*: Original nest identifier before cleanup for standardization.

**data_by_population.txt**
- *pop_id*: Factor indicating what population the row describes.
- *pop_description*: Short description of the population/site.
- *location*: General location of the population (e.g., nearest town).
- *pop_contact*: Contact person for the population/site.
- *start_yr*: Starting year of data collection.
- *end_yr*: Ending year of data collection.
- *total_years*: Total number of years with data.
- *total_nests*: Total number of nests observed.
- *ad_morph_tot*: Total number of adult morphological measurements.
- *ad_morph_uni*: Number of unique adults with morphological measurements.
- *chick_morph_tot*: Total number of chick morphological measurements.
- *chick_morph_uni*: Number of unique chicks with morphological measurements.
- *avg_latitude*: Average latitude of the population.
- *avg_longitude*: Average longitude of the population.
- *in_nestwatch*: Is this population part of NestWatch data?
- *use_in_sensitivity*: Should this population be used in the main sensitivity analysis?
- *sub_pops*: Are there any sub populations within this population?
- *box_descriptions*: Descriptions of nest boxes used in this population.
- *notes*: Any additional notes about the population.
- *pop_coauthors*: Names of coauthors associated with this population.
- *num_authors*: Number of coauthors associated with this population.

**bbs_trends_80_22.txt**
- *state_province*: State or province.
- *trend_80_22*: BBS population trend estimate.
- *low_ci*: Lower confidence interval for trend.
- *high_ci*: Upper confidence interval for trend.
- *center_latitude*: Center latitude of state or province.

**tree_swallow_2024.csv**
- *species*: Indicates species, in this case all are tree swallows.
- *cell*: Cell code used to join to populations.
- *year*: Year of observation.
- *elev_factor*: Factor indicating elevation band (see manuscript).
- *arr_GAM_mean*: Estimated arrival date for cell/year in day of year.
- *arr_GAM_sd*: Standard deviation of estimated arrival date for cell/year.