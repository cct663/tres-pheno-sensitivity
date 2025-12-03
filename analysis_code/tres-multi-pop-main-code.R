## Code for analysis of tree swallow nest and morphology data from many populations
    # See project README for details
    # Code written by Conor Taff; cct663@gmail.com
    # Note that in order to run code you will need to download publicly available
      # datasets and place them in the correct directories as noted below

# Load libraries ----
    pacman::p_load(here, tidyverse, sf, rnaturalearth, ebirdst, terra, daymetr,
                   lubridate, sjPlot, climwin, dbscan, viridis, MODIStsp, dplyr,
                   geosphere, elevatr, exactextractr, MODIS, climateR, nasapower,
                   lme4, dggridR, mgcv, ggeffects)

    # access key for ebird. user will need to get and enter your own free code
      #ebird_key <- "########" # replace with your own key

# Read in population data ----
  # data organized with one nest per row
        d_nest <- read.delim(here::here("raw_data", "data_by_nest.txt"))
      # this year and site has a single late first egg date that is wrong. removed here
        d_nest$first_egg_date[d_nest$pop_id == "tsp_011" & d_nest$sub_pop == "B" & d_nest$year == 2017] <- NA
        
  # data organized with one row per population, for summary plotting only
        d_pop <- read.delim(here::here("raw_data", "data_by_population.txt"))
        
  # read in ebird arrival dates
        grid_match <- read.csv(here::here("ebird_arrival_info", "all_pop_locations_cell_elev.csv"))
        arrival <- read.csv(here::here("ebird_arrival_info", "tree_swallow_2024.csv"))

# Create some summarized datasets ----        
  # make new lat/lon. Use exact locations if possible or approximate locations if exact isnt' available
        d_nest <- d_nest %>%
          dplyr::mutate(
            latitude = as.numeric(ifelse(!is.na(exact_latitude) & exact_latitude != "", exact_latitude, approx_latitude)),
            longitude = as.numeric(ifelse(!is.na(exact_longitude) & exact_longitude != "", exact_longitude, approx_longitude))
          )
        
  # make new egg date column using exact or estimated first egg date
        d_nest <- d_nest %>%
          dplyr::mutate(
            lay_doy = as.numeric(ifelse(!is.na(first_egg_date) & first_egg_date != "", first_egg_date, est_first_egg_date))
          )
        
  # create a numeric number fledged column
        d_nest$number_fledged_num <- as.numeric(d_nest$number_fledged)
        
  # summarize the nest file to unique sub-pops and code by number of years and total nests
        d_nest <- d_nest %>% filter(lay_doy > 70, lay_doy < 225)
        dn_sum <- d_nest %>%
          dplyr::group_by(pop_id, sub_pop) %>%
          dplyr::summarise(tot_nests = n(), n_years = length(unique(year)), 
                    avg_lat = mean(latitude, na.rm = TRUE), avg_lon = mean(longitude, na.rm = TRUE),
                    first_year = min(year), last_year = max(year),
                    mu_lay = mean(lay_doy, na.rm = TRUE), md_lay = median(lay_doy, na.rm = TRUE),
                    mu_fledged = mean(number_fledged_num, na.rm = TRUE))
           
# Load and clean nestwatch data ----
    # this is the project wide data issued in February 2025 through 2024 breeding season
    # In order to run this section you will need to download the project wide data from the 
        # nestwatch website and add it to the nestwatch_dat directory for this project
        
        
    # to reproduce code please download the full nestwatch data from their website and link here
          n_dat <- read.csv(here::here("nestwatch_dat", "nestwatch_data_2024.csv"))
    # filter to tree swallows and drop some columns
          n_dat <- n_dat %>% filter(Species.Code == "treswa") %>%
            dplyr::select(-Substrate.Relationship, -Substrate, -Substrate.Other.Description,
                   -Cavity.Entrance.Diameter.cm, -Habitat.1m, -Habitat.100m, -Location.Entry.Technique,
                   -Predator.Guard, -Predator.Guard.Other, -Attempt.Entry.Technique)
          
    # make better column names
          colnames(n_dat) <- c("attempt_id", "loc_id", "latitude", "longitude", "subnat_code",
                               "height_m", "ent_orientation", "elevation_m", "observer_id",
                               "species_code", "species_name", "year", "ci_date", "ci_date_estimate",
                               "visited_during_lay", "hatch_date", "hatch_date_estimate", "fledge_date",
                               "fledge_date_estimate", "young_fledged", "clutch_size", "young_total",
                               "unhatched_eggs", "outcome")
    
    # wrangle and clean up
          n_dat$lay_dt <- as.POSIXct(n_dat$ci_date, format = "%Y-%m-%d")
          n_dat$hatch_dt <- as.POSIXct(n_dat$hatch_date, format = "%Y-%m-%d")
          n_dat$fledge_dt <- as.POSIXct(n_dat$fledge_date, format = "%Y-%m-%d")
          n_dat$lay_doy <- yday(n_dat$lay_dt)
          n_dat$hatch_doy <- yday(n_dat$hatch_dt)
          n_dat$fledge_doy <- yday(n_dat$fledge_dt)
          n_dat$year <- as.numeric(n_dat$year)
          n_dat <- n_dat %>% filter(year > 1975)
          
    # filter out nests with impossible dates or no date info
        # for lay date analyses I'm not worrying about if the clutch/fledge numbers are impossible
        # as long as we have lay date information
          n_dat <- n_dat %>% 
            filter(lay_doy > 70, lay_doy < 225, is.na(lay_doy) == FALSE)
          n_dat <- n_dat %>% filter(attempt_id != "A1051387") # this is a bad record
          
# Cluster nestwatch to pseudo populations ----
    # define coordinates of tree swallow nests in nestwatch
          coords <- n_dat[, c("latitude", "longitude")]
          coords_matrix <- as.matrix(coords)          

    # using density based spatial clustering to identify groups (dbscan)
        # eps is the radius of clusters, this is about 5-6km diameter
          tres_dbs <- dbscan(coords_matrix, eps = 0.03, minPts = 100)
          n_dat$cluster <- tres_dbs$cluster
          
    # determine stats for the clusters
          cluster_stats <- n_dat %>%
            filter(cluster != 0) %>%      # these ones aren't assigned to a cluster
            group_by(cluster, year) %>%   
            summarise(nyr = n()) %>%
            filter(nyr > 14) %>%           # minimum of 10 nests per year
            group_by(cluster) %>%
            summarise(years = n()) %>%
            filter(years > 7) %>%         # minimum of 8 years at same cluster
            as.data.frame()
          #nrow(cluster_stats)
    # subset to only nests in clusters
          n_dat2 <- subset(n_dat, n_dat$cluster %in% cluster_stats$cluster)
          # something is wrong with monitoring at cluster 91 in 2010. removing
            n_dat2 <- n_dat2 %>% filter(!(cluster == 91 & year == 2010))
    
    # summarize to cluster level
          n_dat3 <- n_dat2 %>%
            group_by(cluster) %>%
            summarise(n_nests = n(), n_years = length(unique(year)),
                      latitude = mean(latitude, na.rm = TRUE), longitude = mean(longitude, na.rm = TRUE),
                      mu_lay = mean(lay_doy), md_lay = median(lay_doy),
                      first_year = min(year), last_year = max(year),
                      mu_fledged = mean(young_fledged))
          
          n_dat4 <- n_dat2 %>%
            group_by(cluster, year) %>%
            summarise(n_nests = n(),
                      mu_lay = mean(lay_doy), md_lay = median(lay_doy))
          
          
    # figure out the nearest professionally monitored population distance
        coords_n <- as.matrix(n_dat3[, c("longitude", "latitude")])
        coords_dn <- as.matrix(dn_sum[, c("avg_lon", "avg_lat")])
        
        # distance matrix
            dist_matrix <- distm(coords_n, coords_dn, fun = distHaversine)
            
        # minimum distance to monitored pop from nestwatch clusters
            n_dat3$min_dist_km <- apply(dist_matrix, 1, function(x){
              if(all(is.na(x))) NA else min(x, na.rm = TRUE)
            }) / 1000
            
    # subset based on distance
            n_dat2 <- plyr::join(n_dat2, n_dat3[, c("cluster", "min_dist_km")], "cluster", "left", "first")
            n_dat2 <- n_dat2 %>% filter(min_dist_km > 15)
            
            n_dat3 <- n_dat3 %>% filter(min_dist_km > 15)
          
    # make the population summary into a spatial object
          sf_n_dat3 <- st_as_sf(n_dat3, coords = c("longitude", "latitude"), crs = 4326)
          
    # convert to Albers equal area projection
          sf_n_dat3 <- st_transform(sf_n_dat3, crs = 5070)
          
# Combined Nestwatch & university dataset ----
  # first make a combined data frame that has all locations in both contributed and nestwatch data
      # with one row per population
        all_pop_locs <- data.frame(
          source = c(rep("contribution", nrow(dn_sum)), rep("nestwatch", nrow(n_dat3))),
          pop_id = c(paste(dn_sum$pop_id, dn_sum$sub_pop, sep = "_"),
                     paste("nwp", n_dat3$cluster, sep = "_")),
          n_years = c(dn_sum$n_years, n_dat3$n_years),
          n_nests = c(dn_sum$tot_nests, n_dat3$n_nests),
          first_year = c(dn_sum$first_year, n_dat3$first_year),
          last_year = c(dn_sum$last_year, n_dat3$last_year),
          latitude = c(dn_sum$avg_lat, n_dat3$latitude),
          longitude = c(dn_sum$avg_lon, n_dat3$longitude),
          mu_first_egg = c(dn_sum$mu_lay, n_dat3$mu_lay),
          md_first_egg = c(dn_sum$md_lay, n_dat3$md_lay),
          mu_fledged = c(dn_sum$mu_fledged, n_dat3$mu_fledged)
        )
        all_pop_locs <- all_pop_locs %>% filter(is.na(latitude) == FALSE)
        
  # make the combined data set with one row per nest (this drops some columns to make the two datasets match up)
        
        all_nests <- data.frame(
          pop_id = c(d_nest$pop_id, paste("nwp", n_dat2$cluster, sep = "_")),
          sub_pop = c(d_nest$sub_pop, rep(NA, nrow(n_dat2))),
          latitude = c(d_nest$latitude, n_dat2$latitude),
          longitude = c(d_nest$longitude, n_dat2$longitude),
          year = c(d_nest$year, n_dat2$year),
          lay_doy = c(d_nest$lay_doy, n_dat2$lay_doy),
          num_fledged = c(d_nest$number_fledged_num, n_dat2$young_fledged),
          source = c(rep("contribution", nrow(d_nest)), rep("nestwatch", nrow(n_dat2))),
          exclude_fitness = c(d_nest$exclude_fitness, rep("no", nrow(n_dat2))),
          clutch_size = c(d_nest$clutch_size, n_dat2$clutch_size)
        )
        all_nests <- all_nests %>% filter(year > 1970)
        
        # summarise this to one row per site-year
          any <- all_nests %>%
            group_by(pop_id, sub_pop, year) %>%
            summarise(n_nests_this_year = n(), 
                      mu_num_fled = mean(num_fledged),
                      first_egg_this_year = if(all(is.na(lay_doy))) NA_real_ else quantile(
                        lay_doy, probs = 0.1, na.rm = TRUE), .groups = "drop")
          
          
        all_nests <- plyr::join(all_nests, any, c("pop_id", "sub_pop", "year"), "left", "first")
        all_nests$day_from_first_egg <- all_nests$lay_doy - all_nests$first_egg_this_year
        
        # count years with 15+ nests within 42 days (7 weeks) of 10th percentile first egg
          any2 <- all_nests %>%
            filter(lay_doy - first_egg_this_year < 42) %>%
            group_by(pop_id, sub_pop, year) %>%
            summarise(n_new = n()) %>%
            filter(n_new > 14) %>%
            group_by(pop_id, sub_pop) %>%
            summarise(used_years = n())
          
        # get the overall grand mean and median lay date for each site pooling across years
          any3 <- all_nests %>%
            filter(lay_doy - first_egg_this_year < 43) %>%
            group_by(pop_id, sub_pop) %>%
            summarise(grand_mu_lay = mean(lay_doy), grand_md_lay = median(lay_doy))
          
        # add mean and median lay dates back into the other dataframes
          all_nests <- plyr::join(all_nests, any2, c("pop_id", "sub_pop"), "left", "first")
          all_nests <- plyr::join(all_nests, any3, c("pop_id", "sub_pop"), "left", "first")
        
  # nests summarised by year and site with lay date
        nests_by_year <- all_nests %>%
          filter(n_nests_this_year > 14, lay_doy - first_egg_this_year < 43) %>%
          group_by(pop_id, sub_pop, year, grand_mu_lay, grand_md_lay) %>%
          summarise(included_nests = n(), mu_lay = mean(lay_doy), md_lay = median(lay_doy), mu_fled = mean(num_fledged, na.rm = TRUE))

  # get elevation from elevatr package     
        all_pop_locs_sf <- st_as_sf(all_pop_locs, coords = c("longitude", "latitude"), crs = 4326)
        elev_data <- get_elev_point(locations = all_pop_locs_sf, units = "meters", src = "aws")
        all_pop_locs$elevation_m <- elev_data$elevation         
        
  # make a description histogram
        nest_hist <- ggplot(all_nests, mapping = aes(x = year, fill = source)) +
          geom_histogram(binwidth = 1, alpha = 0.6, color = "gray25") +
          facet_wrap(~ source, labeller = labeller(source = c("contribution" = "Contributed",
                                                             "nestwatch" = "Nestwatch"))) +
          theme_classic() +
          scale_fill_manual(values = palette.colors(3)[2:3],
                            labels = c("Contributed", "Nestwatch"),
                            name = "Source") +
          theme(panel.grid = element_blank(), axis.title = element_text(size = 14),
                axis.text = element_text(size = 12)) +
          labs(x = "Year", y = "Nests per year") +
          guides(fill = "none")
          
        ggsave(here::here("saved_figures", "nest_by_year_histogram.png"), nest_hist,
               device = "png", width = 6.9, height = 3.7, units = "in")
            
# NALCMS land cover data ----
    # downloaded from their project site on 4/21/2025 for 2010 land cover map
    # in order to run this section download the data and place it in the correct directory
            
    # load the land cover raster
        lc_raster <- raster(here::here("land_cover_dl", "land_cover_2010v3_30m_tif",
                                       "NA_NALCMS_landcover_2010v3_30m", "data", 
                                       "NA_NALCMS_landcover_2010v3_30m.tif"))
            
    # create 5km buffers around the points of each location
            # test <- data.frame(place = c("ohio", "nebraska", "minnesota"),
            #                   latitude = c(40.761286, 42.3721, 47.485), 
            #                   longitude = c(-82.512968, -92.7492, -93.602))
        # turn all pops locations into sf object in wgs projection
            all_pop_locs_sf <- st_as_sf(all_pop_locs, coords = c("longitude", "latitude"), crs = 4326) #wgs
        # transform points of locations into albers equal area
            points_proj <- st_transform(all_pop_locs_sf, crs = 5070) #albers equal area
        # create the buffers with 2.5 km radius buffer around each point
            buffers = st_buffer(points_proj, dist = 1000) # 2.5km buffer around points
        # transform projection to match the raster for clipping
            buffers2 <- st_transform(buffers, crs = crs(lc_raster)) # transform to match raster
        # extract from raster for each buffer
            zonal_stats <- exact_extract(lc_raster, buffers2, progress = TRUE)
        # for each extracted buffer, calculate zonal statistics of percent cover for each category   
            zonal_summary <- lapply(seq_along(zonal_stats), function(i) {
              stats_i <- zonal_stats[[i]]
              
              stats_i %>%
                group_by(value) %>%
                summarise(class_cov = sum(coverage_fraction, na.rm = TRUE)) %>%
                mutate(id = i,
                       prop = class_cov / sum(class_cov)) %>%
                dplyr::select(id, value, prop)
            })
        # Combine into a single data frame
            zonal_summary_df <- bind_rows(zonal_summary)
        # cast land cover into a wide format with a column for each land cover type    
            zonal_wide <- zonal_summary_df %>%
              mutate(value = paste0("lc_", value)) %>%  
              tidyr::pivot_wider(
                names_from = value,
                values_from = prop,
                values_fill = 0  # Fill missing land cover classes with 0%
              )
        # using a new version of all_pop_locs here and adding some columns to join
            a_pops <- all_pop_locs
            a_pops$id <- 1:nrow(a_pops)
            a_pops_lc <- left_join(a_pops, zonal_wide, by = "id")
            
        # converting land cover columns to percentages
            lc_cols <- grep("^lc_", names(a_pops_lc), value = TRUE)
            a_pops_lc[lc_cols] <- round(a_pops_lc[lc_cols] * 100, 2)
           
        # adding in the descriptions and plotting color for land cover categories 
              lc_cats <- data.frame(lc_code = paste0("lc_", seq(1, 19, 1)),
                                    lc_num = seq(1, 19, 1),
                                    lc_cat = c("forest_1", "forest_2", "forest_3", "forest_4", "forest_5", "forest_6",
                                               "shrub_1", "shrub_2", "grass_1", "grass_2",
                                               "shrub_3", "grass_3", "barren_2", "wetland",
                                               "agriculture", "barren_1", "urban", "open_water", "snow_ice"))
        # renaming to match the land cover names I want
              name_map <- setNames(lc_cats$lc_cat, lc_cats$lc_code)
              existing_lc_cols <- intersect(names(name_map), names(a_pops_lc))
              
              a_pops_lc2 <- a_pops_lc %>%
                rename_with(.fn = ~ name_map[.x], .cols = all_of(existing_lc_cols))
              
        # calculating some grouped percentages for similar land categories
              a_pops_lc2$woody <- a_pops_lc2$forest_1 + a_pops_lc2$forest_5 + a_pops_lc2$forest_6 
              a_pops_lc2$grass_shrub <- a_pops_lc2$shrub_2 + a_pops_lc2$shrub_1 + a_pops_lc2$grass_1 + a_pops_lc2$grass_2
              a_pops_lc2$wet <- a_pops_lc2$wetland + a_pops_lc2$open_water
              a_pops_lc2$urb_barren <- a_pops_lc2$urban + a_pops_lc2$barren_1
          
        # pivot longer for plotting    
              a_pops_lc2_long <- a_pops_lc2 %>%
                pivot_longer(cols = c(woody, grass_shrub, agriculture, urb_barren, wet),
                             names_to = "land_cover_type",
                             values_to = "value")
        # make a land cover plot histogram for all sites and save it to file   
              land_cover_summary <- ggplot(a_pops_lc2_long, aes(x = value, fill = land_cover_type)) +
                geom_histogram(binwidth = 5, alpha = 0.8, color = "gray50") +
                facet_wrap(~ land_cover_type, labeller = labeller(land_cover_type = c(agriculture = "Agricultural",
                                                                                    grass_shrub = "Grass or Shrubland",
                                                                                    urb_barren = "Urban or Barren",
                                                                                    wet = "Wetland or Water",
                                                                                    woody = "Forrested"))) +
                theme_classic() +
                theme(panel.grid = element_blank(), axis.text = element_text(size = 12),
                      axis.title = element_text(size = 14), legend.position = c(0.82, 0.23)) +
                labs(x = "Percent of land area", y = "Number of sites") +
                scale_fill_manual(values = c("goldenrod", "lightgreen", "coral3", "lightblue", "darkgreen")) +
                guides(fill = "none")
              
              ggsave(here::here("saved_figures", "land_cover_summary.png"), land_cover_summary,
                     device = "png", units = "in", dpi = 300, width = 7, height = 4.8)
          
          
# Land cover PCA ----
    # make a reduced dataframe for pca
          lc_pca <- a_pops_lc2[, c("pop_id", "woody", "grass_shrub", "agriculture", "wetland", "open_water", "urban", "barren_1")]  
              
      # because these are proportional data, I'm using centered log ratio transformation
          # Extract only the land cover proportions (not pop_id)
            land_cover_matrix <- lc_pca %>%
              dplyr::select(-pop_id) %>%
              as.matrix()
            
          # Replace 0s with small positive number to avoid log(0)
            land_cover_matrix[land_cover_matrix == 0] <- 0.0001
            
          # Apply centered log-ratio (CLR) transformation
            clr_matrix <- compositions::clr(land_cover_matrix)
            
          # Perform PCA on CLR-transformed data
            pca_result <- prcomp(clr_matrix, center = TRUE, scale. = FALSE)
            
          # Add PC scores to your original data
            lc_pca_with_scores <- bind_cols(
              lc_pca %>% dplyr::select(pop_id),
              as.data.frame(pca_result$x[, 1:3])  # First 3 PCs
            )
            
      # make pca plot
            
            # Variance explained for custom facet titles
                var_explained <- summary(pca_result)$importance["Proportion of Variance", 1:3] * 100
                pc_labels <- paste0("PC", 1:3, " (", round(var_explained, 1), "%)")
                names(pc_labels) <- paste0("PC", 1:3)
            
            # Cleaned-up land cover labels (renaming woody to Forested)
                nice_labels <- c(
                  woody = "Forested     ",
                  grass_shrub = "Grass/Shrub     ",
                  agriculture = "Agriculture     ",
                  wetland = "Wetland     ",
                  open_water = "Open Water     ",
                  urban = "Urban     ",
                  barren_1 = "Barren     "
                )
            
            # Prepare data
                loadings_df <- as.data.frame(pca_result$rotation[, 1:2]) %>%
                  rownames_to_column("land_cover") %>%
                  pivot_longer(cols = starts_with("PC"),
                               names_to = "PC",
                               values_to = "loading") %>%
                  mutate(
                    land_cover = factor(land_cover, levels = names(nice_labels), labels = nice_labels),
                    PC = factor(PC, levels = names(pc_labels), labels = pc_labels),
                    pos = loading > 0
                  )
            
            # Define y-position for label placement
                label_y_offset <- min(loadings_df$loading) - .08 * diff(range(loadings_df$loading))
            
            # Plot
                pc_plot <- ggplot(loadings_df, aes(x = land_cover, y = loading, fill = pos)) +
                  geom_col(show.legend = FALSE) +
                  facet_wrap(~ PC, scales = "free_y", labeller = labeller(PC = pc_labels)) +
                  coord_flip(clip = "off") +  # allow labels to overflow into margin
                  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
                  labs(x = NULL, y = "Loading") +
                  theme_minimal(base_size = 13) +
                  theme(
                    strip.text = element_text(face = "bold"),
                    panel.spacing.x = unit(2, "lines"),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    plot.margin = unit(c(1, 1, 1, 7), "lines")  # more space on the left
                  ) +
                  # Manual land cover labels for PC1 only, nudged slightly left
                      geom_text(
                        data = filter(loadings_df, PC == pc_labels[1]),
                        aes(x = land_cover, y = label_y_offset - 0.05, label = land_cover),
                        inherit.aes = FALSE,
                        hjust = 1,
                        size = 3.5
                      ) +
                  geom_hline(yintercept = 0)
                
                ggsave(here::here("saved_figures", "land_cover_pca_loadings.png"), pc_plot,
                       device = "png", width = 7.1, height = 3.5, units = "in", dpi = 300)
            
              
# Make a map of populations ----
  # make the population summary into a spatial object
        sf_dn_sum <- st_as_sf(dn_sum[is.na(dn_sum$avg_lat) == FALSE, ], coords = c("avg_lon", "avg_lat"), crs = 4326)
        sf_all_pop <- st_as_sf(all_pop_locs[is.na(all_pop_locs$latitude) == FALSE &
                                              all_pop_locs$n_nests > 120, ], coords = c("longitude", "latitude"), crs = 4326)
        
  # convert to Albers equal area projection
        sf_dn_sum <- st_transform(sf_dn_sum, crs = 5070)
        sf_all_pop <- st_transform(sf_all_pop, crs = 5070)
        
  # get tree swallow range map from ebird
        #set_ebirdst_access_key(ebird_key)
        range_data <- ebirdst_download_status(species = "treswa")
        range_raster <- load_raster("treswa", period = "seasonal", resolution = "9km")
        breeding_raster <- terra::project(range_raster[[1]], "EPSG:5070")
        
        raster_df <- as.data.frame(breeding_raster, xy = TRUE, na.rm = TRUE)
        colnames(raster_df) <- c("x", "y", "value")
        
        # make bins
          raster_df$bin <- NA
          raster_df$bin[raster_df$value == 0] <- "0"
          raster_df$bin[raster_df$value > 0] <- cut_number(raster_df$value[raster_df$value > 0], n = 8, na.rm = TRUE)
          
        # Generate color palette with transparent color for the zero bin
          n_bins <- length(unique(raster_df$bin[!is.na(raster_df$bin)]))
          blue_palette <- colorRampPalette(c("#dbe9f6", "#878FAB"))(n_bins - 1)  #used to be #000033
          color_values <- c("#00000000", blue_palette)  # Add transparent color for bin "0"  
        
        # from ebirdst vignette: https://ebird.github.io/ebirdst/articles/applications.html
          v <- values(breeding_raster, na.rm = TRUE, mat = FALSE)
          v <- v[v > 0] 
          breaks <- quantile(v, seq(0, 1, by = 0.1))
          #add a bin for 0
            breaks <- c(0, breaks)
            
          pal <- ebirdst_palettes(length(breaks) -2 )
          
          pal <- c("#e6e6e600", pal)
        
  # get north america map
        n_amer <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
        n_amer <- st_transform(n_amer, crs = 5070)
        
  # bounding box
        points_bbox <- st_bbox(sf_dn_sum)

  # create the plot
        # ggplot() +
        #   geom_sf(data = n_amer, fill = "lightgray", color = "white") +
        #   geom_sf(data = sf_dn_sum, aes(size = tot_nests, color = n_years), alpha = 0.7) +
        #   scale_color_viridis_c(option = "plasma") +
        #   geom_raster(raster_df, aes(x = x, y = y, fill = value)) +
        #   theme_bw() +
        #   coord_sf(
        #     xlim = c(points_bbox["xmin"] - 8e5, points_bbox["xmax"] + 5e5),
        #     ylim = c(points_bbox["ymin"] - 1e6, points_bbox["ymax"] + 5e5),
        #     expand = FALSE  # Prevents ggplot from adding extra margins
        #   ) +
        #   labs(
        #     color = "years",
        #     size = "nests"
        #   ) +
        #   xlab("Longitude") +
        #   ylab("Latitude") +
        #   theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
        
        # for contributed data from populations
            # ggplot() +
            #        # Base map
            #        geom_sf(data = n_amer, fill = "lightgray", color = "white") +
            #        
            #        # Add raster layer with equal area bins
            #        geom_tile(data = raster_df, aes(x = x, y = y, fill = bin)) +
            #        scale_fill_manual(
            #              values = color_values,
            #              breaks = levels(raster_df$bin),
            #              labels = function(x) ifelse(x == "0", "None", gsub("[\\[\\]()]", "", x))  # Clean labels for legend
            #          ) +
            #        
            #        # Add points
            #        geom_sf(data = sf_dn_sum, aes(size = tot_nests, color = n_years), alpha = 0.85) +
            #        scale_color_viridis_c(option = "plasma", begin = 0.4) +
            #        
            #        # Map adjustments
            #        theme_bw() +
            #        coord_sf(
            #              xlim = c(points_bbox["xmin"] - 8e5, points_bbox["xmax"] + 5e5),
            #              ylim = c(points_bbox["ymin"] - 1e6, points_bbox["ymax"] + 5e5),
            #              expand = FALSE  # Prevents ggplot from adding extra margins
            #          ) +
            #        
            #        # Labels and styling
            #        labs(
            #              color = "Years",
            #              size = "Nests",
            #              fill = "Abundance"
            #         ) +
            #        xlab("Longitude") +
            #        ylab("Latitude") +
            #        theme(
            #              axis.text = element_text(size = 12),
            #              axis.title = element_text(size = 14),
            #              legend.position = "right"
            #          )
            
          # for pseudo populations from nestwatch
            # ggplot() +
            #        # Base map
            #        geom_sf(data = n_amer, fill = "lightgray", color = "white") +
            #        
            #        # Add raster layer with equal area bins
            #        geom_tile(data = raster_df, aes(x = x, y = y, fill = bin)) +
            #        scale_fill_manual(
            #              values = color_values,
            #              breaks = levels(raster_df$bin),
            #              labels = function(x) ifelse(x == "0", "None", gsub("[\\[\\]()]", "", x))  # Clean labels for legend
            #          ) +
            #        
            #        # Add points
            #        geom_sf(data = sf_n_dat3, aes(size = n_nests, color = n_years), alpha = 0.85) +
            #        scale_color_viridis_c(option = "plasma", begin = 0.4) +
            #        
            #        # Map adjustments
            #        theme_bw() +
            #        coord_sf(
            #              xlim = c(points_bbox["xmin"] - 8e5, points_bbox["xmax"] + 5e5),
            #              ylim = c(points_bbox["ymin"] - 1e6, points_bbox["ymax"] + 5e5),
            #              expand = FALSE  # Prevents ggplot from adding extra margins
            #          ) +
            #        
            #        # Labels and styling
            #        labs(
            #              color = "Years",
            #              size = "Nests",
            #              fill = "Abundance"
            #         ) +
            #        xlab("Longitude") +
            #        ylab("Latitude") +
            #        theme(
            #              axis.text = element_text(size = 12),
            #              axis.title = element_text(size = 14),
            #              legend.position = "right"
            #          )
            
          # for both contributions and nestwatch in one plot
            pop_map <- ggplot() +
                   # Base map
                   geom_sf(data = n_amer, fill = "gray90", color = "white") +
                   
                   # Add raster layer with equal area bins
                   geom_tile(data = raster_df, aes(x = x, y = y, fill = bin)) +
                   scale_fill_manual(
                         values = color_values,
                         breaks = levels(raster_df$bin),
                         labels = function(x) ifelse(x == "0", "None", gsub("[\\[\\]()]", "", x))  # Clean labels for legend
                     ) +
                   
                   # Add points
                   geom_sf(data = sf_all_pop, aes(size = n_nests, color = n_years, shape = source), alpha = 0.9) +
                   scale_color_viridis_c(option = "plasma", begin = 0.4) +
                   scale_size_continuous(breaks = c(500, 1000, 2000, 4000), limits = c(0, 8000)) +
                   
                   # Map adjustments
                   theme_bw() +
                   coord_sf(
                         xlim = c(points_bbox["xmin"] - 8e5, points_bbox["xmax"] + 5e5),
                         ylim = c(points_bbox["ymin"] - 1e6, points_bbox["ymax"] + 5e5),
                         expand = FALSE  # Prevents ggplot from adding extra margins
                     ) +
                   
                   # Labels and styling
                   labs(
                         color = "Years",
                         size = "Nests",
                         fill = "Abundance"
                    ) +
                   xlab("Longitude") +
                   ylab("Latitude") +
                   theme(
                         axis.text = element_text(size = 12),
                         axis.title = element_text(size = 14),
                         legend.position = "right"
                     ) +
              guides(shape = "none")
            
            ggsave(here::here("saved_figures", "tres_pop_plot.png"), device = "png", 
                   units = "in", width = 6.9, height = 4.8, dpi = 300)
            
# Extract relative abundance from eBird s&t raster ----
  # make spatial object of all population locations
      ap84 <- st_as_sf(all_pop_locs, coords = c("longitude", "latitude"), crs = 4326)
            
  # confirm it matches raster projection
      ap84 <- st_transform(ap84, crs = st_crs(breeding_raster))
      
  # extract point values at centroid of each population
      point_values <- extract(breeding_raster, ap84)
      
  # add points to the plot
      all_pop_locs$rel_abundance <- point_values$breeding
            
# Download temperature from Daymet ----
            
            #test_data <- all_pop_locs[1:5, ]
       
    # This is the actual daymet download. It takes a while and makes a bit file, so only run it if needed
      # when populations or dates are updated. Otherwise just load saved data.
        # # function to download single year/location temperature from daymet
              # download_daymet_year <- function(lat, lon, year, pop_id){
              #   tryCatch({
              #     dat <- download_daymet(
              #       lat = lat,
              #       lon = lon,
              #       start = year,
              #       end = year,
              #       internal = TRUE
              #     )
              #     # add ID and year
              #     df <- dat$data %>%
              #       mutate(pop_id = pop_id, year = year)
              #     return(df)
              #   }, error = function(e){
              #     message(glue::glue("Error: {e$message} for {pop_id} in {year}"))
              #     return(NULL)
              #   })
              # }
        # 
        #   # loop through year location combos and download
              # all_temp_data <- purrr::pmap_dfr(all_pop_locs, function(source, pop_id, n_years, n_nests,
              #                                               first_year, last_year, latitude, longitude,
              #                                               mu_first_egg, md_first_egg, mu_fledged,
              #                                               elevation_m, rel_abundance) {
              #   purrr::map_dfr(first_year:last_year, function(y){
              #     download_daymet_year(lat = latitude, lon = longitude, year = y, pop_id = pop_id)
              #   })
              # })

              # colnames(all_temp_data) <- c("year", "yday", "daylength_s", "precip_mm", "solar_rad_watts_m",
              #                              "snow_water_equiv", "tmax", "tmin", "vapor_pressure", "pop_id")

      # # the date format that climwin wants
          # all_temp_data$date <- mapply(function(year, yday){
          #   format(as.Date(yday, origin = paste0(year, "-01-01")) - 1, "%d/%m/%Y")
          # }, all_temp_data$year, all_temp_data$yday)

      # # add an average temperature column
            # all_temp_data$tavg <- (all_temp_data$tmax + all_temp_data$tmin) / 2

      # Load saved temperature data
        #saveRDS(all_temp_data, here::here("saved_data_downloads", "all_temperature.rds"))
        #all_temp_data <- readRDS(here::here("saved_data_downloads", "all_temperature.rds"))
        
# Download temperature from MERRA-2 ----
  # Daymet is limited to only 1980-2023 at present. Only three pops have data before 1980, but a lot
        # have data from 2024 so this shortens most time series by a full year
        # Use nasapower
 
          # for(i in 1:nrow(all_pop_locs)){
          # 
          #   if(all_pop_locs$first_year[i] < 1981){first_yr <- "1981-01-01"}
          #   if(all_pop_locs$first_year[i] > 1980){first_yr <- paste0(all_pop_locs$first_year[i], "-01-01")}
          # 
          #   np <- get_power(
          #     community = "AG",
          #     lonlat = c(all_pop_locs$longitude[i], all_pop_locs$latitude[i]),
          #     pars = c("T2M", "T2M_MAX", "T2M_MIN", "RH2M", "PS", "WS2M", "PRECTOTCORR"),
          #     dates = c(first_yr, paste0(all_pop_locs$last_year[i], "-12-31")),
          #     temporal_api = "DAILY"
          #   )
          # 
          #   np$merra_elevation <- as.numeric(strsplit(attr(np,"POWER.Elevation"),split=" ")[[1]][13])
          #   np$pop_id <- all_pop_locs$pop_id[i]
          #   np$source <- all_pop_locs$source[i]
          #   np$site_elev <- all_pop_locs$elevation_m[i]
          # 
          #   np <- as.data.frame(np)
          # 
          #   if(i == 1){
          #     all_merra_temperature <- np
          #   }
          #   if(i > 1){
          #     all_merra_temperature <- rbind(all_merra_temperature, np)
          #   }
          # 
          #   print(paste("Completed", i, "of", nrow(all_pop_locs), sep = " "))
          # }
      # 
      # 
      # # use lapse rate to correct for elevation differences
        # lapse_rate <- 0.0065
        # all_merra_temperature$T2M_c <- all_merra_temperature$T2M +
        #   lapse_rate * (all_merra_temperature$site_elev - all_merra_temperature$merra_elevation)
        # all_merra_temperature$T2M_MAX_c <- all_merra_temperature$T2M_MAX +
        #   lapse_rate * (all_merra_temperature$site_elev - all_merra_temperature$merra_elevation)
        # all_merra_temperature$T2M_MIN_c <- all_merra_temperature$T2M_MIN +
        #   lapse_rate * (all_merra_temperature$site_elev - all_merra_temperature$merra_elevation)
      # 
      # # change the column names to better ones
        # colnames(all_merra_temperature) <- c("longitude", "latitude", "year", "month", "day", "yday",
        #                                      "yyymmdd", "tavg", "tmax", "tmin", "rh2m", "ps", "ws2m",
        #                                      "precip", "merra_elevation", "pop_id", "source", "site_elev",
        #                                      "tavg_c", "tmax_c", "tmin_c")
      # 
      # # the date format that climwin wants
        # all_merra_temperature$date <- mapply(function(year, yday){
        #   format(as.Date(yday, origin = paste0(year, "-01-01")) - 1, "%d/%m/%Y")
        # }, all_merra_temperature$year, all_merra_temperature$yday)

       #Load saved temperature data
       #saveRDS(all_merra_temperature, here::here("saved_data_downloads", "all_merra_temperature.rds"))
        all_merra_temperature <- readRDS(here::here("saved_data_downloads", "all_merra_temperature.rds"))

# Make maps of land cover for specific example populations ----
      # this is purely for making example circular plots for specific sites
      # I chose three sites with very different land cover
      # the code here produces the parts of the plots, but I put the three
      # together and set up axes in illustrator so reproducing will require
      # running this code a few times with different locations saved and combined
      
      # Note that the figure included in the paper includes additional annotation done manually
      
      # # working with raster from previous section loaded
      #     # Test location data
      #         test <- data.frame(
      #           place = all_pop_locs$pop_id,
      #           latitude = all_pop_locs$latitude, 
      #           longitude = all_pop_locs$longitude
      #         )
      # 
      #         test <- test %>% filter(place == "tsp_001_E" | place == "tsp_012_A" | place == "tsp_018_D")
      #     
      #     # Convert to sf and reproject
      #         test_sf <- st_as_sf(test, coords = c("longitude", "latitude"), crs = 4326)
      #         test_proj <- st_transform(test_sf, crs = crs(lc_raster))  # match raster CRS
      #         
      #     # Pick one location
      #         location <- test_proj %>% filter(place == "tsp_018_D")
      #     
      #     # Create 2.5km circular buffer
      #         buffer <- st_buffer(location, dist = 2500)
      #     
      #     # Convert sf buffer to Spatial for use with raster::crop/mask
      #        buffer_sp <- as(buffer, "Spatial")
      #     
      #     # Crop and mask the raster
      #         lc_crop <- crop(lc_raster, extent(buffer_sp))
      #         lc_masked <- mask(lc_crop, buffer_sp)
      #         #lc_masked_utm <- projectRaster(lc_masked, crs = CRS("+proj=utm +zone=18 +datum=WGS84"))  # montreal
      #         #lc_masked_utm <- projectRaster(lc_masked, crs = CRS("+proj=utm +zone=10 +datum=WGS84"))  # davis
      #         lc_masked_utm <- projectRaster(lc_masked, crs = CRS("+proj=utm +zone=17 +datum=WGS84"))  # davidson
      #         
      #       # Convert the raster to a data frame for plotting
      #         lc_df <- as.data.frame(lc_masked, xy = TRUE)
      #         head(lc_df)  # Check the first few rows of the data
      #         
      #         lc_df2 <- as.data.frame(lc_masked_utm, xy = TRUE)
      #       
      #   # make a plot with the right coordinates to be used for axes  
      #         lc_axes <- ggplot(lc_df2[is.na(lc_df2$Class_EN) == FALSE, ]) +
      #           geom_raster(aes(x = x, y = y, fill = Class_EN)) +
      #           theme_bw() +
      #           theme(panel.grid = element_blank()) +
      #           coord_fixed(ratio = 1) +
      #           guides(fill = "none")
      #         
      #         ggsave("davidson_axes.svg", lc_axes, device = "svg", units = "in", width = 4, height = 4, dpi = 300)
      #         
      #         
      #    # make a plot of the site     
      #         lc_plot <- ggplot(lc_df[is.na(lc_df$Class_EN_Class_EN) == FALSE, ]) +
      #           geom_raster(aes(x = x, y = y, fill = Class_EN_Class_EN)) +
      #           theme_bw() +
      #           theme(panel.grid = element_blank()) +
      #           scale_fill_manual(values = c("Urban and Built-up" = "red",
      #                                        "Temperate or Subpolar Broadleaf Deciduous Forest" = "forestgreen",
      #                                        "Mixed Forest" = "forestgreen",
      #                                        "Cropland" = "darkorange",
      #                                        "Baren Land" = "gray70",
      #                                        "Temperate or Subpolar Grassland" = "gold",
      #                                        "Temperate or Subpolar Shrubland" = "gold",
      #                                        "Wetland" = "slateblue",
      #                                        "Temperate or Subpolar Needleaf Forest" = "forestgreen",
      #                                        "Water" = "lightblue")) +
      #           guides(fill = "none") + coord_fixed(ratio = 1) +
      #           labs(x = "Easting (m)", y = "Northing (m)")
      #         
      #         ggsave("davidson.svg", lc_plot, device = "svg", units = "in", width = 4, height = 4, dpi = 300)
            
# Run climwin analysis to identify windows ----
        nests_by_year$pop_full <- paste(nests_by_year$pop_id, nests_by_year$sub_pop, sep = "_")
        nests_by_year$pop_full <- gsub("_NA", "", nests_by_year$pop_full)
        length_pops <- length(unique(nests_by_year$pop_full))
        
  # add reference day to each pop 14 days after mean egg laying date across years
        rd <- nests_by_year %>%
          dplyr::group_by(pop_full) %>%
          summarise(ref_doy = round(mean(mu_lay), 0) + 14)
        nests_by_year <- plyr::join(nests_by_year, rd, "pop_full", "left", "first")
        
  # make data frame to store climwin results
        climwin_results <- data.frame(
          pop = unique(nests_by_year$pop_full),
          open_w = NA,
          close_w = NA,
          open_b = NA,
          close_b = NA,
          intercept = NA,
          slope = NA,
          r2 = NA,
          mod_p = NA,
          n_years = NA,
          ref_doy = NA,
          clim_low = NA,
          clim_hi = NA,
          open_w_doy = NA,
          open_w_low = NA,
          open_w_hi = NA,
          close_w_doy = NA,
          close_w_low = NA,
          close_w_hi = NA,
          open_b_doy = NA,
          close_b_doy = NA,
          mw_intercept = NA,
          mw_slope = NA
        )
        
    # function to calculate weighted confidence intervals for window open and close based on models
        weighted_quantile <- function(x, w, probs = c(0.05, 0.95)) {
          # Remove missing
          ok <- !(is.na(x) | is.na(w))
          x <- x[ok]
          w <- w[ok]
          
          # Sort x and weights by x
          ord <- order(x)
          x_sorted <- x[ord]
          w_sorted <- w[ord]
          
          # Cumulative weight
          cum_w <- cumsum(w_sorted) / sum(w_sorted)
          
          # Interpolate quantiles
          approx(cum_w, x_sorted, xout = probs, ties = "ordered")$y
        }
        
      
      # add in the number of years with 15 nests to results dataframe
        for(i in 1:nrow(climwin_results)){
          climwin_results$years_with_15_nests[i] <-
            nrow(subset(nests_by_year, nests_by_year$pop_full == climwin_results$pop[i] &
                          nests_by_year$included_nests > 14))
        }
        
  # little function to convert to month and day in the format climwin wants
        doy_to_month_day <- function(doy){
          ref_date <- as.Date("2001-01-01")
          date <- ref_date + (doy - 1)
          data.frame(doy = doy, month = format(date, "%m"), day = format(date, "%d"))
        }
        
  # model weight function to recalculate after subsetting clmwin
        akaike_weights <- function(delta) {
          rel_lik <- exp(-0.5 * delta)
          rel_lik / sum(rel_lik)
        }
        
  # run climwin in a loop through each population and save output
        for(i in 1:nrow(climwin_results)){
          # subset to the right population
              temp_pop <- nests_by_year %>%
                filter(pop_full == climwin_results$pop[i], year > 1980) 
                    # daymet only 1980-2023 right now
                    # merra only 1981+
              temp_pop$date <- mapply(function(year, mu_lay){
                format(as.Date(mu_lay, origin = paste0(year, "-01-01")) - 1, "%d/%m/%Y")
              }, temp_pop$year, temp_pop$mu_lay)
              temp_pop <- temp_pop %>% filter(included_nests > 14)
              
          if(nrow(temp_pop) > 0){
            
            # subset to the temperature data for the matching population and remove duplicates
            # use daymet temperature
            # temp_temp <- all_temp_data %>%
            #   filter(pop_id == climwin_results$pop[i])
            # use merra temperature
                temp_temp <- all_merra_temperature %>%
                  filter(pop_id == climwin_results$pop[i])
            # check for any duplicate downloads and remove  
                temp_temp <- temp_temp[!duplicated(temp_temp$date), ]
            
            # determine the reference date for this population
            # this is 14 days after the overall average mean egg date for the population
            # but note that some nests are still being laid after this point so it
            # makes sense to be later
                ref_doy <- round(mean(temp_pop$mu_lay), 0) + 14
                ref <- doy_to_month_day(ref_doy)
            
            if(nrow(temp_pop) > 1){
              # run the sliding window analysis
                  cwin <- slidingwin(xvar = list(Temp = temp_temp$tavg_c),
                                     cdate = temp_temp$date,
                                     bdate = temp_pop$date,
                                     baseline = lm(mu_lay ~ 1, data = temp_pop),
                                     cinterval = "day",
                                     range = c(56, 0),
                                     type = "absolute", 
                                     refday = c(as.numeric(ref[1, 3]), as.numeric(ref[1, 2])),
                                     stat = "mean",
                                     func = "lin",
                                     cmissing = "method1")
                  cwin_out <- cwin[[1]]$Dataset
                  
                  cwin_out_r <- cwin_out %>% filter(WindowOpen - WindowClose > 13, WindowOpen - WindowClose < 29,
                                                    WindowClose < 28)
                  cwin_out_r$deltaAICc <- cwin_out_r$deltaAICc - min(cwin_out_r$deltaAICc)
                  cwin_out_r$ModWeight <- akaike_weights((cwin_out_r$deltaAICc))
                  
              # get p value for comparison
                  sample.size = nrow(cwin[[1]]$BestModelData)
              
              # collect all the info from this run into the dataframe
                  if(is.na(cwin_out$ModWeight[1]) == FALSE){
                    if(cwin_out$ModWeight[1] < 0.95){
                      mwin <- medwin(cwin_out_r)
                      climwin_results$open_w[i] <- mwin$`Median Window Open`
                      climwin_results$close_w[i] <- mwin$`Median Window Close`
                      open_ci <- weighted_quantile(cwin_out_r$WindowOpen, cwin_out_r$ModWeight)
                      close_ci <- weighted_quantile(cwin_out_r$WindowClose, cwin_out_r$ModWeight)
                      climwin_results$open_w_low[i] <- open_ci[1]
                      climwin_results$open_w_hi[i] <- open_ci[2]
                      climwin_results$close_w_low[i] <- close_ci[1]
                      climwin_results$close_w_hi[i] <- close_ci[2]
                    }
                  }
              
                  if(is.na(cwin_out$ModWeight[1]) == FALSE){
                    if(cwin_out_r$ModWeight[1] > 0.95){
                      climwin_results$open_w[i] <- cwin_out_r$WindowOpen[1]
                      climwin_results$close_w[i] <- cwin_out_r$WindowClose[1]
                    }
                  }
              
                # collect all the output into a dataframe
                  climwin_results$open_b[i] <- cwin_out_r$WindowOpen[1]
                  climwin_results$close_b[i] <- cwin_out_r$WindowClose[1]
                  
                  climwin_results$intercept[i] <- coef(cwin[[1]]$BestModel)[1]
                  climwin_results$slope[i] <- coef(cwin[[1]]$BestModel)[2]
                  climwin_results$r2[i] <- summary(cwin[[1]]$BestModel)$r.squared
                  climwin_results$mod_p[i] <- summary(cwin[[1]]$BestModel)$coefficients[2, 4]
                  
                  climwin_results$n_years[i] <- sample.size
                  
                  climwin_results$ref_doy[i] <- ref_doy
                  
                  climwin_results$open_w_doy[i] <- ref_doy - climwin_results$open_w[i]
                  climwin_results$close_w_doy[i] <- ref_doy - climwin_results$close_w[i]
                  climwin_results$open_b_doy[i] <- ref_doy - climwin_results$open_b[i]
                  climwin_results$close_b_doy[i] <- ref_doy - climwin_results$close_b[i]
                  
                  climwin_results$mw_intercept[i] <- weighted.mean(cwin_out_r$ModelInt, cwin_out_r$ModWeight, na.rm = TRUE)
                  climwin_results$mw_slope[i] <- weighted.mean(cwin_out_r$ModelBeta, cwin_out_r$ModelBeta, na.rm = TRUE)
              
              
              # consensus window across all populations
              # the overall consensus median window is 35 days before to 14 days before the ref day
              # ref day is set at 14 days after the overall mean lay date at this site across all years
              # So this window goes from 21 days before to the day of the overall mean laying date for the pop
              
                # get average temperature for the consensus window
                    c_temp <- temp_temp %>%
                      filter(yday > ref_doy - 36, yday < ref_doy - 13) %>%
                      dplyr::group_by(year) %>%
                      dplyr::summarise(mu_tavg = mean(tavg_c),
                                       mu_tmax = mean(tmax_c),
                                       mu_tmin = mean(tmin_c))
                    temp_pop <- plyr::join(temp_pop, c_temp, "year", "left", "first")
                    
                    temp_m <- lm(mu_lay ~ mu_tavg, data = temp_pop)
                    temp_ms <- summary(temp_m)
                    
                # add consensus model stats to results
                    climwin_results$c_intercept[i] <- coefficients(temp_m)[1]
                    climwin_results$c_slope[i] <- coefficients(temp_m)[2]
                    climwin_results$c_pvalue[i] <- coefficients(temp_ms)[2, 4]
                    climwin_results$c_r2[i] <- temp_ms$r.squared 
                    climwin_results$clim_low[i] <- min(na.omit(c_temp$mu_tavg))
                    climwin_results$clim_hi[i] <- max(na.omit(c_temp$mu_tavg))
                    
                # get temperature from identified best window for this individual population
                    c_temp2 <- temp_temp %>%
                      filter(yday > climwin_results$open_b_doy[i] - 1, yday < climwin_results$close_b_doy[i] + 1) %>%
                      dplyr::group_by(year) %>%
                      dplyr::summarise(mu_tavg_i = mean(tavg_c),
                                       mu_tmax_i = mean(tmax_c),
                                       mu_tmin_i = mean(tmin_c))
                    temp_pop <- plyr::join(temp_pop, c_temp2, "year", "left", "first")
              
                # save output with temperature by year-site
                    if(i == 1){pop_with_climate <- temp_pop}
                    if(i > 1){pop_with_climate <- rbind(pop_with_climate, temp_pop)}
              
            }
            
            if(nrow(temp_pop) > 1){
                # save data frame of deviations for plotting
                    ct <- c_temp %>% filter(year %in% temp_pop$year)
                    mu_climate <- mean(ct$mu_tavg)
                    dev_df <- data.frame(year = temp_pop$year,
                                         pop_id = temp_pop$pop_full,
                                         mu_lay = temp_pop$mu_lay,
                                         climate = ct$mu_tavg,
                                         lay_diff = temp_pop$mu_lay - temp_pop$grand_mu_lay,
                                         clim_diff = ct$mu_tavg - mu_climate)
    
                    if(i == 1){dev_df_out <- dev_df}
                    if(i > 1){dev_df_out <- rbind(dev_df_out, dev_df)}
            }
            
            # progress tracker
            print(paste("Completed", i, "of", nrow(climwin_results), sep = " "))
            
          }
          }
           
  # make some plots from climwin output these can be filtered to what is useful   
        
        # first plot the  actual line centered for each population
            climwin_plotdata <- climwin_results %>%
              mutate(
                x_start = clim_low,
                x_end = clim_hi,
                y_start = mw_intercept + mw_slope * x_start,
                y_end = mw_intercept + mw_slope * x_end
              )
            
            for(i in 1:nrow(climwin_plotdata)){
              climwin_plotdata$c_x_start[i] <- climwin_plotdata$x_start[i] - mean(c(climwin_plotdata$x_start[i], climwin_plotdata$x_end[i]))
              climwin_plotdata$c_x_end[i] <- climwin_plotdata$x_end[i] - mean(c(climwin_plotdata$x_start[i], climwin_plotdata$x_end[i]))
              climwin_plotdata$c_y_start[i] <- climwin_plotdata$y_start[i] - mean(c(climwin_plotdata$y_start[i], climwin_plotdata$y_end[i]))
              climwin_plotdata$c_y_end[i] <- climwin_plotdata$y_end[i] - mean(c(climwin_plotdata$y_start[i], climwin_plotdata$y_end[i]))
            }
            
          # can adjsut filtering to only include years with lon term data
            climwin_plotdata %>%
              filter(n_years > 12) %>% 
            ggplot() +
              geom_segment(aes(x = c_x_start, xend = c_x_end, y = c_y_start, yend = c_y_end), alpha = 0.5) +
              theme_classic() +
              geom_hline(yintercept = 0, color = "coral3") +
              geom_vline(xintercept = 0, color = "coral3") +
              coord_cartesian(ylim = c(-10, 10), xlim = c(-4, 4))
        
      # plot overall opening and closing median windows estimated from each population with at least 8 years
            climwin_results %>%
              filter(years_with_15_nests > 7) %>%
              ggplot(aes(x = (open_w - 14) * -1)) +
              geom_histogram(alpha = 0.6, binwidth = 2) +
              geom_histogram(aes(x = (close_w - 14) * -1), binwidth = 2)
    
    
      # plot lay date anomaly vs. spring temperature anomaly in consensus window with each point as one site-year
          # and comparing contributed vs. nestwatch populations
            dev_df_out <- plyr::join(dev_df_out, all_pop_locs[, c("pop_id", "source")], "pop_id", "left", "first")
            
            p <- ggplot(dev_df_out, aes(x = clim_diff, y = lay_diff)) + #removed fill = source, color = source
              geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
              geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
              geom_point(alpha = 0.4, size = 0.9, color = "steelblue", fill = "steelblue") +
              geom_smooth(method = "lm", color = "black", fill = "gray20") +
              theme_classic() +
              theme(panel.grid = element_blank(), legend.position = c(0.85, 0.85), legend.title = element_blank(),
                    axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
              xlim(-5, 5) + ylim(-12, 12) +
              #scale_color_manual(values = c("#1B9E77", "#D95F02")) + 
              #scale_fill_manual(values = c("#1B9E77", "#D95F02")) +  # "#E66100", "#5D3A9B"
              labs(x = "Temperature anomaly (C)", y = "Lay date anomaly (days)") +
              annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.8, vjust = 1.6, size = 6) +
              annotate("text", x = -3, y = -9, label = "-0.95 days\nper degree", color = "black", size = 4, hjust = 0.5)
            
              anom_plot <- ggExtra::ggMarginal(p, type = "histogram", fill = "gray")
    
        # test for random slope across populations: no difference found
            sens_slp <- lmer(mu_lay ~ climate + (climate|pop_id), data = dev_df_out)
            sens_r <- lmer(mu_lay ~ climate + (1|pop_id), data = dev_df_out)
            
            anova(sens_slp, sens_r)

        
# Add arrival time to nests_by_year ----
  # this is down here because nests_by_year is only complete at this point
  # code from Ben Tonelli
        
      # rearrange a bit to get the right columns
        nby2 <- nests_by_year %>% dplyr::select(-pop_id, -sub_pop)
        colnames(nby2)[1] <- "pop_id"
        nby2 <- plyr::join(nby2, all_pop_locs[, c("pop_id", "source", "latitude", "longitude", "elevation_m")],
                           "pop_id", "left", "first")
        
      # join to object that has cell ids for each population
        nby2 <- plyr::join(nby2, grid_match[, c("pop_id", "cell")], "pop_id", "left", "first")
        
      # join to yearly arrival date
        nby2 <- plyr::join(nby2, arrival[, c("cell", "year", "elev_factor", "arr_GAM_mean", "arr_GAM_sd")],
                           by = c("cell", "year"), "left", "first")
        
      # add some columns
        nby2$arr_to_lay <- nby2$mu_lay - nby2$arr_GAM_mean
        
      # add grand arrival mean across years by site
        g_arr <- nby2 %>%
          group_by(pop_id) %>%
          summarise(grand_arr_mu = mean(arr_GAM_mean, na.rm = TRUE))
        nby2 <- plyr::join(nby2, g_arr, "pop_id")
        
        nby2$arr_anomaly <- nby2$arr_GAM_mean - nby2$grand_arr_mu
        nby2$lay_anomaly <- nby2$mu_lay - nby2$grand_mu_lay
        
# Put all data together ----
  # sections above pull in data from various sources. this puts them all together
    # to make one dataframe with a row per population and one with a row per nest
    # with all relevant info saved together
        
    # start with nby2 which already has a lot of info integrated
        # add ebird abundance
            apl <- all_pop_locs[, c("pop_id", "n_years", "first_year", "last_year", "rel_abundance")]
            nby3 <- plyr::join(nby2, apl, "pop_id", "left", "first")
        # add climwin results
            cwr <- climwin_results[, c("pop", "open_w", "close_w", "ref_doy", "clim_low", "clim_hi", "open_w_doy", "close_w_doy",
                                       "mw_intercept", "mw_slope", "years_with_15_nests", "c_intercept", "c_slope")]
            colnames(cwr)[1] <- "pop_id"
            nby4 <- plyr::join(nby3, cwr, "pop_id", "left", "first")
        # add spring temperature in window for each year
            pwc <- pop_with_climate[, c("pop_full", "year", "mu_tavg", "mu_tmax", "mu_tmin",
                                        "mu_tavg_i", "mu_tmax_i", "mu_tmin_i")]
            colnames(pwc)[1] <- "pop_id"
            nby4 <- plyr::join(nby4, pwc, c("pop_id", "year"), "left", "first")
        
        # add land cover
            lcall <- a_pops_lc2[, c("pop_id", "wetland", "agriculture", "urban", "open_water", "woody", "grass_shrub", "wet", "urb_barren")]
            nby5 <- plyr::join(nby4, lcall, "pop_id", "left", "first")
            nby5$pop_id <- as.factor(nby5$pop_id)
            
        # add land cover PCA
            nby5 <- plyr::join(nby5, lc_pca_with_scores, "pop_id", "left", "first")
            
    # now filter down to pop level information
            nbp <- nby5 %>%
              filter(!duplicated(pop_id)) %>% # next remove all year specific columns
              dplyr::select(-year, -included_nests, -mu_lay, -md_lay, -ref_doy, -arr_GAM_mean, -arr_GAM_sd,
                            -arr_to_lay, -arr_anomaly, -lay_anomaly)
        
# Filter data for selection gradients by population ----
    # filter only to non manipulated nests that have fitness data
         s_all <- all_nests %>% filter(exclude_fitness == "no", num_fledged > -1, num_fledged < 10)  
         s_all <- s_all %>% dplyr::select(pop_id, sub_pop, year, latitude, longitude, lay_doy, num_fledged,
                                    source, exclude_fitness, clutch_size)
        
    # summarise this to one row per site-year
        s_any <- s_all %>%
          group_by(pop_id, sub_pop, year) %>%
          summarise(n_nests_this_year = n(), 
                    mu_num_fled = mean(num_fledged),
                    first_egg_this_year = if(all(is.na(lay_doy))) NA_real_ else quantile(
                      lay_doy, probs = 0.1, na.rm = TRUE), .groups = "drop")
        
        
        s_all <- plyr::join(s_all, s_any, c("pop_id", "sub_pop", "year"), "left", "first")
        s_all$day_from_first_egg <- s_all$lay_doy - s_all$first_egg_this_year  
        
    # add in the number of nests per population and year
        s_temp <- s_all %>%
          group_by(pop_id, sub_pop, year) %>%
          dplyr::summarise(nests_this_year_pop = n())
        
        s_all <- plyr::join(s_all, s_temp, c("pop_id", "sub_pop", "year"), "left", "first")
        
    # subset to only years with 30 or more nests
        s_all <- s_all %>% filter(nests_this_year_pop > 29)
        
    # add subpop label to pop_id
        for(i in 1:nrow(s_all)){
          if(s_all$source[i] == "contribution"){
            s_all$pop_id[i] <- paste(s_all$pop_id[i], s_all$sub_pop[i], sep = "_")
          }
        }
        
    # add in the population year level data to each nest row
        s_all <- plyr::join(s_all, nby5, c("pop_id", "year"), "left", "first")
        s_all <- s_all[, !duplicated(as.list(s_all))]
        
    # set up dataframe for selection gradients
        sel_grads <- s_all[, c("pop_id", "year")]
        sel_grads <- sel_grads[!duplicated(sel_grads), ]
        
    # calculate selection separately for each population
        s_all$clutch_size <- as.numeric(s_all$clutch_size)
        for(i in 1:nrow(sel_grads)){
          temp_dat <- s_all %>% filter(pop_id == sel_grads$pop_id[i], year == sel_grads$year[i])
          if(sum(temp_dat$num_fledged > 0)){
            temp_dat$rel_fitness <- temp_dat$num_fledged / mean(temp_dat$num_fledged)
            
            m <- lm(rel_fitness ~ scale(lay_doy), data = temp_dat)
            m2 <- lm(rel_fitness ~ scale(lay_doy) + scale(clutch_size), data = temp_dat)
            
            sel_grads$gradient[i] <- coefficients(m)[2]
            sel_grads$p_value[i] <- summary(m)$coefficients[2, 4]
            
            sel_grads$gradient_clutch[i] <- coefficients(m2)[2]
            sel_grads$p_value_clutch[i] <- summary(m2)$coefficients[2, 4]
          }
        }
        
    # make a second frame with arrival times
        # summarise this to one row per site-year
          s_any2 <- s_all %>%
            group_by(pop_id, sub_pop, year) %>%
            summarise(n_nests_this_year = n(), 
                      mu_num_fled = mean(num_fledged),
                      mu_arr_doy = mean(arr_GAM_mean),
                      mu_lay_doy = mean(lay_doy, na.rm = TRUE),
                      first_egg_this_year = if(all(is.na(lay_doy))) NA_real_ else quantile(
                        lay_doy, probs = 0.1, na.rm = TRUE), .groups = "drop")
        
    # convert p to -1*log(p) for plotting
        sel_grads$neg_p <- -1 * log(sel_grads$p_value)
        sel_grads$neg_p_clutch <- -1 * log(sel_grads$p_value_clutch)
        
    # Adjust the p-values for false discovery rate
        sel_grads$p_adj <- p.adjust(sel_grads$p_value, method = "BH")
        
    # Mark significant results at FDR < 0.05
        sel_grads$significant_fdr <- sel_grads$p_adj < 0.05     
        
        pop_gradients_plot <- ggplot(sel_grads, aes(x = gradient, y = neg_p)) +
          geom_point(alpha = 0.7, aes(color = significant_fdr), size = 0.9) +
          scale_color_manual(values = c("gray70", "coral3"), labels = c("Not Significant", "Significant after FDR"), name = "") +
          geom_hline(yintercept = c(-1*log(0.05), -1*log(0.001)), linetype = "dotted", color = "gray50") +
          geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
          xlim(c(-0.5, 0.5)) +
          ylim(c(0, 30)) +
          theme_classic() +
          theme(panel.grid = element_blank(), 
                axis.title = element_text(size = 14), 
                axis.text = element_text(size = 12)) +
          guides(color = "none") +
          xlab("Selection gradient") +
          ylab("P-value (-log)") +
          annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.8, vjust = 1.6, size = 6) +
          annotate("text", x = -0.455, y = -log(0.05), label = "P = 0.05", color = "gray40", vjust = -.4) +
          #annotate("text", x = -0.47, y = -log(0.01), label = "P = 0.01", color = "gray40", vjust = -.4) +
          annotate("text", x = -0.441, y = -log(0.001), label = "P = 0.001", color = "gray40", vjust = -.4) +
          #annotate("text", x = -0.446, y = -log(0.0001), label = "P = 0.0001", color = "gray40", vjust = -.4) +
          annotate("text", x = -0.26, y = 26, label = "Selection for \n earlier breeding", color = "black") +
          annotate("text", x = 0.26, y = 26, label = "Selection for \n later breeding", color = "black")
        
    # add back in population data
        sel_grads <- plyr::join(sel_grads, nby5[, c("pop_id", "year", "latitude", "longitude",
                                                    "source", "elevation_m", "arr_GAM_mean", "arr_to_lay",
                                                    "lay_anomaly", "rel_abundance", "mu_tavg")], 
                                c("pop_id", "year"), "left", "first")
        
    # add in relative fitness and relative lay date
        s_all <- s_all %>%
          group_by(pop_id, year) %>%
          mutate(
            lay_doy_z = scale(lay_doy)[,1],  # z-score (mean = 0, sd = 1)
            rel_fitness = num_fledged / mean(num_fledged, na.rm = TRUE)
          ) %>%
          ungroup() %>%
          as.data.frame()
        
    # make a binned version of relative fitness
        # Set your bin width
            bin_width <- 0.4
        
        # Create binned variable
            s_binned <- s_all %>%
              mutate(lay_bin = cut(lay_doy_z, 
                                   breaks = seq(floor(min(lay_doy_z, na.rm = TRUE)),
                                                ceiling(max(lay_doy_z, na.rm = TRUE)),
                                                by = bin_width),
                                   include.lowest = TRUE)) %>%
              group_by(lay_bin) %>%
              summarise(
                n = sum(!is.na(rel_fitness)),
                mean_fitness = mean(rel_fitness, na.rm = TRUE),
                sd_fitness = sd(rel_fitness, na.rm = TRUE),
                sem_fitness = sd_fitness / sqrt(n),
                .groups = "drop"
              )
            
            #get midpoint of each bin
            s_binned <- s_binned %>%
              mutate(lay_bin_mid = as.numeric(sub("\\((.+),(.+)\\]", "\\1", lay_bin)) +
                       bin_width / 2) %>%
              filter(lay_bin_mid > -3 & lay_bin_mid < 4.3) %>%
              as.data.frame()
            
        # add land use PC
            s_all <- plyr::join(s_all, lc_pca_with_scores[, c("pop_id", "PC1", "PC2", "PC3")], "pop_id", "left", "first")
            
            s_all$pc_f <- cut(s_all$PC1, 3)
            
            s_all$pop_year <- interaction(s_all$pop_id, s_all$year)
            s_all$lat_s <- as.numeric(scale(s_all$latitude))
            s_all$pc1_s <- as.numeric(scale(s_all$PC1))
            s_all$pc2_s <- as.numeric(scale(s_all$PC2))
            
        # and land use to selection gradient
            gradients <- plyr::join(sel_grads, s_all[, c("pop_id", "lat_s", "PC1", "PC2", "pc1_s", "pc2_s")], "pop_id", "left", "first")
            
            grad_pc1 <- lmer(gradient ~ pc1_s + (1|pop_id), data = gradients)
            grad_pc2 <- lmer(gradient ~ pc2_s + (1|pop_id), data = gradients)
            
            gradients$pop_id <- as.factor(gradients$pop_id)
            grad_pc1_gamm <- gam(gradient ~ pc1_s + s(pop_id, bs = "re") +
                                   s(longitude, latitude), data = gradients)
            grad_pc2_gamm <- gam(gradient ~ pc2_s + s(pop_id, bs = "re") +
                                   s(longitude, latitude), data = gradients)
            
        # fit a gamm for relative fitness
            s_all$pop_id <- as.factor(s_all$pop_id)
            s_all$clutch_size <- as.numeric(s_all$clutch_size)
            rel_fit1 <- gam(rel_fitness ~ s(lay_doy_z, k = 5) +
                              s(pop_id, bs = "re") + te(longitude, latitude, k = c(5, 5)),
                            data = s_all)  
            rel_fit1b <- gam(rel_fitness ~ lay_doy_z +
                              s(pop_id, bs = "re") + te(longitude, latitude, k = c(5, 5)),
                            data = s_all)  
            rel_fit2 <- gam(rel_fitness ~ s(lay_doy_z, k = 5) + s(clutch_size, k = 5) +
                              s(pop_id, bs = "re") + te(longitude, latitude),
                            data = s_all)  
            preds <- ggeffects::ggpredict(rel_fit1, terms = "lay_doy_z [-3.5:4 by=0.1]")
            preds_c <- ggeffects::ggpredict(rel_fit2, terms = "lay_doy_z [-3.5:4 by=0.1]")
            
        # make a plot
            global_sel_plot <- ggplot(preds, aes(x = x, y = predicted)) +
              coord_cartesian(xlim = c(-2.6, 3), ylim = c(0.5, 1.2)) +
              geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
              geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
              #geom_line(data = preds_c, color = "gray50") +
              #geom_ribbon(data = preds_c, aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "gray70") +
              geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = "coral3") +
              geom_line(color = "coral3") +
              theme_classic() +
              theme(panel.grid = element_blank(), axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14)) +
              geom_point(data = s_binned, aes(x = lay_bin_mid, y = mean_fitness), 
                             size = 2, shape = 21, fill = "coral3") +
              geom_segment(data = s_binned, aes(x = lay_bin_mid, 
                                              xend = lay_bin_mid, 
                                              y = mean_fitness + sem_fitness, yend = mean_fitness-sem_fitness), 
                         color = "black") +
              ylab("Relative reproductive success") +
              xlab("Standardized lay date (z-score)") +
              annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.8, vjust = 1.6, size = 6)
            
          #select_plot <- ggpubr::ggarrange(pop_gradients_plot, global_sel_plot, nrow = 1)
           # ggsave(here::here("saved_figures", "selection_gradients.png"), select_plot, device = "png",
            #       units = "in", dpi = 300, width = 8.4, height = 4.25)
            
          #select_plot2 <- ggpubr::ggarrange(anom_plot, global_sel_plot, nrow = 1) #this one has marginal histograms
          select_plot2 <- ggpubr::ggarrange(p, global_sel_plot, nrow = 1)
          ggsave(here::here("saved_figures", "select_plot2.png"), select_plot2, device = "png",
                 units = "in", dpi = 300, width = 8.4, height = 4.25)
          ggsave(here::here("saved_figures", "select_plot2.svg"), select_plot2, device = "svg",
                 units = "in", dpi = 300, width = 8.4, height = 4.25)
          
          
        # Alternate version of this plot where temperature anomaly is stratified into thirds
                  avg_window_temp <- nby5 %>%
                    dplyr::group_by(pop_id) %>%
                    dplyr::summarise(mean_window_temp = mean(mu_tavg, na.rm = TRUE))
                  
                  s_all2 <- plyr::join(s_all, avg_window_temp, "pop_id", "left", "first")
                  s_all2$temp_anomaly <- s_all2$mu_tavg - s_all2$mean_window_temp
                  
                  s_all2$temp_cat <- cut(s_all2$temp_anomaly, 3)
                  cats <- data.frame(temp_cat = levels(s_all2$temp_cat),
                                     temp_cat2 = c("Cooler than average", "Average", "Warmer than average"))
                  s_all2 <- plyr::join(s_all2, cats, "temp_cat", "left", "first")
                  s_all2$temp_cat2 <- as.factor(s_all2$temp_cat2)
              
              # fit the model again with smooth split by temperature category
                rf_split <- gam(rel_fitness ~ s(lay_doy_z, k = 5, by = temp_cat2) + temp_cat2 +
                                  s(pop_id, bs = "re") + te(longitude, latitude, k = c(5, 5)),
                                data = s_all2)
                
              # now predict and get the values for each level of temperature anomaly
                  taf_levs <- levels(s_all2$temp_cat2)
                  newdat <- expand.grid(lay_doy_z = seq(min(s_all2$lay_doy_z, na.rm = TRUE),
                                                        max(s_all2$lay_doy_z, na.rm = TRUE),
                                                        length.out = 200),
                                        temp_cat2 = taf_levs)
                  newdat$pop_id <- s_all2$pop_id[1]
                  newdat$longitude <- mean(s_all2$longitude, na.rm = TRUE)
                  newdat$latitude <- mean(s_all2$latitude, na.rm = TRUE)
                  pr <- predict(rf_split, newdata = newdat, se.fit = TRUE, exclude = c("s(pop_id)", "te(longitude,latitude"))            
                  newdat$fit <- pr$fit
                  newdat$se <- pr$se.fit
                  newdat$lower <- newdat$fit - 2 * newdat$se
                  newdat$upper <- newdat$fit + 2 * newdat$se
                  
              # Create binned variable
                  bin_width <- 0.6
                  s_binned2 <- s_all2 %>%
                    mutate(lay_bin = cut(lay_doy_z, 
                                         breaks = seq(floor(min(lay_doy_z, na.rm = TRUE)),
                                                      ceiling(max(lay_doy_z, na.rm = TRUE)),
                                                      by = bin_width),
                                         include.lowest = TRUE)) %>%
                    group_by(lay_bin, temp_cat2) %>%
                    summarise(
                      n = sum(!is.na(rel_fitness)),
                      mean_fitness = mean(rel_fitness, na.rm = TRUE),
                      sd_fitness = sd(rel_fitness, na.rm = TRUE),
                      sem_fitness = sd_fitness / sqrt(n),
                      .groups = "drop"
                    )
                  
                #get midpoint of each bin
                  s_binned2 <- s_binned2 %>%
                    mutate(lay_bin_mid = as.numeric(sub("\\((.+),(.+)\\]", "\\1", lay_bin)) +
                             bin_width / 2) %>%
                    filter(lay_bin_mid > -3 & lay_bin_mid < 4.3) %>%
                    as.data.frame()   
                  s_binned2 <- s_binned2 %>% filter(is.na(temp_cat2) == FALSE)
                  
                  
                  # make the plot
                  
                      global_sel_plot2 <- ggplot(preds, aes(x = x, y = predicted)) +
                        coord_cartesian(xlim = c(-2.6, 3), ylim = c(0.5, 1.3)) +
                        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
                        geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
                        geom_line(data = newdat, aes(x = lay_doy_z, y = fit, color = temp_cat2), size = .8) +
                        #geom_line(data = preds_c, color = "gray50") +
                        #geom_ribbon(data = preds_c, aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "gray70") +
                        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = "gray50") +
                        geom_line(color = "black") +
                        theme_classic() +
                        theme(panel.grid = element_blank(), axis.text = element_text(size = 12),
                              axis.title = element_text(size = 14)) +
                        # geom_segment(data = s_binned, aes(x = lay_bin_mid,
                        #                                  xend = lay_bin_mid,
                        #                                  y = mean_fitness + sem_fitness, yend = mean_fitness-sem_fitness),
                        #             color = "gray50") +
                        # geom_point(data = s_binned, aes(x = lay_bin_mid, y = mean_fitness),
                        #           size = 1.5, shape = 21, fill = "gray50") +
                        ylab("Relative reproductive success") +
                        xlab("Standardized lay date (z-score)") +
                        annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.8, vjust = 1.6, size = 6) +
                        guides(color = "none") +
                        scale_color_manual(values = c("#FDAE61", "#2B83BA", "#D7191C"))
                      ggsave(here::here("saved_figures", "sel_plot_by_temp.svg"), global_sel_plot2, device = "svg",
                             units = "in", dpi = 300, width = 4.2, height = 4.25)
                  
            
 
          
                   
  # Analyze whether selection in stronger in warmer years and at higher latitudes      
      # add in temperature anomaly in window
          avg_window_temp <- nby5 %>%
            dplyr::group_by(pop_id) %>%
            dplyr::summarise(mean_window_temp = mean(mu_tavg, na.rm = TRUE))
          
          sel_grads2 <- plyr::join(sel_grads, avg_window_temp, "pop_id", "left", "first")
          sel_grads2$temp_anomaly <- sel_grads2$mu_tavg - sel_grads2$mean_window_temp
          
      # fit a model for temperature anomaly and latitude
          sel_grads2$pop_id <- as.factor(sel_grads2$pop_id)
          sel_temp <- gam(gradient ~ temp_anomaly + latitude + s(pop_id, bs = "re") + s(longitude, latitude),
                          data = sel_grads2)
          
          
      # make a plot of model effects for each of temperature anomaly and latitude
          # Hold-at values for the "other" predictors
              lat_bar <- mean(sel_grads2$latitude,  na.rm = TRUE)
              lon_bar <- mean(sel_grads2$longitude, na.rm = TRUE)
              tmp_bar <- mean(sel_grads2$temp_anomaly, na.rm = TRUE)
              
          # A valid level for pop_id (placeholder; excluded in prediction)
              pop0 <- if (is.factor(sel_grads2$pop_id)) levels(sel_grads2$pop_id)[1] else sel_grads2$pop_id[1]
              
          # Panel 1: parametric effect of temp_anomaly (latitude held at mean)
              grid_temp <- tibble(
                temp_anomaly = seq(-5,
                                   5, length.out = 200),
                latitude  = lat_bar,
                longitude = lon_bar,
                pop_id    = pop0
              )
              
              pred_temp <- predict(
                sel_temp2, newdata = grid_temp, se.fit = TRUE,
                exclude = c("s(pop_id)", "s(longitude,latitude)")
              )
              
              panel_temp <- grid_temp %>%
                transmute(
                  x = temp_anomaly,
                  predicted = pred_temp$fit,
                  conf.low  = pred_temp$fit - 1.96 * pred_temp$se.fit,
                  conf.high = pred_temp$fit + 1.96 * pred_temp$se.fit,
                  panel = "Selection gradient ~ temp anomaly (parametric)"
                )
              
          # Panel 2: parametric effect of latitude (temp_anomaly held at mean)
                  grid_lat <- tibble(
                    latitude = seq(min(sel_grads2$latitude, na.rm = TRUE),
                                   max(sel_grads2$latitude, na.rm = TRUE), length.out = 200),
                    temp_anomaly = tmp_bar,
                    longitude    = lon_bar,
                    pop_id       = pop0
                  )
                  
                  pred_lat <- predict(
                    sel_temp2, newdata = grid_lat, se.fit = TRUE,
                    exclude = c("s(pop_id)", "s(longitude,latitude)")
                  )
                  
                  panel_lat <- grid_lat %>%
                    transmute(
                      x = latitude,
                      predicted = pred_lat$fit,
                      conf.low  = pred_lat$fit - 1.96 * pred_lat$se.fit,
                      conf.high = pred_lat$fit + 1.96 * pred_lat$se.fit,
                      panel = "Selection gradient ~ latitude (parametric)"
                    )
                  
              # make the plots
                  panel_temp_p <- panel_temp %>%
                    ggplot(aes(x = x, y = predicted)) +
                    #geom_vline(xintercept = 0, color = "gray70", linetype = "dashed") +
                    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "coral3") +
                    geom_line(linewidth = 1, color = "coral3") +
                    theme_classic() +
                    theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
                    labs(x = "Temperature anomaly (C)", y = "Standardized linear \n selection gradient") +
                    coord_cartesian(xlim = c(-5, 5)) +
                    annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.8, vjust = 1.6, size = 6)
                  
                  panel_lat_p <- panel_lat %>%
                    ggplot(aes(x = x, y = predicted)) +
                    #geom_vline(xintercept = 0, color = "gray70", linetype = "dashed") +
                    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "coral3") +
                    geom_line(linewidth = 1, color = "coral3") +
                    theme_classic() +
                    theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
                    labs(x = "Latitude", y = "Standardized linear \n selection gradient") +
                    #coord_cartesian(xlim = c(-5, 5))
                    annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.8, vjust = 1.6, size = 6)
                  
                  sgrad_p <- ggpubr::ggarrange(panel_temp_p, panel_lat_p, nrow = 2)
                  ggsave(here::here("saved_figures", "sgrad_plot.png"), sgrad_p, device = "png",
                         units = "in", dpi = 300, width = 2.9, height = 4.6)
                  ggsave(here::here("saved_figures", "sgrad_plot.svg"), sgrad_p, device = "svg",
                         units = "in", dpi = 300, width = 2.9, height = 4.6)
                  
              
          
         
          
        
          
            
# Analyze & plot sliding window results ----
     
      # make plot of sensitive windows
          window_plot <- nbp %>% filter(n_years > 9) %>%
            ggplot(aes(x = (open_w - 14) * -1)) +
            geom_histogram(fill = "#B12A90FF", alpha = 0.5, color = "gray25", binwidth = 3) +
            geom_histogram(aes(x = (close_w - 14) * -1), fill = "#FCA636FF", color = "gray25", alpha = 0.5, binwidth = 3) +
            geom_histogram(aes(x = grand_arr_mu - grand_mu_lay), fill = "#0D0887FF", color = "gray25", alpha = 0.7, binwidth = 3) +
            theme_classic() +
            theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
            labs(x = "Offset from site mean lay date (days)", y = "Count of populations") +
            annotate("text", x = -60, y = 12, label = "arrival", color = "#0D0887FF", size = 5) +
            annotate("segment", x = -21, xend = 0, y = 40, yend = 40, color = "black", linewidth = 1.5) +
            annotate("point", x = c(-21, 0), y = c(40, 40), color = "black", size = 2.5) +
            annotate("text", x = -30, y = 28, label = "window \n open", color = "#B12A90FF", size = 5) +
            annotate("text", x = -9, y = 28, label = "window \n close", color = "#FCA636FF", size = 5) +
            annotate("text", x = -10.5, y = 45, label = "consensus \n window", color = "black", size = 5) +
            coord_cartesian(ylim = c(0, 48)) +
            annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.8, vjust = 1.6, size = 6) 
            
        #windows <- ggpubr::ggarrange(window_plot, anom_plot, nrow = 1)
        
        #ggsave(here::here("saved_figures", "windows.png"), windows, device = "png", units = "in",
        #       width = 9, height = 4)
        
    # alternative window plot
        window_plot2 <- nbp %>% filter(n_years > 9) %>%
          ggplot() +
          geom_segment(aes(x = latitude, xend = latitude, y = open_w_doy, yend = close_w_doy), color = "gray50", size = 0.2) +
          geom_segment(aes(x = latitude, xend = latitude, y = open_w_doy, yend = grand_arr_mu), linetype = "11", color = "gray80", size = 0.2) +
          geom_point(aes(x = latitude, y = open_w_doy), shape = 21, fill = "#B12A90FF") +
          geom_point(aes(x = latitude, y = close_w_doy), shape = 21, fill = "#FCA636FF") +
          geom_point(aes(x = latitude, y = grand_arr_mu), shape = 21, fill = "#0D0887FF") +
          theme_classic() +
          theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
          annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.8, vjust = 1.6, size = 6) +
          labs(x = "Latitude", y = "Day of year")
        
        windows2 <- ggpubr::ggarrange(window_plot2, window_plot, nrow = 1)
        ggsave(here::here("saved_figures", "windows2.png"), windows2, device = "png", units = "in",
               width = 9, height = 4)
        
    # check for any evidence of variation in window open, close, slope by latitude or by land use
    # these models are only reported very briefly
        
        win_m1 <- gam(open_w ~ s(latitude, k = 5) + te(longitude, latitude), data = nbp)
        win_m2 <- gam(open_w ~ s(PC1, k = 5) + te(longitude, latitude), data = nbp)
        win_m3 <- gam(open_w ~ s(PC2, k = 5) + te(longitude, latitude), data = nbp)
        
        cls_m1 <- gam(close_w ~ s(latitude, k = 5) + te(longitude, latitude), data = nbp)
        cls_m2 <- gam(close_w ~ s(PC1, k = 5) + te(longitude, latitude), data = nbp)
        cls_m3 <- gam(close_w ~ s(PC2, k = 5) + te(longitude, latitude), data = nbp)
        
        slp_m1 <- gam(mw_slope ~ s(latitude, k = 5) + te(longitude, latitude), data = nbp)
        slp_m2 <- gam(mw_slope ~ s(PC1, k = 5) + te(longitude, latitude), data = nbp)
        slp_m3 <- gam(mw_slope ~ s(PC2, k = 5) + te(longitude, latitude), data = nbp)
      
 
# Plot arrival date and waiting time results ----
        
    # latitude vs. waiting time between arrival and egg laying  
        ggplot(nbp, aes(x = latitude, y = grand_mu_lay - grand_arr_mu)) +
          geom_point(alpha = 0.5) +
          geom_smooth(method = "lm", color = "slateblue", fill = "slateblue") +
          theme_classic() +
          theme(panel.grid = element_blank(), axis.title = element_text(size = 14),
                axis.text = element_text(size = 12)) +
          labs(x = "latitude", y = "average arrival to \n average first egg (days)")
        
    # Arrival anomaly vs. lay anomaly
        # fit a model
            arr_lay_m <- gam(lay_anomaly ~ arr_anomaly + s(pop_id, bs = "re") + te(longitude, latitude), data = nby5)
            alm_pred <- ggeffects::ggpredict(model = arr_lay_m, terms = "arr_anomaly [-12:12 by=1]")
        
        # make a figure
            arr_lay_p <- ggplot(nby5, aes(x = arr_anomaly)) +
              geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
              geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
              geom_point(color = "slateblue", alpha = 0.4, aes(y = lay_anomaly)) +
              theme_classic() +
              theme(axis.title = element_text(size = 14), 
                    axis.text = element_text(size = 12), panel.grid = element_blank()) +
              xlab("Arrival date anomaly (days)") +
              ylab("Lay date anomaly (days)") +
              coord_cartesian(ylim = c(-11, 11), xlim = c(-11, 11)) +
              geom_line(data = alm_pred, aes(x = x, y = predicted), color = "black") +
              geom_ribbon(data = alm_pred, aes(x = x, ymin = conf.low, ymax = conf.high), 
                            alpha = 0.2, fill = "gray30") +
              annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.8, vjust = 1.6, size = 6)
        
      # arrival, window, lay in same plot
            nby_p <- nby5 #%>% filter(elevation_m < 900)
            nbp_p <- nbp #%>% filter(elevation_m < 900)
            
            # fit three models
                arr_gm <- gam(arr_GAM_mean ~ s(latitude, k = 5) + s(elevation_m, k = 5) + s(pop_id, bs = "re"), data = nby_p)
                win_gm <- gam(open_w_doy ~ s(latitude, k = 5) + s(elevation_m, k = 5), data = nbp_p)
                lay_gm <- gam(mu_lay ~ s(latitude, k = 5) + s(elevation_m, k = 5) + s(pop_id, bs = "re"), data = nby_p)
                
            # get predictions for each model
                arr_pred <- ggeffects::ggpredict(arr_gm, terms = "latitude [32:66 by=1]")
                win_pred <- ggeffects::ggpredict(win_gm, terms = "latitude [32:66 by=1]")
                lay_pred <- ggeffects::ggpredict(lay_gm, terms = "latitude [32:66 by=1]")
                
            # make a plot
                waiting_p <- ggplot(nbp_p, aes(x = latitude)) +
                  geom_point(aes(y = grand_arr_mu), color = "#0D0887FF", alpha = 0.7, size = 0.9, shape = 17) +
                  geom_point(aes(y = open_w_doy), color = "#B12A90FF", alpha = 0.7, size = 0.9, shape = 16) +
                  geom_point(aes(y = grand_mu_lay), color = "forestgreen", alpha = 0.7, size = 0.9, shape = 15) +
                  theme_classic() +
                  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
                  xlab("Latitude") +
                  ylab("Day of year \n delete") +
                  coord_cartesian(xlim = c(34, 64)) +
                  geom_line(data = arr_pred, aes(x = x, y = predicted), color = "#0D0887FF") +
                  geom_line(data = win_pred, aes(x = x, y = predicted), color = "#B12A90FF") +
                  geom_line(data = lay_pred, aes(x = x, y = predicted), color = "forestgreen") +
                  geom_ribbon(data = arr_pred, aes(x = x, ymin = conf.low, ymax = conf.high), 
                              alpha = 0.2, fill = "#0D0887FF") +
                  geom_ribbon(data = win_pred, aes(x = x, ymin = conf.low, ymax = conf.high), 
                              alpha = 0.2, fill = "#B12A90FF") +
                  geom_ribbon(data = lay_pred, aes(x = x, ymin = conf.low, ymax = conf.high), 
                              alpha = 0.2, fill = "forestgreen") +
                  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.8, vjust = 1.6, size = 6) +
                  annotate("text", x = 58, y = 50, label = "Arrival date", color = "#0D0887FF", size = 5) +
                  annotate("text", x = 58, y = 60, label = "Window open", color = "#B12A90FF", size = 5) +
                  annotate("text", x = 58, y = 70, label = "Lay date", color = "forestgreen", size = 5)
                
              arrival_plot <- ggpubr::ggarrange(arr_lay_p, waiting_p, nrow = 1)
              
              ggsave(here::here("saved_figures", "arrival_lay.png"), arrival_plot, device = "png",
                     units = "in", dpi = 300, width = 8.4, height = 4.25)
              
              
          # reduction in waiting time with latitude
              time_gm <- gam(mu_lay - arr_GAM_mean ~ s(latitude, k = 5) + s(elevation_m, k = 5) + s(pop_id, bs = "re"), data = nby_p)
              time_pred <- ggeffects::ggpredict(time_gm, terms = "latitude [32:66 by=1]")
              
              
              time_p <- ggplot(nbp_p, aes(x = latitude)) +
                geom_point(aes(y = grand_mu_lay - grand_arr_mu, color = elevation_m, size = n_years), alpha = 0.6) +
                scale_color_viridis_c(option = "D", end = 0.9, name = "Elevation (m)") +
                theme_classic() +
                theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
                xlab("Latitude") +
                ylab("Arrival to egg laying (days) \n delete") +
                labs(size = "Years of data") +
                #guides(color = "none", size = "none") + 
                coord_cartesian(xlim = c(34, 64)) +
                geom_line(data = time_pred, aes(x = x, y = predicted), color = "gray10") +
                geom_ribbon(data = time_pred, aes(x = x, ymin = conf.low, ymax = conf.high), 
                            alpha = 0.3, fill = "gray30") 
              
              ggsave(here::here("saved_figures", "time.png"), time_p, device = "png",
                     units = "in", dpi = 300, width = 5.6, height = 4.4)
              
              fig3ab <- ggpubr::ggarrange(waiting_p, time_p, nrow = 1)
              # ggsave(here::here("saved_figures", "fig3ab.svg"), fig3ab, device = "svg",
              #       units = "in", dpi = 300, width = 8.4, height = 4.25)
              # ggsave(here::here("saved_figures", "fig3blegend.svg"), time_p, device = "svg",
              #                   units = "in", dpi = 300, width = 4.2, height = 4.25)
        
# Change in egg laying time range wide ----
    # model of overall change in egg laying date across all sites
          
        # add centered year
              s_all$year_c <- s_all$year - mean(s_all$year)
              nby5$year_c <- nby5$year - mean(nby5$year)    
              
          m_rs <- gam(mu_lay ~ year_c * latitude + s(longitude, latitude) +
                                 s(pop_id, bs = "re") + s(pop_id, by = year_c, bs = "re"), 
                                data = nby5, method = "REML")
          
    # make a figure from the model
          # ---- settings 
              h    <- 1        # 1-year finite difference (±0.5)
              half <- h/2
              lon0 <- mean(nby5$longitude, na.rm = TRUE)
          
          # Identify the pop_id smooth labels to exclude when making the population-level curve
              smooth_labels <- sapply(m_rs$smooth, function(sm) sm$label)
              re_labs <- grep("pop_id", smooth_labels, value = TRUE)  # catches s(pop_id) and s(pop_id):year_c
          
          # Pull coefs & vcov 
              b <- coef(m_rs)
              V <- vcov(m_rs)
          
          #  Population-level slope vs latitude (exclude site REs) 
              lat_seq <- seq(min(nby5$latitude), max(nby5$latitude), length.out = 500)
              
              nd_plus  <- data.frame(year_c = +half, latitude = lat_seq, longitude = lon0,
                                     pop_id  = levels(nby5$pop_id)[1])
              nd_minus <- nd_plus; nd_minus$year_c <- -half
              
              Xp <- predict(m_rs, newdata = nd_plus,  type = "lpmatrix", exclude = re_labs)
              Xm <- predict(m_rs, newdata = nd_minus, type = "lpmatrix", exclude = re_labs)
              
              Xd <- (Xp - Xm) / h
              mu <- as.numeric(Xd %*% b)
              se <- sqrt(rowSums((Xd %*% V) * Xd))
              
              cont <- tibble(
                latitude = lat_seq,
                slope_decade = mu * 10,
                lwr_decade   = (mu - 1.96*se) * 10,
                upr_decade   = (mu + 1.96*se) * 10
              )
          
          # Site-specific random slopes at site mean latitude 
              site_xy <- nby5 %>%
                group_by(pop_id) %>%
                summarise(latitude = mean(latitude, na.rm = TRUE),
                          longitude = mean(longitude, na.rm = TRUE),
                          .groups = "drop") %>%
                filter(is.finite(latitude), is.finite(longitude))
              
              ndp_site <- transform(site_xy, year_c = +half)
              ndm_site <- transform(site_xy, year_c = -half)
              
              Xp_site <- predict(m_rs, newdata = ndp_site, type = "lpmatrix")  # includes REs
              Xm_site <- predict(m_rs, newdata = ndm_site, type = "lpmatrix")
              Xd_site <- (Xp_site - Xm_site) / h
              
              mu_site <- as.numeric(Xd_site %*% b)
              se_site <- sqrt(rowSums((Xd_site %*% V) * Xd_site))
          
              site_pts <- site_xy %>%
                mutate(slope_decade = mu_site * 10,
                       lwr_decade   = (mu_site - 1.96*se_site) * 10,
                       upr_decade   = (mu_site + 1.96*se_site) * 10)
          
          # Plot: continuous line + CI, plus per-site points 
              panel_a <- ggplot() +
                
                geom_errorbar(data = site_pts, color = "gray55",
                              aes(x = latitude, ymin = lwr_decade, ymax = upr_decade),
                              width = 0, alpha = 0.25, inherit.aes = FALSE, show.legend = FALSE) +
                geom_point(data = site_pts, color = "gray35", 
                           aes(x = latitude, y = slope_decade),
                           size = 1.4, inherit.aes = FALSE, show.legend = FALSE, alpha = 0.4) +
                geom_ribbon(data = cont, fill = "coral3",
                            aes(x = latitude, ymin = lwr_decade, ymax = upr_decade),
                            alpha = 0.3, inherit.aes = FALSE) +
                geom_line(data = cont, color = "coral3",
                          aes(x = latitude, y = slope_decade),
                          linewidth = 1.1, inherit.aes = FALSE) +
                geom_hline(yintercept = 0, linetype = "dashed") +
                labs(x = "Latitude",
                     y = "Change in laying date \n (days per decade 1975-2024)") +
                theme_classic() +
                theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
                guides(color = "none", fill = "none") +
                annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.8, vjust = 1.6, size = 6)
              
          
              
# Population trend over time ----
  # using data from BBS taken from USGS model website for each state/province
      bbs_dat <- read.delim(here::here("raw_data", "bbs_trends_80_22.txt"))  
      bbs_dat$state_province[7] <- "Québec" #just fixing a naming issue
  
  # add state/province to each population
        # get locations from natural earth
          st_prov <- rnaturalearth::ne_states(
            country     = c("United States of America", "Canada"),
            returnclass = "sf"
          ) %>%
            transmute(
              state_province = name,            # e.g., "New York", "Ontario"
              iso_3166_2     = iso_3166_2,      # e.g., "US-NY", "CA-ON"
              country        = admin            # "United States of America" or "Canada"
            )
        
        # make shape file of population points
          pts <- st_as_sf(
            all_pop_locs,
            coords = c("longitude", "latitude"),
            crs = 4326,
            remove = FALSE
          )
        
        # spatial join states and points
          joined_state <- st_join(pts, st_prov, join = st_within)
          
        # deal with any just out of bounds
          if (anyNA(joined_state$state_province)) {
            idx <- which(is.na(joined_state$state_province))
            nn  <- st_nearest_feature(joined_state[idx, ], st_prov)
            joined_state$state_province[idx] <- st_prov$state_province[nn]
            joined_state$iso_3166_2[idx]     <- st_prov$iso_3166_2[nn]
            joined_state$country[idx]        <- st_prov$country[nn]
            joined_state$assigned_by         <- "nearest"
            joined_state$assigned_by <- ifelse(is.na(joined_state$assigned_by), "within", joined_state$assigned_by)
          }
          
          js2 <- st_drop_geometry(joined_state[, c("pop_id", "state_province")])
          
    # join trend back to all pop locations
          all_pop_decline <- plyr::join(all_pop_locs, js2, "pop_id", "left", "first")
          all_pop_decline2 <- plyr::join(all_pop_decline, bbs_dat, "state_province", "left", "first")
          all_pop_decline2 <- plyr::join(all_pop_decline2, nbp[, c("pop_id", "grand_arr_mu", "grand_mu_lay")],
                                         "pop_id", "left", "first")
          all_pop_decline2$arr_to_lay <- all_pop_decline2$grand_mu_lay - all_pop_decline2$grand_arr_mu
          
    # fit a model
          # all_pop_locs$stprov <- as.factor(all_pop_locs$state_province)
          # trend_m <- gam(trend_80_22 ~ latitude + s(longitude, latitude) + s(stprov, bs = "re"), data = all_pop_locs)
          
          
    # aglommerate to state level
          ap_state <- all_pop_decline2 %>%
            dplyr::group_by(state_province) %>%
            dplyr::summarise(n_nests = sum(n_nests), trend_80_22 = mean(trend_80_22), n = n(),
                             latitude = mean(latitude), longitude = mean(longitude),
                             low_bbs = mean(low_ci), high_bbs = mean(high_ci),
                             arr_to_lay = mean(arr_to_lay))
          
          trend_m1 <- lm(trend_80_22 ~ latitude, data = ap_state)
          
          t_preds <- ggpredict(trend_m1, terms = "latitude [32:66 by=1]")
          
          
    # make a plot
        panel_b <- ggplot(ap_state, aes(x = latitude, y = trend_80_22)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
          
          geom_ribbon(data = t_preds, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), 
                      alpha = 0.3, fill = "coral3") +
          geom_segment(aes(x = latitude, xend = latitude,
                           y = low_bbs, yend = high_bbs), color = "gray65") +
          #geom_point(aes(size = n_nests), shape = 21, fill = "gray70", color = "gray20") +
          geom_point(shape = 21, fill = "gray70", color = "gray20") +
          geom_line(data = t_preds, aes(x = x, y = predicted), color = "coral3", linewidth = 1.1) +
          scale_size_continuous(range = c(.8, 4)) +
          #geom_smooth(method = "lm", color = "steelblue", fill = "steelblue") +
          theme_classic() +
          theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12),
                legend.position = c(0.85, 0.8)) +
          #guides(size = guide_legend(title = "nests")) +
          labs(x = "Latitude", y = "Breeding Bird Survey trend \n (yearly % change 1966-2022)") +
          annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.8, vjust = 1.6, size = 6)
        
        change_bbs <- ggpubr::ggarrange(panel_a, panel_b, nrow = 1)
        
        ggsave(here::here("saved_figures", "change_bbs.png"), change_bbs, device = "png",
               units = "in", dpi = 300, width = 8.4, height = 4.25)
        
        ggsave(here::here("saved_figures", "fig3CD.svg"), change_bbs, device = "svg",
               units = "in", dpi = 300, width = 8.4, height = 4.25)

                  
# Change in temperature over time Berkeley Earth data ----
    # Will need to download Berkeley earth data and place it in the directory here    
        
        
    library(ncdf4)
    library(geosphere)
    
  # location of the downloaded data
        best_nc <- here::here("berkeley", "Complete_TAVG_Daily_EqualArea.nc")
        months_keep <- 3:7
        
  # settings
        k_neighbors <- 1
        idw_power <- 2
            
  # functions
      # dates from nc
        read_var_cell <- function(nc, var, cell_index, ncell_total = NULL){
          v <- nc$var[[var]]
          if (is.null(v)) stop(sprintf("Variable '%s' not found in file.", var))
          dim_names <- sapply(v$dim, `[[`, "name")
          dim_lens  <- sapply(v$dim, `[[`, "len")
          
          # Guess which dimension is "cell"
          if (is.null(ncell_total)) {
            # try to infer from file-wide lon/lat vectors
            lon <- try(ncvar_get(nc, "longitude"), silent = TRUE)
            ncell_total <- if (inherits(lon, "try-error")) NA_integer_ else length(lon)
          }
          cell_dim <- which(
            grepl("cell", dim_names, ignore.case = TRUE) |
              grepl("grid", dim_names, ignore.case = TRUE) |
              (!is.na(ncell_total) & dim_lens == ncell_total)
          )
          if (length(cell_dim) != 1L) {
            stop(sprintf("Couldn't identify cell dimension for '%s'. Dims: %s",
                         var, paste(sprintf("%s=%d", dim_names, dim_lens), collapse=", ")))
          }
          
          # Build 1-based start/count vectors
          start <- rep(1L, length(dim_lens))
          count <- dim_lens
          start[cell_dim] <- as.integer(cell_index)
          count[cell_dim] <- 1L
          
          # Read and drop singleton dims -> vector
          as.numeric(drop(ncvar_get(nc, var, start = start, count = count)))
        }
        
        get_dates_from_nc <- function(nc){
          has <- \(nm) nm %in% names(nc$var) || nm %in% names(nc$dim)
          if (has("year") && has("day_of_year")){
            yr  <- ncvar_get(nc, "year")
            doy <- ncvar_get(nc, "day_of_year")
            as.Date(sprintf("%04d-01-01", yr)) + (doy - 1L)
          } else if (has("time")){
            tt    <- ncvar_get(nc, "time")
            units <- ncatt_get(nc, "time", "units")$value  # e.g., "days since 1870-01-01"
            origin <- as.Date(sub(".*since\\s+", "", units))
            origin + tt
          } else stop("Couldn't find time in BEST file.")
        }
        
        read_cells_tbl <- function(nc){
          lon <- ncvar_get(nc, "longitude")
          lat <- ncvar_get(nc, "latitude")
          if (max(lon, na.rm = TRUE) > 180) lon <- ifelse(lon > 180, lon - 360, lon)
          tibble(cell = seq_along(lon), lon = lon, lat = lat)
        }
            
    # open file and read it in
        nc <- nc_open(best_nc); on.exit(nc_close(nc))
        
        dates_all <- get_dates_from_nc(nc)
        keep_time <- dates_all >= as.Date("1950-01-01") & dates_all <= as.Date("2024-12-31")
        if (!is.null(months_keep)) keep_time <- keep_time & (month(dates_all) %in% months_keep)
        
        dates   <- dates_all[keep_time]
        doy     <- yday(dates)
        doy_fix <- pmin(doy, 365)
        
    # map each population to equal area cells
        cells_tbl <- read_cells_tbl(nc)
        
        nearest_k <- function(lat, lon, k = 1){
          d <- distHaversine(cbind(cells_tbl$lon, cells_tbl$lat), c(lon, lat)) # meters
          ord <- order(d)[seq_len(k)]
          list(idx = cells_tbl$cell[ord], dist_m = as.numeric(d[ord]))
        }
        
        nbp_map <- nbp %>%
          mutate(neigh = pmap(list(latitude, longitude), nearest_k, k = k_neighbors))
        
        unique_cells <- sort(unique(unlist(lapply(nbp_map$neigh, `[[`, "idx"))))
        
    # read anomalies and climatology for each needed cell
        cell_tas <- setNames(vector("list", length(unique_cells)), unique_cells)
        ncell_total <- nrow(read_cells_tbl(nc))  # pass to reader for robustness
        
        clim_cache <- NULL  # cache the 365-day climatology per cell if needed (kept simple here)
        
        for (cc in unique_cells) {
          anom_full <- read_var_cell(nc, "temperature", cc, ncell_total)  # full time
          clim_365  <- read_var_cell(nc, "climatology", cc, ncell_total)  # length 365
          tas <- anom_full[keep_time] + clim_365[doy_fix]                 # absolute °C
          cell_tas[[as.character(cc)]] <- tas
        }
        
    # per population daily series
        pop_daily <- map_dfr(seq_len(nrow(nbp_map)), function(i){
          pid <- nbp_map$pop_id[i]
          nn  <- nbp_map$neigh[[i]]$idx
          dm  <- nbp_map$neigh[[i]]$dist_m
          wts <- if (k_neighbors == 1) 1 else 1 / (dm^idw_power)
          wts <- wts / sum(wts)
          
          ts_mat <- do.call(cbind, lapply(nn, \(cidx) cell_tas[[as.character(cidx)]]))
          tas_C  <- if (is.null(dim(ts_mat))) ts_mat else as.numeric(ts_mat %*% wts)
          
          tibble(pop_id = pid, date = dates, tas_C = tas_C)
        })
        
    # summarize average within sensitive window
        pop_year <- pop_daily %>%
          mutate(year = year(date), doy = yday(date)) %>%
          left_join(nbp %>% dplyr::select(pop_id, open_w_doy, close_w_doy), by = "pop_id") %>%
          filter(
            (close_w_doy >= open_w_doy & doy >= open_w_doy & doy <= close_w_doy) |
              (close_w_doy <  open_w_doy & (doy >= open_w_doy | doy <= close_w_doy))
          ) %>%
          group_by(pop_id, year) %>%
          summarise(tavg_C = mean(tas_C, na.rm = TRUE),
                    n_days = dplyr::n(), .groups = "drop")
        
        nbp$common_open <- nbp$grand_mu_lay - 21.1
        nbp$common_close <- nbp$grand_mu_lay + 0.1
        pop_year <- pop_daily %>%
          mutate(year = year(date), doy = yday(date)) %>%
          left_join(nbp %>% dplyr::select(pop_id, common_open, common_close), by = "pop_id") %>%
          filter(
            (common_close >= common_open & doy >= common_open & doy <= common_close) |
              (common_close <  common_open & (doy >= common_open | doy <= common_close))
          ) %>%
          group_by(pop_id, year) %>%
          summarise(tavg_C = mean(tas_C, na.rm = TRUE),
                    n_days = dplyr::n(), .groups = "drop")
        
      # calculate on per 25 year
        pop_era <- pop_year %>%
          mutate(era = case_when(
            year <= 1974 ~ "1950-1974",
            year <= 1999 ~ "1975-1999",
            TRUE         ~ "2000-2024"
          )) %>%
          group_by(pop_id, era) %>%
          summarise(tavg_C = mean(tavg_C, na.rm = TRUE),
                    .groups = "drop")
        
      # pop_era wide
        pop_era_w <- pivot_wider(pop_era, names_from = era, values_from = c(tavg_C))
        colnames(pop_era_w)[2:4] <- c("a_1950", "b_1975", "c_2000")
        pop_era_w <- as.data.frame(pop_era_w)
        
        bnbp <- plyr::join(nbp, pop_era_w, "pop_id", "left")
        
        t_delta <- gam(c_2000 - a_1950 ~ latitude + s(longitude, latitude), data = bnbp)

# Change in temp variability over time ----
    py2 <- plyr::join(pop_year, nbp[, c("pop_id", "latitude", "longitude")], "pop_id", "left")    
        
    py2$decade <- floor((py2$year - 1950) / 15) * 15 + 1950
    
    py3 <- py2 %>%
      group_by(pop_id, latitude, longitude) %>%
      summarise(n = n(), sd = sd(tavg_C, na.rm = TRUE), mu = mean(tavg_C, na.rm = TRUE))
    
    py3$lat_band <- cut(py3$latitude, 3)
    
    sd_delta <- gam(sd ~ scale(latitude) + s(longitude, latitude), data = py3)
        
        
              
        