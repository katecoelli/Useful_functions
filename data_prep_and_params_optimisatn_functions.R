
######################################################################################
# File Overview
#######################################################################################
# Author: Sabastine 

# This script contains helper functions for:

# 1. Processing of input data required for running RothC. These datasets
# include: SILO climate data, Landsat NDVI, MODIS ET, soil moisture output
# from a water balance model, and mesaured soil fractions. While the data
# are for the Muttama catchment, the code, with slight modification, can
# be used for sites with similar data.

# 2. Running RothC and parameter optimisation (Both for the conventional and Bayesian
# approaches of calibaration).

#######################################################################################
# Data processing code
#######################################################################################

extract_clim_vars <- function(climate_files, clim_var, site_sp) {
  # Extract climate information from SILO data

  # Parameters
  # ----------
  # climate_files: Vector
  #   Names of climate files

  # clim_var: String
  #   Name of climate variable to extract

  # site_sp: Shapefile
  #   Sites locations

  # Returns
  # -------
  # Dataframe
  #   Outputs monthly time series of the given climate variable

  var_files <- climate_files[grep(clim_var, climate_files)]
  var_rasters <- var_files %>%
    map(brick)
  ras <- stack(var_rasters)

  # Extract monthly rainfall values
  if (clim_var == "daily_rain") {
    return(
      # Extract and aggregate daily climate (min or maximum temperatures or) observations
      raster::extract(ras, site_sp, fun = mean, df = TRUE, sp = TRUE)@data %>%
        pivot_longer(starts_with("X"),
          names_to = c(".value", "date"),
          names_pattern = "(.)(.+)"
        ) %>%
        mutate_at("date", ~ format(ymd(gsub("X", "", .)), "%Y-%m")) %>%
        dplyr::group_by(site_id, date) %>%
        summarise(X = sum(X, na.rm = TRUE)) %>%
        ungroup() %>%
        rename_at(vars("X"), ~ as_string(clim_var))
    )
  } else {
    # Extract and aggregate daily climate (min or maximum temperatures) observations
    raster::extract(ras, site_sp, fun = mean, df = TRUE, sp = TRUE)@data %>%
      pivot_longer(starts_with("X"),
        names_to = c(".value", "date"),
        names_pattern = "(.)(.+)"
      ) %>%
      mutate_at("date", ~ format(ymd(gsub("X", "", .)), "%Y-%m")) %>%
      dplyr::group_by(site_id, date) %>%
      summarise(X = mean(X, na.rm = TRUE)) %>%
      ungroup() %>%
      rename_at(vars("X"), ~ as_string(clim_var)) %>%
      dplyr::select(as_string(clim_var))
  }
}


get_site_data <- function(site, monthly_data, soil_data) {
  # Extract data for the a given site

  # Parameters
  # ----------
  # site: String
  #   Site id

  # monthly_data: Dataframe
  #   Monthly timeseries of vegation, climate and soil moisture for study sites

  # soil_data: Dataframe
  #   sampled soil information for the study sites

  # Returns
  # -------
  # List
  #   Vegation, climate, soil moisture monthly timeseries, and soil information of the site

  # get site's bucket size
  bucket_size <- soil_data %>%
    dplyr::filter(site_id == site) %>%
    dplyr::filter(year == min(year)) %>%
    pull(bucket_size)

  # get site's initial C pools
  Cpools_initial <- soil_data %>%
    dplyr::filter(site_id == site) %>%
    dplyr::filter(year == min(year)) %>%
    dplyr::select(TOC, RPM, HUM, IOM)

  # get site's final C pools
  Cpools_final <- soil_data %>%
    dplyr::filter(site_id == site) %>%
    dplyr::filter(year == max(year)) %>%
    dplyr::select(TOC, RPM, HUM, IOM)

  # get site's clay content
  pclay <- soil_data %>%
    dplyr::filter(site_id == site) %>%
    dplyr::filter(year == min(year)) %>%
    pull(clay)

  # GET TEMPORAL VARIABLES
  # get site's sm series
  soil_moisture <- monthly_data %>%
    dplyr::filter(site_id == site) %>%
    pull(SW_buckets123)

  # get NDVI and cover info
  NDVI_cover <- monthly_data %>%
    dplyr::filter(site_id == site) %>%
    dplyr::select(NDVI, cover)

  # get MODIS ET
  ET <- monthly_data %>%
    dplyr::filter(site_id == site) %>%
    pull(ET)

  # get temperature
  temperature <- monthly_data %>%
    dplyr::filter(site_id == site) %>%
    pull(mean_temp)

  # get rain
  rain <- monthly_data %>%
    dplyr::filter(site_id == site) %>%
    pull(rain)

  # get evap
  evap <- monthly_data %>%
    dplyr::filter(site_id == site) %>%
    pull(evap)

  # get time
  time <- monthly_data %>%
    dplyr::filter(site_id == site) %>%
    pull(yr_mon)

  return(list(
    "bucket_size" = bucket_size,
    "Cpools_initial" = Cpools_initial,
    "Cpools_final" = Cpools_final,
    "pclay" = pclay,
    "soil_moisture" = soil_moisture,
    "NDVI_cover" = NDVI_cover,
    "ET" = ET,
    "temperature" = temperature,
    "rain" = rain,
    "evap" = evap,
    "time" = time
  ))
}


##############################################################################
# Helper functions for running, evaluating and optimising parameters for RothC
# results across study sites
##############################################################################

# The function  is designed to run with the parameters placed in a vector.
# Using a vector to hold the parameters (rather than say a named list) is because
# the optimization functions we are using- "optim" (R base package) and routines
# from the BayesianTools package- are designed to use vector
calc_site_C_turnover <- function(data_list,
                                 params) {
  # Compute rothC for a site

  # Parameters
  # ----------
  # data_list: List
  #   Monthly input data and soil information to use to run RothC

  # params: Vector
  #   Numeric values of RothC parameters (e.g. the C fractions decomposition rates)

  # Returns
  # -------
  # Vector
  #   C fraction and total C "DPM", "RPM", "BIO", "HUM", "IOM", "CO2", TOC

  # Calculate temperature modifying rate
  fT <- calc_fT_RothC(data_list$temperature)

  # Calculate moisture modifying rate using data froom soil moisture
  fM <- calc_fM_WB(
    data_list$soil_moisture,
    data_list$bucket_size,
    data_list$NDVI_cover$cover
  )$b

  # Calculate cover modifying rate
  fC <- calc_fC_RothC(data_list$NDVI_cover$cover)

  # Calculate abc
  abc <- fT * fM * fC

  # Compute plant residue
  plant_residue <- calc_Cinput(
    data_list$NDVI_cover$NDVI,
    data_list$ET,
    params[1]
  )

  # FYM
  fym <- rep(0, length(plant_residue))

  # Put plant residue and fmy in a tibble
  C_in <- tibble(plant_residue, fym)

  # Run RothC
  model <- calc_RothC_model(
    data_list$time,
    initial_C_pools = list(
      DPM = 0,
      RPM = data_list$Cpools_initial$RPM,
      BIO = 0,
      HUM = data_list$Cpools_initial$HUM,
      IOM = data_list$Cpools_initial$IOM,
      CO2 = 0
    ),
    ks = list(
      DPM = params[2],
      RPM = params[3],
      BIO = params[4],
      HUM = params[5],
      IOM = 0,
      CO2 = Inf
    ),
    C_in,
    data_list$pclay,
    abc,
    DPM_RPM_ratio = params[6]
  )

  # Return the C fractions for the last month of simulation
  as.vector(unlist((tail(model, n = 1)[,1:7]))) # "DPM", "RPM", "BIO", "HUM", "IOM", "CO2", TOC
}


get_sites_C_turnover <- function(params,
                                 monthly_data,
                                 soil_data) {
  # Compute RothC for sites

  # Parameters
  # ----------
  # params: Vector
  #   Numeric values of RothC parameters (e.g. the C fractions decomposition rates)

  # monthly_data: Dataframe
  #   Monthly timeseries of vegation, climate and soil moisture for study sites

  # soil_data: Dataframe
  #   Sampled soil information for the study sites

  # Returns
  # -------
  # Dataframe
  #   C fractions and total C ("DPM", "RPM", "BIO", "HUM", "IOM", "CO2", TOC)

  # Get sites' names from soil data
  sites <- unique(soil_data$site_id)
  # Instantiate result list
  pred_results <- vector("list", length(sites))

  # Get observed and simulated C
  for (site_no in seq_along(pred_results)) {
    ## Get site data
    data_list <- get_site_data(sites[site_no], monthly_data, soil_data)

    pred_results[[site_no]] <- c(sites[site_no], calc_site_C_turnover(data_list, params))
  }
  df <- do.call(rbind.data.frame, pred_results)
  names(df) <- c(
    "site",
    "DPM",
    "RPM",
    "BIO",
    "HUM",
    "IOM",
    "CO2",
    "TOC"
  )
  df %>%
    mutate_at(vars(!site), as.numeric)
}


#############################################################################
#  Helper functions specific to Optimisation with 'optim'
#############################################################################

get_sites_obs_pred_C <- function(params,
                                 sites,
                                 monthly_data,
                                 soil_data) {
  # Put together observed and simulated C fractions for all sites

  # Parameters
  # ----------
  # params: Vector
  #   Numeric values of RothC parameters (e.g. the C fractions decomposition rates)

  # sites: Vector
  #   Site ids

  # monthly_data: Dataframe
  #   Monthly timeseries of vegation, climate and soil moisture for study sites

  # soil_data: Dataframe
  #   Sampled soil information for the study sites

  # Returns
  # -------
  # Dataframe
  #   Simulated and predicted C fractions and total C

  # Initiate result list
  obs_pred_results <- vector("list", length(sites))

  # Get observed and simulated C
  for (site_no in seq_along(obs_pred_results)) {
    ## Get site data
    data_list <- get_site_data(sites[site_no], monthly_data, soil_data)

    # Get observed and simulated fractions
    observed <- c(
      data_list$Cpools_final$RPM,
      data_list$Cpools_final$HUM,
      data_list$Cpools_final$IOM,
      data_list$Cpools_final$TOC
    )

    obs_pred_results[[site_no]] <- c(sites[site_no], observed, calc_site_C_turnover(data_list, params))
  }
  df <- do.call(rbind.data.frame, obs_pred_results)
  names(df) <- c(
    "site",
    "RPM_obs",
    "HUM_obs",
    "IOM_obs",
    "TOC_obs",
    "DPM_pred",
    "RPM_pred",
    "BIO_pred",
    "HUM_pred",
    "IOM_pred",
    "CO2_pred",
    "TOC_pred"
  )
  df
}


# This function is specifically designed to run within
# some functions for parameter optimization with "optim"
compute_stat_obs_pred <- function(obs_pred_df) {

  # Compute the RMSE of observed and simulated C fractions across the sites

  # Parameters
  # ----------
  # obs_pred_df: Dataframe
  #   Observed and simulated C fractions for all sites

  # Returns
  # -------
  # Scalar
  #   Numeric RMSE

  round(RMSE(obs_pred_df$obs, obs_pred_df$pred), 3)
}


# This function is specifically designed to run with
# the "optim" package for parameter optimization
get_obs_pred_stat <- function(
                              params,
                              sites,
                              monthly_data,
                              soil_data,
                              C_fractions_to_use = c("RPM", "HUM")) {

  # Get the RMSE of observed and simulated C fractions for all the sites
  # Also allows for the specification of what C fractions should be used in
  # Optimisation

  # Parameters
  # ----------
  # params: Vector
  #   Numeric values of RothC parameters (e.g. the C fractions decomposition rates)

  # sites: Vector
  #   Site ids

  # monthly_data: Dataframe
  #   Monthly timeseries of vegation, climate and soil moisture for study sites

  # soil_data: Dataframe
  #   Sampled soil information for the study sites

  # C_fractions_to_use: Scalar Vector
  #   C fractions to used for RMSE computation

  # Returns
  # -------
  # Scalar
  #   Returns RMSE

  # Get observed and predicted C fractions across all sites
  obs_pred_df <- get_sites_obs_pred_C(
    params,
    sites,
    monthly_data,
    soil_data
  )

  # Make a two column df (one for observed and the other
  # for predicted values) for the chosen fractions
  selected_fractions_df <- obs_pred_df %>%
    dplyr::select(contains(C_fractions_to_use)) %>%
    pivot_longer(everything(),
      names_to = c("fraction", ".value"),
      names_pattern = "(.*)_(.*)"
    ) %>%
    dplyr::select(!"fraction") %>%
    mutate_if(is.character, as.numeric)

  # Compute RMSE
  stat <- compute_stat_obs_pred(selected_fractions_df)
  stat
}


# TODO: May delete this function
run_rothC_on_parameter_grid <- function(site,
                                        param_grid,
                                        monthly_data,
                                        soil_data) {
  # Given a grid of parameters, runs RothC on the grid

  # Parameters
  # ----------
  # site: String
  #   Site id

  # param_grid: Dataframe
  #   Parameters of RothC and that the C model

  # soil_data: Dataframe
  #   sampled soil information for the study sites

  # Returns
  # -------
  # List
  #   Returns parameter values, predicted and observed RPM and HUMIC fractions,
  #   and Total C

  # 1. Get required data sets for site
  data_list <- get_site_data(site, monthly_data, soil_data)

  # 2. Run RothC for site

  # Calculate temperature modifying rate
  fT <- calc_fC_RothC(data_list$temperature)

  # Calculate moisture modifying rate using soil moisture results from water balance model
  # We made this modification to RothC
  fM <- calc_fM_WB(data_list$soil_moisture, data_list$bucket_size, data_list$NDVI_cover$cover)$b

  # Below is the default RothC routine for calculating moisture modifying rate. I left this here as a
  # reminder but may integrate it into the code by adding a condition to switch between the two
  # fM  <- calc_fM_RothC(rain, evap, pan_ET = TRUE, pclay, soil_thick = 30, bare = FALSE)$b

  # Calculate cover modifying rate
  fC <- calc_fC_RothC(data_list$NDVI_cover$cover)

  # Calculate abc
  abc <- fT * fM * fC

  # Iterate through the parameter grid and use the different
  # combinations of the parameters to compute Cinput and then
  # run RothC

  # - Initialise the vector to save the results of each parameter combination
  results <- vector("list", nrow(param_grid))
  for (l in seq_along(results)) {

    #--- Compute C inputs ---
    # Get coefficient, b, that integrates harvest index, root-shoot ratio, and transpiration rate
    b <- param_grid$b[l]

    # Compute plant residue input
    plant_residue <- calc_Cinput(data_list$NDVI_cover$NDVI, data_list$ET, b)

    # Save row values
    row_values <- c(site, b, l)

    # FYM
    fym <- rep(0, length(plant_residue))

    # Put plant residue and fmy in a table
    C_in <- tibble(plant_residue, fym)

    #--- Run RothC ---
    model <- calc_RothC_model(
      data_list$time,
      initial_C_pools = list(
        DPM = 0,
        RPM = data_list$Cpools_initial$RPM,
        BIO = 0,
        HUM = data_list$Cpools_initial$HUM,
        IOM = data_list$Cpools_initial$IOM,
        CO2 = 0
      ),
      ks = list(
        DPM = param_grid$k_DPM[l],
        RPM = param_grid$k_RPM[l],
        BIO = param_grid$k_BIO[l],
        HUM = param_grid$k_HUM[l],
        IOM = 0,
        CO2 = Inf
      ),
      C_in,
      data_list$pclay,
      abc,
      DPM_RPM_ratio = param_grid$DPM_RPM_ratio[l]
    )

    # Save row values to list
    results[[l]] <- append(row_values, c(
      param_grid$k_DPM[l],
      param_grid$k_RPM[l],
      param_grid$k_BIO[l],
      param_grid$k_HUM[l],
      param_grid$DPM_RPM_ratio[l],
      data_list$Cpools_final$RPM,
      tail(model, n = 1)$RPM,
      data_list$Cpools_final$HUM,
      tail(model, n = 1)$HUM,
      data_list$Cpools_final$TOC,
      tail(model, n = 1)$Total_C
    ))
  }

  df <- do.call(rbind.data.frame, results)

  names(df) <- c(
    "site",
    "b",
    "model",
    "k_dpm",
    "k_rpm",
    "k_bio",
    "k_hum",
    "DPM_RPM_ratio",
    "RPM_observed",
    "RPM_predicted",
    "HUM_observed",
    "HUM_predicted",
    "TOC_obsereved",
    "TOC_predicted"
  )
  df
}
