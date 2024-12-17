#
# PBT Data Preparation Function Library
#
# A collection of functions to support statistical analysis of Dee's pygmy
# bluetongue fieldwork data 2020/2021.
# Written by David Mansueto, Aug 2021 - 2024.

# 2024 Update: Now works with behaviour scoring undertaken in Solomon Coder

library(dplyr) # used to concatenate data-frames

# A note on type prefixes
# Detailed backstory: David uses type prefixes to reduce the confusion
# inherent in "loosely typed" languages.
# Loosely typed languages are those where variables are not tied to a
# specific data type.
# In 'strictly typed' languages, such as c, variables must have a defined
# data type, for example:
# > float my_var = 3.14;
# > my_var = "some string";
# will give an error: a float cannot hold a string.
# However, in 'loosely typed' languages, variables can change type
# seemingly at whim, for example:
# > my_var <- 3.14
# > my_var <- "some string"
# R says "No worries!" The confusion really sets in when functions return
# a datatype you're not expecting.
#
#| # | Datatype                                                           |
#|---|--------------------------------------------------------------------|
#| n | Numeric (integer, floating point, or complex?)                     |
#| s | String ("character" in R vernacular)                               |
#| l | Logical (boolean)                                                  |
#|---|--------------------------------------------------------------------|
#| i | Integer                                                            |
#| j | iterator (an integer that increments inside of loops)              |
#| f | Factor                                                             |
#|---|--------------------------------------------------------------------|
#| v | Vector (1-dim, single datatype)                                    |
#| m | Matrix (2-dim, single datatype)                                    |
#| a | Array (n-dim, single datatype)                                     |
#| l | List (1-dim, multiple datatypes)                                   |
#| d | Data.frame (2-dim, multiple datatypes)                             |
#| t | dateTime (POSIXct or POSIXlt; there is no isolated time datatype)  |
#| td| TimeDiff (a difference between two datetimes)                      |




#' Import and Prepare Dee's PBT Honours Data.
#'
#' \code{d_parse_dee_pbt_data} is a utility that imports and prepares for
#' statistical analysis, raw experimental data from Dee's 2020/2021 Honours
#' project field work on \emph{Tiliqua Adelaidensis}—Pygmy Blue-Tongue
#' (PBT) skinks.
#'
#' @section{Import}
#' 1. Loads behavioural and thermal data from .xlsx
#'
#' 1.1 Behavioural Data
#'
#' @param s_behavioural_workbook filename of xlsx with behavioural data.
#' Must include absolute path—if not in current working directory.
#'
#' @param s_behavioural_sheet name of 'tab' in workbook.
#' Default [\code{""}]: use first sheet.
#'
#' @param i_behavioural_row_skip number of rows to skip in behavioural xlsx
#' (e.g. combined headings, etc).
#' Default [{\code0}]: skip zero rows.
#'
#' @param s_behavioural_na_match e.g. \code{"#DIV/0!"}
#'
#' 1.2 Thermal Data
#'
#' @param s_thermal_workbook filename of xlsx with thermal data. Must
#' include absolute path—if not in current working directory.
#'
#' @param s_thermal_sheet name of 'tab' in workbook.
#' Default [\code{""}]: use first sheet.
#'
#' @param i_thermal_row_skip number of rows to skip in thermal xlsx (e.g.
#' combined headings, etc).
#' Default [{\code0}]: skip zero rows.
#'
#' @param s_behavioural_na_match e.g. \code{"#DIV/0!"}
#'
#'
#' @section{Summary}
#' 2. Does some clean up and sanity checking, then performs preparations
#' & produces output data. Used for internal processes.
#'
#' 2.1 Summary Data
#' This is a summary of conditions, counts, etc. for each unique
#' date/site/burrow in the dataset.
#'
#' @param s_summary_data_filename name of the file to write to, e.g.
#' \code{"2020-2022_NoFeb-Summary.rda"}.
#' Default [\code{""}]: do not save to disk.
#'
#'
#' @section{Output}
#' 3. Preparation of output data
#'
#' This is the primary purpose of this function; a combination of the input
#' behavioural and thermal data, ready for statistical analysis.
#'
#' 3.1 Temporal-Regulation (& "Binary" activity)
#' This stage ensures each (statistically speaking) 'obersvation' has the
#' same temporal weighting, in code terms by ensuring each row has the same
#' time duration.
#'
#' @param i_target_row_mins this is the number of minutes that each row
#' in the resulting dataframe covers.
#' Set to fewer minutes for more granular data, or more minutes to have
#' greater simplification.
#' If setting more than thirty minutes per row, consider setting
#' \code{i_activity_test = "majority"}
#'
#' @param s_activity_test tells the utility how to decide if an obersvation
#' is 'on' in a given row. For example, if a PBT was basking for 18 minutes
#' out of a 30 minute row, \code{"any"} and \code{"majority"} would return
#' \code{TRUE}, while \code{"entirity"} would return \code{FALSE}.
#' Default [\code{"any"}]: returns \code{TRUE} if there is any activity
#' within the row.
#'
#' DEPRECATED @param l_binarize if \code{TRUE}, activity observations
#' such as "basking" or "rotate" will be reduced to a
#' single "active" column: \code{TRUE} if the lizard is out
#' of the burrow and false otherwise.
#' Default [\code{FALSE}]: do not binarize.
#' DEPRECATED: parameter removed; always created.
#'
#' 3.2 Temporal-Trimming
#' This stage 'crops' the observations to an evenly-weighted total time of
#' overservations for each date/site/burrow.
#'
#' @param i_target_total_mins the total amount of time for a burrow.
#' Set to \code{0} or \code{NA} to disable trimming.
#' Default [\code{6 * 60}]: 6 hours.
#'
#' @param s_target_median_time the time to try and centralise the trimmed
#' data about.
#' Note, in cases of uneven distribution the code will keep as much of the
#' 'underrepresented' side as possible.
#' In cases where the median is not present, footage is trimmed from the
#' end. For example: a median of \code{"12:00"} but footage runs
#' from \code{"14:33"} to \code{"19:23"}, the cropped footage would also
#' start at \code{"14:33"}.
#' Default [\code{12:00}]: noon.
#'
#' 3.3 Optional, save to disk.
#'
#' @param s_pbt_data_filename the name of the file to save the resulting
#' dataframe too, e.g. \code{"2020-2022_NoFeb.rda"}.
#'
#'
#' @return a dataframe containing both behavioural and thermal data,
#' rebased to even time step per row, with a binary activity column
#' \code{vl_active}, optionally trimmed to a specified total duration.
#'
#' Miscellaneous
#'
#' @param l_verbose_output whether to do \emph{lots} of \code{View}ing, of
#' non-necessary (but perhaps still interesting) things.
#' Default [\code{FALSE}]: only print necessary messages.
#'
#' @usage
#' First, ensure source code is available:
#' \code{source("dee_pbt_utils.r")}
#' Then this function can be called via:
#' \code{d_pbt <- d_parse_dee_pbt_data(...)}

d_parse_behave_vidObs_data <- function(
  d_vidObs_raw,
  d_enviro_raw,
  vs_months_to_omit = "",
  i_target_row_mins = 30,
  s_activity_test = "any",
  s_target_median_time = "12:00",
  i_target_total_mins = 6 * 60,
  s_session_name = "",
  s_behave_data_filename = "",
  s_summary_data_filename = "",
  l_verbose_output = FALSE
) {

  if (s_activity_test != "any") {
    stop(
      'Sorry, so far only s_activity_test = "any" has been implemented.'
      )
  }

  if ((i_target_total_mins %% 2) != 0) {
    stop(
      '"i_target_total_mins" must be an even integer.'
      )
  }


  #* ----- PREAMBLE -------------------------------------------------------

  s_function_version <- "0.2: 25 Mar 2022"

  s_logdir <- "Dee PBT Utils Logs"
  if (!file.exists(s_logdir)) dir.create(s_logdir)
  s_logfile <- paste(
    s_logdir,
    "/Log_",
    if (length(s_session_name) > 0) {
      paste(
        s_session_name,
        "_",
        sep = ""
      )
    },
    strftime(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
    ".txt",
    sep = ""
  )

  # Convenient (succinct) version of print_message
  print_msg <- function(s_msg,
    l_urgent = FALSE,
    l_verbose = l_verbose_output,
    s_file = s_logfile
  ) {
    print_message(s_msg, s_file, l_verbose, l_urgent)
  }

  print_session_details(
    sprintf(
        "dee_pbt_utils::d_parse_behave_vidObs_data
      Version %s
      Logfile for session %s, %s on %s",
        s_function_version,
        s_session_name,
        s_from_t(Sys.time(), "%H:%M"),
        s_from_t(Sys.time(), "%a %d %b %Y")
    ),
    s_function_version,
    s_logfile,
    data.frame(
      "vs_months_to_omit" = vs_months_to_omit,
      "i_target_row_mins" = i_target_row_mins,
      "s_activity_test" = s_activity_test,
      "s_target_median_time" = s_target_median_time,
      "i_target_total_mins" = i_target_total_mins,
      "s_session_name" = s_session_name,
      "s_behave_data_filename" = s_behave_data_filename,
      "s_summary_data_filename" = s_summary_data_filename,
      "l_verbose_output" = l_verbose_output
    ),
    print_msg
  )


  #* ----- PROCESSING START -----------------------------------------------

  print_msg("Section 1. Tidy")

  d_vidObs_tidy <- d_tidy_behave_vidObs_data(d_vidObs_raw, vs_months_to_omit)

  d_enviro_tidy <- d_tidy_enviro_data(d_enviro_raw, vs_months_to_omit)


  print_msg("Section 2. Transform")

  d_pbt <- d_combine_pbt_data(
    d_vidObs_tidy,
    d_enviro_tidy
  )

  d_summary <- d_summarise_pbt_data(
    d_pbt,
    d_enviro_tidy,
    print_msg
  )

  d_pbt <- d_regulate_pbt_data(
    d_pbt = d_vidObs_tidy,
    d_summary,
    i_target_row_mins,
    i_target_total_mins,
    print_msg
  )

  d_pbt <- d_binarize_pbt_data(
    d_regulated = d_pbt,
    d_behave = d_vidObs_tidy,
    print_msg
  )

  d_pbt <- d_combine_pbt_data(
    d_behave = d_pbt,
    d_enviro_tidy
  )

  d_pbt <- d_trim_pbt_data(
    d_pbt,
    d_summary,
    i_target_total_mins,
    i_target_row_mins,
    s_target_median_time,
    print_msg
  )

  d_summary <- d_trim_then_summarise_pbt_data(
    d_behave = d_vidObs_tidy,
    d_thermo = d_enviro_tidy,
    d_summary,
    i_target_total_mins,
    s_target_median_time,
    print_msg
  )

  print_msg("Section 4. Save to Disk")

  save_data_to_file(d_pbt, s_behave_data_filename, print_msg)

  save_data_to_file(d_summary, s_summary_data_filename, print_msg)


  return(d_pbt)
}


###
d_parse_approach_data <- function(
    d_apprch_raw,
    d_enviro_raw,
    vs_months_to_omit = "",
    s_session_name = "",
    s_approach_data_filename = ""
) {
  
  
  #* ----- PREAMBLE -------------------------------------------------------
  
  s_function_version <- "0.2: 30 Aug 2021"
  
  s_logdir <- "Dee PBT Utils Logs"
  if (!file.exists(s_logdir)) dir.create(s_logdir)
  s_logfile <- paste(
    s_logdir,
    "/Log_",
    if (length(s_session_name) > 0) {
      paste(
        s_session_name,
        "_",
        sep = ""
      )
    },
    strftime(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
    ".txt",
    sep = ""
  )
  
  # Convenient (succinct) version of print_message
  print_msg <- function(s_msg,
                        l_urgent = FALSE,
                        l_verbose = FALSE,
                        s_file = s_logfile
  ) {
    print_message(s_msg, s_file, l_verbose, l_urgent)
  }
  
  print_session_details(
    sprintf(
      "dee_pbt_utils::d_parse_approach_data
      Version %s
      Logfile for session %s, %s on %s",
      s_function_version,
      s_session_name,
      s_from_t(Sys.time(), "%H:%M:%S"),
      s_from_t(Sys.time(), "%a %d %b %Y")
    ),
    s_function_version,
    s_logfile,
    data.frame(
      "s_session_name" = s_session_name,
      "s_approach_data_filename" = s_approach_data_filename
    ),
    print_msg
  )
  
  
  #* ----- PROCESSING START -----------------------------------------------
  
  print_msg("Section 2. Tidy")
  
  d_apprch_tidy <- d_tidy_apprch_data(d_apprch_raw, vs_months_to_omit)
  
  d_enviro_tidy <- d_tidy_enviro_data(d_enviro_raw, vs_months_to_omit)
  
  
  print_msg("Section 3. Transform")
  
  d_pbt <- d_combine_pbt_data(
    d_apprch_tidy,
    d_enviro_tidy
  )
  
  
  print_msg("Section 4. Save to Disk")
  
  save_data_to_file(d_pbt, s_approach_data_filename, print_msg)
  
  
  return(d_pbt)
}


# * SECTION FUNCTIONS =====================================================

# No longer used
d_import_data_from_xlsx <- function(
  s_workbook,
  s_sheet,
  i_row_skip,
  s_na_match,
  print_msg
) {

  d_data <- readxl::read_excel(
    path = s_workbook,
    sheet = s_sheet,
    skip = i_row_skip,
    na = s_na_match
  )
  print_msg(sprintf(
    "Loaded data from %s; %i observations of %i variables.",
    s_workbook,
    nrow(d_data),
    ncol(d_data)
  ))
  return(d_data)
}

#*=========================================================================

d_tidy_enviro_data <- function(d_enviro, vs_months_to_omit=NULL) {

  # Rename time and date columns
  colnames(d_enviro)[colnames(d_enviro) == "Date Text"] <- "vs_dates"
  colnames(d_enviro)[colnames(d_enviro) == "Time Text"] <- "vs_times"
  colnames(d_enviro)[colnames(d_enviro) == "Site"]      <- "vs_sites"

  # Source data time is a bit wonky; sometimes out by ± 1 min, sometimes 15
  # Use floor to ensure we still have consistent time steps but all times
  # are precisely thirty minutes.
  # Note: #HARDCODE for 30-minute incremental data.
  vt_times <- t_floor(
    t_from_s(d_enviro$vs_times, "%H:%M"),
    as.difftime(30, units = "mins", tz = "UTC")
  )
  # Update in dataframe
  d_enviro$vs_times <- s_from_t(vt_times, "%H:%M")

  # Create compound ID column
  d_enviro$vs_date_time30_site <- paste(
    d_enviro$vs_dates,
    d_enviro$vs_times,
    d_enviro$vs_sites,
    sep = "_"
  )

  # Remove any data in months slated for omission
  if(!is.null(vs_months_to_omit)){
    d_enviro <- d_omitting_months(d_enviro, vs_months_to_omit)
  }
  
  return(d_enviro)
}

#*-------------------------------------------------------------------------

d_tidy_vidObs_data <- function(d_vidObs, vs_months_to_omit) {

  # Tidy up video-observed behavioural data

  # Set up columns for internal function covenience

  # Date string
  colnames(d_vidObs)[colnames(d_vidObs) == "Date"] <- "vs_dates"
  d_vidObs$vs_dates <- s_from_t(d_vidObs$vs_dates, "%Y-%m-%d")
  # First stab at generic time column; this one will change a bit
  d_vidObs$vs_times <- d_vidObs$"Actual Start Time"
  d_vidObs$vt_times <- t_from_s(d_vidObs$vs_times, "%H:%M:%S")
  # Other time columns
  d_vidObs$vt_vidObss <- t_from_s(d_vidObs$"Actual Start Time", "%H:%M:%S")
  d_vidObs$vt_activities <- t_from_s(d_vidObs$"Actual Time", "%H:%M:%S")
  # Site is complex in Behave data:
  #  1. it encodes "enclosure" at Kelly Hill, so we need to strip out the
  #     enclosure number (i.e. "5" or "7")
  #  2. Tiliqua is any of "Til", "TilW", or "Til(W)", etc.), so we make it
  #     exactly "Til"
  #  3. We do TC, too, just 'cos.
  colnames(d_vidObs)[colnames(d_vidObs) == "Site"] <- "vs_sites"
  vs_sites <- d_vidObs$vs_sites # for code clarity
  vs_sites[grep("k", vs_sites, ignore.case = TRUE)] <- "KH"
  vs_sites[grep("til", vs_sites, ignore.case = TRUE)] <- "Til"
  vs_sites[grep("tc", vs_sites, ignore.case = TRUE)] <- "TC"
  d_vidObs$vs_sites <- vs_sites
  colnames(d_vidObs)[colnames(d_vidObs) == "Burrow ID"] <- "vs_burrows"
  # 2021-08-17: Adding column for site and lineage
  # Need to rebuild enclosure from cleaned site
  vs_first_char_of_burrows <- substr(
    d_vidObs$vs_burrows,
    start = 1,
    stop = 1
  )
  # Delete chars that are 't' = Til or TC
  vs_first_char_of_burrows[grep("t", vs_sites, ignore.case = TRUE)] <- ""
  d_vidObs$vs_enclosures <- paste(
    vs_sites,
    vs_first_char_of_burrows,
    sep = ""
  )
  colnames(d_vidObs)[colnames(d_vidObs) == "Lineage"] <- "vs_lineages"
  # Now we can combine enclosure and lineage as "Treatment Group"
  d_vidObs$vs_treatment_groups <- paste(
    d_vidObs$vs_enclosures, # "KH5", "KH7", "Til" or "TC"
    d_vidObs$vs_lineages,
    sep = "-"
  )

  # Make compound columns

  # Thermal alignment column
  d_vidObs$vs_date_time30_site <- paste(
    d_vidObs$vs_dates,
    s_from_t(
      t_floor(
        t_from_s(d_vidObs$vs_times, "%H:%M:%S"),
        as.difftime(30, units = "mins", tz = "UTC")
      ),
      "%H:%M"
    ),
    d_vidObs$vs_sites,
    sep = "_"
  )

  # Per Burrow column
  d_vidObs$vs_date_site_burrow <- paste(
    d_vidObs$vs_dates,
    d_vidObs$vs_sites,
    d_vidObs$vs_burrows,
    sep = "_"
  )

  # Remove any data in months slated for omission
  d_vidObs <- d_omitting_months(d_vidObs, vs_months_to_omit)

  return(d_vidObs)
}

#*-------------------------------------------------------------------------

d_tidy_apprch_data <- function(d_apprch, vs_months_to_omit) {

  # Parse columns
  vs_dates          <- v_get_col(d_apprch, "Date")
  vs_times          <- v_get_col(d_apprch, "Time")
  vs_sites          <- v_get_col(d_apprch, "Site")
  vs_burrows        <- v_get_col(d_apprch, "Burrow ID")
  vs_enclosures     <- v_get_col(d_apprch, "Enclosure")
  vs_weather_codes  <- v_get_col(d_apprch, "Weather")
  vn_retreat_metres <- v_get_col(d_apprch, "Distance (m)")

  # Tidy up site and enclosure coding - can have different versions, we
  # need to have exactly the same coding everywhere

  # Kelly Hill is easy: always starts with "K"
  vi_kelly_hill <- grep("k",  vs_sites, ignore.case = TRUE)
  # Burra is always "Tiliqua" or "Til", so "ti" works
  vi_tiliqua    <- grep("ti", vs_sites, ignore.case = TRUE)
  # Twin Creek can be "TC" or "Twin Creek", so need some logic
  vi_twin_creek <- if (sum(grep("tw", vs_sites, ignore.case = TRUE))) {
    grep("tw", vs_sites, ignore.case = TRUE)
  } else {
    grep("tc", vs_sites, ignore.case = TRUE)
  }

  if ((length(vi_kelly_hill) + length(vi_tiliqua) + length(vi_twin_creek) !=
    nrow(d_apprch))) {
    stop(paste("Unexpected site encoding in approach data! Must be:\n",
      '  Kelly Hill: "KH", "KH5" or "KH7"\n',
      '  Twin Creek: "TC" or "Twin Creek"\n',
      '     Tiliqua: "Til" or "Tiliqua"\n',
      "Note: upper or lower case does not matter."))
  }

  # For enclosure, just retain '5' or '7' for now
  vs_enclosures[vi_kelly_hill] <- s_get_last_char(
    vs_enclosures[vi_kelly_hill]
  )
  vs_enclosures[vi_tiliqua]    <- ""
  vs_enclosures[vi_twin_creek] <- ""

  # Tidy up sites
  vs_sites[vi_kelly_hill] <- "KH"
  vs_sites[vi_twin_creek] <- "TC"
  vs_sites[vi_tiliqua]    <- "Til"

  # Use clean site as basis for enclosure (plus enc number from before)
  vs_enclosures <- paste(vs_sites, vs_enclosures, sep = "")

  # Lineage is HARD-CODE from enclosure
  vs_lineages <- rep(NA, nrow(d_apprch))
  vs_lineages[vs_sites == "Til" | vs_enclosures == "KH5"] <- "Tiliqua"
  vs_lineages[vs_sites == "TC"  | vs_enclosures == "KH7"] <- "Twin Creek"

  # Treatment Group = Enclosure + Lineage
  vs_treatment_groups <- paste(
    vs_enclosures,
    vs_lineages,
    sep = "_"
  )

  vl_actives <- rep(FALSE, nrow(d_apprch))
  vl_actives[!is.na(vn_retreat_metres)] <- TRUE

  vl_translocated <- rep(FALSE, nrow(d_apprch))
  vl_translocated[vs_sites == "KH"] <- TRUE

  vt_raw <- t_from_s(
    paste(vs_dates, vs_times, sep = "_"),
    "%Y-%m-%d_%H:%M"
  )

  vs_date_time30_site <- paste(
    s_from_t(vt_raw, "%Y-%m-%d"),
    s_from_t(t_floor(vt_raw, as.difftime(30, units = "mins")), "%H:%M"),
    vs_sites,
    sep = "_"
  )

  d_apprch <- data.frame(
    vt_raw,
    vs_dates,
    vs_times,
    vs_sites,
    vs_burrows,
    vs_enclosures,
    vs_lineages,
    vs_treatment_groups,
    vs_weather_codes,
    vn_retreat_metres,
    vl_actives,
    vl_translocated,
    vs_date_time30_site
  )

  d_apprch <- d_omitting_months(d_apprch, vs_months_to_omit)

  return(d_apprch)
}

#*-------------------------------------------------------------------------

d_omitting_months <- function(d_data, vs_months_to_omit) {
  vs_months_to_omit <- tolower(vs_months_to_omit)
  vt_months <- t_from_s(d_data$vs_dates, "%Y-%m-%d")
  vl_omit <- vector(mode = "logical", length = nrow(d_data))
  for (s_month in vs_months_to_omit) {
    vl_omit <- vl_omit |
      s_month == tolower(s_from_t(vt_months, "%b")) |
      s_month == tolower(s_from_t(vt_months, "%B"))
  }
  return(d_data[!vl_omit, ])
}

#' Gets a column or stops trying.
v_get_col <- function(d_data, s_col_name) {
  if (is.na(match(s_col_name, names(d_data)))) {
    stop(sprintf(
      paste('Dataframe does not contain required column "%s";',
        'does include: "%s"'
      ),
      s_col_name,
      paste(names(d_data), collapse = '", "')
      ))
  }
  return(eval(substitute(d_data$s, list(s = s_col_name))))
}

#' Gets the last character from a string of characters
s_get_last_char <- function(s) {
  substr(s, start = nchar(s), stop = nchar(s))
}

#*=========================================================================

#' Adds thermal data to behavioural, expanding time base of thermal data
#' to i_target_row_mins such that a match should probably occur...
d_combine_pbt_data <- function(d_behave, d_enviro) {

  # Remove redundant columns from thermal data (these are present in
  # behavioural anyway, having both when merge() leads to weirdness...)
  vs_thermo_headings <- names(d_enviro)
  vs_headings_to_drop <- c("vs_dates", "vs_times")
  d_enviro <- d_enviro[, vs_thermo_headings != vs_headings_to_drop]

  # Combine behavioural and thermal data. This is done by "going through"
  # each row in behavioural data and finding the corresponding thermal row,
  # the columns of which are stuck to the end of the behavioural row.
  # Note that unused thermal rows are discarded.
  d_combined <- merge(
    d_behave,
    d_enviro,
    by = "vs_date_time30_site",
    all.x = TRUE,
    sort = FALSE
  )

  return(d_combined)
}


#*=========================================================================

#' Summarise Dee's pygmy blue-tongue data at date/site/burrow level
#'
#' @return a dataframe with summary data, one row for each date/site/burrow
d_summarise_pbt_data <- function(
  d_Behave,
  d_enviro,
  print_msg
) {

  # Split df into a list of df by date/site/burrow
  ld_Behave_by_burrows <- split(d_Behave, d_Behave$vs_date_site_burrow)
  # Sort each burrow by time
  ld_Behave_by_burrows <- lapply(
    ld_Behave_by_burrows,
    function(d_burrow) d_burrow[order(d_burrow$vs_times), ]
  )

  # Create dataframe, preallocate with addressable ID column
  d_summary <- data.frame(
    "vs_date_site_burrow" = unique(d_Behave$vs_date_site_burrow)
  )

  # Loop through burrows, synthesising summary data
  for (j_burrow in seq_len(nrow(d_summary))) {

    # Note that the row of a given burrow, say, "2021-03-17_KH_5U5" in the
    # imported data is by no means necessarily located in the same row number
    # in this new summary dataframe. Thus, we index by that which does always
    # match: our ID string, precisely "2021-03-17_KH_5U5" in this example.
    # Get the burrow ID
    s_burrow <- d_summary$vs_date_site_burrow[[j_burrow]]
    # Fetch the actual data for this burrow
    d_burrow <- ld_Behave_by_burrows[[s_burrow]]
    # # Apparently time isn't always increasing, so be sure that it is
    # d_burrow <- d_burrow[order(d_burrow$vs_times), ]

    # Make a copy of the data.frame with only one row per Behave
    # Note: "Actual Start Time" is a corrected-to-local-time version of
    # "Behave Start Time" as camera time may have systematic error.
    d_burrow_unique_vids <- d_trimmed_to_unique_rows(
      d_burrow,
      s_col = "vs_times"
    )

    # Parse Behave Start Times & Calculate Durations
    lt_Behave_starts <- t_from_s(
      d_burrow_unique_vids$vs_times,
      "%H:%M:%S"
    )

#! Note: using "u" prefix here as I don't know what these containers are!
#! They're various containers of multiple values, but they're not a list,
#! vector, array, etc. Quite mystifying.
#! > u_td_Behaves
#! Time differences in mins
#!  [1] 2  3 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25
#! > u_td_Behaves[1]
#! Time difference of 2 mins
#! > u_td_Behaves[[1]]
#! [1] 2
#! > class(u_td_Behaves)
#! [1] "difftime"
#! > typeof(u_td_Behaves)
#! [1] "double"
#! > str(u_td_Behaves)
#!   difftime' num [1:23] 2 3 25 25 ...
#!   - attr(*, "units")= chr "mins"

    u_td_Behaves <- diff(lt_Behave_starts, lag = 1)
    vn_Behave_mins <- n_datetime_as_mins(u_td_Behaves)

    # Calculate Total Duration of Behave Data for this Burrow
    # Estimate how long the last Behave was (impossible to include in diff())
    # R's mode() returns class() not "most-frequently-occuring-item-in-list"...
    # Returns a table of numbers - could be more than one, so just take first.
    n_mode_Behave_mins <- tab_statistical_mode(vn_Behave_mins)[[1]]
    n_last_Behave_mins <- n_mode_Behave_mins
    # n_total_Behave_mins <- sum(as.numeric(vn_Behave_mins)) + n_last_Behave_mins
    t_Behave_start <- lt_Behave_starts[1]
    t_Behave_end <- lt_Behave_starts[length(lt_Behave_starts)] +
      as.difftime(n_last_Behave_mins, units = "mins", tz = "UTC")
    n_total_Behave_mins <- as.numeric(t_Behave_end - t_Behave_start, "mins")

    # Doing some sanity checking here
    # Looking for impossible time things in the source data
    if (any(as.logical(lapply(
      vn_Behave_mins,
      function(n_mins) n_mins < 0
    )))) {
      print_msg(sprintf(
        "Warning: negative Behave mins in %s; durations are: %s",
        s_burrow,
        paste(vn_Behave_mins, collapse = " ")
      ))
    }

    if (any(as.logical(lapply(
      vn_Behave_mins,
      function(n_mins) n_mins > (n_mode_Behave_mins * 1.1)
    )))) {
      print_msg(sprintf(
        paste(
          "Warning: Behave longer than modal duration (%f mins) in %s;",
          "durations are: %s"
        ),
        n_mode_Behave_mins,
        s_burrow,
        paste(vn_Behave_mins, collapse = " ")
      ))
    }

    
    # Calculate Number of Observations
    # Check if this is a Solomon behaviour data-frame, and adapt accordingly
    if("Non.basking.beh"%in%names(d_burrow)){
      # Solomon data has "Basking" column; any text in here means PBT is visible,
      # so we count any non-blank rows as 'active'
      n_activity_count<-sum(d_burrow$Basking!="")
    }else{
      # Tease out just columns with observations
      ls_activity_cols <- c(
        "Basking",
        "Rotate",
        "Emergence",
        "Retreat",
        "Dispersal",
        "Novel Lizard",
        "NL Retreat",
        "NL Disperse",
        "NL Return",
        "Mating",
        "Feeding"
      )
      # Make a temporary dataframe with just these columns
      d_activity <- d_burrow[, ls_activity_cols]
      # Add up rows that have any observations in them
      n_activity_count <- sum(sum(d_activity))
    }

    # Prepare dates for summary environmental data
    s_today <- d_burrow$vs_dates[1]
    s_yesterday <- s_from_t(
      t_from_s(d_burrow$vs_dates[1], "%Y-%m-%d") -
      as.difftime("24", format = "%H", tz = "UTC"),
      "%Y-%m-%d"
    )
    s_tomorrow <- s_from_t(
      t_from_s(d_burrow$vs_dates[1], "%Y-%m-%d") +
      as.difftime("24", format = "%H", tz = "UTC"),
      "%Y-%m-%d"
    )

    #* Collate data into output dataframe

    # Rename for brevity
    r <- j_burrow

    # Add meta data
    d_summary$vs_dates[r]             <- d_burrow$vs_dates[[1]]
    d_summary$vs_sites[r]             <- d_burrow$vs_sites[[1]]
    d_summary$vs_burrows[r]           <- d_burrow$vs_burrows[[1]]
    d_summary$vs_lineages[r]          <- d_burrow$vs_lineages[[1]]
    d_summary$vs_treatment_groups[r]  <- d_burrow$vs_treatment_groups[[1]]

    # Push summaries into the Summary Data data.frame
    d_summary$vn_entries[r]           <- nrow(d_burrow)
    d_summary$vn_Behaves[r]           <- nrow(d_burrow_unique_vids)
    d_summary$vn_Behave_minutes[r]    <- n_total_Behave_mins
    d_summary$vs_Behave_starts[r]     <- s_from_t(t_Behave_start, "%H:%M:%S")
    d_summary$vs_Behave_ends[r]       <- s_from_t(t_Behave_end, "%H:%M:%S")
    d_summary$vn_modal_Behave_mins[r] <- n_mode_Behave_mins

    d_summary$vn_activity_counts[r]   <- n_activity_count

    d_summary$vn_temp_amb[r] <- mean(d_burrow$"MeanAmbT")
    d_summary$vn_temp_top[r] <- mean(d_burrow$"MeanTopT")
    d_summary$vn_temp_mid[r] <- mean(d_burrow$"MeanMidT")
    d_summary$vn_temp_bot[r] <- mean(d_burrow$"MeanBaseT")

    d_summary$vn_humid_amb[r] <- mean(d_burrow$"MeanAmbH")
    d_summary$vn_humid_top[r] <- mean(d_burrow$"MeanTopH")
    d_summary$vn_humid_mid[r] <- mean(d_burrow$"MeanMidH")
    d_summary$vn_humid_bot[r] <- mean(d_burrow$"MeanBaseH")

    d_summary$"Yesterday_Max_Celsius_Top" <- max(
      d_enviro$vn_temp_top[d_enviro$vs_dates == s_yesterday])
    d_summary$"Yesterday_Min_Celsius_Top" <- min(
      d_enviro$vn_temp_top[d_enviro$vs_dates == s_yesterday])
    d_summary$"Yesterday_Mean_Celsius_Top" <- mean(
      d_enviro$vn_temp_top[d_enviro$vs_dates == s_yesterday])
    d_summary$"Yesterday_Max_RH_Amb" <- max(
      d_enviro$vn_humid_amb[d_enviro$vs_dates == s_yesterday])
    d_summary$"Yesterday_Min_RH_Amb" <- min(
      d_enviro$vn_humid_amb[d_enviro$vs_dates == s_yesterday])
    d_summary$"Yesterday_Mean_RH_Amb" <- mean(
      d_enviro$vn_humid_amb[d_enviro$vs_dates == s_yesterday])

    d_summary$"Today_Max_Celsius_Top" <- max(
      d_enviro$vn_temp_top[d_enviro$vs_dates == s_today])
    d_summary$"Today_Min_Celsius_Top" <- min(
      d_enviro$vn_temp_top[d_enviro$vs_dates == s_today])
    d_summary$"Today_Mean_Celsius_Top" <- mean(
      d_enviro$vn_temp_top[d_enviro$vs_dates == s_today])
    d_summary$"Today_Max_RH_Amb" <- max(
      d_enviro$vn_humid_amb[d_enviro$vs_dates == s_today])
    d_summary$"Today_Min_RH_Amb" <- min(
      d_enviro$vn_humid_amb[d_enviro$vs_dates == s_today])
    d_summary$"Today_Mean_RH_Amb" <- mean(
      d_enviro$vn_humid_amb[d_enviro$vs_dates == s_today])

    d_summary$"Tomorrow_Max_Celsius_Top" <- max(
      d_enviro$vn_temp_top[d_enviro$vs_dates == s_tomorrow])
    d_summary$"Tomorrow_Min_Celsius_Top" <- min(
      d_enviro$vn_temp_top[d_enviro$vs_dates == s_tomorrow])
    d_summary$"Tomorrow_Mean_Celsius_Top" <- mean(
      d_enviro$vn_temp_top[d_enviro$vs_dates == s_tomorrow])
    d_summary$"Tomorrow_Max_RH_Amb" <- max(
      d_enviro$vn_humid_amb[d_enviro$vs_dates == s_tomorrow])
    d_summary$"Tomorrow_Min_RH_Amb" <- min(
      d_enviro$vn_humid_amb[d_enviro$vs_dates == s_tomorrow])
    d_summary$"Tomorrow_Mean_RH_Amb" <- mean(
      d_enviro$vn_humid_amb[d_enviro$vs_dates == s_tomorrow])
  }

  print_msg(sprintf(
    "Summary data complete: %i columns for %i burrows.",
    ncol(d_summary), nrow(d_summary)
  ))

  return(d_summary)
}


#*=========================================================================

#' Apply Regulation to Dee's PBT Data
#' Creates a new dataframe with continuous, regular timescale (column)
d_regulate_pbt_data <- function(
  d_pbt,
  d_summary,
  i_target_row_mins,
  i_target_total_mins,
  print_msg
) {

  dt_target_row_duration <- as.difftime(i_target_row_mins,
    units = "mins", tz = "UTC")

  # Split df into a list of df by date/site/burrow
  ld_pbt_by_burrows <- split(d_pbt, f = d_pbt$vs_date_site_burrow)
  # Sort each burrow by time
  ld_pbt_by_burrows <- lapply(
    ld_pbt_by_burrows,
    function(d_burrow) d_burrow[order(d_burrow$vs_times), ]
  )

  # Here we make a new dataframe in which all rows are of the
  # same time duration, as specified by @param i_target_row_mins.
  # This makes the data evenly-weighted temporally.

  # We also add single column activity observation, thus 'binary'.
  # This is 'true' if the lizard is 'active' during that row's time,
  # based on:
  # a) if the lizard is out of the burrow (defined as from
  # 'basking' until 'retreat')
  # b) what proportion of the row's duration the lizard is
  # active, either "any", "majority" or "entirety", specified
  # in @param s_activity_test.

  # Create placeholder dataframe for time-regulated data
  d_regulated <- data.frame()

  for (s_burrow in names(ld_pbt_by_burrows)) {

    d_burrow <- ld_pbt_by_burrows[[s_burrow]]
    # Find the row in d_summary that corresponds to this burrow
    r_summary <- match(s_burrow, d_summary$vs_date_site_burrow)

    #TODO: Check col names!
    # We have to round outwards in these bounds lest we lose data
    # Hopefully trimming later removes any "invented" data
    t_new_start <- t_floor(
      t_from_s(d_summary$vs_vidObs_starts[r_summary], "%H:%M:%S"),
      dt_target_row_duration
    )
    t_new_end   <- t_ceiling(
      t_from_s(d_summary$vs_vidObs_ends[r_summary], "%H:%M:%S"),
      dt_target_row_duration
    )

    # Create a vector of start times for each row
    n_total_mins <- floor(
      as.numeric(t_new_end - t_new_start, units = "mins"))
    n_new_rows <- n_total_mins / i_target_row_mins

    vt_row_starts <- 1:n_new_rows * dt_target_row_duration + t_new_start

    # Build the new dataframe for this burrow
    d_burrow_regulated <- data.frame(
      vs_date_site_burrow = s_burrow,
      vs_dates = d_burrow$vs_dates[[1]],
      vs_times = s_from_t(vt_row_starts, "%H:%M:%S"),
      vs_sites = d_burrow$vs_sites[[1]],
      vs_burrows = d_burrow$vs_burrows[[1]],
      vs_lineages = d_burrow$vs_lineages[[1]],
      vs_enclosures = d_burrow$vs_enclosures[[1]],
      vs_treatment_groups = d_burrow$vs_treatment_groups[[1]],
      vs_date_time30_site = vapply(
        vt_row_starts,
        function(t) paste(
          d_burrow$vs_dates[[1]],
          s_from_t(
            t_floor(t, as.difftime(30, units = "mins", tz = "UTC")),
            "%H:%M"
          ),
          d_burrow$vs_sites[[1]],
          sep = "_"
        ),
        ""
      )
    )

    d_regulated <- rbind(d_regulated, d_burrow_regulated)
  }

  return(d_regulated)
}

#*=========================================================================

d_binarize_pbt_data <- function(
  d_regulated,
  d_vidObs,
  print_msg
) {

  ld_regulated_by_burrows <- split(
    d_regulated,
    d_regulated$vs_date_site_burrow
  )
  # Sort each burrow by time
  ld_regulated_by_burrows <- lapply(
    ld_regulated_by_burrows,
    function(d_burrow) d_burrow[order(d_burrow$vs_times), ]
  )

  ld_vidObs_by_burrows <- split(
    d_vidObs,
    d_vidObs$vs_date_site_burrow
  )
  # Sort each burrow by time
  ld_vidObs_by_burrows <- lapply(
    ld_vidObs_by_burrows,
    function(d_burrow) d_burrow[order(d_burrow$vs_times), ]
  )

  d_binarized <- data.frame()

  for (s_burrow in names(ld_regulated_by_burrows)) {
    d_burrow_regulated <- ld_regulated_by_burrows[[s_burrow]]
    d_burrow_behave    <- ld_vidObs_by_burrows[[s_burrow]]

    d_burrow_regulated$vl_actives <- vl_binarize_activity(
      s_burrow,
      d_burrow_regulated,
      d_burrow_behave,
      print_msg
    )

    d_binarized <- rbind(d_binarized, d_burrow_regulated)
  }

  return(d_binarized)
}

#*-------------------------------------------------------------------------

vl_binarize_activity <- function(
  s_burrow,
  d_burrow_regulated,
  d_burrow_behave,
  print_msg
) {

  # Create vector for activity state ("binary column")
  vl_active <- vector(
    mode = "logical", # default value is FALSE
    length = nrow(d_burrow_regulated)
  )
  # Time step of rows in regulated data
  dt_row_duration <- t_from_s(d_burrow_regulated$vs_times[2], "%H:%M:%S") -
                     t_from_s(d_burrow_regulated$vs_times[1], "%H:%M:%S")

  # Activity = out of burrow: [theoretically] always initiated by "Basking"
  vl_active_starts <- as.logical(d_burrow_behave$"Basking")
  vl_active_ends <- as.logical(d_burrow_behave$"Retreat") |
    as.logical(d_burrow_behave$"Return")

  if (sum(vl_active_starts) == 0 | sum(vl_active_ends) == 0) {
    print_msg(sprintf(
      "Notice: No activity flagged in %s.",
      s_burrow
    ))
    return(vl_active)
  }

  # Try and fix mismatched activity bounds.
  # As of 18 Aug 2021, protocol is:
  # If the earliest observation is an "end", assume active prior
  # If the latest observation is a "start", assume active until end
  vr_active_starts <- which(vl_active_starts) # indexes of active rows
  vr_active_ends <- which(vl_active_ends)
  # Does the first end occur before the first start?
  if (vr_active_ends[[1]] < vr_active_starts[[1]]) {
    # then assume active from beginning
    vl_active_starts[[1]] <- TRUE
    vr_active_starts <- c(1, vr_active_starts)
  }
  # Does the last start occur after the last end?
  if (vr_active_starts[[length(vr_active_starts)]] >
  vr_active_ends[[length(vr_active_ends)]]) {
    # then assume active until finish
    vl_active_ends[[length(vl_active_ends)]] <- TRUE
    vr_active_ends <- c(vr_active_ends, length(vl_active_ends))
  }

  # The above code will fix mismatched activity at either end, but not
  # in the middle. Add a warning so user knows to check the raw data.
  # Using if/else if since the min calc will cause an error if the sums do
  # not match.
  # print_msg(sprintf(
  #   "%s: %i s, %i e; (%s) s, (%s) e",
  #   s_burrow,
  #   sum(vl_active_starts),
  #   sum(vl_active_ends),
  #   paste(vr_active_starts, collapse = ", "),
  #   paste(vr_active_ends, collapse = ", ")
  # ))
  l_match_fault <- FALSE
  if (sum(vl_active_starts) != sum(vl_active_ends)) {
    l_match_fault <- TRUE
  } else if (min(vr_active_ends - vr_active_starts) < 0) {
    l_match_fault <- TRUE
  }
  if (l_match_fault) {
    print_msg(sprintf(
      paste(
        "Warning: IMPORTANT! Activity not binarized for %s due to",
        "mismatched activity bounds:\n",
        "  %i starts, rows (%s)\n",
        "  %i ends,   rows (%s)."
      ),
      s_burrow,
      sum(vl_active_starts), paste(vr_active_starts, collapse = ", "),
      sum(vl_active_ends), paste(vr_active_ends, collapse = ", ")
    ))
    return(vl_active)
  }

  # Convert bounds to times
  vt_active_starts <- d_burrow_behave$vt_activities[vl_active_starts]
  vt_active_ends   <- d_burrow_behave$vt_activities[vl_active_ends]

  # Loop through each activity from behave data
  for (j_activity in seq_len(length(vt_active_starts))) {
    t_active_start <- vt_active_starts[[j_activity]]
    t_active_end <- vt_active_ends[[j_activity]]
    l_active <- FALSE # Flag whether lizard is already active
    # Loop through each row in regulated data
    for (j_row in seq_len(nrow(d_burrow_regulated))) {
      t_row_start <- t_from_s(
        d_burrow_regulated$vs_times[j_row],
        "%H:%M:%S"
      )
      t_row_end <- t_row_start + dt_row_duration
      # If not active, check if this is the start of activity
      if (
        !l_active &
        t_active_start >= t_row_start &
        t_active_start <= t_row_end) {
          # This is the start of an activity period, so
          # raise the flag
          l_active <- TRUE
      }
      # If active in this row
      if (l_active) {
        # Set boolean activity value
        vl_active[[j_row]] <- TRUE
        # Check if this activity period ends this row
        if (t_active_end < t_row_end) {
          # It does end this row, so we're done for this activity
          break # End for loop here
        }
      }
    }
  }

  return(vl_active)
}

#*-------------------------------------------------------------------------

#' Reattaches behavioural activity to a temporally-regulated dataset
d_reattach_activity <- function(
  d_regulated,
  d_vidObs,
  print_msg
) {

  # Make new ID column that is specific to the minute
  d_vidObs$vs_date_time1_site_burrow <-
    vs_get_date_time1_site_burrow(d_vidObs)

  d_regulated$vs_date_time1_site_burrow <-
    vs_get_date_time1_site_burrow(d_regulated)

  vs_activity_cols <- c(
    "Basking",
    "Rotate",
    "Emergence",
    "Retreat",
    "Dispersal",
    "Novel Lizard",
    "NL Retreat",
    "NL Disperse",
    "NL Return",
    "Mating",
    "Feeding"
  )

  vs_cols_to_keep <- c(
    vs_activity_cols,
    "vs_date_time1_site_burrow"
  )

  d_activity <- d_vidObs[, vs_cols_to_keep]

  d_reattached <- merge(
    d_regulated,
    d_activity,
    by = "vs_date_time1_site_burrow",
    all.x = TRUE,
    sort = FALSE
  )

  d_reattached[is.na(d_reattached)] <- 0

  return(d_reattached)
}

#*-----------------

vs_get_date_time1_site_burrow <- function(d_data) {
    d_data$vs_date_time1_site_burrow <- paste(
    d_data$vs_dates,
    s_from_t(
      if (any(names(d_data) == "vt_activities")) {
        t_round(
          d_data$vt_activities,
          as.difftime(1, units = "mins", tz = "UTC")
        )
      } else {
        t_from_s(d_data$vs_times, "%H:%M:%S")
      },
      "%H:%M"
    ),
    d_data$vs_sites,
    d_data$vs_burrows,
    sep = "_"
  )
}


#*=========================================================================

#' Trim Dee's PBT Honours Data to equal temporal length
d_trim_pbt_data <- function(
  d_vidObs,
  d_summary,
  i_target_total_mins,
  i_target_row_mins,
  s_target_median_time,
  print_msg
) {

  # Check if trimming is enabled
  if (is.na(i_target_total_mins) | i_target_total_mins <= 0) {
    s_msg <- sprintf(
      "Attempted to trim pbt data by time, but targeted %i mins",
      i_target_total_mins
    )
    print_msg(paste("Error:", s_msg))
    stop(s_msg)
  }

  # Split long dataframe into a list of dataframes by burrow
  ld_burrows <- split(d_vidObs, d_vidObs$vs_date_site_burrow)
  # Ensure each burrow sorted by time
  ld_burrows <- lapply(
    ld_burrows,
    function(d_burrow) d_burrow[order(d_burrow$vs_time), ]
  )

  # Parse input arguments
  dt_target_total_duration <- as.difftime(
    i_target_total_mins,
    units = "mins",
    tz = "UTC"
  )
  dt_target_row_duration <- as.difftime(
    i_target_row_mins,
    unit = "mins",
    tz = "UTC"
  )

  # Calculate number of rows of equal duration required
  i_target_rows_per_burrow <- as.integer(
    n_datetime_as_mins(dt_target_total_duration) /
    n_datetime_as_mins(dt_target_row_duration)
  )

  print_msg(sprintf(
    "Trimming each burrow to %i mins (%i rows)",
    i_target_total_mins,
    i_target_row_mins
  ))

  # Create placeholder dataframe for temporally-trimmed data
  d_trimmed <- data.frame()


  # * Temporal Trimming Loop ----------------------------------------------

  # Loop through each burrow, trim, add to output dataframe
  for (s_burrow in names(ld_burrows)) {

    # Fetch the data.frame for this burrow from the list
    d_burrow <- ld_burrows[[s_burrow]]

    # Fetch duration from Summary data
    r_summary <- match(s_burrow, d_summary$vs_date_site_burrow)
    n_vidObs_mins <- d_summary$vn_vidObs_minutes[[r_summary]]

    # Check if burrow has enough footage to warrant trimming
    if (
      nrow(d_burrow) < i_target_rows_per_burrow |
      n_vidObs_mins < i_target_total_mins
    ) {
      print_msg(sprintf(
        paste(
          "Notice: Didn't trim '%s': only %i rows (need >= %i) /",
          " only %i mins (need >= %i)."
        ),
        s_burrow,
        nrow(d_burrow),
        i_target_rows_per_burrow,
        as.integer(n_vidObs_mins),
        i_target_total_mins
      ))
      # Add to dataframe anyway
      d_trimmed <- rbind(d_trimmed, d_burrow)
      next # ignore rest of this loop cycle and start next burrow
    }

    # Check for exact number of rows
    if (nrow(d_burrow) == i_target_rows_per_burrow) {
      # No need to crop: exactly the right number of blocks
      # Naeively assuming each row is a block, no weirds...
      # Add to dataframe
      d_trimmed <- rbind(d_trimmed, d_burrow)
      next # ignore rest of this loop cycle and start next
    }

    # "Implicit else" we do need to trim.

    # Step 1: fact finding. Find row that is median time,
    # then see how many rows we have before and after it...
    vr_medians <- which(d_burrow$vs_time %in% s_target_median_time)

    # Define trimming bounding rows based on median placement
    if (length(vr_medians) > 0) { # Check if we found target median

      # Since which() can return a vector, ensure we have have just one
      # row defined by throwing away any extra rows
      r_median <- vr_medians[[1]]

      # Found target median time in dataset; set crop range accordingly
      i_rows_after_median <- nrow(d_burrow) - r_median
      i_rows_before_median <- r_median - 1
      # Step 2: work out where to do the cropping from.
      # If either side has less than half the target rows, keep
      # all of them, then crop other side.
      # We can be sure only one side can be less than target, as
      # code wouldn't let us get here otherwise.
      i_target_rows_either_side <- floor((i_target_rows_per_burrow - 1) / 2)
      if (i_rows_before_median <= i_target_rows_either_side) {
        r_start <- 1
        r_end <- i_target_rows_per_burrow
      } else if (i_rows_after_median < i_target_rows_either_side) {
        r_end <- nrow(d_burrow)
        r_start <- r_end - i_target_rows_per_burrow + 1
      } else {
        r_start <- r_median - i_target_rows_either_side - 1
        r_end <- r_median + i_target_rows_either_side
      }
    } else {
      print_msg(sprintf(
        "Notice: Target median time %s was not found in %s",
        s_target_median_time,
        s_burrow
      ))
      # Target median time not in dataset; make a choice
      # Choosing to do simple crop from end.
      r_start <- 1
      r_end <- i_target_rows_per_burrow
    }

    # Perform the trimming
    d_burrow_trimmed <- d_burrow[r_start:r_end, ]

    if (nrow(d_burrow_trimmed) != i_target_rows_per_burrow) {
      r_summary <- match(s_burrow, d_summary$vs_date_site_burrow)
      print_msg(sprintf(
        paste(
          "ERROR! Trimming algorithm failed for '%s': %i rows /",
          " %i; %f minutes of vidObs."
        ),
        s_burrow,
        nrow(d_burrow_trimmed),
        i_target_rows_per_burrow,
        d_summary$vn_vidObs_minutes[[r_summary]]
      ))
    }

    # Add this dataframe to the new master dataframe
    # N.B. this is NOT a list!
    d_trimmed <- rbind(d_trimmed, d_burrow_trimmed)
  }

  print_msg("Temporal-trimming complete.")

  return(d_trimmed)
}

#*=========================================================================

d_trim_then_summarise_pbt_data <- function(
  d_vidObs,
  d_enviro,
  d_summary,
  i_target_total_mins,
  s_target_median_time,
  print_msg
) {

  #* Parse input arguments

  # Calculate start and end times of target period
  dt_target_duration <- as.difftime(
    i_target_total_mins,
    units = "mins",
    tz = "UTC")
  t_target_median <- t_from_s(s_target_median_time, "%H:%M")

  # Regulate rows in behavioural data to be 1 minute each
  d_pbt <- d_regulate_pbt_data(
    d_vidObs,
    d_summary,
    i_target_row_mins = 1,
    print_msg
  )

  # Re-attach behavioural activity to regulated data
  d_pbt <- d_reattach_activity(
    d_pbt,
    d_vidObs,
    print_msg
  )

  # Re-attach thermal data
  d_pbt <- d_combine_pbt_data(
    d_vidObs = d_pbt,
    d_enviro
  )

  # Split df into a list of df by date/site/burrow
  ld_pbt_by_burrows <- split(d_pbt, d_pbt$vs_date_site_burrow)
  # Sort each burrow by time
  ld_pbt_by_burrows <- lapply(
    ld_pbt_by_burrows,
    function(d_burrow) d_burrow[order(d_burrow$vs_times), ]
  )

  vs_activity_cols <- c(
    "Basking",
    "Rotate",
    "Emergence",
    "Retreat",
    "Dispersal",
    "Novel Lizard",
    "NL Retreat",
    "NL Disperse",
    "NL Return",
    "Mating",
    "Feeding"
  )

  #* Loop through burrows, synthesizing summary data

  for (j_burrow in seq_len(nrow(d_summary))) {

    # Note that the row of a given burrow, say, "2021-03-17_HK_5U5" in the
    # imported data is by no means necessarily located in the same row number
    # in this new summary dataframe. Thus, we index by that which does always
    # match: our ID string, precisely "2021-03-17_HK_5U5" in this example.
    # Get the burrow ID
    s_burrow <- d_summary$vs_date_site_burrow[[j_burrow]]
    # Fetch the combined & trimmed dataframes for this burrow
    d_burrow <- ld_pbt_by_burrows[[s_burrow]]

    dt_burrow_duration <-
      t_from_s(d_summary$vs_vidObs_ends[j_burrow], "%H:%M:%S") -
      t_from_s(d_summary$vs_vidObs_starts[j_burrow], "%H:%M:%S")

    if (dt_burrow_duration <= dt_target_duration) {

      r_start <- 1
      r_end   <- nrow(d_burrow)

    } else {

      lr_bounds <- lr_duration_bounds(
        d_burrow,
        t_target_median,
        dt_target_duration
      )

      r_start <- lr_bounds[[1]]
      r_end   <- lr_bounds[[2]]

    }

    # Trim burrow to these bounds
    d_burrow <- d_burrow[r_start:r_end, ]

    # Isolate behavioural activity data
    d_activity <- d_burrow[vs_activity_cols]

    #* Push summary data into dataframe

    # Rebind this burrow's row index in the summary dataframe to be succinct
    r <- j_burrow

    # Activity Summaries
    d_summary$vn_trimmed_activity_count[r] <- sum(sum(d_activity))
    d_summary$vn_trimmed_basking_count[r] <- sum(d_activity$"Basking")
    # Positional temperature
    d_summary$vn_trimmed_temp_amb[r] <- mean(d_burrow$"MeanAmbT")
    d_summary$vn_trimmed_temp_top[r] <- mean(d_burrow$"MeanTopT")
    d_summary$vn_trimmed_temp_mid[r] <- mean(d_burrow$"MeanMidT")
    d_summary$vn_trimmed_temp_bot[r] <- mean(d_burrow$"MeanBaseT")
    # Positional humidity
    d_summary$vn_trimmed_humid_amb[r] <- mean(d_burrow$"MeanAmbH")
    d_summary$vn_trimmed_humid_top[r] <- mean(d_burrow$"MeanTopH")
    d_summary$vn_trimmed_humid_mid[r] <- mean(d_burrow$"MeanMidH")
    d_summary$vn_trimmed_humid_bot[r] <- mean(d_burrow$"MeanBaseH")
    # Mean of all temperatures
    d_summary$vn_trimmed_temp[r] <- mean(c(
      mean(d_burrow$"MeanAmbT"),
      mean(d_burrow$"MeanTopT"),
      mean(d_burrow$"MeanMidT"),
      mean(d_burrow$"MeanBaseT")
    ))
    # Mean of all humidity
    d_summary$vn_trimmed_humid[r] <- mean(c(
      mean(d_burrow$"MeanAmbH"),
      mean(d_burrow$"MeanTopH"),
      mean(d_burrow$"MeanMidH"),
      mean(d_burrow$"MeanBaseH")
    ))
  }

  return(d_summary)
}

#*-------------------------------------------------------------------------

#' Returns the \emph{r}ow containing the time nearest to the median time.
r_nearest_time <- function(
  vt_contenders,
  t_target
) {
  which.min(
    abs(
      as.numeric(
        vt_contenders - t_target,
        unit = "secs"
  )))
}

#*-------------------------------------------------------------------------

lr_duration_bounds <- function(
  d_burrow,
  t_target_median,
  dt_target_duration
) {

  # Parse time from string to POSIX
  vt_times <- t_from_s(d_burrow$vs_times, "%H:%M:%S")

  # Find the index of the row closest to the target median time
  r_nearest_median <- r_nearest_time(vt_times, t_target_median)

  # Find index of rows bounding target duration

  if (r_nearest_median == 1) {
    # Footage is entirely after target median
    t_start <- vt_times[1]
    r_start <- r_nearest_median
    t_end <- vt_times[r_start] + dt_target_duration
    r_end <- r_nearest_time(vt_times, t_end)

  } else if (r_nearest_median == nrow(d_burrow)) {
    # Footage is entirely before target median
    r_end <- r_nearest_median
    t_end <- vt_times[length(vt_times)]
    t_start <- vt_times[r_start] - dt_target_duration
    r_start <- r_nearest_time(vt_times, t_start)

  } else {
    # Median is somewhere within footage
    t_target_start <- t_target_median - dt_target_duration / 2
    t_target_end   <- t_target_median + dt_target_duration / 2

    r_start <- r_nearest_time(vt_times, t_target_start)
    r_end   <- r_nearest_time(vt_times, t_target_end)

    t_start <- vt_times[[r_start]]
    t_end   <- vt_times[[r_end]]

    td_duration <- t_end - t_start

    # Ensure we've got enough time, trying to alternate which end is grown
    l_previously_before <- FALSE
    while (td_duration < dt_target_duration) {

      if (r_start == 1 & r_end == nrow(d_burrow)) {
        # already as big as can be!
        break
      }

      if (l_previously_before) {

        if (r_end < nrow(d_burrow)) {
          r_end <- r_end + 1
          l_previously_before <- FALSE
        } else {
          r_start <- r_start - 1
          l_previously_before <- TRUE
        }

      } else {

        if (r_start > 1) {
          r_start <- r_start - 1
          l_previously_before <- TRUE
        } else {
          r_end <- r_end + 1
          l_previously_before <- FALSE
        }

      }

      t_start <- vt_times[[r_start]]
      t_end   <- vt_times[[r_end]]
      td_duration <- t_end - t_start

    }

  }

  return(list(r_start, r_end, t_start, t_end))
}

#==============================================================================#
# ===== SOLOMON BEHAVIOUR FUNCTIONS (2024) ====================================
#==============================================================================#

#' Create a CSV file for manual entry of local start times for each video
#' recording session.
#'
#'@param pSolomonCsvs [string] The absolute path to directory containing
#' Solomon CSVs; function will include CSVs from sub-directories. 
#' 
generate_start_times_csv <- function( pSolomonCsvs) {
  
  # find all the CSV inside folders inside folder
  vpFilesRelative = get_list_of_files( pSolomonCsvs)
  # vpFilesRelative will be col of str: "Xxxx/SSBBB_YYYY-MM-DD.csv", where
  #		Xxxx is sub directory name (eg "Dee_scor_JT", "Wendi")
  #		SSBBB_YYYY-MM-DD.csv is a Solomon score file
  
  # Prepare header line of output CSV
  sHeader = "Filename, Start Time as HH:mm"
  
  # Create the CSV
  fOutput <- file( "StartTimes.csv", "w")
  writeLines( sHeader, fOutput)
  
  for ( iFile in 1:length( vpFilesRelative)) {
    vpFile = vpFilesRelative[ iFile]
    sLine = ""
    for ( iPart in 1:length( vpFile)) {
      sLine = paste( sLine, vpFile[ iPart], ", ", sep = "", collapse = '')
    }
    writeLines( sLine, fOutput)
  }
  close( fOutput)
  
}

#' Import a CSV with manually-entered start times for each recording session
#' 
#' @param pStartTimesCsv [string] absolute path to CSV
#' 
#' @return dStartTimes [data-frame] containing Session ID and Starts, where Starts
#' is local, absolute time in POSIXct format.
#' 
d_from_start_times_csv <- function( pStartTimesCsv) {
	
	dStartTimes = read.csv( pStartTimesCsv, header=TRUE)
	# dStartTimes has two columns:
	#   1. Relative path & filename, eg "Dee_Scored_Bur/BU022_2022-12-09.csv"
	#   2. Manually-entered start time, eg "07:46"
	# To create a column with just the unique session ID, we need to:
	#   1. Split up the path (then only keep the bit right of the "/")
	#   2. Remove the ".csv"
	
	vpFilesRelative = dStartTimes[[ 1]]
	lvsFileParts = strsplit( vpFilesRelative, "/")
	# lvsFileParts is a list, one element for each row in dStartTimes
	# Each list item contains a vector of strings: first, any sub-directories,
	# and in the last vector position, the filename string
	lsFileNames = lapply( lvsFileParts, FUN=function( x) x[ length( x)])
	# lapply is the 'list apply' function: it applies a specified function to each
	# element of a list. Here the function is defined inline (making it a so
	# -called "anonymous function") to return the last part of the vector.
	# So, lsFileNames is a list of strings taken from the end of each string
	# vector in lvsFileParts, ideally being the '"SSBBB_YYYY-MM-DD.csv"' bit.
	lsSessionIds = lapply( lsFileNames, FUN=function( x) substring( x, 1, nchar( x) - nchar( ".csv")))
	# Now with the last 4 characters removed, we should just have "SSBBB_YYYY-MM-DD"
	dStartTimes$Session = unlist( lsSessionIds)
	# Store the sessions back in the data frame; unlist() converts list of strings
	# into a vector of strings.
	
	# We also want to get the time into a format R understands. Currently we have
	# the date as part of the session ID, and time in another column, both strings.
	# First, split out the date from the session ID: we can split after the "_"
	lvsSessionParts = strsplit( unlist( lsSessionIds), "_")
	lsDates = lapply( lvsSessionParts, FUN=function( x) x[ length( x)])
	# Now we combine the dates and times into a 'datetime'
	vtStartTimes = t_from_s( paste( unlist( lsDates), trimws( dStartTimes$Start.Time.as.HH.mm), sep=" "), "%Y-%m-%d %H:%M")
	# And stick this back into the data frame
	dStartTimes$StartTime = vtStartTimes
	
	return( dStartTimes)
}

#' Import a batch of CSV files containing PBT behaviour scored using Solomon
#' 
#' @param pSolomonCsvs [string] The absolute path to directory containing
#' Solomon CSVs; function will include CSVs from sub-directories.
#' @param dStartTimes [data frame] The start time for each CSV found in 
#' pSolomonCsvs and its sub-directories, uniquely indexed by Session string.
#' @param nRealSecondsPerFilmSecond [numeric] Scalar that converts time-lapse
#' footage into real time durations. Default = 60.
#' 
#' @return [data frame] containing a row for each film second from the Solomon
#' CSV files, each marked with Site/Burrow, Session & (local, absolute) Time.
#' 
d_from_solomon_csvs_and_start_times <- function( pSolomonCsvs, dStartTimes, nRealSecondPerFilmSecond=60) {
  
  # Check for overloaded variable: passing path to CSV of start times, not a
  # pre-imported data-frame of start times
  if(is.character(dStartTimes)&length(dStartTimes)==1){
    pStartTimesCsv=dStartTimes
    dStartTimes=d_from_start_times_csv(pStartTimesCsv)
  }
  
  # find all the CSV inside folders inside folder
  vpFilesRelative = get_list_of_files( pSolomonCsvs)
  vpFilesAbsolute = paste( pSolomonCsvs, "/", vpFilesRelative, sep="")
  
  # Create dataframe to store data; pre-set column names (avoids issue where
  # some dFile might not have a Comments column, causing issues in rbind())
  dBehaveSolomon = data.frame()
  lInitialPass = TRUE
  
  for ( iFile in 1:length( vpFilesAbsolute)) {
    
    # Isolate the Session ID (eg "SSBBB_YYYY-MM-DD") from the file-name
    vsFilenameParts = unlist( strsplit( vpFilesRelative[ iFile], "/"))
    sSessionWithExt = vsFilenameParts[ length( vsFilenameParts)]
    sSession = substring( sSessionWithExt, 1, nchar( sSessionWithExt) - nchar( ".csv"))
    
    # Get the site & burrow
    vsSessionParts = unlist( strsplit( sSession, "_"))
    sSiteBurrow = vsSessionParts[ 1]
    # sSite   = substr( vsSessionParts[ 1], 1, 2)
    # sBurrow = substr( vsSessionParts[ 1], 3, nchar( vsSessionParts[ 1]))
    
    # Get the local, absolute start time for this session
    tStartTime = dStartTimes$StartTime[ which( dStartTimes$Session == sSession)]
    
    # Read in the CSV
    dFile = read.csv( vpFilesAbsolute[ iFile], header=TRUE, stringsAsFactors=FALSE)
    # dFile will be a data-frame with x obs. of 4 variables:
    # "Time", "Basking", "Non.basking.beh", "Novel.lizards"
    
    # Sanitising the CSV: there's some strange things that happen, so we try and
    # fix the derps automatically...
    nMaxColumns = length(c("Time",
                           "Basking",
                           "Non.basking.beh",
                           "Novel.lizards",
                           "Comments"))
    if(dim(dFile)[2]>nMaxColumns){
      warning(sprintf("Detected more than %i columns in %s.csv: truncated.",
                      nMaxColumns,
                      sSession))
      dFile=dFile[,1:nMaxColumns]
    }
    
    # Rename "Time" column name for clarity: this is (time-lapsed) film time
    names( dFile)[ names( dFile) == "Time"] = "TimeFilm"
    
    # Convert time-lapse film time into real time scale
    anTimeReal_s = dFile$TimeFilm * nRealSecondPerFilmSecond
    
    # Convert real, relative time into local, absolute time
    atLocal = anTimeReal_s + tStartTime
    
    # Save time into data-frame
    dFile$Time = atLocal
    
    # Add Session, Site & Burrow to data-frame
    dFile$Session = sSession
    dFile$SiteBurrow = sSiteBurrow
    # dFile$Site = sSite
    # dFile$Burrow = sBurrow
    
    # There is a good chance the first n CSVs we parse won't have a "Comments"
    # column, then when rbind does hit one it gets confused by the additional
    # column. So we check each time if we have more columns now, and if so, we
    # try to add a Comments column to the main data frame and hope that fixes it
    if ( dim( dBehaveSolomon)[ 2] < dim( dFile)[ 2]) {
      dBehaveSolomon$Comments = rep( "", dim( dBehaveSolomon)[ 1])
    }
    
    # Bind the rows from this CSV (dFile) onto the main dataframe (dBehaveSolomon)
    dBehaveSolomon = dplyr::bind_rows( dBehaveSolomon, dFile)
    
  }
  
  return( dBehaveSolomon)
}


#' Creates a list of all files inside a directory and it's subdirectories that
#' match the specified file extension.
#' 
#' @param pFiles [string] absolute path of directory
#' @param file_extension [string] of matching files; default = "csv"
#' 
#' @return [char vec] of file-names with relative folder paths, eg
#' "All csv/BU001_2022-12-09.csv" "All csv/BU001_2022-12-10.csv"
#'
get_list_of_files <- function( pFiles, file_extension="csv") {
	return(
		list.files(
			path = pFiles, # char vec of full path names
			pattern = paste( "\\.", file_extension, "$", sep=""), # optional regular expression matching; $ means 'at end'
			ignore.case = TRUE, # pattern matching case insensitive?
			#all.files = FALSE, # include hidden files
			#full.names = FALSE, # pre-pend path to file names
			recursive = TRUE, # recurse into sub-directories?
			#include.dirs = FALSE, # include sub-directory names in recursive listings?
			#no.. = FALSE # exclude `"."` and `".."` from non-recursive listings?
		)
	)
}

#' d_from_solomon_behave_data()
#' 
#' 2024 adaptation of d_parse_behave_vidObs_data() to suit Solomon Coder data
#' 
d_from_solomon_behave_data<-function(
    dBehaveInput,
    dEnviroInput,
    sSessionName="",
    pBehaveOutput="",
    pSummaryOutput="",
    lVerboseOutput=FALSE
){
  
  #
  # Preparation
  #
  
  sFunctionVersion <- "0.01: 2024-11-26"
  
  pLogs<-"Dee PBT Utils Logs"
  if(!file.exists(pLogs))dir.create(pLogs)
  pLogfile<-paste(
    pLogs,
    "/Log_",
    if(length(sSessionName)>0){
      paste(
        sSessionName,
        "_",
        sep=""
      )
    },
    strftime(Sys.time(),"%Y-%m-%d_%H-%M-%S"),
    ".txt",
    sep = ""
  )
  
  # Convenient (succinct) version of print_message
  print_msg<-function(s_msg,
                      l_urgent = FALSE,
                      l_verbose = lVerboseOutput,
                      s_file = pLogfile
  ) {
    print_message(s_msg, s_file, l_verbose, l_urgent)
  }
  
  print_session_details(
    sprintf(
      "dee_pbt_utils::d_parse_behave_solomon_data
      Version %s
      Logfile for session %s, %s on %s",
      sFunctionVersion,
      sSessionName,
      s_from_t(Sys.time(), "%H:%M"),
      s_from_t(Sys.time(), "%a %d %b %Y")
    ),
    sFunctionVersion,
    pLogfile,
    data.frame(
      "sSessionName" = sSessionName,
      "pBehaveOutput" = pBehaveOutput,
      "pSummaryOutput" = pSummaryOutput,
      "lVerboseOutput" = lVerboseOutput
    ),
    print_msg
  )
  
  #
  # Processing
  #
  
  print_msg("Section 1. Tidy")
  
  #TODO: Done
  dBehaveTidy <- d_tidy_behave_solomon_data(dBehaveInput)
  
  #TODO: Done
  dEnviroTidy <- d_tidy_enviro_data(dEnviroInput)
  
  print_msg("Section 2. Transform")
  
  #TODO: Done
  dPbt <- d_combine_pbt_data(
    dBehaveTidy,
    dEnviroTidy
  )
  
  #TODO: Done
  dSummary <- d_summarise_solomon_pbt_data(
    dPbt,
    dEnviroTidy,
    print_msg
  )
  
  # #TODO: Delete?
  # dBehave <- d_regulate_pbt_data(
  #   dBehave = dBehaveTidy,
  #   dSummary,
  #   i_target_row_mins,
  #   i_target_total_mins,
  #   print_msg
  # )
  
  #TODO: Check
  dPbt <- d_binarize_solomon_pbt_data(
    dPbt,
    print_msg
  )
  
  # #TODO: Check
  # dBehave <- d_combine_pbt_data(
  #   d_behave = dBehave,
  #   dEnviroTidy
  # )
  
  #TODO: Delete?
  # dBehave <- d_trim_pbt_data(
  #   dBehave,
  #   dSummary,
  #   i_target_total_mins,
  #   i_target_row_mins,
  #   s_target_median_time,
  #   print_msg
  # )
  
  #TODO: Recreate. Actually, trying delete.
  # dSummary <- d_trim_then_summarise_pbt_data(
  #   d_behave = dBehaveTidy,
  #   d_thermo = dEnviroTidy,
  #   dSummary,
  #   i_target_total_mins,
  #   s_target_median_time,
  #   print_msg
  # )
  
  print_msg("Section 4. Save to Disk")
  
  save(dPbt,file=pBehaveOutput)
  save(dSummary,file=pSummaryOutput)
  
  # save_data_to_file(dPbt, pBehaveOutput, print_msg)
  # 
  # save_data_to_file(dSummary, pSummaryOutput, print_msg)
  
  
  return(dPbt)
}

d_tidy_behave_solomon_data <- function(dBehave) {
  
  # # Tidy up behavioural data
  # dBehave$Time = dBehave$StartTime
  # # Set up columns for internal function convenience
  
  # Time & Date
  dBehave$vs_dates<-s_from_t(dBehave$Time,"%Y-%m-%d")
  dBehave$vs_times<-s_from_t(dBehave$Time,"%H:%M:%S")
  dBehave$vt_times<-t_from_s(dBehave$vs_times, "%H:%M:%S") # because old fun had this
  
  # Other time columns
  dBehave$vt_vidObss    <- dBehave$vs_times #t_from_s(dBehave$"Actual Start Time", "%H:%M:%S")
  dBehave$vt_activities <- dBehave$vs_times #t_from_s(dBehave$"Actual Time", "%H:%M:%S")
  
  # SiteBurrow should follow the convention:
  # SSBBB
  # where:
  #		SS = Site, two-letter form, eg JT, KP, BU
  #		BBB = Burrow ID, eg 012, 6J22
  #TODO WARNING Parsing code for enclosures at Tarlee are not done!
  vsSitesShort<-substr(dBehave$SiteBurrow,1,2)
  dBehave$vs_sites<-lapply(
    vsSitesShort,
    function(sSiteShort)
      switch(sSiteShort,BU={"Burra"},JT={"Jamestown"},KA={"Kapunda"},TL={"Tarlee"})
    )
  dBehave$vs_burrows<-substr(dBehave$SiteBurrow,3,length(dBehave$SiteBurrow)) #TODO TEST
  
  #TODO WARNING Lineage not coded!
  # colnames(dBehave)[colnames(dBehave) == "Lineage"] <- "vs_lineages"
  # # Now we can combine enclosure and lineage as "Treatment Group"
  # dBehave$vs_treatment_groups <- paste(
  #   dBehave$vs_enclosures, # "KH5", "KH7", "Til" or "TC"
  #   dBehave$vs_lineages,
  #   sep = "-"
  # )
  dBehave$vs_enclosures      <-rep("",dim(dBehave)[1])
  dBehave$vs_lineages        <-dBehave$vs_sites
  dBehave$vs_treatment_groups<-rep("",dim(dBehave)[1])
  
  #
  # Make compound columns
  #
  
  # Thermal alignment column
  dBehave$vs_date_time30_site <- paste(
    dBehave$vs_dates,
    s_from_t(
      t_floor(
        t_from_s(dBehave$vs_times, "%H:%M:%S"),
        as.difftime(30, units = "mins", tz = "UTC")
      ),
      "%H:%M"
    ),
    dBehave$vs_sites,
    sep = "_"
  )
  
  # Per Burrow column
  dBehave$vs_date_site_burrow <- paste(
    dBehave$vs_dates,
    dBehave$vs_sites,
    dBehave$vs_burrows,
    sep = "_"
  )
  
  return(dBehave)
}

#' Summarise Dee's pygmy blue-tongue data at date/site/burrow level
#'
#' @return a dataframe with summary data, one row for each date/site/burrow
d_summarise_solomon_pbt_data <- function(
    d_Behave,
    d_enviro,
    print_msg
) {
  
  # Split df into a list of df by date/site/burrow
  ld_Behave_by_burrows <- split(d_Behave, d_Behave$vs_date_site_burrow)
  # Sort each burrow by time
  ld_Behave_by_burrows <- lapply(
    ld_Behave_by_burrows,
    function(d_burrow) d_burrow[order(d_burrow$vs_times), ]
  )
  
  # Create dataframe, preallocate with addressable ID column
  d_summary <- data.frame(
    "vs_date_site_burrow" = unique(d_Behave$vs_date_site_burrow)
  )
  
  # Loop through burrows, synthesising summary data
  for (j_burrow in seq_len(nrow(d_summary))) {
    
    # Note that the row of a given burrow, say, "2021-03-17_KH_5U5" in the
    # imported data is by no means necessarily located in the same row number
    # in this new summary dataframe. Thus, we index by that which does always
    # match: our ID string, precisely "2021-03-17_KH_5U5" in this example.
    # Get the burrow ID
    s_burrow <- d_summary$vs_date_site_burrow[[j_burrow]]
    # Fetch the actual data for this burrow
    d_burrow <- ld_Behave_by_burrows[[s_burrow]]
    # # Apparently time isn't always increasing, so be sure that it is
    # d_burrow <- d_burrow[order(d_burrow$vs_times), ]
    
    # Calculate Number of Observations
    # Solomon data has "Basking" column; any text in here means PBT is visible,
    # so we count any non-blank rows as 'active'
    vlActiveInFrame=d_burrow$Basking!="" # any observation -> visible
    nActiveFrames=sum(vlActiveInFrame)
    vlActiveDiff=vlActiveInFrame[2:length(vlActiveInFrame)]-
                 vlActiveInFrame[1:length(vlActiveInFrame)-1]
    vlBlockStarts=vlActiveDiff==1
    nBlockStarts=sum(vlBlockStarts)
    nActiveRatio=nActiveFrames/dim(d_burrow)[1]
    
    
    
    # Prepare dates for summary environmental data
    s_today <- d_burrow$vs_dates[1]
    s_yesterday <- s_from_t(
      t_from_s(d_burrow$vs_dates[1], "%Y-%m-%d") -
        as.difftime("24", format = "%H", tz = "UTC"),
      "%Y-%m-%d"
    )
    s_tomorrow <- s_from_t(
      t_from_s(d_burrow$vs_dates[1], "%Y-%m-%d") +
        as.difftime("24", format = "%H", tz = "UTC"),
      "%Y-%m-%d"
    )
    
    #* Collate data into output dataframe
    
    # Rename for brevity
    r <- j_burrow
    
    # Add meta data
    d_summary$vs_dates[r]             <- d_burrow$vs_dates[[1]]
    d_summary$vs_sites[r]             <- d_burrow$vs_sites[[1]]
    d_summary$vs_burrows[r]           <- d_burrow$vs_burrows[[1]]
    d_summary$vs_lineages[r]          <- d_burrow$vs_lineages[[1]]
    d_summary$vs_treatment_groups[r]  <- d_burrow$vs_treatment_groups[[1]]
    
    # # Push summaries into the Summary Data data.frame
    # d_summary$vn_entries[r]           <- nrow(d_burrow)
    # d_summary$vn_Behaves[r]           <- nrow(d_burrow_unique_vids)
    # d_summary$vn_Behave_minutes[r]    <- n_total_Behave_mins
    # d_summary$vs_Behave_starts[r]     <- s_from_t(t_Behave_start, "%H:%M:%S")
    # d_summary$vs_Behave_ends[r]       <- s_from_t(t_Behave_end, "%H:%M:%S")
    # d_summary$vn_modal_Behave_mins[r] <- n_mode_Behave_mins
    
    d_summary$Sum_Active_Frames[r]   <- nActiveFrames
    d_summary$Sum_Active_Periods[r]  <- nBlockStarts
    d_summary$Ratio_Active_to_Inactive_Frames[r] <- nActiveRatio
    
    d_summary$vn_temp_amb[r] <- mean(d_burrow$"MeanAmbT")
    d_summary$vn_temp_top[r] <- mean(d_burrow$"MeanTopT")
    d_summary$vn_temp_mid[r] <- mean(d_burrow$"MeanMidT")
    d_summary$vn_temp_bot[r] <- mean(d_burrow$"MeanBaseT")
    
    d_summary$vn_humid_amb[r] <- mean(d_burrow$"MeanAmbH")
    d_summary$vn_humid_top[r] <- mean(d_burrow$"MeanTopH")
    d_summary$vn_humid_mid[r] <- mean(d_burrow$"MeanMidH")
    d_summary$vn_humid_bot[r] <- mean(d_burrow$"MeanBaseH")
    
  #   d_summary$"Yesterday_Max_Celsius_Top" <- max(
  #     d_enviro$vn_temp_top[d_enviro$vs_dates == s_yesterday])
  #   d_summary$"Yesterday_Min_Celsius_Top" <- min(
  #     d_enviro$vn_temp_top[d_enviro$vs_dates == s_yesterday])
  #   d_summary$"Yesterday_Mean_Celsius_Top" <- mean(
  #     d_enviro$vn_temp_top[d_enviro$vs_dates == s_yesterday])
  #   d_summary$"Yesterday_Max_RH_Amb" <- max(
  #     d_enviro$vn_humid_amb[d_enviro$vs_dates == s_yesterday])
  #   d_summary$"Yesterday_Min_RH_Amb" <- min(
  #     d_enviro$vn_humid_amb[d_enviro$vs_dates == s_yesterday])
  #   d_summary$"Yesterday_Mean_RH_Amb" <- mean(
  #     d_enviro$vn_humid_amb[d_enviro$vs_dates == s_yesterday])
  #   
  #   d_summary$"Today_Max_Celsius_Top" <- max(
  #     d_enviro$vn_temp_top[d_enviro$vs_dates == s_today])
  #   d_summary$"Today_Min_Celsius_Top" <- min(
  #     d_enviro$vn_temp_top[d_enviro$vs_dates == s_today])
  #   d_summary$"Today_Mean_Celsius_Top" <- mean(
  #     d_enviro$vn_temp_top[d_enviro$vs_dates == s_today])
  #   d_summary$"Today_Max_RH_Amb" <- max(
  #     d_enviro$vn_humid_amb[d_enviro$vs_dates == s_today])
  #   d_summary$"Today_Min_RH_Amb" <- min(
  #     d_enviro$vn_humid_amb[d_enviro$vs_dates == s_today])
  #   d_summary$"Today_Mean_RH_Amb" <- mean(
  #     d_enviro$vn_humid_amb[d_enviro$vs_dates == s_today])
  #   
  #   d_summary$"Tomorrow_Max_Celsius_Top" <- max(
  #     d_enviro$vn_temp_top[d_enviro$vs_dates == s_tomorrow])
  #   d_summary$"Tomorrow_Min_Celsius_Top" <- min(
  #     d_enviro$vn_temp_top[d_enviro$vs_dates == s_tomorrow])
  #   d_summary$"Tomorrow_Mean_Celsius_Top" <- mean(
  #     d_enviro$vn_temp_top[d_enviro$vs_dates == s_tomorrow])
  #   d_summary$"Tomorrow_Max_RH_Amb" <- max(
  #     d_enviro$vn_humid_amb[d_enviro$vs_dates == s_tomorrow])
  #   d_summary$"Tomorrow_Min_RH_Amb" <- min(
  #     d_enviro$vn_humid_amb[d_enviro$vs_dates == s_tomorrow])
  #   d_summary$"Tomorrow_Mean_RH_Amb" <- mean(
  #     d_enviro$vn_humid_amb[d_enviro$vs_dates == s_tomorrow])
  
    }
  
  print_msg(sprintf(
    "Summary data complete: %i columns for %i burrows.",
    ncol(d_summary), nrow(d_summary)
  ))
  
  return(d_summary)
}

d_binarize_solomon_pbt_data <- function(
    dBehave,
    print_msg
) {
  
  dBehave$vl_actives = lapply(
    dBehave$Basking,
    function(sBasking)
      sBasking != ""
  )
  
  return(dBehave)
}

#==============================================================================#
#===== UTILITY FUNCTIONS ======================================================
#==============================================================================#
#
# These are custom functions, basically a bunch of steps lumped together.
# The advantage is that code is (hopefully) easier to follow.
#
#==============================================================================#

#' Trims a \link(data.frame) to exclude duplicate rows
#'
#' Returns a version of data.frame \code{d} with only the first occurence
#' of each unique value in column \code{s_col}
#'
#' @param d the dataframe to trim
#' @param s_col the name of the column to find unique values within
d_trimmed_to_unique_rows <- function(d, s_col) {
  d <- d[order(d[[s_col]]), ] # Sort d by s_col
  v <- d[[s_col]] # Extract d$s_col
  n <- length(v)
  last <- v[1:(n - 1)]
  this <- v[2:n]
  v_unique <- c(TRUE, this != last) # N.B. defines first as unique
  return(d[v_unique, ])
}

#' Calculates the statistical mode: most often occuring value.
#'
#' R's \link{mode} seemingly wraps \link{class}, not what one expects
#' in a statisics focused language...
#'
#' @param v a vector of values
#' @return a tablulation of values being the most often occuring value(s)
tab_statistical_mode <- function(v) {
  v_unique <- unique(v)
  tab <- tabulate(match(v, v_unique))
  return(v_unique[tab == max(tab)])
}

#' Wraps \code{strftime} so I don't have to keep working out if I want
#' \emph{s} or \emph{p}...
s_from_t <- function(t_time, s_format, tz = "UTC") {
  return(strftime(t_time, s_format, tz))
}

#' Wraps \code{strptime} so I don't have to keep working out if I want
#' \emph{s} or \emph{p}...
t_from_s <- function(t_time, s_format, tz = "UTC") {
  return(strptime(t_time, s_format, tz))
}

#' Round datetime \code{t} to timediff \code{d}.
#' @param t datetime to convert
#' @param d timeDiff target
t_round <- function(t, d) {
  n_seconds_in_base <- n_datetime_as_secs(d)
  s_date_origin <- "1970-01-01"
  return(
    as.POSIXct(
      round(as.numeric(t) / n_seconds_in_base) * n_seconds_in_base,
      origin = s_date_origin
    )
  )
}

#' Floor datetime \code{t} to timediff \code{d}.
#' @param t datetime to convert
#' @param d timeDiff target
t_floor <- function(t, d) {
  n_seconds_in_base <- n_datetime_as_secs(d)
  s_date_origin <- "1970-01-01"
  return(
    as.POSIXct(
      floor(as.numeric(t) / n_seconds_in_base) * n_seconds_in_base,
      origin = s_date_origin
    )
  )
}

#' Ceiling datetime \code{t} to timediff \code{d}.
#' @param t datetime to convert
#' @param d timeDiff target
t_ceiling <- function(t, d) {
  n_seconds_in_base <- n_datetime_as_secs(d)
  s_date_origin <- "1970-01-01"
  return(
    as.POSIXct(
      ceiling(as.numeric(t) / n_seconds_in_base) * n_seconds_in_base,
      origin = s_date_origin
    )
  )
}

#' Converts a datetime into numeric minutes.
#' @param t datetime to convert
n_datetime_as_mins <- function(t) {
  as.numeric(t) * switch(
    units(t),
    "secs" = 0.01666667, #>          1 / 60
    "mins" = 1, #>                   1
    "hours" = 60, #>            60 * 1
    "days" = 1440, #>      24 * 60 * 1
    "weeks" = 10080 #> 7 * 24 * 60 * 1
  )
}

#' Converts a datetime into numeric seconds.
#' @param t datetime to convert
n_datetime_as_secs <- function(t) {
  as.numeric(t) * switch(
    units(t),
    "secs"  = 1, #>                        1
    "mins"  = 60, #>                  60 * 1
    "hours" = 3600, #>           60 * 60 * 1
    "days"  = 86400, #>     24 * 60 * 60 * 1
    "weeks" = 604800 #> 7 * 24 * 60 * 60 * 1
  )
}

#' Prints a message string to file and, conditionally, logfile
print_message <- function(s_msg, s_file, l_verbose, l_urgent = FALSE) {
  print_to_file(s_msg, s_file)
  if (l_verbose | l_urgent) print_to_console(s_msg)
}

#' Print message string to file
print_to_file <- function(s_msg, s_file) {
  sink(
    file = s_file,
    append = TRUE
  )
  cat(s_msg, "\n", sep = "")
  sink()
}

#' Print message string to console
print_to_console <- function(s_msg) {
  print(s_msg)
}

#' Prints session introduction and details
print_session_details <- function(
  s_intro_msg,
  s_function_version,
  s_logfile,
  d_args,
  print_msg
) {

  print_msg(
    paste(
      s_intro_msg,
      "\n",
      "Session parameters:",
      sep = ""
  ))

  unwanted_return <- sapply(
    seq_len(ncol(d_args)),
    function(r)
    print_msg(
      sprintf(
        paste("  ",
        names(d_args)[r],
        if (is.numeric(d_args[[r]])) " = %i" else " = \"%s\"",
        sep = ""),
        d_args[[r]]
      )
    )
  )
}

#' Save to Disk
#' Writes a dataframe out to the working directory.
#' @param d_data \code{data.frame} to be saved
#' @param s_file \code{character} name the output is to be called;
#' will save as extension specified or default to \code{Rda}:
#' \code{".csv"}, for character separated values, or
#' \code{".Rda"}, for R data storage (can be \code{readRDS()}'d back in)
save_data_to_file <- function(d_data, s_file, print_msg) {

  # If no filename provided, don't do anything!
  if (s_file == "") {
    return(FALSE)
  }

  # Parse filename and split out the filetype/extension
  # '[' and ']' escape '.', which is otherwise 'wildcard for character'
  ls_file_parts <- strsplit(s_file, split = "[.]")[[1]]
  # assume last part is extension
  s_extension <- ls_file_parts[[length(ls_file_parts)]]

  # Perform the write, using appropriate function for filetype
  switch(
    tolower(s_extension),
    csv = write.csv(d_data, file = s_file),
    rda = saveRDS(d_data, file = s_file),
    {
      s_file <- paste(s_file, "Rda", sep = ".")
      saveRDS(d_data, file = s_file)
    }
  )

  print_msg(sprintf(
    "Data saved to %s.",
    s_file
  ))

  return(TRUE)

}

get_version <- function() {
  return("0.139: 23 Aug 2021")
}
