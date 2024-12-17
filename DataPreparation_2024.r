# PBT Data Preparation - November 2024
#
# This is a single-run 'script' with edit-inline variables.
# Edit strings to reflect your session ID, folder structure and file names.
#
# Pulls in behaviour observations scored from video using Solomon Coder
# https://solomon.andraspeter.com/

# Load functions from library file
source("pbtutils.r")

# Session name: how you identify this data preparation batch
sSessionName<-"Solomon2022_20241207i"

# Set relative folder root (so this code works from any computer with Dee's
# shared OneDrive folder locally synchronized)
#
# Expected folder layout:
#
# [Data analysis] (OneDrive folder)(<- pAnalysis Root)
#   ├─[Code]
#       ├─ ...
#       ├─[DataPreparation_2024] (this project)(<- working directory)
#       └─ ...
#   └─[Data]
#       ├─[Behaviour]
#           ├─ ...
#           ├─[Some folder containing Solomon scoring CSVs] (<- pSolomonCsvs)
#           └─ ...
#       └─[Environment]
#           ├─ ...
#
# Expect the root folder to be two levels up
pAnalysisRoot<-dirname(dirname(getwd())) # dirname means 'go up one level'

# Path to folder containing Solomon Coder scoring CSV files. Can be a parent
# folder containing CSVs in sub-folders.
# Note: R doesn't like Microsoft Windows' use of backslash ("\") as path separator char,
# manually change to either forward slash ("/") or two consecutive backslash ("\\")
pSolomonCsvs  <-file.path(pAnalysisRoot,"Data/Behaviour/SolomonScoring_2024/SolomonCsvs")

# CSV file containing manually-entered start times for each recording session.
# NB! Start times CSVs must not be put inside the Solomon scoring CSVs folder,
# as this program expects all CSVs in there to be scoring CSVs.
pStartTimesCsv<-file.path(pAnalysisRoot,"Data/Behaviour/SolomonScoring_2024/StartTimes_2022Season_2024-11-27.csv")

#==============================================================================#
# 1. Import & Parse Input Data ================================================
#==============================================================================#

# Generate Filming Session Start Times CSV -------------------------------------

# Scans through the file path provided looking for CSV files and adds their
# relative file path & name to a CSV that it saves to the current working
# directory as "StartTimes.csv".
# This shouldn't really be needed more than once.
# Un-comment to run: -->
#generate_start_times_csv(pSolomonCsvs)
# <-- end of code block (yes, it's just one function)

# 1.1 Import Solomon scored behaviour Data -------------------------------------

cat(sprintf("Reading Solomon scoring data from all CSVs in:\n%s",
               pSolomonCsvs))

# Raw Solomon scoring data-frame save location
pSolomonScoring<-file.path(getwd(),"OutputData/dSolomonScoring.Rda")

# Once Solomon data has been imported, this program saves the data-frame to disk
# So we can just load that previous file, if it exists
if(file.exists(pSolomonScoring)){
  load(pSolomonScoring)
  cat(sprintf("Loaded Solomon scoring data-frame:\n%s",
                pSolomonScoring))
}else{
  # If there was no existing data-frame saved to disk, do the import, parsing
  # Solomon scoring CSVs into a data-frame
  dSolomonScoring<-d_from_solomon_csvs_and_start_times(pSolomonCsvs,pStartTimesCsv)
  cat(sprintf("Imported Solomon Scoring CSVs"))
  # Save data-frame to disk
  save(dSolomonScoring,file=pSolomonScoring)
  cat(sprintf("Saved Solomon scoring data-frame to disk at:\n%s",
                pSolomonScoring))
}

# 1.2 Import Environment Data --------------------------------------------------

pEnviroInput <-file.path(pAnalysisRoot,"Data/Environment/2022 to 2023 season/ENV_22 to 23 Master.xlsx")
pEnviroOutput<-file.path(getwd(),paste("OutputData/dEnviro_",sSessionName,".Rda",sep=""))

# Environmental data: air temperature and relative humidity
dEnviro <- readxl::read_excel(
  pEnviroInput,
  path = ,
  sheet = "Collated",
  skip = 0,
  na = "#DIV/0!"
)
# Save data-frame to disk
save(dEnviro,file=pEnviroOutput)


# 1.3 Parse Behaviour and Environment Data -------------------------------------

pPbtOutput    <-file.path(getwd(),paste("OutputData/dPbt_",sSessionName,".Rda",sep=""))
pSummaryOutput<-file.path(getwd(),paste("OutputData/dSummary_",sSessionName,".Rda",sep=""))

dBehave<-d_from_solomon_behave_data(
  dBehaveInput=dSolomonScoring,
  dEnviroInput=dEnviro,
  sSessionName=sSessionName,
  pBehaveOutput=pPbtOutput,
  pSummaryOutput=pSummaryOutput,
  lVerboseOutput=FALSE
)

# 1.4 Parse Approach Data ------------------------------------------------------

pApproachInput <-file.path(pAnalysisRoot,"Data/Behaviour/PBT Dee Filming log_approach_2022 to 2024.xlsx")
pApproachOutput<-file.path(getwd(),paste("OutputData/dApproach_",sSessionName,".Rda",sep=""))

# Import the raw Approach Distance data from excel spreadsheet
dApproachImported = readxl::read_excel(
  path = pApproachInput,
  sheet = "Sheet1",
  skip = 0,
  na = "#DIV/0!"
)

# Parse Approach Distance data and add Environmental data to it, save to disk
dApproach <- d_parse_approach_data(
    d_apprch_raw = dApproachImported,
    d_enviro_raw = dEnviro,
    vs_months_to_omit = "",
    s_session_name = sSessionName,
    s_approach_data_filename = pApproachOutput
)