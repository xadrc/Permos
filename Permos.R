# ---- SETTINGS ----

cat("\014") ; options(warn = -1)

benchmain_1 <- Sys.time()

directory_main    <- getwd()
message(paste0("working directory:"))
message(paste0("  ~ ", directory_main))
directory_data    <- paste0(directory_main, "/inputs")    # Dir. data
directory_outputs <- paste0(directory_main, "/outputs")   # Dir. outputs
directory_scripts <- paste0(directory_main, "/scripts")   # Dir. scripts


# ---- PACKAGES ----

source(file.path(paste0(directory_scripts, "/PCKG.R")), echo = FALSE)

# ---- FUNCTIONS ----

source(file.path(paste0(directory_scripts, "/FUNC.R")), echo = FALSE)

# ---- DATA ----

# load data
source(file.path(paste0(directory_scripts, "/DATA.R")), echo = FALSE)

# ---- FILLING METHODS ----

# compute different gap-filling methods
source(file.path(paste0(directory_scripts, "/FNAS.R")),   echo = FALSE)
# plots
source(file.path(paste0(directory_scripts, "/FNAS_P.R")), echo = FALSE)

# ---- MODELLING & BACKTESTING ----

# model & test
source(file.path(paste0(directory_scripts, "/MFTS.R")),   echo = FALSE)
# plot results
source(file.path(paste0(directory_scripts, "/MFTS_P.R")), echo = FALSE)

# ---- BENCHMARKING ----

benchmain_2 <- Sys.time()

message("__________________________")
message(paste0("TOTAL RUN TIME: ", round(difftime(benchmain_2, benchmain_1,units = "mins"), 3), " MIN"))
rm(benchmain_1, benchmain_2)

 