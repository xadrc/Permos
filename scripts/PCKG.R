###

message("installing required packages...")

req <- readLines(paste0(directory_main, "/", "requirements.txt"), warn = FALSE)

for (package in req){
    if (!require(package, character.only = TRUE, quietly = TRUE)){
        install.packages(package)
        library(package, character.only = TRUE)
    }
}

rm(req, package)

message("Required packages installed")

###