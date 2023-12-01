# Packages
#BiocManager::install(version = "3.14")
#BiocManager::install(c("GenomicRanges", "GenomeInfoDb", "TailRank", "IRanges", "Gviz"))
#install.packages("JuliaCall")
library(JuliaCall)

#install.packages("RTIGER")
library(RTIGER)


### SETUP
# Done once
setupJulia(JULIA_HOME="/Applications/Julia-1.9.app/Contents/Resources/julia/bin")
# Needs to be run everytime we load RTIGER
sourceJulia()
