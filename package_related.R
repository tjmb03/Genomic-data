# Replace "pkg" with the name of the package you're interested in
pkg <- "sva"

# Use available.packages() to get information about all available packages
av <- available.packages(filters = list())

# Subset the data frame to get information about the specific package
package_info <- av[av[, "Package"] == pkg, ]

# Print or inspect the information
print(package_info)
print(av)

BiocManager::valid()
help("repositories", package =
       "BiocManager")

brew install gcc
Sys.setenv(PATH = paste("/opt/homebrew/bin", Sys.getenv("PATH"), sep=":"))

# Check if the PATH is updated
Sys.getenv("PATH")
BiocManager::install(c("genefilter", "sva", "edge"))

# Get the current working directory
current_directory <- getwd()

# Print the current working directory
print(current_directory)



install.packages("genefilter_1.84.0.tgz", repos = NULL, type = "source")
library(genefilter)
install.packages("sva_3.50.0.tgz", repos = NULL, type = "source")
library(sva)
install.packages("edge_2.34.0.tgz", repos = NULL, type = "source")
library(edge)

