# get tidyverse docker from rocker
FROM rocker/r-ver:4.3.1

# system libraries
# Try to only install system libraries you actually need
# Package Manager is a good resource to help discover system deps
RUN apt-get update --yes \
 && apt-get upgrade --yes \
 && apt-get install --yes \
libglpk-dev libxml2-dev

# install R packages required 
# Change the packages list to suit your needs
RUN R -e 'install.packages("https://packagemanager.posit.co/cran/latest/src/contrib/Archive/renv/renv_1.0.7.tar.gz")'

# Copy renv files 
WORKDIR /csf_immune_atlas
COPY renv.lock renv.lock

ENV RENV_PATHS_LIBRARY renv/library

# Restore the R environment
RUN R -e "renv::restore()"


