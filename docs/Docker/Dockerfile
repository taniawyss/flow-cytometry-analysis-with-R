FROM bioconductor/bioconductor_docker:RELEASE_3_17

RUN sudo apt-get update

RUN sudo apt-get install -y \
libxt6 \
libgl1-mesa-glx \
libcurl4-openssl-dev \
libglpk-dev \
libxml2-dev \
libproj-dev \
libudunits2-dev \
libgdal-dev \
tree \
libcairo2-dev \
libxt-dev \
less

RUN cd /

RUN wget https://raw.githubusercontent.com/taniawyss/flow-cytometry-analysis-with-R/master/docs/assets/scripts/install_packages_TW.r

RUN Rscript /install_packages_TW.r

# download data on S3 with wget