FROM rocker/rstudio:4.0.5

RUN apt-get update && apt-get install -y \
  tcl \
  tk \
  libxml2 \
&& rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('remotes'); \
  remotes::install_version('fitdistrplus', '1.0.9'); \
  remotes::install_version('tidyr', '1.1.0'); \
  remotes::install_version('dplyr', '1.0.0'); \
  remotes::install_version('ggplot2', '3.3.3'); \
  remotes::install_version('lme4', '1.1.27'); \
  remotes::install_version('AICcmodavg', '2.1.1'); \
  remotes::install_version('magrittr', '1.5.0'); \
  remotes::install_version('rootSolve', '1.7.0'); \
  remotes::install_version('numDeriv', '2016.8.1'); \
  remotes::install_version('rappdirs', '0.3.1'); \
  remotes::install_version('RColorBrewer', '1.1.2'); \
  remotes::install_version('scales', '0.5.0'); \
  remotes::install_version('plot3D', '1.1.0'); \
  remotes::install_version('gtable', '0.2.0'); \
  remotes::install_version('gridExtra', '2.2.1'); \
  remotes::install_version('latex2exp', '0.4.0'); \
  remotes::install_version('doParallel', '1.0.10'); \
  remotes::install_version('DHARMa', '0.1.5'); \
  remotes::install_version('flextable', '0.6.5'); \
  remotes::install_version('rmarkdown', '2.8'); \
  remotes::install_version('vegan', '2.5.7'); \
  remotes::install_version('webshot', '0.5.2'); \
  webshot::install_phantomjs()"
