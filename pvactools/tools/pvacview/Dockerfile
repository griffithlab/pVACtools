FROM rocker/shiny-verse

RUN R -e "install.packages(c('shiny',\
 'shinyjs',\
  'shinythemes',\
  'plotly', \
  'RMySQL',\
  'ggplot2',\
  'DT',\
  'reshape2',\
  'jsonlite',\
  'tibble',\
  'tidyr',\
  'plyr',\
  'dplyr',\
  'stringr', \
  'shinydashboard', \
  'shinydashboardPlus', \
  'shinycssloaders', \
  'fresh', \
  'RCurl' ), repos='https://cran.rstudio.com/')"

RUN R -e "devtools::install_github('eclarke/ggbeeswarm', ref='v0.6.1')"

RUN rm -rf /srv/shiny-server/*

COPY . /srv/shiny-server/

EXPOSE 3333

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', port = 3333, host='0.0.0.0')"]