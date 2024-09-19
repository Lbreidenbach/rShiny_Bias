FROM rocker/shiny-verse:4

COPY app.R plot_dag_ui2.R /srv/shiny-server/McBias/
COPY install.R /app/

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    jags

RUN Rscript /app/install.R

