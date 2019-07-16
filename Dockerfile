FROM rocker/rstudio

# install packrat
RUN R -e 'install.packages("packrat", repos="http://cran.rstudio.com", dependencies=TRUE, lib="/usr/local/lib/R/site-library");'

USER rstudio

# copy lock file & install deps
RUN R -e 'packrat::init("~/Documents/ProteoSign-2.0/")';
COPY --chown=rstudio:rstudio ~/Documents/ProteoSign-2.0/packrat/packrat.* /home/rstudio/project/packrat/
RUN R -e 'packrat::restore(project="/home/rstudio/project");'

# copy the rest of the directory
# .dockerignore can ignore some files/folders if desirable
COPY --chown=rstudio:rstudio ~/Documents/ProteoSign-2.0/ /home/rstudio/project

USER root
