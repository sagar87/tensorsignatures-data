FROM bioconductor/bioconductor_docker:devel
RUN R -e 'BiocManager::install("VariantAnnotation")'
RUN R -e 'BiocManager::install("rhdf5")'

WORKDIR /usr/src/app
COPY . /usr/src/app/

ENTRYPOINT ["/usr/src/app/entrypoint.sh"]