# DockerFile for PCR-GLOB model. The ini-file should be mounted as config.ini,
# the input data root directory should be mounted as /data
FROM continuumio/anaconda3
LABEL maintainer="Willem van Verseveld <Willem.vanVerseveld@deltares.nl>"

# build
COPY . /opt/wflow/
WORKDIR /opt/wflow
RUN conda env create -f environment.yml
RUN conda run -n wflow python setup.py install
VOLUME /data
WORKDIR /
ENTRYPOINT ["conda", "run", "-n", "wflow", "python3", "-m", "wflow.wflow_sbm"]
