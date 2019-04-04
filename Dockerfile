# DockerFile for PCR-GLOB model. The ini-file should be mounted as config.ini,
# the input data root directory should be mounted as /data
FROM ewatercycle/pcraster-container:421
MAINTAINER Gijs van den Oord <g.vandenoord@esciencecenter.nl>
RUN apt-get update -y

# INSTALL compilers and build toold
RUN apt-get install -y python3-setuptools python3-pip python3-gdal

# INSTALL pip packages
RUN pip3 install numba

# build
COPY . /opt/wflow/
WORKDIR /opt/wflow
RUN python3 setup.py install
VOLUME /data
ENV PYTHONPATH /opt/pcraster/python/
WORKDIR /
ENTRYPOINT ["python3","/usr/local/bin/wflow_sbm.py","-C","/data"]
