FROM continuumio/miniconda3:24.7.1-0

COPY environment.yml /tmp/environment.yml

RUN conda env update -n base -f /tmp/environment.yml && rm /tmp/environment.yml

