FROM python:2.7

RUN pip --no-cache-dir install \
        numpy \
        scipy \
        astropy \
        h5py \
        healpy

WORKDIR /pygdsm
COPY . .
RUN python setup.py install

CMD ["/bin/bash"]
