FROM continuumio/miniconda
FROM centos:7

# Install scl/epel
RUN yum install --assumeyes \
    centos-release-scl \
    epel-release \
    scl-utils
    
# gcc 7
RUN yum install --assumeyes \
    devtoolset-7

# Copy over conda from miniconda base image
COPY --from=0 /opt/conda /opt/conda

# Install crystalplan conda dependencies
RUN /opt/conda/bin/conda install -c conda-forge --yes \
    vtk \
    xtensor \
    xtensor-python \
    pybind11 \
    Mayavi

# Build wxpython
RUN yum install --assumeyes \
    bzip2 \
    wxGTK-devel
COPY ./build-and-install-wxpython.sh /tmp
RUN bash /tmp/build-and-install-wxpython.sh
RUN rm -f /tmp/build-and-install.sh

