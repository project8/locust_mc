ARG IMG_USER=project8
ARG IMG_REPO=p8compute_dependencies
ARG IMG_TAG=v1.0.0

FROM ${IMG_USER}/${IMG_REPO}:${IMG_TAG} as locust_common

ARG build_type=Release
ENV LOCUST_BUILD_TYPE=$build_type

ARG LOCUST_TAG=beta
ENV LOCUST_TAG=${LOCUST_TAG}
ENV LOCUST_BUILD_PREFIX=/usr/local/p8/locust/$LOCUST_TAG

ARG CC_VAL=gcc
ENV CC=${CC_VAL}
ARG CXX_VAL=g++
ENV CXX=${CXX_VAL}

SHELL ["/bin/bash", "-c"]

RUN mkdir -p $LOCUST_BUILD_PREFIX &&\
    chmod -R 777 $LOCUST_BUILD_PREFIX/.. &&\
    cd $LOCUST_BUILD_PREFIX &&\
    echo "source ${COMMON_BUILD_PREFIX}/setup.sh" > setup.sh &&\
    echo "export LOCUST_TAG=${LOCUST_TAG}" >> setup.sh &&\
    echo "export LOCUST_BUILD_PREFIX=${LOCUST_BUILD_PREFIX}" >> setup.sh &&\
    echo 'ln -sfT $LOCUST_BUILD_PREFIX $LOCUST_BUILD_PREFIX/../current' >> setup.sh &&\
    echo 'export PATH=$LOCUST_BUILD_PREFIX/bin:$PATH' >> setup.sh &&\
    echo 'export LD_LIBRARY_PATH=$LOCUST_BUILD_PREFIX/lib:$LD_LIBRARY_PATH' >> setup.sh &&\
    echo 'export LD_LIBRARY_PATH=$LOCUST_BUILD_PREFIX/lib64:$LD_LIBRARY_PATH' >> setup.sh &&\
    /bin/true

########################
FROM locust_common as locust_done

COPY Config /tmp_source/Config
RUN true
COPY Data /tmp_source/Data
RUN true
COPY kassiopeia /tmp_source/kassiopeia
RUN true
COPY monarch /tmp_source/monarch
RUN true
COPY Scarab /tmp_source/Scarab
RUN true
COPY Source /tmp_source/Source
RUN true
COPY Config /tmp_source/Config
RUN true
COPY CMakeLists.txt /tmp_source/CMakeLists.txt
RUN true
COPY .git /tmp_source/.git

# repeat the cmake command to get the change of install prefix to set correctly (a package_builder known issue)
RUN source $LOCUST_BUILD_PREFIX/setup.sh &&\
    cd /tmp_source &&\
    mkdir build &&\
    cd build &&\
    cmake -D CMAKE_BUILD_TYPE=$LOCUST_BUILD_TYPE \
          -D CMAKE_INSTALL_PREFIX:PATH=$LOCUST_BUILD_PREFIX \
          -D DATA_INSTALL_DIR=$LOCUST_BUILD_PREFIX/data \
          -D SET_INSTALL_PREFIX_TO_DEFAULT=FALSE \
          -D locust_mc_BUILD_WITH_KASSIOPEIA=TRUE .. &&\
    cmake -D CMAKE_BUILD_TYPE=$LOCUST_BUILD_TYPE \
          -D CMAKE_INSTALL_PREFIX:PATH=$LOCUST_BUILD_PREFIX \
          -D DATA_INSTALL_DIR=$LOCUST_BUILD_PREFIX/data \
          -D SET_INSTALL_PREFIX_TO_DEFAULT=FALSE \
          -D locust_mc_BUILD_WITH_KASSIOPEIA=TRUE .. &&\
    make -j3 install &&\
    /bin/true


########################
FROM locust_common

COPY --from=locust_done $LOCUST_BUILD_PREFIX $LOCUST_BUILD_PREFIX
