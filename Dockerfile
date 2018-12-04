FROM project8/p8compute_dependencies:v0.2.0 as locust_common

ARG build_type=Release
ENV LOCUST_BUILD_TYPE=$build_type

ENV LOCUST_TAG=v1.8.2
ENV LOCUST_BUILD_PREFIX=/usr/local/p8/locust/$LOCUST_TAG

RUN mkdir -p $LOCUST_BUILD_PREFIX &&\
    cd $LOCUST_BUILD_PREFIX &&\
    echo "source ${COMMON_BUILD_PREFIX}/setup.sh" > setup.sh &&\
    echo "export LOCUST_TAG=${LOCUST_TAG}" >> setup.sh &&\
    echo "export LOCUST_BUILD_PREFIX=${LOCUST_BUILD_PREFIX}" >> setup.sh &&\
    echo 'ln -sf $LOCUST_BUILD_PREFIX $LOCUST_BUILD_PREFIX/../current || /bin/true' >> setup.sh &&\
    echo 'export PATH=$LOCUST_BUILD_PREFIX/bin:$PATH' >> setup.sh &&\
    echo 'export LD_LIBRARY_PATH=$LOCUST_BUILD_PREFIX/lib:$LD_LIBRARY_PATH' >> setup.sh &&\
    /bin/true

########################
FROM locust_common as locust_done

# repeat the cmake command to get the change of install prefix to set correctly (a package_builder known issue)
RUN source $LOCUST_BUILD_PREFIX/setup.sh &&\
    mkdir /tmp_install &&\
    cd /tmp_install &&\
    git clone https://github.com/project8/locust_mc &&\
    cd locust_mc &&\
    git fetch && git fetch --tags &&\
    git checkout $LOCUST_TAG &&\
    git submodule update --init --recursive &&\
    mkdir build &&\
    cd build &&\
    cmake -D CMAKE_BUILD_TYPE=$LOCUST_BUILD_TYPE \
          -D CMAKE_INSTALL_PREFIX:PATH=$LOCUST_BUILD_PREFIX \
          -D locust_mc_BUILD_WITH_KASSIOPEIA=TRUE .. &&\
    cmake -D CMAKE_BUILD_TYPE=$LOCUST_BUILD_TYPE \
          -D CMAKE_INSTALL_PREFIX:PATH=$LOCUST_BUILD_PREFIX \
          -D locust_mc_BUILD_WITH_KASSIOPEIA=TRUE .. &&\
    make -j3 install &&\
    /bin/true

########################
FROM locust_common

COPY --from=locust_done $LOCUST_BUILD_PREFIX $LOCUST_BUILD_PREFIX
