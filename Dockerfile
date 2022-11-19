ARG final_img_repo=ghcr.io/project8/luna_base
ARG final_img_tag=v1.3.0

ARG img_repo=ghcr.io/project8/luna_base
ARG img_tag=v1.3.0-dev

########################
FROM ${final_img_repo}:${final_img_tag} AS final_base

########################
FROM ${img_repo}:${img_tag} AS base

ARG build_type=Release
ENV LOCUST_BUILD_TYPE=$build_type
ARG build_tests_exe=FALSE
ENV LOCUST_BUILD_TESTS_EXE=$build_tests_exe

ARG locust_tag=beta
ENV LOCUST_TAG=${locust_tag}
ARG build_with_kassiopeia=TRUE
ENV LOCUST_BUILD_WITH_KASSIOPEIA=$build_with_kassiopeia
ENV LOCUST_BUILD_PREFIX=/usr/local/p8/locust/${LOCUST_TAG}

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
FROM base AS build

COPY Config /tmp_source/Config
COPY Data /tmp_source/Data
COPY kassiopeia /tmp_source/kassiopeia
COPY monarch /tmp_source/monarch
COPY Scarab /tmp_source/Scarab
COPY Source /tmp_source/Source
COPY Config /tmp_source/Config
COPY CMakeLists.txt /tmp_source/CMakeLists.txt
COPY .git /tmp_source/.git

# repeat the cmake command to get the change of install prefix to set correctly (a package_builder known issue)
RUN source $LOCUST_BUILD_PREFIX/setup.sh &&\
    cd /tmp_source &&\
    mkdir build &&\
    cd build &&\
    cmake -D CMAKE_BUILD_TYPE=$LOCUST_BUILD_TYPE \
          -D CMAKE_INSTALL_PREFIX:PATH=$LOCUST_BUILD_PREFIX \
          -D DATA_INSTALL_DIR=$LOCUST_BUILD_PREFIX/data \
          -D locust_mc_ENABLE_TESTING:BOOL=$LOCUST_BUILD_TESTS_EXE \
          -D locust_mc_BUILD_WITH_KASSIOPEIA:BOOL=$LOCUST_BUILD_WITH_KASSIOPEIA .. &&\
    cmake .. &&\
    make -j install &&\
    /bin/true

########################
FROM final_base

COPY --from=build $LOCUST_BUILD_PREFIX $LOCUST_BUILD_PREFIX
