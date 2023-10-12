ARG final_img_repo=ghcr.io/project8/luna_base
ARG final_img_tag=v1.3.4

ARG build_img_repo=ghcr.io/project8/kassiopeia_builder
ARG build_img_tag=v4.0.0-dev

########################
FROM ${build_img_repo}:${build_img_tag} AS build

ARG build_type=Release
ARG build_tests_exe=TRUE
ARG locust_tag=beta
ARG locust_subdir=locust
ARG build_with_kassiopeia=TRUE
ARG nproc=4

ENV LOCUST_PREFIX=${P8_ROOT}/${locust_subdir}/${locust_tag}

RUN source ${P8_ROOT}/kassiopeia/current/setup.sh &&\
    mkdir -p $LOCUST_PREFIX &&\
    chmod -R 777 $LOCUST_PREFIX/.. &&\
    cd $LOCUST_PREFIX &&\
    echo "source ${KASS_PREFIX}/setup.sh" > setup.sh &&\
    echo "export LOCUST_TAG=${locust_tag}" >> setup.sh &&\
    echo "export LOCUST_PREFIX=${LOCUST_PREFIX}" >> setup.sh &&\
    echo 'ln -sfT $LOCUST_PREFIX $LOCUST_PREFIX/../current' >> setup.sh &&\
    echo 'export PATH=$LOCUST_PREFIX/bin:$PATH' >> setup.sh &&\
    echo 'export LD_LIBRARY_PATH=$LOCUST_PREFIX/lib:$LD_LIBRARY_PATH' >> setup.sh &&\
    echo 'export LD_LIBRARY_PATH=$LOCUST_PREFIX/lib64:$LD_LIBRARY_PATH' >> setup.sh &&\
    /bin/true

COPY Config /tmp_loc_source/Config
COPY Data /tmp_loc_source/Data
COPY kassiopeia /tmp_loc_source/kassiopeia
COPY monarch /tmp_loc_source/monarch
COPY Scarab /tmp_loc_source/Scarab
COPY Source /tmp_loc_source/Source
COPY Config /tmp_loc_source/Config
COPY CMakeLists.txt /tmp_loc_source/CMakeLists.txt
COPY .git /tmp_loc_source/.git

# repeat the cmake command to get the change of install prefix to set correctly (a package_builder known issue)
RUN source $LOCUST_PREFIX/setup.sh &&\
    cd /tmp_loc_source &&\
    mkdir build &&\
    cd build &&\
    cmake .. &&\
    cmake -D CMAKE_BUILD_TYPE=$build_type \
          -D CMAKE_INSTALL_PREFIX:PATH=$LOCUST_PREFIX \
          -D DATA_INSTALL_DIR=$LOCUST_PREFIX/data \
          -D locust_mc_ENABLE_TESTING:BOOL=$build_tests_exe \
          -D locust_mc_BUILD_WITH_KASSIOPEIA:BOOL=$build_with_kassiopeia \
          -D locust_mc_PREBUILT_KASS_PREFIX:PATH=$KASS_PREFIX \
          -D locust_mc_KASS_NPROC=$nproc \
          .. &&\
    make -j$nproc install &&\
    /bin/true

########################
#FROM final_base
# Status:  Below lines are commented out for now.  If commented back in, the
# executables fail due to missing libraries. 
FROM ${final_img_repo}:${final_img_tag}

COPY --from=build $P8_ROOT $P8_ROOT
