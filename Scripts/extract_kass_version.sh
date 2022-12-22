#! /bin/bash

# Extracts the version of Kassiopeia from its CMakeLists.txt file
# Usage:
#   > extract_kass_version.sh /path/to/kassiopeia/CMakeLists.txt

kass_cmake_file=${1}

ver_major=$(grep "set(KASPER_VERSION_MAJOR" ${kass_cmake_file} | tr -dc '0-9')
ver_minor=$(grep "set(KASPER_VERSION_MINOR" ${kass_cmake_file} | tr -dc '0-9')
ver_patch=$(grep "set(KASPER_VERSION_PATCH" ${kass_cmake_file} | tr -dc '0-9')
echo ${ver_major}.${ver_minor}.${ver_patch}
