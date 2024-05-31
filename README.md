# locust_mc

locust_mc is a monte carlo simulation package for simulating RF signals.

### Instructions For Use:

After installation, Locust runs on the command line as
```
  $ LocustSim config=/path/to/config/LocustConfigFile.json
```
The output of the simulation is a file in egg format.  It is recommended to process the egg files with Katydid with an appropriate config file, for example as
```
  $ Katydid -c /path/to/config/katydid.json
```
### Tutorial

Tutorial documents with examples are located in the locust-tutorial repository, [here](https://github.com/project8/locust-tutorial). 

Issues should be posted via [GitHub](https://github.com/project8/locust_mc/issues).

Installation
---------------

### Developing in the Project 8 luna computing environment (recommended):
*  Pull the luna development container:  ```ghcr.io/project8/luna:latest-dev``` .
*  Clone Locust:  ```git clone git@github.com:project8/locust_mc.git``` .
*  Start the luna container, replacing ```/path/to/locust_mc``` with the path to your local clone of Locust:  ```docker run -it --rm -v /path/to/locust_mc:/locust ghcr.io/project8/luna:latest-dev /bin/bash```.
*  Compile your local installation:
    *  source /usr/local/p8/compute/v1.1.1/setup.sh 
    *  cd /locust
    *  mkdir cbuild
    *  cd cbuild
    *  cmake -Dlocust_mc_PREBUILT_KASS_PREFIX=${KASS_PREFIX} -Dlocust_mc_BUILD_WITH_KASSIOPEIA=ON ../
    *  make install
    *  source /Kassiopeia/install/bin/kasperenv.sh
    *  export LD_LIBRARY_PATH=/locust/cbuild/lib:$LD_LIBRARY_PATH
* Test the executable:  ```bin/LocustSim```.

### Running the software in the pre-built Docker container:

* A Docker container with Locust is available at ```ghcr.io/project8/locust_mc:latest``` .  
* This Docker container can be used with Docker (with superuser privileges).
* Pull the container to your system with
```
> sudo docker pull ghcr.io/project8/locust_mc:latest
```
* Open a shell inside the container with
```
> sudo docker run --rm -it ghcr.io/project8/locust_mc:latest /bin/bash
```
* The Docker container can also be used with Singularity (as a regular user without sudo).  Pull the container to your system with
```
> singularity pull docker://ghcr.io/project8/locust_mc:latest
```

* Open a shell inside the container with
```
> singularity shell /path/to/local/sif/[nameOfContainer].sif
```
* Optionally, mount a local directory for e.g. I/O to the container
```
> singularity shell --no-home --bind /path/to/local/directory:/usr/local/p8/locust/current/output ./[nameOfContainer].sif 
```



Directory Structure
-------------------

*  Config - Example configuration files
*  Documentation - Doxygen-based code documentation.
*  kassiopeia - Submodule for simulation of charged particle trajectories
*  monarch - Submodule library for file I/O in Project 8
*  Scarab - Submodule
*  Source - Locust source code
