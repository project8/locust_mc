Directory Structure
-------------------

*  Tutorial - Tutorial outlining usage of locust in Phase II

File Conventions
-------------------
*  .xml files are kassiopeia configuration files, and are compatible with the current release of Kassiopeia
   Typically, to improve the modularity of the XML files, users create separate xml files specifically specifying the geometry of the electromagnetic trap (ie. Project8Geometry.xml) which can be included in other files.

    https://github.com/KATRIN-Experiment/Kassiopeia

*  .json files are generally locust configuration files. Locust uses an ordered list of generators (ie. kass-signal, lpf-fft), which either generate signals, or apply transformations, such as filtering, to these signals. These files specify both the list of generators used and the parameters of these generators, in standard JSON format.
  In the Tutorial subdirectory, there are katydid.json files, which are used as configuration for katydid. These are clearly marked and explained in the tutorial.


*  .and files are HFSS input files. If using the hfss-signal generator, this file is used as input in locust to generate an .nfd file, which is another HFSS configuration file, describing the electron fields. Both the .and and .nfd files are compatible with up to HFSS version 18.1. An issue should be raised on github if future versions are not back-compatible to the current file format.


