Contents
-------------------

*  This directory contains a selection of config files for running short examples with locust.

File Conventions
-------------------
*  .xml files are Kassiopeia configuration files, and are compatible with the current release of Kassiopeia.
   Typically, to improve the modularity of the XML files, users create separate xml files for specifying the geometry of the electromagnetic trap (i.e. WaveguideGeometry.xml) which can then be included by other xml files.  See https://github.com/KATRIN-Experiment/Kassiopeia .

*  .json files are generally locust configuration files. Locust uses an ordered list of generators (i.e. kass-signal, lpf-fft), which either generate signals, or apply transformations, such as filtering, to these signals. These files specify both the list of generators used and the parameters of these generators, in standard JSON format.  This directory also contains one JSON file for use with Katydid.

*  .and files are HFSS input files. If using the hfss-signal generator, this file is used as input in locust to generate an .nfd file, which is another HFSS configuration file, describing the electron fields. Both the .and and .nfd files are compatible with up to HFSS version 18.1. An issue should be raised on github if future versions are not back-compatible to the current file format.

Examples
--------------------
* Test signal:  A sinusoidal voltage signal drives the LMCTestSignalGenerator.  Configurable parameters are "rf-frequency", "lo-frequency", and [voltage] "amplitude".  Katydid signal processing results in a signal at 50 MHz above DC, with total power=amplitude^2.  A system impedance of 50 ohms is applied explicitly in Locust and in Katydid.  To run this example, enter
```LocustSim config=[/path/to/config]/LocustTestSignal.json```

* Fake Track:  Tracks with predefined characteristics are generated in LMCFakeTrackSignalGenerator.  Configurable parameters are described in the header file ```LMCFakeTrackSignalGenerator.hh```.  To run this example, enter ```LocustSim config=[/path/to/config]/LocustFakeTrack.json```

* Waveguide simulation:  A simulated Kassiopeia electron is trapped magnetically inside a waveguide.  Two waveguide geometries are presently implemented in Locust:  WR42 and circular 0.396'' diameter, where WR42 is the default case selected in ```LMCKassSignalGenerator.cc``` and in ```LMCCyclotronRadiationExtractor.cc``` .  Electron kinematics and magnetic fields are defined in ```LocustKass_Waveguide_Template.xml``` which includes the general trap geometry file ```WaveguideGeometry.xml```.  Electron duration is limited to 1.e-6 s by the field `term_max_time` in `LocustKass_Waveguide_Template.xml`.  The number of electrons to simulate serially is set with the field `events`.  Time between electrons is set to 15000 samples with the parameter `event-spacing-samples` in `LocustWaveguideTemplate.json`.  To run this example, enter ```LocustSim config=[/path/to/config]/LocustWaveguideTemplate.json```

* Free Space simulation:  A simulated Kassiopeia electron is trapped magnetically in free space.  Electron kinematics and magnetic fields are defined in `LocustKass_FreeSpace_Template.xml` which includes the general trap geometry file `FreeSpaceGeometry.xml`.  Electron duration is limited to 1.e-6 s by the field `term_max_time` in `LocustKass_FreeSpace_Template.xml`.  The number of electrons to simulate serially is set with the field `events`.  Time between electrons is set to 15000 samples with the parameter `event-spacing-samples` in `LocustFreeSpaceTemplate.json`.  To run this example, enter ```LocustSim config=[/path/to/config]/LocustFreeSpaceTemplate.json```

* Kassiopeia application:  Kassiopeia is compiled as a standalone application within the Locust installation.  A general configuration file for simulating trapped electrons with keV energies is `JustKass.xml`.  It includes the general magnetic trap geometry file `WaveguideGeometry.xml` .  To run this application, enter `LMCKassiopeia [/path/to/config]/JustKass.xml`

* Kassiopeia field map:  A config file `JustKassFieldMap.xml` is provided to check the magnetic field solution in Kassiopeia.  First enter `LMCKassiopeia [/path/to/config]/JustKassFieldMap.xml` .  The output Root file will contain magnetic field information.  After generating the Root file, the longitudinal magnetic field component can be inspected in interactive Root by entering `.L plotfieldmap.c; fieldmap();`

