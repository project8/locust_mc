Contents
-------------------

*  This directory contains a selection of config files for running short examples with locust.

File Conventions
-------------------
*  .xml files are Kassiopeia configuration files, and are compatible with the current release of Kassiopeia.
   Typically, to improve the modularity of the XML files, users create separate xml files for specifying the geometry of the electromagnetic trap (i.e. WaveguideGeometry.xml) which can then be included by other xml files.  See https://github.com/KATRIN-Experiment/Kassiopeia .

*  .json files are generally locust configuration files. Locust uses an ordered list of generators (i.e. kass-signal, lpf-fft), which either generate signals, or apply transformations, such as filtering, to these signals. These files specify both the list of generators used and the parameters of these generators, in standard JSON format.  This directory also contains one JSON file for use with Katydid.

Examples
--------------------
* Test signal:  A sinusoidal voltage signal drives the LMCTestSignalGenerator.  Configurable parameters are "rf-frequency", "lo-frequency", and [voltage] "amplitude".  Katydid signal processing results in a signal at 50 MHz above DC, with total power=amplitude^2.  A system impedance of 50 ohms is applied explicitly in Locust and in Katydid.  To run this example, enter
```LocustSim config=[/path/to/config]/LocustTestSignal.json```

* Fake Track:  Tracks with predefined characteristics are generated in LMCFakeTrackSignalGenerator.  Configurable parameters are described in the header file ```LMCFakeTrackSignalGenerator.hh```.  To run this example, enter ```LocustSim config=[/path/to/config]/LocustFakeTrack.json```

* Plane wave:  A plane wave is incident on a strip of passively-combined patch antenna elements, each with transfer function PatchTFLocust.txt.  Configurable parameters of the plane wave are frequency (transmitter-frequency), planewave-amplitude in V/m, and angle of incidence (AOI).  

* Magnetic dipole antenna transmitter:  A magnetic dipole antenna is driven with a voltage signal with amplitude antenna-voltage-amplitude.  Its emitted electric field is incident on a single patch antenna with transfer function tf-transmitter-filename.  Configurable parameters are frequency (transmitter-frequency), amplitude, and its transfer function tf-receiver-filename.  Additional configurable parameters such as position and orientation are described in classes LMCTransmitterHardware and in LMCDipoleAntenna.

* Waveguide simulation:  A simulated Kassiopeia electron is trapped magnetically inside a waveguide.  Two waveguide geometries are presently implemented in Locust:  WR42 and circular 0.396'' diameter, where WR42 is the default case selected in ```LMCKassSignalGenerator.cc``` and in ```LMCCyclotronRadiationExtractor.cc``` .  Electron kinematics and magnetic fields are defined in ```LocustKass_Waveguide_Template.xml``` which includes the general trap geometry file ```WaveguideGeometry.xml```.  Electron duration is limited to 1.e-6 s by the field `term_max_time` in `LocustKass_Waveguide_Template.xml`.  The number of electrons to simulate serially is set with the field `events`.  Time between electrons is set to 15000 samples with the parameter `event-spacing-samples` in `LocustWaveguideTemplate.json`.  To run this example, enter ```LocustSim config=[/path/to/config]/LocustWaveguideTemplate.json```

* Free Space simulation:  A simulated Kassiopeia electron is trapped magnetically in free space.  Electron kinematics and magnetic fields are defined in `LocustKass_FreeSpace_Template.xml` which includes the general trap geometry file `FreeSpaceGeometry.xml`.  Electron duration is limited to 1.e-6 s by the field `term_max_time` in `LocustKass_FreeSpace_Template.xml`.  The number of electrons to simulate serially is set with the field `events`.  Time between electrons is set to 15000 samples with the parameter `event-spacing-samples` in `LocustFreeSpaceTemplate.json`.  To run this example, enter ```LocustSim config=[/path/to/config]/LocustFreeSpaceTemplate.json```

* Kassiopeia application:  Kassiopeia is compiled as a standalone application within the Locust installation.  A general configuration file for simulating trapped electrons with keV energies is `JustKass.xml`.  It includes the general magnetic trap geometry file `WaveguideGeometry.xml` .  To run this application, enter `LMCKassiopeia [/path/to/config]/JustKass.xml`

* Kassiopeia field map:  A config file `JustKassFieldMap.xml` is provided to check the magnetic field solution in Kassiopeia.  First enter `LMCKassiopeia [/path/to/config]/JustKassFieldMap.xml` .  The output Root file will contain magnetic field information.  After generating the Root file, the longitudinal magnetic field component can be inspected in interactive Root by entering `.L plotfieldmap.c; fieldmap();`

