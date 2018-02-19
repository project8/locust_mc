# Configuration - Number of channels and physical dimensions
Number of amplifier channels "nchannels" is defined in LMCSimulationController.
CENTER_TO_SHORT:  distance, defined in LMCGlobals.
CENTER_TO_ANTENNA:  distance, defined in LMCGlobals.

nominal values:
Test Signal:  nchannels = 5
Phase 1:  nchannels = 1, CENTER_TO_SHORT = 0.045 m, CENTER_TO_ANTENNA = 0.045 m 
Phase 2:  nchannels = 1, CENTER_TO_SHORT = 0.075 m, CENTER_TO_ANTENNA = 0.075 m
Phase 3:  nchannels = 10+

# Tutorial -- Simulate a Test Signal.
The LMCTestSignalGenerator defines a sinusoidal signal that generates a time series of complex voltages.  In the generator it is mixed with a sinusoidal local oscillator signal and is sampled in the resulting intermediate band.  After low-pass filtering, decimation, digitization (all in subsequent Locust generators) and fft (in Katydid) it results in a single test pulse at the intermediate frequency.

First, compile Locust without Kassiopeia by switching the ccmake option locust_mc_BUILD_WITH_KASSIOPEIA to "OFF".  After compiling, edit the "generators" field of LocustTemplate.json to contain:  "test-signal", "lpf-fft", "decimate-signal", "digitizer".  Run the simulation as
    ```
         /path/to/LocustSim config=/path/to/LocustTemplate.json
    ```
The result will be an egg file containing a time series of complex voltages.  Process the egg file with Katydid as:
    ```
         /path/to/Katydid config=/path/to/katydid_new.json
    ```
which will generate a file called basic.root.  Open the file and verify that the test pulse in each channel has the correct frequency and power.


# Tutorial -- Phase II Simulation
In this tutorial, we will be simulating a 89.96 degree electron in a Phase II Project 8 setup. This tutorial assumes that the user has already successfully installed locust and added the executables to their $PATH. Installation instructions and prerequisites can be found in locust_mc/README.md. In addition, this tutorial uses the optional dependency VTK, for visualization.

## Visualize the Electron Trajectory (requires VTK)
1. Open LocustTemplate.json with your preferred text editor. Edit the "xml-filename" and "egg-filename" paths to point to the locust_mc directory on your machine.
2. In Project8Phase2_WithRoot_Template.xml, modify the path of the Phase II geometry file to that appropriate for your machine.

    ```
        <geometry>
            <include name="/home/XXX/locust_mc/Config/Tutorial/Project8Phase2Geometry.xml"/>
        </geometry>
    ```
3. Edit the maximum time in the Project8Phase2_WithRoot_Template xml file to 3e-5 s.
    ```
        <ksterm_max_time name="term_max_time" time="5.e-4"/>
    ```
4. Uncomment the kswrite_vtk and vtk_window blocks in the above xml file.

5. In the terminal, run locust with the command:

    ```
         /path/to/LocustSim config=/path/to/LocustTemplate.json
    ```
A VTK window should appear showing the particle track.

## Visualize the Field Map (requires VTK)
1. Edit LocustTemplate.json so that the "xml-filename" points to the "Project8Phase2_FieldMap.xml" file in the current directory.  Alternatively, make the following edits in "Project8Phase2_WithRoot_Template.xml":  define the "generator" field to be gen_bfieldlines, define the trajectory field to be "traj_magnetic", and set P8Phase=0 in the step modifier "cycl_rad_extr".
2. Look at the parameters in the field "gen_bfieldlines":

    ```
            <r_set value_start="0.0" value_stop="0.07" value_count="1"/>
            <phi_set value_start="0.0" value_stop="360.0" value_count="11"/>
            <z_list add_value="-0.12"/>
    ```
3. Run the simulation with:

    ```
         /path/to/LocustSim config=/path/to/LocustTemplate.json
    ```
4. Look at the resulting field map in the VTK window. Now there will be a ROOT file with the magnetic field information.
5. Plot the field map in ROOT. In the root terminal:
```
 >.L plotoutput.C
 > fieldmap()
```


## Create a simulated Waterfall Plot with Locust+Katydid
1. Edit LocustTemplate.json so that "xml-filename" points to "Project8Phase2_WithRoot_Template.xml"
2. Change the max_time field in the xml file to 5e-4.
3. Comment the kswrite_vtk and vtk_window blocks that you have uncommented in the first exercise. Visualization becomes unwieldy for large simulations.
4. Run the simulation as before. The expected run time is about half an hour!  Make sure the following generators are listed in the "generators" field of LocustTemplate.json:  "kass-signal", "lpf-fft", "decimate-signal", "digitizer".
5. Process the resulting egg file of IQ voltages using Katydid and the config file katydid_new.json with:

    ```
         /path/to/Katydid config=/path/to/katydid_new.json
    ```
This will result in 2 root files:  katydidwaterfall.root and basic.root .
6. Open the file plotoutput.C and edit any paths to point to the correct ones on your machine. 
7. Start root:
```
 >.L plotoutput.C
 > katydidwaterfall()
```

# Tutorial -- Phase I Simulation
Follow the instructions for the Phase II simulation above.  Filenames should be changed to Project8Phase1_WithRoot_Template.xml and Project8Phase1Geometry.xml .  The LO will need to be retuned.

