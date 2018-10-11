# Running Locust
----------------

## Phase 1

Run the	simulation like this:
```
/path/to/LocustSim config=~/locust_mc/Config/Tutorial/LocustPhase1Template.json
```

In the LocustPhase1Template.json file, the "generators" field should contain these generators (in order):  "kass-signal", "lpf-fft", "decimate-signal", "gaussian-noise", "digitizer".  The "gaussian-noise" generator can be omitted as needed.


Check or configure the following Kassiopeia fields in Project8Phase1_WithRoot_Template.xml:

- Select a kinematic generator (gen_krypton, gen_uniform, etc.) near the top of the xml file.
- Check that the path to the *Geometry*.xml file is correct, also near the top of the xml file.
- Trap currents are specified near the top of the xml file.
- Set the range of starting positions and pitch angles in the generator field.
- "max_time" (max time per electron)
- The number of electrons that will be generated is defined near the bottom of the xml file, in the "Simulation" field.


Check these fields in LocustPhase1Template.json:
- The generator we are using for this application is "kass-signal".  Paths to I/O files for this generator need to be edited.
- If Gaussian noise is to be included, then the "gaussian-noise" generator should be listed in the "generators" field.
- "n-records" is usually 1.  We have not been actively testing multi-record runs.
- "record-size" can be shortened to fit around a single electron track, or can be as long as 4194304.
- In the "digitizer" field, check the voltage range and offset.  They should be optimized for your signal + noise power.

There are two tuning distances that are hard-coded into the Phase 1 simulation:  CENTER_TO_SHORT and CENTER_TO_ANTENNA.  They are global in scope and are defined in LMCCyclotronRadiationExtractor.cc .

## Phase 2
Run the	simulation like this:
```
/path/to/LocustSim config=~/locust_mc/Config/Tutorial/LocustPhase2Template.json
```
The details in the Phase 1 section apply similarly to Phase 2.  The two tuning distances CENTER_TO_SHORT and CENTER_TO_ANTENNA have definitions unique to Phase 2, also defined in LMCCyclotronRadiationExtractor.cc .

## Phase 3
Run the	simulation like this:
```
/path/to/LocustSim config=~/locust_mc/Config/Tutorial/LocustPhase3Template.json
```

In the LocustPhase3Template.json file, the "generators" field should contain these generators (in order):  "patch-signal", "lpf-fft", "decimate-signal", "gaussian-noise", "digitizer".  The "gaussian-noise" generator can be omitted as needed.  Configuration of the Kassiopeia xml file is similar to Phases 1-2, except that a path to a Phase 3 Geometry xml file should be specified near the top of the xml file.

The patch antenna array is configured from the LocustPhase3Template.json file.  The patch array populates itself based on the parameters in the json file.  The fields to check are:
- The generator we are using for this application is "patch-signal".  Paths to I/O files for this generator need to be edited.
- nchannels in the "simulation" field is usually 30.
- "array-radius": 0.05,
- "npatches-per-strip": 21,
- "patch-spacing": 0.0054,
- "feed": "corporate",
- If Gaussian noise is to be included, then the "gaussian-noise" generator should be listed in the "generators" field.
- "n-records" is usually 1.  We have not been actively testing multi-record runs.
- "record-size" can be shortened to fit around a single electron track, or can be as long as 4194304.
- In the "digitizer" field, check the voltage range and offset.  They should be optimized for your signal + noise power.

## Phase 4
Run the	simulation like this:
```
/path/to/LocustSim config=~/locust_mc/Config/Tutorial/LocustPhase4Template.json
```
Configuration options are similar to those in Phase 3, except that for these fields in the LocustPhase4Template.json file:
- nchannels in the "simulation" field has been tested up to 240.
- "array-radius": 0.2.  
The file Project8Phase4_WithRoot_Template.xml presently references a large placeholder trap in Project8Phase4Geometry.xml .

## Fake Tracks
Run like this:
```
/path/to/LocustSim config=~/locust_mc/Config/Tutorial/LocustFakeTrack.json
```
There is a pdf to describe this generator here:  https://github.com/project8/locust_mc/blob/develop/Config/Tutorial/locust_faketrack_tutorial.pdf


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


## Plot a waterfall spectrogram
Process the simulated file of IQ voltages using Katydid and the config file katydid_locust.json :
    ```
         /path/to/Katydid config=/path/to/katydid_locust.json
    ```
This will result in 2 root files:  katydidwaterfall.root and basic.root .
Open the file plotoutput.C and edit any paths to point to the correct ones on your machine. 
Start root:
```
 >.L plotoutput.C
 > katydidwaterfall()
```

