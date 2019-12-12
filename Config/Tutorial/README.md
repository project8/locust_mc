# Running Locust
----------------

## Phase 1

Run the	simulation like this:
```
/path/to/LocustSim config=~/locust_mc/Config/Tutorial/LocustPhase1Template.json
```

In the LocustPhase1Template.json file, the "generators" field should contain these generators (in order):  "kass-signal", "lpf-fft", "decimate-signal", "gaussian-noise", "digitizer".  The "gaussian-noise" generator can be omitted as needed.  For Locust 1.14 or greater, each generator listed in the "generators" field should have a code block to define its parameters, regardless of whether there are any parameters.  For example, generators "lpf-fft" and "decimate-signal" do not have parameters, but a section of the config file should look like this:
```
    "generators":
    [          
       "fake-track",
       "lpf-fft",
       "decimate-signal",
       "gaussian-noise",
       "digitizer"
    ],

    "lpf-fft":
    {

    },

    "decimate-signal":
    {

    },

```


Check or configure the following Kassiopeia fields in Project8Phase1_WithRoot_Template.xml:

- <cycl_rad_extr name="my_rad_extr" P8Phase="1" /> 
- Select a kinematic generator (gen_krypton, gen_uniform, etc.) near the top of the xml file.
- Check that the path to the *Geometry*.xml file is correct, also near the top of the xml file.
- Trap currents are specified near the top of the xml file.
- Set the range of starting positions and pitch angles in the generator field.
- "max_time" (max time per electron)
- The number of electrons that will be generated in Kassiopeia is defined near the bottom of the xml file, in the "Simulation:events" field.  The electrons are generated sequentially in time and are spaced by a parameter presently defined inside the Locust generators.


Check these fields in LocustPhase1Template.json:
- The generator we are using for this application is "kass-signal".  Paths to I/O files for this generator need to be edited.
- If Gaussian noise is to be included, then the "gaussian-noise" generator should be listed in the "generators" field.
- "n-records" is usually 1.  We have not been actively testing multi-record runs.
- "record-size" can be shortened to fit around a single electron track, or can be as long as 4194304.
- In the "digitizer" field, check the voltage range and offset.  They should be optimized for your signal + noise power.

There are two tuning distances that are hard-coded into the Phase 1 simulation:  CENTER_TO_SHORT and CENTER_TO_ANTENNA.  They are global in scope and are defined in LMCCyclotronRadiationExtractor.cc .  A git gist describing details of the field calculations from Section 3.2.2 of the Locust paper is here: https://gist.github.com/pslocum/3ee7b440a947974caa1b881cd89c5071

## Phase 2
Run the	simulation like this:
```
/path/to/LocustSim config=~/locust_mc/Config/Tutorial/LocustPhase2Template.json
```
The details in the Phase 1 section apply similarly to Phase 2.  The two tuning distances CENTER_TO_SHORT and CENTER_TO_ANTENNA have definitions unique to Phase 2, also defined in LMCCyclotronRadiationExtractor.cc .  Make sure that you have specified Phase 2 in the file Project8Phase2_WithRoot_Template.xml, like this:
```
<cycl_rad_extr name="my_rad_extr" P8Phase="2" /> 
```

## Phase 3
Run the	simulation like this:
```
/path/to/LocustSim config=~/locust_mc/Config/Tutorial/LocustPhase3Template.json
```

In the LocustPhase3Template.json file, the "generators" field should contain these generators (in order):  "patch-signal", "lpf-fft", "decimate-signal", "gaussian-noise", "digitizer".  The "gaussian-noise" generator can be omitted as needed.  Configuration of the Kassiopeia xml file is similar to Phases 1-2, except that a path to a Phase 3 Geometry xml file should be specified near the top of the xml file.  Make sure that you have specified Phase 3 in the file Project8Phase3_WithRoot_Template.xml, like this:
```
<cycl_rad_extr name="my_rad_extr" P8Phase="3" /> 
```


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

The file Project8Phase4_WithRoot_Template.xml presently references a large placeholder trap in Project8Phase4Geometry.xml .  Make sure that you have specified Phase 4 in the file Project8Phase4_WithRoot_Template.xml, like this:
```
<cycl_rad_extr name="my_rad_extr" P8Phase="4" /> 
```
## Transfer Functions

Config file options:

- "tf-transmitter-filename": "../Data/TransferFunctions/CoupledeDipoleTF.txt", [required]
- "tf-transmitter-bin-width": 0.01e9, [default = 10 MHz, needs to match bin width in tf-transmitter-filename]
- "tf-transmitter-nskips": 1, [default = 1]
- "window-function-type": 1, [default = 1 which is a Tukey window]
- "zero-padding-size": 97199, [default = 100000], use default unless debugging something specific.
- "shift-n-bins": 5000, [default 5000]
- "fir-receiver-filename": "../out/Dipole/firReceiverFile.txt", [required]
- "fir-receiver-dt": 1.0e-12, [default 1.e-12, needs to be specified to match fir-receiver-dt time resolution]
- "fir-receiver-nskips": 1, [default = 1]



## Fake Events/Tracks

There are two ways to generate fake events/tracks: 1. One event/track at a time using Locust directly or 2. Multiple events/tracks using the LocustFakeEvent.py script (available in project8 scripts repo under scripts/FakeEventGenerator). There is a tutorial pdf to describe these options in detail here: https://github.com/project8/scripts/blob/master/FakeEventGenerator/locust_fakeevent_tutorial.pdf

Note that the fake track parameters are defined in the LMCFakeTrackSignalGenerator class through probability distribution functions (PDF) whose parameters may be set through the LocustFakeTrack.json file. The following parameters use configurable PDFs:
- Start frequency: Uniform PDF
    - "start-frequency-min" in Hz
    - "start-frequency-max" in Hz
- Track slope: Gaussian PDF
    - "slope-mean" in MHz/ms
    - "slope-std" in MHz/ms
- Track start sime: Uniform PDF
    - "start-time-min" in s
    - "start-time-max" in s
- Track length: Exponential PDF
    - "track-length-mean" in s
- Number of tracks in event: Exponential PDF
    - "ntracks-mean" unitless double
- Random generator seed
    - "random-seed" unitless integer

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

## Visualize the Field Map

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
4. Now there will be a ROOT file with the magnetic field information.
5. Plot the field map in ROOT. In the root terminal:
```
 >.L plotoutput.C
 > fieldmap()
```


## Plot a waterfall spectrogram

Process the simulated egg file using Katydid and the config file katydid_locust.json :
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

# Instructions for running on Yale cluster

Log in to the grace cluster either like this:
```ssh netID@grace.hpc.yale.edu```
or, if you want to be able to pull up Xwindows for plotting, like this:
```ssh -Y netID@grace.hpc.yale.edu```
From a Mac OS laptop, you may have to install XQuarts to pull up XWindows.  Test your XWindow by logging into Grace and then typing ```xclock```.  You should see a clock.

Copy the setup script into your home directory on Grace:
```
cp ~ps48/setup_script .
```

It should only need to be run one time, like this:
```
./setup_script
```
After you have run the setup_script, the files in your directories are yours to edit.  The next step will have to be run every time you log into a session on Grace, so you may want to put this command into your .bashrc file:
```
module load Tools/project8
```

To run a Locust simulation as a batch job with the Slurm queue https://research.computing.yale.edu/support/hpc/user-guide/slurm , first cd into the directory ~/project8/managePhaseN where 1<N<4.  For Phase 1,
```
cd ~/project8/managePhase1
```
Edit the file SimulateSeed to set the range of Monte Carlo seeds to be run.  Each seed will have a separate parallel batch job and will take 3 hours.  30 seeds at a time has been manageable; 1 seed is a good place to start.  Start the job like this:
```
./SimulateSeed
```
Check your job status like this:
```
squeue -u netID
```
Also, while the job with e.g. Seed = 55 is running, the file ~/project8/managePhase1/locust_jobSeed55 should be growing as it collects the terminal output.  When the job is finished, check that the egg file has been written in ~/data/Simulation/Phase1 directory, as
```
ls -l ~/data/Simulation/Phase1/*.egg
```
Next, process the egg files with Katydid using the script ProcessEggFiles (edit it to check the seed range), either interactively as
```
./ProcessEggFiles
```
or in batch mode as
```sbatch ProcessEggFilesBatch```
There should now be processed root tree files in your directory ~/data/Simulation/Phase1.  You will have to delete or remove the raw egg files and any root spectrogram files, as they are bulky and we do not yet have space allocated for data storage on Grace.  

There are data processing macros in ~/project8/managePhase(N) that can be used to generate rough figures for inspection.  For Phases 1-2, first get onto an interactive compute node (with x11 forwarding include --x11) like this:
```
srun --pty --x11 -p interactive -c4 bash
```
Then start root and run the macro: ```
```
root -l
.L PlotSpectrum.c
PlotKrypton()
```
If you would rather look at the processed data files locally on your laptop, then use scp to download them
```
scp netID@grace.hpc.yale.edu:~/data/Simulation/PhaseN/*filename*.root .
```
