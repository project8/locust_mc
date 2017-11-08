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
         LocustSim config=LocustTemplate.json
    ```
A VTK window should appear showing the particle track.

## Visualize the Field Map (requires VTK)
1. Edit LocustTemplate.json so that the "xml-filename" points to the "Project8Phase2_FieldMap.xml" file in the current directory.
2. Look at the parameters in the field "gen_bfieldlines":

    ```
            <r_set value_start="0.0" value_stop="0.07" value_count="1"/>
            <phi_set value_start="0.0" value_stop="360.0" value_count="11"/>
            <z_list add_value="-0.12"/>
    ```
3. Run the simulation with:

    ```
         LocustSim config=LocustTemplate.json
    ```
4. Look at the resulting field map in the VTK window. Now there will be a ROOT file with the magnetic field information.
5. Plot the field map in ROOT. In the root terminal:
```
 >.L plotoutput.C
 > fieldmap()

```


##Create a simulated Waterfall Plot with Locust+Katydid
1. Edit LocustTemplate.json so that "xml-filename" points to "Project8Phase2_WithRoot_Template.xml"
2. Change the max_time field in the xml file to 5e-4.
3. Comment the kswrite_vtk and vtk_window blocks that you have uncommented in the first exercise. Visualization becomes unwieldy for large simulations.
4. Run the simulation as before. The expected run time is about half an hour!
5. Open the file plotoutput.C and edit any paths to point to the correct ones on your machine. 
6. Start root:
```
 >.L plotoutput.C
 > katydidwaterfall()

```
