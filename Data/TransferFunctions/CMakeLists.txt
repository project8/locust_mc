# CMakeLists.txt for Locust_mc/Data/TransferFunctions

set( LOCUST_MC_DATA_FILES
    PatchTFLocust.txt
    UncoupledHalfWaveeDipoleTF.txt
    CoupledeDipoleTF.txt                 
    PatchTFLocust_wide.txt
    UncoupledMagDipoleTF.txt
    PatchTFLocust_cropped.txt            
    TF_Template.txt
)

pbuilder_install_data( ${LOCUST_MC_DATA_FILES} )
