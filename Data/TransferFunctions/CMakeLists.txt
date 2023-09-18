# CMakeLists.txt for Locust_mc/Data/TransferFunctions

set( LOCUST_MC_DATA_FILES
    PatchTFLocust.txt
    UncoupledHalfWaveeDipoleTF.txt
    CoupledeDipoleTF.txt                 
    PatchTFLocust_wide.txt
    UncoupledMagDipoleTF.txt
    PatchTFLocust_cropped.txt
    cavityTF.txt            
    WEGA_Impedance_Center.txt
    CCA_TE011_sParams_V01.txt
    TF_Template.txt
)

pbuilder_install_data( ${LOCUST_MC_DATA_FILES} )
