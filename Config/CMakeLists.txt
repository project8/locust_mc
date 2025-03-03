# CMakeLists.txt for locust_mc/Config
# Author: N.S. Oblath

configure_file( JustKass.xml.in JustKass.xml )
configure_file( JustKassFieldMap.xml.in JustKassFieldMap.xml )
configure_file( LocustKass_FreeSpace_Template.xml.in LocustKass_FreeSpace_Template.xml )
configure_file( LocustKass_Cavity_CCA.xml.in LocustKass_Cavity_CCA.xml )
configure_file( LocustKass_Cavity_CCA_Pileup.xml.in LocustKass_Cavity_CCA_Pileup.xml )
configure_file( LocustKass_Cavity_1GHz.xml.in LocustKass_Cavity_1GHz.xml )
configure_file( LocustWaveguideTemplate.json.in LocustWaveguideTemplate.json )
configure_file( LocustCavityCCA.json.in LocustCavityCCA.json )
configure_file( LocustCavityCCA_Pileup.json.in LocustCavityCCA_Pileup.json )
configure_file( LocustCavityCCA_ModeMapTest.json.in LocustCavityCCA_ModeMapTest.json )
configure_file( LocustCavityRectangular.json.in LocustCavityRectangular.json )
configure_file( LocustCavity1GHz.json.in LocustCavity1GHz.json )
configure_file( LocustSignalPlusNoise.json.in LocustSignalPlusNoise.json )
configure_file( LocustFreeSpaceTemplate.json.in LocustFreeSpaceTemplate.json )
configure_file( LocustPlaneWaveTemplate.json.in LocustPlaneWaveTemplate.json )
configure_file( LocustFakeTrack.json.in LocustFakeTrack.json )
configure_file( LocustTestSignal.json.in LocustTestSignal.json )
configure_file( plotfieldmap.c.in plotfieldmap.c )
configure_file( katydid.json.in katydid.json )
configure_file( ChannelAggregatorConfig.yaml.in ChannelAggregatorConfig.yaml )
configure_file( LocustBlankTemplate.json.in LocustBlankTemplate.json )


set (EXAMPLE_CONFIGFILES
    ${CMAKE_CURRENT_BINARY_DIR}/JustKass.xml
    ${CMAKE_CURRENT_BINARY_DIR}/JustKassFieldMap.xml
    ${CMAKE_CURRENT_BINARY_DIR}/LocustKass_FreeSpace_Template.xml 
    ${CMAKE_CURRENT_BINARY_DIR}/LocustKass_Cavity_CCA.xml
    ${CMAKE_CURRENT_BINARY_DIR}/LocustKass_Cavity_CCA_Pileup.xml
    ${CMAKE_CURRENT_BINARY_DIR}/LocustKass_Cavity_1GHz.xml
    FreeSpaceGeometry.xml
	CavityGeometry_VTinyCoil.xml
	CavityGeometry_V00_00_00.xml
    ${CMAKE_CURRENT_BINARY_DIR}/LocustWaveguideTemplate.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustCavityCCA.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustCavityCCA_ModeMapTest.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustCavityCCA_Pileup.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustCavityRectangular.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustCavity1GHz.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustFreeSpaceTemplate.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustPlaneWaveTemplate.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustFakeTrack.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustTestSignal.json
    ${CMAKE_CURRENT_BINARY_DIR}/LocustSignalPlusNoise.json    
    ${CMAKE_CURRENT_BINARY_DIR}/katydid.json
    ${CMAKE_CURRENT_BINARY_DIR}/ChannelAggregatorConfig.yaml
    ${CMAKE_CURRENT_BINARY_DIR}/LocustBlankTemplate.json
    ${CMAKE_CURRENT_BINARY_DIR}/plotfieldmap.c
)

pbuilder_install_config (${EXAMPLE_CONFIGFILES})
