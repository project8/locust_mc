# CMakeLists.txt for locust_mc/Config
# Author: N.S. Oblath

configure_file( Project8_JustKass.xml.in Project8_JustKass.xml )
configure_file( Project8_JustKassFieldMap.xml.in Project8_JustKassFieldMap.xml )
configure_file( Project8Phase1_WithRoot_Template.xml.in Project8Phase1_WithRoot_Template.xml )
configure_file( Project8Phase2_WithRoot_Template.xml.in Project8Phase2_WithRoot_Template.xml )
configure_file( Project8Phase3_WithRoot_Template.xml.in Project8Phase3_WithRoot_Template.xml )
configure_file( LocustPhase1Template.json.in LocustPhase1Template.json )
configure_file( LocustPhase2Template.json.in LocustPhase2Template.json )
configure_file( LocustPhase3Template.json.in LocustPhase3Template.json )
configure_file( LocustTestSignal.json.in LocustTestSignal.json )
configure_file( katydid.json.in katydid.json )


set (EXAMPLE_CONFIGFILES
    ${CMAKE_CURRENT_BINARY_DIR}/Project8_JustKass.xml
    ${CMAKE_CURRENT_BINARY_DIR}/Project8_JustKassFieldMap.xml
    ${CMAKE_CURRENT_BINARY_DIR}/Project8Phase1_WithRoot_Template.xml
    ${CMAKE_CURRENT_BINARY_DIR}/Project8Phase2_WithRoot_Template.xml
    ${CMAKE_CURRENT_BINARY_DIR}/Project8Phase3_WithRoot_Template.xml 
    Project8Phase1Geometry.xml
    Project8Phase2Geometry.xml
    Project8Phase3Geometry.xml
    ${CMAKE_CURRENT_BINARY_DIR}/LocustPhase1Template.json    
    ${CMAKE_CURRENT_BINARY_DIR}/LocustPhase2Template.json    
    ${CMAKE_CURRENT_BINARY_DIR}/LocustPhase3Template.json    
    ${CMAKE_CURRENT_BINARY_DIR}/LocustTestSignal.json    
    ${CMAKE_CURRENT_BINARY_DIR}/katydid.json    
)

pbuilder_install_config (${EXAMPLE_CONFIGFILES})
