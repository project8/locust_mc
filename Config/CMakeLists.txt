# CMakeLists.txt for locust_mc/Config
# Author: N.S. Oblath

configure_file( Project8_JustKass.xml.in Project8_JustKass.xml )
configure_file( Project8Phase1_WithRoot_Template.xml.in Project8Phase1_WithRoot_Template.xml )
configure_file( Project8Phase2_WithRoot_Template.xml.in Project8Phase2_WithRoot_Template.xml )
configure_file( Project8Phase3_WithRoot_Template.xml.in Project8Phase3_WithRoot_Template.xml )


set (EXAMPLE_CONFIGFILES
    ${CMAKE_CURRENT_BINARY_DIR}/Project8_JustKass.xml
    ${CMAKE_CURRENT_BINARY_DIR}/Project8Phase1_WithRoot_Template.xml
    ${CMAKE_CURRENT_BINARY_DIR}/Project8Phase2_WithRoot_Template.xml
    ${CMAKE_CURRENT_BINARY_DIR}/Project8Phase3_WithRoot_Template.xml 
    Project8Phase1Geometry.xml
    Project8Phase2Geometry.xml
    Project8Phase3Geometry.xml
)

pbuilder_install_config (${EXAMPLE_CONFIGFILES})
