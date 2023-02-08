# CMakeLists for locust_mc/Source/IO
# Author: P. Slocum

include_directories( BEFORE 
    ${CMAKE_CURRENT_SOURCE_DIR}
)    


if (ROOT_FOUND)

    set (IO_HEADERFILES
        LMCFileWriter.hh
        LMCRootGraphWriter.hh
        LMCRootHistoWriter.hh
        LMCRootTreeWriter.hh
    )   
         
    set (IO_SOURCEFILES
        LMCTrack.cc
        LMCEvent.cc
        LMCRunParameters.cc
        LMCFileWriter.cc
        LMCRootGraphWriter.cc
        LMCRootHistoWriter.cc
        LMCRootTreeWriter.cc
    )

    set (IO_DICT_HEADERFILES
        LMCTrack.hh
        LMCEvent.hh
        LMCRunParameters.hh
    )
    
    set (IO_HEADERFILES
        ${IO_HEADERFILES}
        ${IO_DICT_HEADERFILES}
    )

    set (IO_SOURCEFILES
        ${IO_SOURCEFILES}
    )

    set (IO_LINKDEF_HEADERFILE LinkDef/LMCLinkDef.hh)
    set (IO_DICT_OUTFILE ${CMAKE_CURRENT_BINARY_DIR}/IODict.cxx)
    set (IO_DICT_PCMFILE ${CMAKE_CURRENT_BINARY_DIR}/IODict_rdict.pcm)

endif (ROOT_FOUND)



##################################################

if (ROOT_FOUND)
    ROOT_GENERATE_DICTIONARY (IODict ${IO_DICT_HEADERFILES} LINKDEF ${IO_LINKDEF_HEADERFILE} OPTIONS -inlineInputHeader)
    pbuilder_install_files (${LIB_INSTALL_DIR} ${IO_DICT_PCMFILE})
endif (ROOT_FOUND)

add_library (LocustIO ${IO_SOURCEFILES} ${IO_DICT_OUTFILE})
target_link_libraries (LocustIO ${EXTERNAL_LIBRARIES})

pbuilder_install_headers (${IO_HEADERFILES})
pbuilder_install_libraries (LocustIO)


