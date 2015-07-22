#include(UseLATEX)

#ADD_LATEX_DOCUMENT(Manual_v2.0.tex
#    IMAGE_DIRS Figures
#    FORCE_PDF
#)

#ADD_CUSTOM_COMMAND(TARGET Manual_v2.0 POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#        "${CMAKE_CURRENT_BINARY_DIR}/Manual_v2.0.pdf"
#        "${CMAKE_BINARY_DIR}/bin")
