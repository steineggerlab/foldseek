add_library(3di
        structureto3di.cpp
        structureto3di.h
        structureto3diseqdist.h
        structureto3diseqdist.cpp
        )
mmseqs_setup_derived_target(3di)

target_include_directories(3di PRIVATE ..)  # needed for kerasify/keras_model.h
