add_library(3di
        structureto3di.cpp
        structureto3di.h
        structureto3diseqdist.h
        structureto3diseqdist.cpp
        structureto12st.cpp
        structureto12st.h
        )
mmseqs_setup_derived_target(3di)
target_include_directories(3di
        PRIVATE 
        .. # kerasify/keras_model.h
        ${PROJECT_BINARY_DIR}/generated # encoder_weights_3di.kerasify.h
)
add_dependencies(3di local-generated)
