include(MMseqsResourceCompiler)

set(COMPILED_RESOURCES
        easystructuresearch.sh
        structurecluster.sh
        structuresearch.sh
        structdatabases.sh
        mat3di.out
        3di_encoder_weights.kerasify
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(local-generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})
