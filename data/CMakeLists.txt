include(MMseqsResourceCompiler)

set(COMPILED_RESOURCES
        easystructuresearch.sh
        structurecluster.sh
        structureindex.sh
        structuresearch.sh
        structureiterativesearch.sh
        structurerbh.sh
        structdatabases.sh
        mat3di.out
        encoder_weights_3di.kerasify
        evalue_nn.kerasify
        main.js
        vendor.js.zst
        complexsearch.sh
        easycomplexsearch.sh
        easycomplexcluster.sh
        complexcluster.sh
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(local-generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})
