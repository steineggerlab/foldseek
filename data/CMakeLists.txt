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
<<<<<<< HEAD
        complexsearch.sh
        easycomplexsearch.sh
        complexcluster.sh
        easycomplexcluster.sh
=======
        multimersearch.sh
        easymultimersearch.sh
>>>>>>> 25812ffa585248b146fb0217b981b507dc92e851
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(local-generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})
