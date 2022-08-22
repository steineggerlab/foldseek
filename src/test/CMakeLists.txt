include(FoldseekSetupTest)

set(TESTS
	TestGenerateTree.cpp
        )


FOREACH (TEST ${TESTS})
    foldseek_setup_test(${TEST})
ENDFOREACH ()
