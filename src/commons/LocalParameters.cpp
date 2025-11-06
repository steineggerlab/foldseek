#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include "mat3di.out.h"

const int LocalParameters::DBTYPE_CA_ALPHA = 101;
const int LocalParameters::DBTYPE_TMSCORE = 102;

LocalParameters::LocalParameters() :
        Parameters(),
        PARAM_TMSCORE_THRESHOLD(PARAM_TMSCORE_THRESHOLD_ID,"--tmscore-threshold", "TMscore threshold", "accept alignments with a tmsore > thr [0.0,1.0]",typeid(float), (void *) &tmScoreThr, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_TMSCORE_THRESHOLD_MODE(PARAM_TMSCORE_THRESHOLD_MODE_ID,"--tmscore-threshold-mode", "TMscore threshold mode", "0: alignment, 1: query 2: target length",typeid(int), (void *) &tmScoreThrMode, "^[0-2]{1}$"),
        PARAM_TMALIGN_HIT_ORDER(PARAM_TMALIGN_HIT_ORDER_ID,"--tmalign-hit-order", "TMalign hit order", "order hits by 0: (qTM+tTM)/2, 1: qTM, 2: tTM, 3: min(qTM,tTM) 4: max(qTM,tTM)",typeid(int), (void *) &tmAlignHitOrder, "^[0-4]{1}$"),
        PARAM_LDDT_THRESHOLD(PARAM_LDDT_THRESHOLD_ID,"--lddt-threshold", "LDDT threshold", "accept alignments with a lddt > thr [0.0,1.0]",typeid(float), (void *) &lddtThr, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_SORT_BY_STRUCTURE_BITS(PARAM_SORT_BY_STRUCTURE_BITS_ID,"--sort-by-structure-bits", "Sort by structure bit score", "sort by bits*sqrt(alnlddt*alntmscore)",typeid(int), (void *) &sortByStructureBits, "^[0-1]{1}$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MASK_BFACTOR_THRESHOLD(PARAM_MASK_BFACTOR_THRESHOLD_ID,"--mask-bfactor-threshold", "Mask b-factor threshold", "mask residues for seeding if b-factor < thr [0,100]",typeid(float), (void *) &maskBfactorThreshold, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_ALIGNMENT_TYPE(PARAM_ALIGNMENT_TYPE_ID,"--alignment-type", "Alignment type", "How to compute the alignment:\n0: 3di alignment\n1: TM alignment\n2: 3Di+AA\n3: LoL alignmnet",typeid(int), (void *) &alignmentType, "^[0-3]{1}$"),
        PARAM_CHAIN_NAME_MODE(PARAM_CHAIN_NAME_MODE_ID,"--chain-name-mode", "Chain name mode", "Add chain to name:\n0: auto\n1: always add\n",typeid(int), (void *) &chainNameMode, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
        PARAM_MODEL_NAME_MODE(PARAM_MODEL_NAME_MODE_ID,"--model-name-mode", "Model name mode", "Add model to name:\n0: auto\n1: always add\n",typeid(int), (void *) &modelNameMode, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
        PARAM_WRITE_MAPPING(PARAM_WRITE_MAPPING_ID, "--write-mapping", "Write mapping file", "write _mapping file containing mapping from internal id to taxonomic identifier", typeid(int), (void *) &writeMapping, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
        PARAM_WRITE_FOLDCOMP(PARAM_WRITE_FOLDCOMP_ID, "--write-foldcomp", "Write Foldcomp", "write _fcz Foldcomp database", typeid(int), (void *) &writeFoldcomp, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
        PARAM_TMALIGN_FAST(PARAM_TMALIGN_FAST_ID,"--tmalign-fast", "TMalign fast","turn on fast search in TM-align" ,typeid(int), (void *) &tmAlignFast, "^[0-1]{1}$"),
        PARAM_EXACT_TMSCORE(PARAM_EXACT_TMSCORE_ID,"--exact-tmscore", "Exact TMscore","turn on fast exact TMscore (slow), default is approximate" ,typeid(int), (void *) &exactTMscore, "^[0-1]{1}$"),
        PARAM_N_SAMPLE(PARAM_N_SAMPLE_ID, "--n-sample", "Sample size","pick N random sample" ,typeid(int), (void *) &nsample, "^[0-9]{1}[0-9]*$"),
        PARAM_COORD_STORE_MODE(PARAM_COORD_STORE_MODE_ID, "--coord-store-mode", "Coord store mode", "Coordinate storage mode: \n1: C-alpha as float\n2: C-alpha as difference (uint16_t)", typeid(int), (void *) &coordStoreMode, "^[1-2]{1}$",MMseqsParameter::COMMAND_EXPERT),
        PARAM_MIN_ASSIGNED_CHAINS_THRESHOLD(PARAM_MIN_ASSIGNED_CHAINS_THRESHOLD_ID, "--min-assigned-chains-ratio", "Minimum assigned chains percentage Threshold", "Minimum ratio of assigned chains out of all query chains > thr [0.0,1.0]", typeid(float), (void *) & minAssignedChainsThreshold, "^[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_MONOMER_INCLUDE_MODE(PARAM_MONOMER_INCLUDE_MODE_ID, "--monomer-include-mode", "Monomer inclusion Mode for MultimerSerch", "Monomer Complex Inclusion 0: include monomers, 1: NOT include monomers", typeid(int), (void *) & monomerIncludeMode, "^[0-1]{1}$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_CLUSTER_SEARCH(PARAM_CLUSTER_SEARCH_ID, "--cluster-search", "Cluster search", "first find representative then align all cluster members", typeid(int), (void *) &clusterSearch, "^[0-1]{1}$",MMseqsParameter::COMMAND_MISC),
        PARAM_FILE_INCLUDE(PARAM_FILE_INCLUDE_ID, "--file-include", "File Inclusion Regex", "Include file names based on this regex", typeid(std::string), (void *) &fileInclude, "^.*$"),
        PARAM_FILE_EXCLUDE(PARAM_FILE_EXCLUDE_ID, "--file-exclude", "File Exclusion Regex", "Exclude file names based on this regex", typeid(std::string), (void *) &fileExclude, "^.*$"),
        PARAM_INDEX_EXCLUDE(PARAM_INDEX_EXCLUDE_ID, "--index-exclude", "Index Exclusion", "Exclude parts of the index:\n0: Full index\n1: Exclude k-mer index (for use with --prefilter-mode 1)\n2: Exclude C-alpha coordinates (for use with --sort-by-structure-bits 0)\nFlags can be combined bit wise", typeid(int), (void *) &indexExclude, "^[0-3]{1}$", MMseqsParameter::COMMAND_EXPERT),
        PARAM_MULTIMER_REPORT_MODE(PARAM_MULTIMER_REPORT_MODE_ID, "--multimer-report-mode", "Complex report mode", "Complex report mode:\n0: No report\n1: Write complex report", typeid(int), (void *) &multimerReportMode, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
        PARAM_MULTIMER_REPORT_MODE_BC_COMPAT(PARAM_MULTIMER_REPORT_MODE_BC_COMPAT_ID, "--complex-report-mode", "", "", typeid(int), (void *) &multimerReportMode, "^[0-1]{1}$", MMseqsParameter::COMMAND_HIDDEN),
        PARAM_EXPAND_MULTIMER_EVALUE(PARAM_EXPAND_MULTIMER_EVALUE_ID, "--expand-multimer-evalue", "Multimer E-value", "E-value threshold for multimer chain expansion (range 0.0-inf)", typeid(double), (void *) &eValueThrExpandMultimer, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_ALIGN),
        PARAM_EXPAND_MULTIMER_EVALUE_BC_COMPAT(PARAM_EXPAND_MULTIMER_EVALUE_BC_COMPAT_ID, "--expand-complex-evalue", "", "", typeid(double), (void *) &eValueThrExpandMultimer, "^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)|[0-9]*(\\.[0-9]+)?$", MMseqsParameter::COMMAND_HIDDEN),
        PARAM_INPUT_FORMAT(PARAM_INPUT_FORMAT_ID, "--input-format", "Input format", "Format of input structures:\n0: Auto-detect by extension\n1: PDB\n2: mmCIF\n3: mmJSON\n4: ChemComp\n5: Foldcomp", typeid(int), (void *) &inputFormat, "^[0-5]{1}$"),
        PARAM_PDB_OUTPUT_MODE(PARAM_PDB_OUTPUT_MODE_ID, "--pdb-output-mode", "PDB output mode", "PDB output mode:\n0: Single multi-model PDB file\n1: One PDB file per chain\n2: One PDB file per complex", typeid(int), (void *) &pdbOutputMode, "^[0-2]{1}$", MMseqsParameter::COMMAND_MISC),
        PARAM_PROSTT5_MODEL(PARAM_PROSTT5_MODEL_ID, "--prostt5-model", "Path to ProstT5", "Path to ProstT5 model", typeid(std::string), (void *) &prostt5Model, "^.*$", MMseqsParameter::COMMAND_COMMON),
        PARAM_DB_EXTRACTION_MODE(PARAM_DB_EXTRACTION_MODE_ID, "--db-extraction-mode", "Createdb extraction mode", "createdb extraction mode: 0: chain 1: interface", typeid(int), (void *) &dbExtractionMode, "^[0-1]{1}$"),
        PARAM_DISTANCE_THRESHOLD(PARAM_DISTANCE_THRESHOLD_ID, "--distance-threshold", "Interface distance threshold", "Residues with C-beta below this threshold will be part of interface", typeid(float), (void *) &distanceThreshold, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_CHAIN_TM_THRESHOLD(PARAM_CHAIN_TM_THRESHOLD_ID,"--chain-tm-threshold", "chain TMscore threshold", "accept alignments with a minimum chain tmsore > thr [0.0,1.0]",typeid(float), (void *) &filtChainTmThr, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_INTERFACE_LDDT_THRESHOLD(PARAM_INTERFACE_LDDT_THRESHOLD_ID,"--interface-lddt-threshold", "Interface LDDT threshold", "accept alignments with a lddt > thr [0.0,1.0]",typeid(float), (void *) &filtInterfaceLddtThr, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_MIN_ALIGNED_CHAINS(PARAM_MIN_ALIGNED_CHAINS_ID, "--min-aligned-chains", "Minimum threshold of aligned chains","save alignments with at least n chain aligned between query and target" ,typeid(int), (void *) &minAlignedChains, "^[0-9]{1}[0-9]*$"),
        PARAM_MULTIDOMAIN(PARAM_MULTIDOMAIN_ID, "--lolalign-multidomain", "MultiDomain Mode", "MultiDomain Mode LoLalign", typeid(int), (void *) &multiDomain, "^[0-1]{1}$")
        {
    PARAM_ALIGNMENT_MODE.description = "How to compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also start_pos and cov\n3: also seq.id";
    PARAM_ALIGNMENT_MODE.regex = "^[0-3]{1}$";
    PARAM_ALIGNMENT_MODE.category = MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT;
    PARAM_NUM_ITERATIONS.description = "Number of iterative profile search iterations (0: auto (select optimal), 1: default, 1-n), N≥2 exactly N)";
    PARAM_NUM_ITERATIONS.regex = "^[0-9]{1}[0-9]*$";
    PARAM_EXHAUSTIVE_SEARCH.description = "Turns on an exhaustive all vs all search by by passing the prefilter step";
    PARAM_EXHAUSTIVE_SEARCH.category = MMseqsParameter::COMMAND_PREFILTER;
    PARAM_MIN_ALN_LEN.category = MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT;
    PARAM_ALIGNMENT_OUTPUT_MODE.category = MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT;
    PARAM_TAR_EXCLUDE.category = MMseqsParameter::COMMAND_EXPERT;
    PARAM_TAR_INCLUDE.category = MMseqsParameter::COMMAND_EXPERT;
    PARAM_TAXON_LIST.category = MMseqsParameter::COMMAND_EXPERT;
    PARAM_ZDROP.category = MMseqsParameter::COMMAND_HIDDEN;

    PARAM_FORMAT_MODE.description = "Output format:\n0: BLAST-TAB\n"
                                    "1: SAM\n2: BLAST-TAB + query/db length\n"
                                    "3: Pretty HTML\n4: BLAST-TAB + column headers\n"
                                    "5: Calpha only PDB super-posed to query\n"
                                    "BLAST-TAB (0) and BLAST-TAB + column headers (4)"
                                    "support custom output formats (--format-output)\n"
                                    "(5) Superposed PDB files (Calpha only)";

    // TODO
    PARAM_FORMAT_MODE.regex = "^[0-5]{1}$";
    PARAM_SEARCH_TYPE.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_TRANSLATION_TABLE.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_TRANSLATION_TABLE.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_PCA.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_PCB.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_ZDROP.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_CORR_SCORE_WEIGHT.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_REALIGN.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_REALIGN_SCORE_BIAS.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_REALIGN_MAX_SEQS.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_SCORE_BIAS.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_WRAPPED_SCORING.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_ALPH_SIZE.category = MMseqsParameter::COMMAND_HIDDEN;
    PARAM_INCLUDE_IDENTITY.category = MMseqsParameter::COMMAND_HIDDEN;

    scoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>("3di.out", "3di.out"));
    seedScoringMatrixFile = MultiParam<NuclAA<std::string>>(NuclAA<std::string>("3di.out", "3di.out"));
    substitutionMatrices.emplace_back("3di.out", mat3di_out, mat3di_out_len);

    // structurecreatedb
    structurecreatedb.push_back(&PARAM_GPU);
    structurecreatedb.push_back(&PARAM_PROSTT5_MODEL);
    structurecreatedb.push_back(&PARAM_CHAIN_NAME_MODE);
    structurecreatedb.push_back(&PARAM_MODEL_NAME_MODE);
    structurecreatedb.push_back(&PARAM_DB_EXTRACTION_MODE);
    structurecreatedb.push_back(&PARAM_DISTANCE_THRESHOLD);
    structurecreatedb.push_back(&PARAM_WRITE_MAPPING);
    structurecreatedb.push_back(&PARAM_WRITE_FOLDCOMP);
    structurecreatedb.push_back(&PARAM_MASK_BFACTOR_THRESHOLD);
    structurecreatedb.push_back(&PARAM_COORD_STORE_MODE);
    structurecreatedb.push_back(&PARAM_WRITE_LOOKUP);
    structurecreatedb.push_back(&PARAM_INPUT_FORMAT);
    // protein chain only
    structurecreatedb.push_back(&PARAM_FILE_INCLUDE);
    structurecreatedb.push_back(&PARAM_FILE_EXCLUDE);
    structurecreatedb.push_back(&PARAM_THREADS);
    structurecreatedb.push_back(&PARAM_V);

    convertalignments.push_back(&PARAM_EXACT_TMSCORE);

    createindex.push_back(&PARAM_INDEX_EXCLUDE);

    // tmalign
    tmalign.push_back(&PARAM_MIN_SEQ_ID);
    tmalign.push_back(&PARAM_C);
    tmalign.push_back(&PARAM_COV_MODE);
    tmalign.push_back(&PARAM_MAX_REJECTED);
    tmalign.push_back(&PARAM_MAX_ACCEPT);
    tmalign.push_back(&PARAM_ADD_BACKTRACE);
    tmalign.push_back(&PARAM_INCLUDE_IDENTITY);
    tmalign.push_back(&PARAM_TMSCORE_THRESHOLD);
    tmalign.push_back(&PARAM_TMSCORE_THRESHOLD_MODE);
    tmalign.push_back(&PARAM_TMALIGN_HIT_ORDER);
    tmalign.push_back(&PARAM_TMALIGN_FAST);
    tmalign.push_back(&PARAM_PRELOAD_MODE);
    tmalign.push_back(&PARAM_THREADS);
    tmalign.push_back(&PARAM_V);
//    tmalign.push_back(&PARAM_GAP_OPEN);
//    tmalign.push_back(&PARAM_GAP_EXTEND);

    //LoLalign
    lolalign.push_back(&PARAM_MULTIDOMAIN);
    lolalign.push_back(&PARAM_MIN_SEQ_ID);
    lolalign.push_back(&PARAM_PRELOAD_MODE);
    lolalign.push_back(&PARAM_MAX_REJECTED);
    lolalign.push_back(&PARAM_MAX_ACCEPT);
    lolalign.push_back(&PARAM_C);
    lolalign.push_back(&PARAM_COV_MODE);
    lolalign.push_back(&PARAM_ADD_BACKTRACE);
    lolalign.push_back(&PARAM_THREADS);
    lolalign.push_back(&PARAM_V);



    structurerescorediagonal.push_back(&PARAM_EXACT_TMSCORE);
    structurerescorediagonal.push_back(&PARAM_TMSCORE_THRESHOLD);
    structurerescorediagonal.push_back(&PARAM_TMSCORE_THRESHOLD_MODE);
    structurerescorediagonal.push_back(&PARAM_LDDT_THRESHOLD);
    structurerescorediagonal.push_back(&PARAM_ALIGNMENT_TYPE);
    structurerescorediagonal = combineList(structurerescorediagonal, align);

    structurealign.push_back(&PARAM_TMSCORE_THRESHOLD);
    structurealign.push_back(&PARAM_TMSCORE_THRESHOLD_MODE);
    structurealign.push_back(&PARAM_LDDT_THRESHOLD);
    structurealign.push_back(&PARAM_SORT_BY_STRUCTURE_BITS);
    structurealign.push_back(&PARAM_ALIGNMENT_TYPE);
    structurealign.push_back(&PARAM_EXACT_TMSCORE);
    structurealign = combineList(structurealign, align);

    // strucclust
    strucclust = combineList(clust, structurealign);
    strucclust = combineList(strucclust, structurerescorediagonal);
    strucclust = combineList(strucclust, kmermatcher);
    strucclust.push_back(&PARAM_REMOVE_TMP_FILES);
    strucclust.push_back(&PARAM_RUNNER);

    databases.push_back(&PARAM_HELP);
    databases.push_back(&PARAM_HELP_LONG);
    databases.push_back(&PARAM_TSV);
    databases.push_back(&PARAM_REUSELATEST);
    databases.push_back(&PARAM_REMOVE_TMP_FILES);
    databases.push_back(&PARAM_COMPRESSED);
    databases.push_back(&PARAM_THREADS);
    databases.push_back(&PARAM_V);

    samplemulambda.push_back(&PARAM_N_SAMPLE);
    samplemulambda.push_back(&PARAM_THREADS);
    samplemulambda.push_back(&PARAM_V);

    //compressca
    compressca.push_back(&PARAM_COORD_STORE_MODE);
    compressca.push_back(&PARAM_THREADS);
    compressca.push_back(&PARAM_V);

    //scorecmultimer
    scoremultimer.push_back(&PARAM_MIN_ASSIGNED_CHAINS_THRESHOLD);
    scoremultimer.push_back(&PARAM_MONOMER_INCLUDE_MODE);
    scoremultimer.push_back(&PARAM_THREADS);
    scoremultimer.push_back(&PARAM_V);
    scoremultimer.push_back(&PARAM_C);
    scoremultimer.push_back(&PARAM_COV_MODE);
    scoremultimer.push_back(&PARAM_INTERFACE_LDDT_THRESHOLD);
    scoremultimer.push_back(&PARAM_CHAIN_TM_THRESHOLD);
    scoremultimer.push_back(&PARAM_TMSCORE_THRESHOLD);
    scoremultimer.push_back(&PARAM_MIN_ALIGNED_CHAINS);
    
    //makepaddeddb
    makepaddeddb.push_back(&PARAM_SUB_MAT);
    makepaddeddb.push_back(&PARAM_SCORE_BIAS);
    makepaddeddb.push_back(&PARAM_MASK_RESIDUES);
    makepaddeddb.push_back(&PARAM_MASK_PROBABILTY);
    makepaddeddb.push_back(&PARAM_WRITE_LOOKUP);
    makepaddeddb.push_back(&PARAM_THREADS);
    makepaddeddb.push_back(&PARAM_V);
    makepaddeddb.push_back(&PARAM_CLUSTER_SEARCH);

    //result2structprofile
    result2structprofile.push_back(&PARAM_SUB_MAT);
    result2structprofile.push_back(&PARAM_E);
    result2structprofile.push_back(&PARAM_MASK_PROFILE);
    result2structprofile.push_back(&PARAM_E_PROFILE);
    result2structprofile.push_back(&PARAM_NO_COMP_BIAS_CORR);
    result2structprofile.push_back(&PARAM_NO_COMP_BIAS_CORR_SCALE);
    result2structprofile.push_back(&PARAM_WG);
    result2structprofile.push_back(&PARAM_ALLOW_DELETION);
    result2structprofile.push_back(&PARAM_FILTER_MSA);
    result2structprofile.push_back(&PARAM_FILTER_MIN_ENABLE);
    result2structprofile.push_back(&PARAM_FILTER_MAX_SEQ_ID);
    result2structprofile.push_back(&PARAM_FILTER_QID);
    result2structprofile.push_back(&PARAM_FILTER_QSC);
    result2structprofile.push_back(&PARAM_FILTER_COV);
    result2structprofile.push_back(&PARAM_FILTER_NDIFF);
    result2structprofile.push_back(&PARAM_PC_MODE);
    result2structprofile.push_back(&PARAM_PCA);
    result2structprofile.push_back(&PARAM_PCB);
    result2structprofile.push_back(&PARAM_PRELOAD_MODE);
    result2structprofile.push_back(&PARAM_GAP_OPEN);
    result2structprofile.push_back(&PARAM_GAP_EXTEND);
#ifdef GAP_POS_SCORING
    result2structprofile.push_back(&PARAM_GAP_PSEUDOCOUNT);
#endif
    result2structprofile.push_back(&PARAM_THREADS);
    result2structprofile.push_back(&PARAM_COMPRESSED);
    result2structprofile.push_back(&PARAM_V);
    result2structprofile.push_back(&PARAM_PROFILE_OUTPUT_MODE);

    //createstructsubdb
    createstructsubdb.push_back(&PARAM_SUBDB_MODE);
    createstructsubdb.push_back(&PARAM_ID_MODE);
    createstructsubdb.push_back(&PARAM_V);
    
    // createmultimerreport
    createmultimerreport.push_back(&PARAM_DB_OUTPUT);
    createmultimerreport.push_back(&PARAM_THREADS);
    createmultimerreport.push_back(&PARAM_V);

    // expandmultimer
    expandmultimer.push_back(&PARAM_THREADS);
    expandmultimer.push_back(&PARAM_V);

    // convert2pdb
    convert2pdb.push_back(&PARAM_PDB_OUTPUT_MODE);
    convert2pdb.push_back(&PARAM_THREADS);
    convert2pdb.push_back(&PARAM_V);

    // structuresearchworkflow
    structuresearchworkflow = combineList(structurealign, prefilter);
    structuresearchworkflow = combineList(structuresearchworkflow, ungappedprefilter);
    structuresearchworkflow = combineList(structuresearchworkflow, tmalign);
    structuresearchworkflow = combineList(structuresearchworkflow, lolalign);
    structuresearchworkflow = combineList(structuresearchworkflow, result2structprofile);
    structuresearchworkflow.push_back(&PARAM_CLUSTER_SEARCH);
    structuresearchworkflow.push_back(&PARAM_EXHAUSTIVE_SEARCH);
    structuresearchworkflow.push_back(&PARAM_NUM_ITERATIONS);
    structuresearchworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
    structuresearchworkflow.push_back(&PARAM_REUSELATEST);
    structuresearchworkflow.push_back(&PARAM_RUNNER);

    easystructuresearchworkflow = combineList(structuresearchworkflow, structurecreatedb);
    easystructuresearchworkflow = combineList(easystructuresearchworkflow, convertalignments);
    easystructuresearchworkflow = combineList(easystructuresearchworkflow, taxonomyreport);
    easystructuresearchworkflow.push_back(&PARAM_GREEDY_BEST_HITS);

    structureclusterworkflow = combineList(prefilter, structurealign);
    structureclusterworkflow = combineList(structureclusterworkflow, structurerescorediagonal);
    structureclusterworkflow = combineList(structureclusterworkflow, tmalign);
    structureclusterworkflow = combineList(structureclusterworkflow, clust);
    structureclusterworkflow.push_back(&PARAM_CASCADED);
    structureclusterworkflow.push_back(&PARAM_CLUSTER_STEPS);
    structureclusterworkflow.push_back(&PARAM_CLUSTER_REASSIGN);
    structureclusterworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
    structureclusterworkflow.push_back(&PARAM_REUSELATEST);
    structureclusterworkflow.push_back(&PARAM_RUNNER);
    structureclusterworkflow = combineList(structureclusterworkflow, linclustworkflow);

    easystructureclusterworkflow = combineList(structureclusterworkflow, structurecreatedb);
    easystructureclusterworkflow = combineList(easystructureclusterworkflow, result2repseq);

    // multimersearchworkflow
    multimersearchworkflow = combineList(structuresearchworkflow, scoremultimer);
    multimersearchworkflow = combineList(multimersearchworkflow, expandmultimer);
    multimersearchworkflow.push_back(&PARAM_EXPAND_MULTIMER_EVALUE);
    multimersearchworkflow.push_back(&PARAM_EXPAND_MULTIMER_EVALUE_BC_COMPAT);
    multimersearchworkflow.push_back(&PARAM_MULTIMER_REPORT_MODE);
    multimersearchworkflow.push_back(&PARAM_MULTIMER_REPORT_MODE_BC_COMPAT);

    // easymultimersearchworkflow
    easymultimersearchworkflow = combineList(structurecreatedb, multimersearchworkflow);
    easymultimersearchworkflow = combineList(easymultimersearchworkflow, convertalignments);
    easymultimersearchworkflow = combineList(easymultimersearchworkflow, createmultimerreport);
    easymultimersearchworkflow = removeParameter(easymultimersearchworkflow, PARAM_PROSTT5_MODEL);

    // multimerclusterworkflow
    multimerclusterworkflow  = combineList(multimersearchworkflow, clust);

    //easymultimerclusterworkflow
    easymultimerclusterworkflow = combineList(structurecreatedb, multimerclusterworkflow);
    easymultimerclusterworkflow = combineList(easymultimerclusterworkflow, result2repseq);

    // set masking
    maskMode = 0;
    maskNrepeats = 6;
    maskProb = 0.999995;
    maskLowerCaseMode = 1;

    // createdb
    maskBfactorThreshold = 0;
    chainNameMode = 0;
    modelNameMode = 0;
    writeMapping = 0;
    writeFoldcomp = 0;
    coordStoreMode = COORD_STORE_MODE_CA_DIFF;
    inputFormat = 0; // auto detect
    fileInclude = ".*";
    fileExclude = "^$";
    prostt5SplitLength = 1024;
    prostt5Model = "";

    // search parameter
    alignmentType = ALIGNMENT_TYPE_3DI_AA;
    tmScoreThr = 0.0;
    tmScoreThrMode = TMSCORE_THRESHOLD_MODE_ALIGNMENT;
    tmAlignHitOrder = TMALIGN_HIT_ORDER_AVG;
    lddtThr = 0.0;
    evalThr = 10;
    sortByStructureBits = 1;
    clusterSearch = 0;
    minDiagScoreThr = 30;
    minAssignedChainsThreshold = 0.0;
    monomerIncludeMode = 0;
    tmAlignFast = 1;
    exactTMscore = 0;
    gapOpen = 10;
    gapExtend = 1;
    nsample = 5000;
    dbSuffixList = "_h,_ss,_ca";
    indexExclude = 0;

    // profiles
    evalProfile = 0.1;
    // multimer
    eValueThrExpandMultimer = 10000.0;
    multimerReportMode = 1;
    dbExtractionMode = DB_EXTRACT_MODE_CHAIN;
    distanceThreshold = 10.0;
    // filtMultimerTmThr = 0.0;
    filtChainTmThr = 0.3;
    filtInterfaceLddtThr = 0;
    minAlignedChains = 2;

    // LoLalign
    multiDomain = 1;

    citations.emplace(CITATION_FOLDSEEK, "van Kempen, M., Kim, S.S., Tumescheit, C., Mirdita, M., Lee, J., Gilchrist, C.L.M., Söding, J., and Steinegger, M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology, doi:10.1038/s41587-023-01773-0 (2023)");
    citations.emplace(CITATION_FOLDSEEK_MULTIMER, "Kim, W., Mirdita, M., Levy Karin, E., Gilchrist, C.L.M., Schweke, H., Söding, J., Levy, E., and Steinegger, M. Rapid and sensitive protein complex alignment with Foldseek-Multimer. Nature Methods, doi:10.1038/s41592-025-02593-7 (2025)");
    citations.emplace(CITATION_PROSTT5, "Heinzinger, M., Weissenow, K., Gomez Sanchez, J., Henkel, A., Mirdita, M., Steinegger, M., and Burkhard, R. Bilingual Language Model for Protein Sequence and Structure. NAR Genomics and Bioinformatics, doi:10.1093/nargab/lqae150 (2024)");
    
    //rewrite param vals.
    PARAM_FORMAT_OUTPUT.description = "Choose comma separated list of output columns from: query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen\ntstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,q3di,t3di,qheader,theader,qaln,taln,q3dialn,t3dialn,mismatch,qcov,tcov\nqset,qsetid,tset,tsetid,taxid,taxname,taxlineage,\nlddt,lddtfull,qca,tca,t,u,qtmscore,ttmscore,alntmscore,rmsd,prob\ncomplexqtmscore,complexttmscore,complexu,complext,qcomplexcoverage,tcomplexcoverage,qchaintms,tchaintms,qchains,tchains,interfacelddt,complexassignid\n";

    // allow higher values for GGML debug trace
    PARAM_V.regex = "^[0-4]{1}$";
}



std::vector<int> LocalParameters::getOutputFormat(
    int formatMode, const std::string &outformat, bool &needSequences, bool &need3Di, bool &needBacktrace, bool &needFullHeaders,
    bool &needLookup, bool &needSource, bool &needTaxonomyMapping, bool &needTaxonomy, bool &needQCa, bool &needTCa, bool &needTMaligner,
    bool &needLDDT) {
    std::vector<int> formatCodes;
    if (formatMode == Parameters::FORMAT_ALIGNMENT_SAM || formatMode == Parameters::FORMAT_ALIGNMENT_HTML) {
        needSequences = true;
        needBacktrace = true;
	    needQCa = true;
        needTCa = true;
        return formatCodes;
    }
    std::vector<std::string> outformatSplit = Util::split(outformat, ",");
    int code = 0;
    for (size_t i = 0; i < outformatSplit.size(); ++i) {
        if(outformatSplit[i].compare("query") == 0){ code = Parameters::OUTFMT_QUERY;}
        else if (outformatSplit[i].compare("target") == 0){ code = Parameters::OUTFMT_TARGET;}
        else if (outformatSplit[i].compare("evalue") == 0){ code = Parameters::OUTFMT_EVALUE;}
        else if (outformatSplit[i].compare("gapopen") == 0){ code = Parameters::OUTFMT_GAPOPEN;}
        else if (outformatSplit[i].compare("pident") == 0){ code = Parameters::OUTFMT_PIDENT;}
        else if (outformatSplit[i].compare("fident") == 0){ code = Parameters::OUTFMT_FIDENT;}
        else if (outformatSplit[i].compare("nident") == 0){ code = Parameters::OUTFMT_NIDENT;}
        else if (outformatSplit[i].compare("qstart") == 0){ code = Parameters::OUTFMT_QSTART;}
        else if (outformatSplit[i].compare("qend") == 0){ code = Parameters::OUTFMT_QEND;}
        else if (outformatSplit[i].compare("qlen") == 0){ code = Parameters::OUTFMT_QLEN;}
        else if (outformatSplit[i].compare("tstart") == 0){ code = Parameters::OUTFMT_TSTART;}
        else if (outformatSplit[i].compare("tend") == 0){ code = Parameters::OUTFMT_TEND;}
        else if (outformatSplit[i].compare("tlen") == 0){ code = Parameters::OUTFMT_TLEN;}
        else if (outformatSplit[i].compare("alnlen") == 0){ code = Parameters::OUTFMT_ALNLEN;}
        else if (outformatSplit[i].compare("bits") == 0){ code = Parameters::OUTFMT_BITS;}
        else if (outformatSplit[i].compare("cigar") == 0){ needBacktrace = true; code = Parameters::OUTFMT_CIGAR;}
        else if (outformatSplit[i].compare("qseq") == 0){ needSequences = true; code = Parameters::OUTFMT_QSEQ;}
        else if (outformatSplit[i].compare("tseq") == 0){ needSequences = true; code = Parameters::OUTFMT_TSEQ;}
        else if (outformatSplit[i].compare("q3di") == 0) { need3Di = true; code = LocalParameters::OUTFMT_Q3DI; }
        else if (outformatSplit[i].compare("t3di") == 0) { need3Di = true; code = LocalParameters::OUTFMT_T3DI; }
        else if (outformatSplit[i].compare("qheader") == 0){ needFullHeaders = true; code = Parameters::OUTFMT_QHEADER;}
        else if (outformatSplit[i].compare("theader") == 0){ needFullHeaders = true; code = Parameters::OUTFMT_THEADER;}
        else if (outformatSplit[i].compare("qaln") == 0){ needBacktrace = true; needSequences = true; code = Parameters::OUTFMT_QALN;}
        else if (outformatSplit[i].compare("taln") == 0){ needBacktrace = true; needSequences = true; code = Parameters::OUTFMT_TALN;}
        else if (outformatSplit[i].compare("q3dialn") == 0){ needBacktrace = true; need3Di = true; code = LocalParameters::OUTFMT_Q3DIALN;}
        else if (outformatSplit[i].compare("t3dialn") == 0){ needBacktrace = true; need3Di = true; code = LocalParameters::OUTFMT_T3DIALN;}
        else if (outformatSplit[i].compare("mismatch") == 0){ code = Parameters::OUTFMT_MISMATCH;}
        else if (outformatSplit[i].compare("qcov") == 0){ code = Parameters::OUTFMT_QCOV;}
        else if (outformatSplit[i].compare("tcov") == 0){ code = Parameters::OUTFMT_TCOV;}
        else if (outformatSplit[i].compare("qca") == 0){ needQCa = true; code = LocalParameters::OUTFMT_QCA;}
        else if (outformatSplit[i].compare("tca") == 0){ needTCa = true; code = LocalParameters::OUTFMT_TCA;}
        else if (outformatSplit[i].compare("u") == 0){ needQCa = true; needTCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_U;}
        else if (outformatSplit[i].compare("t") == 0){ needQCa = true; needTCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_T;}
        else if (outformatSplit[i].compare("alntmscore") == 0){ needQCa = true; needTCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_ALNTMSCORE;}
        else if (outformatSplit[i].compare("qtmscore") == 0){ needQCa = true; needTCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_QTMSCORE;}
        else if (outformatSplit[i].compare("ttmscore") == 0){ needQCa = true; needTCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_TTMSCORE;}
        else if (outformatSplit[i].compare("rmsd") == 0){ needQCa = true; needTCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_RMSD;}
        else if (outformatSplit[i].compare("qset") == 0){ needLookup = true; needSource = true; code = Parameters::OUTFMT_QSET;}
        else if (outformatSplit[i].compare("qsetid") == 0){ needLookup = true; needSource = true; code = Parameters::OUTFMT_QSETID;}
        else if (outformatSplit[i].compare("tset") == 0){ needLookup = true; code = Parameters::OUTFMT_TSET;}
        else if (outformatSplit[i].compare("tsetid") == 0){ needLookup = true; needSource = true; code = Parameters::OUTFMT_TSETID;}
        else if (outformatSplit[i].compare("taxid") == 0){ needTaxonomyMapping = true; code = Parameters::OUTFMT_TAXID;}
        else if (outformatSplit[i].compare("taxname") == 0){ needTaxonomyMapping = true; needTaxonomy = true; code = Parameters::OUTFMT_TAXNAME;}
        else if (outformatSplit[i].compare("taxlineage") == 0){ needTaxonomyMapping = true; needTaxonomy = true; code = Parameters::OUTFMT_TAXLIN;}
        else if (outformatSplit[i].compare("empty") == 0){ code = Parameters::OUTFMT_EMPTY;}
        else if (outformatSplit[i].compare("lddt") == 0) { needQCa = true; needTCa = true; needLDDT = true; needBacktrace = true; code = LocalParameters::OUTFMT_LDDT; }
        else if (outformatSplit[i].compare("lddtfull") == 0) { needQCa = true; needTCa = true; needLDDT = true; needBacktrace = true; code = LocalParameters::OUTFMT_LDDT_FULL; }
        else if (outformatSplit[i].compare("prob") == 0) { needQCa = true; needTCa = true; needLDDT = true; needBacktrace = true; needTMaligner = true; code = LocalParameters::OUTFMT_PROBTP; }
        // TODO
        else if (outformatSplit[i].compare("complexqtmscore")==0 || outformatSplit[i].compare("multimerqtmscore")==0){code=LocalParameters::OUTFMT_Q_COMPLEX_TMSCORE; }
        else if (outformatSplit[i].compare("complexttmscore")==0 || outformatSplit[i].compare("multimerttmscore")==0){code=LocalParameters::OUTFMT_T_COMPLEX_TMSCORE;}
        else if (outformatSplit[i].compare("complexassignid")==0 || outformatSplit[i].compare("multimerassignid")==0){code=LocalParameters::OUTFMT_ASSIGN_ID;}
        else if (outformatSplit[i].compare("complexu")==0 || outformatSplit[i].compare("multimeru")==0){code=LocalParameters::OUTFMT_COMPLEX_U;}
        else if (outformatSplit[i].compare("complext")==0 || outformatSplit[i].compare("multimert")==0){code=LocalParameters::OUTFMT_COMPLEX_T;}
        else if (outformatSplit[i].compare("qcomplexcoverage")==0 || outformatSplit[i].compare("qmultimercoverage")==0){code=LocalParameters::OUTFMT_Q_COMPLEX_COV;}
        else if (outformatSplit[i].compare("tcomplexcoverage")==0 || outformatSplit[i].compare("tmultimercoverage")==0){code=LocalParameters::OUTFMT_T_COMPLEX_COV;}
        else if (outformatSplit[i].compare("qchaintms")==0 ){code=LocalParameters::OUTFMT_COMPLEX_QCHAINTMS;}
        else if (outformatSplit[i].compare("tchaintms")==0 ){code=LocalParameters::OUTFMT_COMPLEX_TCHAINTMS;}
        else if (outformatSplit[i].compare("qchains")==0 ){code=LocalParameters::OUTFMT_COMPLEX_QNAME;}
        else if (outformatSplit[i].compare("tchains")==0 ){code=LocalParameters::OUTFMT_COMPLEX_TNAME;}
        else if (outformatSplit[i].compare("interfacelddt")==0 ){code=LocalParameters::OUTFMT_INTERFACE_LDDT;}
        else {
            Debug(Debug::ERROR) << "Format code " << outformatSplit[i] << " does not exist.";
            EXIT(EXIT_FAILURE);
        }
        formatCodes.push_back(code);
    }
    return formatCodes;
}


std::vector<int> FoldSeekDbValidator::tmscore = {LocalParameters::DBTYPE_TMSCORE};
std::vector<int> FoldSeekDbValidator::cadb = {LocalParameters::DBTYPE_CA_ALPHA};
std::vector<int> FoldSeekDbValidator::flatfileStdinAndFolder = {LocalParameters::DBTYPE_FLATFILE, LocalParameters::DBTYPE_STDIN,LocalParameters::DBTYPE_DIRECTORY};
std::vector<int> FoldSeekDbValidator::flatfileAndFolder = {LocalParameters::DBTYPE_FLATFILE, LocalParameters::DBTYPE_DIRECTORY};
