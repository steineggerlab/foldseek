#include "LocalParameters.h"
#include "Command.h"
#include "Debug.h"
#include "mat3di.out.h"


const int LocalParameters::DBTYPE_CA_ALPHA = 101;
const int LocalParameters::DBTYPE_TMSCORE = 102;

LocalParameters::LocalParameters() :
        Parameters(),
        PARAM_PREF_MODE(PARAM_PREF_MODE_ID,"--prefilter-mode", "Prefilter mode", "prefilter mode: 0: kmer/ungapped 1: ungapped, 2: nofilter",typeid(int), (void *) &prefMode, "^[0-2]{1}$"),
        PARAM_TMSCORE_THRESHOLD(PARAM_TMSCORE_THRESHOLD_ID,"--tmscore-threshold", "TMscore threshold", "accept alignments with a tmsore > thr [0.0,1.0]",typeid(float), (void *) &tmScoreThr, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_TMALIGN_HIT_ORDER(PARAM_TMALIGN_HIT_ORDER_ID,"--tmalign-hit-order", "TMalign hit order", "order hits by 0: (qTM+tTM)/2, 1: qTM, 2: tTM, 3: min(qTM,tTM) 4: max(qTM,tTM)",typeid(int), (void *) &tmAlignHitOrder, "^[0-4]{1}$"),
        PARAM_LDDT_THRESHOLD(PARAM_LDDT_THRESHOLD_ID,"--lddt-threshold", "LDDT threshold", "accept alignments with a lddt > thr [0.0,1.0]",typeid(float), (void *) &lddtThr, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_SORT_BY_STRUCTURE_BITS(PARAM_SORT_BY_STRUCTURE_BITS_ID,"--sort-by-structure-bits", "Sort by structure bit score", "sort by bits*sqrt(alnlddt*alntmscore)",typeid(int), (void *) &sortByStructureBits, "^[0-1]{1}$", MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT),
        PARAM_MASK_BFACTOR_THRESHOLD(PARAM_MASK_BFACTOR_THRESHOLD_ID,"--mask-bfactor-threshold", "Mask b-factor threshold", "mask residues for seeding if b-factor < thr [0,100]",typeid(float), (void *) &maskBfactorThreshold, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_ALIGNMENT_TYPE(PARAM_ALIGNMENT_TYPE_ID,"--alignment-type", "Alignment type", "How to compute the alignment:\n0: 3di alignment\n1: TM alignment\n2: 3Di+AA",typeid(int), (void *) &alignmentType, "^[0-2]{1}$"),
        PARAM_CHAIN_NAME_MODE(PARAM_CHAIN_NAME_MODE_ID,"--chain-name-mode", "Chain name mode", "Add chain to name:\n0: auto\n1: always add\n",typeid(int), (void *) &chainNameMode, "^[0-1]{1}$", MMseqsParameter::COMMAND_EXPERT),
        PARAM_WRITE_MAPPING(PARAM_WRITE_MAPPING_ID, "--write-mapping", "Write mapping file", "write _mapping file containing mapping from internal id to taxonomic identifier", typeid(int), (void *) &writeMapping, "^[0-1]{1}", MMseqsParameter::COMMAND_EXPERT),
        PARAM_TMALIGN_FAST(PARAM_TMALIGN_FAST_ID,"--tmalign-fast", "TMalign fast","turn on fast search in TM-align" ,typeid(int), (void *) &tmAlignFast, "^[0-1]{1}$"),
        PARAM_N_SAMPLE(PARAM_N_SAMPLE_ID, "--n-sample", "Sample size","pick N random sample" ,typeid(int), (void *) &nsample, "^[0-9]{1}[0-9]*$"),
        PARAM_COORD_STORE_MODE(PARAM_COORD_STORE_MODE_ID, "--coord-store-mode", "Coord store mode", "Coordinate storage mode: \n1: C-alpha as float\n2: C-alpha as difference (uint16_t)", typeid(int), (void *) &coordStoreMode, "^[1-2]{1}$",MMseqsParameter::COMMAND_EXPERT),
        PARAM_PROTEIN_CHAINS_ONLY(PARAM_PROTEIN_CHAINS_ONLY_ID, "--protein-chains-only", "store protein chains only", "store protein chains only: \n0: store all types of chains\n1: store protein chains only", typeid(int), (void *) &storeProteinChainsOnly, "^[0-1]{1}$"),
        PARAM_MIN_ASSIGNED_CHAINS_THRESHOLD(PARAM_MIN_ASSIGNED_CHAINS_THRESHOLD_ID, "--min-assigned-chains-ratio", "Minimum assigned chains percentage Threshold", "minimum percentage of assigned chains out of all query chains > thr [0,100] %", typeid(float), (void *) & minAssignedChainsThreshold, "^[0-9]*(\\.[0-9]+)?$"),
        PARAM_CLUSTER_SEARCH(PARAM_CLUSTER_SEARCH_ID, "--cluster-search", "Cluster search", "first find representative then align all cluster members", typeid(int), (void *) &clusterSearch, "^[0-1]{1}$",MMseqsParameter::COMMAND_MISC),
        PARAM_FILE_INCLUDE(PARAM_FILE_INCLUDE_ID, "--file-include", "File Inclusion Regex", "Include file names based on this regex", typeid(std::string), (void *) &fileInclude, "^.*$"),
        PARAM_FILE_EXCLUDE(PARAM_FILE_EXCLUDE_ID, "--file-exclude", "File Exclusion Regex", "Exclude file names based on this regex", typeid(std::string), (void *) &fileExclude, "^.*$")

{
    PARAM_ALIGNMENT_MODE.description = "How to compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also start_pos and cov\n3: also seq.id";
    PARAM_ALIGNMENT_MODE.regex = "^[0-3]{1}$";
    PARAM_ALIGNMENT_MODE.category = MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT;
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
                                    "6: BLAST-TAB + q/db complex tm (only for scorecomplex result)\n"
                                    "BLAST-TAB (0) and BLAST-TAB + column headers (4)"
                                    "support custom output formats (--format-output)\n"
                                    "(5) Superposed PDB files (Calpha only)";

    // TODO
    PARAM_FORMAT_MODE.regex = "^[0-6]{1}$";
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

    scoringMatrixFile = "3di.out";
    seedScoringMatrixFile = "3di.out";
    substitutionMatrices.emplace_back("3di.out", mat3di_out, mat3di_out_len);

    // structurecreatedb
    structurecreatedb.push_back(&PARAM_CHAIN_NAME_MODE);
    structurecreatedb.push_back(&PARAM_WRITE_MAPPING);
    structurecreatedb.push_back(&PARAM_MASK_BFACTOR_THRESHOLD);
    structurecreatedb.push_back(&PARAM_COORD_STORE_MODE);
    structurecreatedb.push_back(&PARAM_WRITE_LOOKUP);
    structurecreatedb.push_back(&PARAM_TAR_INCLUDE);
    structurecreatedb.push_back(&PARAM_TAR_EXCLUDE);
    // protein chain only
    structurecreatedb.push_back(&PARAM_PROTEIN_CHAINS_ONLY);
    structurecreatedb.push_back(&PARAM_FILE_INCLUDE);
    structurecreatedb.push_back(&PARAM_FILE_EXCLUDE);
    structurecreatedb.push_back(&PARAM_THREADS);
    structurecreatedb.push_back(&PARAM_V);

    // tmalign
    tmalign.push_back(&PARAM_MIN_SEQ_ID);
    tmalign.push_back(&PARAM_C);
    tmalign.push_back(&PARAM_COV_MODE);
    tmalign.push_back(&PARAM_MAX_REJECTED);
    tmalign.push_back(&PARAM_MAX_ACCEPT);
    tmalign.push_back(&PARAM_ADD_BACKTRACE);
    tmalign.push_back(&PARAM_INCLUDE_IDENTITY);
    tmalign.push_back(&PARAM_TMSCORE_THRESHOLD);
    tmalign.push_back(&PARAM_TMALIGN_HIT_ORDER);
    tmalign.push_back(&PARAM_TMALIGN_FAST);
    tmalign.push_back(&PARAM_PRELOAD_MODE);
    tmalign.push_back(&PARAM_THREADS);
    tmalign.push_back(&PARAM_V);

    structurerescorediagonal.push_back(&PARAM_TMSCORE_THRESHOLD);
    structurerescorediagonal.push_back(&PARAM_LDDT_THRESHOLD);
    structurerescorediagonal.push_back(&PARAM_ALIGNMENT_TYPE);
    structurerescorediagonal = combineList(structurerescorediagonal, align);

    structurealign.push_back(&PARAM_TMSCORE_THRESHOLD);
    structurealign.push_back(&PARAM_LDDT_THRESHOLD);
    structurealign.push_back(&PARAM_SORT_BY_STRUCTURE_BITS);
    structurealign.push_back(&PARAM_ALIGNMENT_TYPE);
    structurealign = combineList(structurealign, align);
//    tmalign.push_back(&PARAM_GAP_OPEN);
//    tmalign.push_back(&PARAM_GAP_EXTEND);
    // strucclust
    strucclust = combineList(clust, structurealign);
    strucclust = combineList(strucclust, structurerescorediagonal);
    strucclust = combineList(strucclust, kmermatcher);
    strucclust.push_back(&PARAM_REMOVE_TMP_FILES);
    strucclust.push_back(&PARAM_RUNNER);
    // structuresearchworkflow
    structuresearchworkflow = combineList(structurealign, prefilter);
    structuresearchworkflow = combineList(tmalign, structuresearchworkflow);
    structuresearchworkflow.push_back(&PARAM_EXHAUSTIVE_SEARCH);
    structuresearchworkflow.push_back(&PARAM_PREF_MODE);
    structuresearchworkflow.push_back(&PARAM_NUM_ITERATIONS);
    structuresearchworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
    structuresearchworkflow.push_back(&PARAM_RUNNER);
    structuresearchworkflow.push_back(&PARAM_REUSELATEST);
    structuresearchworkflow.push_back(&PARAM_CLUSTER_SEARCH);

    easystructuresearchworkflow = combineList(structuresearchworkflow, structurecreatedb);
    easystructuresearchworkflow = combineList(easystructuresearchworkflow, convertalignments);
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


    databases.push_back(&PARAM_HELP);
    databases.push_back(&PARAM_HELP_LONG);
    databases.push_back(&PARAM_TSV);
    databases.push_back(&PARAM_REUSELATEST);
    databases.push_back(&PARAM_REMOVE_TMP_FILES);
    databases.push_back(&PARAM_COMPRESSED);
    databases.push_back(&PARAM_THREADS);
    databases.push_back(&PARAM_V);
    //easystructureclusterworkflow = combineList(structuresearchworkflow, structurecreatedb);
    samplemulambda.push_back(&PARAM_N_SAMPLE);
    samplemulambda.push_back(&PARAM_THREADS);
    samplemulambda.push_back(&PARAM_V);
    //scorecomplex
    scorecomplex.push_back(&PARAM_THREADS);
    scorecomplex.push_back(&PARAM_V);
    scorecomplex.push_back(&PARAM_MIN_ASSIGNED_CHAINS_THRESHOLD);
    easyscorecomplexworkflow = combineList(easyscorecomplexworkflow, structurecreatedb);
    easyscorecomplexworkflow = combineList(easyscorecomplexworkflow, structuresearchworkflow);
    easyscorecomplexworkflow = combineList(easyscorecomplexworkflow, scorecomplex);
    easyscorecomplexworkflow = combineList(easyscorecomplexworkflow, convertalignments);
    createcomplexreport.push_back(&PARAM_THREADS);
    createcomplexreport.push_back(&PARAM_V);
    prefMode = PREF_MODE_KMER;
    alignmentType = ALIGNMENT_TYPE_3DI_AA;
    tmScoreThr = 0.0;
    tmAlignHitOrder = TMALIGN_HIT_ORDER_AVG;
    lddtThr = 0.0;
    evalThr = 10;
    sortByStructureBits = 1;
    minDiagScoreThr = 30;
    maskBfactorThreshold = 0;
    chainNameMode = 0;
    writeMapping = 0;
    tmAlignFast = 1;
    gapOpen = 10;
    gapExtend = 1;
    nsample = 5000;
    maskLowerCaseMode = 1;
    coordStoreMode = COORD_STORE_MODE_CA_DIFF;
    clusterSearch = 0;
    fileInclude = ".*";
    fileExclude = "^$";

    citations.emplace(CITATION_FOLDSEEK, "van Kempen, M., Kim, S.S., Tumescheit, C., Mirdita, M., Lee, J., Gilchrist, C.L.M., Söding, J., and Steinegger, M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology, doi:10.1038/s41587-023-01773-0 (2023)");

    //rewrite param vals.
//    PARAM_FORMAT_OUTPUT.description = "Choose comma separated list of output columns from: query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen\ntstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,mismatch,qcov,tcov\nqset,qsetid,tset,tsetid,taxid,taxname,taxlineage,\nlddt,lddtfull,qca,tca,t,u,qtmscore,ttmscore,alntmscore,rmsd,prob\n";
    PARAM_FORMAT_OUTPUT.description = "Choose comma separated list of output columns from: query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen\ntstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,mismatch,qcov,tcov\nqset,qsetid,tset,tsetid,taxid,taxname,taxlineage,\nlddt,lddtfull,qca,tca,t,u,qtmscore,ttmscore,alntmscore,rmsd,prob\nqcomplextmscore,tcomplextmscore,assignid\n";
}



std::vector<int> LocalParameters::getOutputFormat(int formatMode, const std::string &outformat, bool &needSequences, bool &needBacktrace, bool &needFullHeaders,
                                             bool &needLookup, bool &needSource, bool &needTaxonomyMapping, bool &needTaxonomy, bool &needCa, bool &needTMaligner,
                                             bool &needLDDT) {
    std::vector<int> formatCodes;
    if (formatMode == Parameters::FORMAT_ALIGNMENT_SAM || formatMode == Parameters::FORMAT_ALIGNMENT_HTML) {
        needSequences = true;
        needBacktrace = true;
	needCa = true;
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
        else if (outformatSplit[i].compare("qheader") == 0){ needFullHeaders = true; code = Parameters::OUTFMT_QHEADER;}
        else if (outformatSplit[i].compare("theader") == 0){ needFullHeaders = true; code = Parameters::OUTFMT_THEADER;}
        else if (outformatSplit[i].compare("qaln") == 0){ needBacktrace = true; needSequences = true; code = Parameters::OUTFMT_QALN;}
        else if (outformatSplit[i].compare("taln") == 0){ needBacktrace = true; needSequences = true; code = Parameters::OUTFMT_TALN;}
        else if (outformatSplit[i].compare("mismatch") == 0){ code = Parameters::OUTFMT_MISMATCH;}
        else if (outformatSplit[i].compare("qcov") == 0){ code = Parameters::OUTFMT_QCOV;}
        else if (outformatSplit[i].compare("tcov") == 0){ code = Parameters::OUTFMT_TCOV;}
        else if (outformatSplit[i].compare("qca") == 0){ needCa = true; code = LocalParameters::OUTFMT_QCA;}
        else if (outformatSplit[i].compare("tca") == 0){ needCa = true; code = LocalParameters::OUTFMT_TCA;}
        else if (outformatSplit[i].compare("u") == 0){ needCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_U;}
        else if (outformatSplit[i].compare("t") == 0){ needCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_T;}
        else if (outformatSplit[i].compare("alntmscore") == 0){ needCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_ALNTMSCORE;}
        else if (outformatSplit[i].compare("qtmscore") == 0){ needCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_QTMSCORE;}
        else if (outformatSplit[i].compare("ttmscore") == 0){ needCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_TTMSCORE;}
        else if (outformatSplit[i].compare("rmsd") == 0){ needCa = true; needTMaligner = true; needBacktrace=true; code = LocalParameters::OUTFMT_RMSD;}
        else if (outformatSplit[i].compare("qset") == 0){ needLookup = true; needSource = true; code = Parameters::OUTFMT_QSET;}
        else if (outformatSplit[i].compare("qsetid") == 0){ needLookup = true; needSource = true; code = Parameters::OUTFMT_QSETID;}
        else if (outformatSplit[i].compare("tset") == 0){ needLookup = true; code = Parameters::OUTFMT_TSET;}
        else if (outformatSplit[i].compare("tsetid") == 0){ needLookup = true; needSource = true; code = Parameters::OUTFMT_TSETID;}
        else if (outformatSplit[i].compare("taxid") == 0){ needTaxonomyMapping = true; code = Parameters::OUTFMT_TAXID;}
        else if (outformatSplit[i].compare("taxname") == 0){ needTaxonomyMapping = true; needTaxonomy = true; code = Parameters::OUTFMT_TAXNAME;}
        else if (outformatSplit[i].compare("taxlineage") == 0){ needTaxonomyMapping = true; needTaxonomy = true; code = Parameters::OUTFMT_TAXLIN;}
        else if (outformatSplit[i].compare("empty") == 0){ code = Parameters::OUTFMT_EMPTY;}
        else if (outformatSplit[i].compare("lddt") == 0) { needCa = true; needLDDT = true; needBacktrace = true; code = LocalParameters::OUTFMT_LDDT; }
        else if (outformatSplit[i].compare("lddtfull") == 0) { needCa = true; needLDDT = true; needBacktrace = true; code = LocalParameters::OUTFMT_LDDT_FULL; }
        else if (outformatSplit[i].compare("prob") == 0) { needCa = true; needLDDT = true; needBacktrace = true; needTMaligner = true; code = LocalParameters::OUTFMT_PROBTP; }
        // TODO
        else if (outformatSplit[i].compare("qcomplextmscore")==0){code=LocalParameters::OUTFMT_Q_COMPLEX_TMSCORE; }
        else if (outformatSplit[i].compare("tcomplextmscore")==0){code=LocalParameters::OUTFMT_T_COMPLEX_TMSCORE;}
        else if (outformatSplit[i].compare("assignid")==0){code=LocalParameters::OUTFMT_ASSIGN_ID;}
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
