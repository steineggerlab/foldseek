#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "TranslateNucl.h"
#include "Sequence.h"
#include "Orf.h"
#include "MemoryMapped.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"
#include "Coordinate16.h"
#include "createcomplexreport.h"

#define ZSTD_STATIC_LINKING_ONLY


#include <zstd.h>

#include "main.js.h"
#include "vendor.js.zst.h"

#include "TMaligner.h"
#include "LDDT.h"
#include "CalcProbTP.h"
#include <map>

#ifdef OPENMP
#include <omp.h>
#endif


const char *singleLetterToThree(char singleLetterAminoacid) {
    switch (singleLetterAminoacid) {
        case 'A':
            return "ALA";
        case 'R':
            return "ARG";
        case 'N':
            return "ASN";
        case 'D':
            return "ASP";
        case 'C':
            return "CYS";
        case 'Q':
            return "GLN";
        case 'E':
            return "GLU";
        case 'G':
            return "GLY";
        case 'H':
            return "HIS";
        case 'I':
            return "ILE";
        case 'L':
            return "LEU";
        case 'K':
            return "LYS";
        case 'M':
            return "MET";
        case 'F':
            return "PHE";
        case 'P':
            return "PRO";
        case 'S':
            return "SER";
        case 'T':
            return "THR";
        case 'W':
            return "TRP";
        case 'Y':
            return "TYR";
        case 'V':
            return "VAL";
        default:
            return "UNK";
    }

}




void caToStr(float *ca, size_t len, std::string & ret) {
    for (size_t i = 0; i < len; i++) {
        ret.append(SSTR(ca[i]));
        ret.push_back(',');
        ret.append(SSTR(ca[len+i]));
        ret.push_back(',');
        ret.append(SSTR(ca[2*len+i]));
        ret.push_back(',');
    }
}

void structurePrintSeqBasedOnAln(std::string &out, const char *seq, unsigned int offset,
                        const std::string &bt, bool reverse, bool isReverseStrand,
                        bool translateSequence, const TranslateNucl &translateNucl) {
    unsigned int seqPos = 0;
    char codon[3];
    for (uint32_t i = 0; i < bt.size(); ++i) {
        char seqChar = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos]) : seq[offset + seqPos];
        if (translateSequence) {
            codon[0] = (isReverseStrand == true) ? Orf::complement(seq[offset - seqPos])     : seq[offset + seqPos];
            codon[1] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+1)]) : seq[offset + (seqPos+1)];
            codon[2] = (isReverseStrand == true) ? Orf::complement(seq[offset - (seqPos+2)]) : seq[offset + (seqPos+2)];
            seqChar = translateNucl.translateSingleCodon(codon);
        }
        switch (bt[i]) {
            case 'M':
                out.append(1, seqChar);
                seqPos += (translateSequence) ?  3 : 1;
                break;
            case 'I':
                if (reverse == true) {
                    out.append(1, '-');
                } else {
                    out.append(1, seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                }
                break;
            case 'D':
                if (reverse == true) {
                    out.append(1, seqChar);
                    seqPos += (translateSequence) ?  3 : 1;
                } else {
                    out.append(1, '-');
                }
                break;
        }
    }
}

/*
query       Query sequence label
target      Target sequence label
evalue      E-value
gapopen     Number of gap opens
pident      Percentage of identical matches
nident      Number of identical matches
qstart      1-based start position of alignment in query sequence
qend        1-based end position of alignment in query sequence
qlen        Query sequence length
tstart      1-based start position of alignment in target sequence
tend        1-based end position of alignment in target sequence
tlen        Target sequence length
alnlen      Number of alignment columns
raw         Raw alignment score
bits        Bit score
cigar       Alignment as string M=letter pair, D=delete (gap in query), I=insert (gap in target)
qseq        Full-length query sequence
tseq        Full-length target sequence
qheader     Header of Query sequence
theader     Header of Target sequence
qaln        Aligned query sequence with gaps
taln        Aligned target sequence with gaps
qframe      Query frame (-3 to +3)
tframe      Target frame (-3 to +3)
mismatch    Number of mismatches
qcov        Fraction of query sequence covered by alignment
tcov        Fraction of target sequence covered by alignment
qset        Query set
tset        Target set
qca        Query ca
tca        Target ca
 */

std::map<unsigned int, unsigned int> structureReadKeyToSet(const std::string& file) {
    std::map<unsigned int, unsigned int> mapping;
    if (file.length() == 0) {
        return mapping;
    }

    MemoryMapped lookup(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char* data = (char *) lookup.getData();
    const char* entry[255];
    while (*data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        mapping.emplace(Util::fast_atoi<unsigned int>(entry[0]), Util::fast_atoi<unsigned int>(entry[2]));
        data = Util::skipLine(data);
    }
    lookup.close();
    return mapping;
}


std::map<unsigned int, std::string> structureReadSetToSource(const std::string& file) {
    std::map<unsigned int, std::string> mapping;
    if (file.length() == 0) {
        return mapping;
    }

    MemoryMapped source(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char* data = (char *) source.getData();
    const char* entry[255];
    while (*data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 2) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        data = Util::skipLine(data);
        std::string source(entry[1], data - entry[1] - 1);
        mapping.emplace(Util::fast_atoi<unsigned int>(entry[0]), source);
    }
    source.close();
    return mapping;
}

int structureconvertalis(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    int format = par.formatAlignmentMode;
    bool addColumnHeaders = false;
    if (format == Parameters::FORMAT_ALIGNMENT_BLAST_TAB_WITH_HEADERS) {
        format = Parameters::FORMAT_ALIGNMENT_BLAST_TAB;
        addColumnHeaders = true;
    }
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    bool needSequenceDB = false;
    bool needBacktrace = false;
    bool needFullHeaders = false;
    bool needLookup = false;
    bool needSource = false;
    bool needTaxonomy = false;
    bool needQCA = false;
    bool needTCA = false;
    bool needTaxonomyMapping = false;
    bool needTMaligner = false;
    bool needLDDT = false;
    std::vector<int> outcodes = LocalParameters::getOutputFormat(format, par.outfmt, needSequenceDB, needBacktrace, needFullHeaders,
                                                                  needLookup, needSource, needTaxonomyMapping, needTaxonomy, needQCA, needTCA, needTMaligner, needLDDT);


    if(LocalParameters::FORMAT_ALIGNMENT_PDB_SUPERPOSED == format){
        needTMaligner = true;
        needQCA = true;
        needTCA = true;
        needSequenceDB = true;
    }
    NcbiTaxonomy* t = NULL;
    if(needTaxonomy){
        std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
        t = NcbiTaxonomy::openTaxonomy(db2NoIndexName);
    }
    MappingReader* mapping = NULL;
    if (needTaxonomy || needTaxonomyMapping) {
        std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
        mapping = new MappingReader(db2NoIndexName);
    }

    bool isTranslatedSearch = false;

    int dbaccessMode = needSequenceDB ? (DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA) : (DBReader<unsigned int>::USE_INDEX);

    std::map<unsigned int, unsigned int> qKeyToSet;
    std::map<unsigned int, unsigned int> tKeyToSet;
    if (needLookup) {
        std::string file1 = par.db1 + ".lookup";
        std::string file2 = par.db2 + ".lookup";
        qKeyToSet = structureReadKeyToSet(file1);
        tKeyToSet = structureReadKeyToSet(file2);
    }

    std::map<unsigned int, std::string> qSetToSource;
    std::map<unsigned int, std::string> tSetToSource;
    if (needSource) {
        std::string file1 = par.db1 + ".source";
        std::string file2 = par.db2 + ".source";
        qSetToSource = structureReadSetToSource(file1);
        tSetToSource = structureReadSetToSource(file2);
    }

    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(FileUtil::parseDbType(par.db3.c_str()));
    bool isExtendedAlignment = extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC;

    IndexReader *tDbr;
    IndexReader *tDbrHeader;
    if (sameDB) {
        tDbr = &qDbr;
        tDbrHeader= &qDbrHeader;
    } else {
        tDbr = new IndexReader(par.db2, par.threads,
                               isExtendedAlignment ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
                               (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tDbrHeader = new IndexReader(par.db2, par.threads,
                                     isExtendedAlignment ? IndexReader::SRC_HEADERS : IndexReader::HEADERS,
                                     (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }
    IndexReader *qcadbr = NULL;
    IndexReader *tcadbr = NULL;

    if(needQCA) {
        qcadbr = new IndexReader(
                par.db1,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca");
    }

    if(needTCA) {
        if (sameDB) {
            tcadbr = qcadbr;
        } else {
            tcadbr = new IndexReader(
                    par.db2,
                    par.threads,
                    isExtendedAlignment ? IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB2) :
                                           IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                    touch ? IndexReader::PRELOAD_INDEX : 0,
                    DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                    isExtendedAlignment ? "_seq_ca" : "_ca"
            );
        }
    }


    bool queryNucs = Parameters::isEqualDbtype(qDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    bool targetNucs = Parameters::isEqualDbtype(tDbr->sequenceReader->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    if (needSequenceDB) {
        // try to figure out if search was translated. This is can not be solved perfectly.
        bool seqtargetAA = false;
        if(Parameters::isEqualDbtype(tDbr->getDbtype(), Parameters::DBTYPE_INDEX_DB)){
            IndexReader tseqDbr(par.db2, par.threads, IndexReader::SEQUENCES, 0, IndexReader::PRELOAD_INDEX);
            seqtargetAA = Parameters::isEqualDbtype(tseqDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_AMINO_ACIDS);
        } else if(targetNucs == true && queryNucs == true && par.searchType == Parameters::SEARCH_TYPE_AUTO){
            Debug(Debug::WARNING) << "It is unclear from the input if a translated or nucleotide search was performed\n "
                                     "Please provide the parameter --search-type 2 (translated) or 3 (nucleotide)\n";
            EXIT(EXIT_FAILURE);
        } else if(par.searchType == Parameters::SEARCH_TYPE_TRANSLATED){
            seqtargetAA = true;
        }

        if((targetNucs == true && queryNucs == false )  || (targetNucs == false && queryNucs == true ) || (targetNucs == true && seqtargetAA == true && queryNucs == true )  ){
            isTranslatedSearch = true;
        }
    }

    int gapOpen, gapExtend;
    SubstitutionMatrix * subMat= NULL;
    if (targetNucs == true && queryNucs == true && isTranslatedSearch == false) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
        gapOpen = par.gapOpen.values.nucleotide();
        gapExtend =  par.gapExtend.values.nucleotide();
    }else{
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
        gapOpen = par.gapOpen.values.aminoacid();
        gapExtend = par.gapExtend.values.aminoacid();
    }
    EvalueComputation *evaluer = NULL;
    bool queryProfile = false;
    bool targetProfile = false;
    if (needSequenceDB) {
        queryProfile = Parameters::isEqualDbtype(qDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_HMM_PROFILE);
        targetProfile = Parameters::isEqualDbtype(tDbr->sequenceReader->getDbtype(), Parameters::DBTYPE_HMM_PROFILE);
        evaluer = new EvalueComputation(tDbr->sequenceReader->getAminoAcidDBSize(), subMat, gapOpen, gapExtend);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads, shouldCompress, dbType);
    resultWriter.open();

    const bool isDb = par.dbOut;
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));

    if (format == Parameters::FORMAT_ALIGNMENT_SAM) {
        char buffer[1024];
        unsigned int lastKey = tDbr->sequenceReader->getLastKey();
        bool *headerWritten = new bool[lastKey + 1];
        memset(headerWritten, 0, sizeof(bool) * (lastKey + 1));
        resultWriter.writeStart(0);
        std::string header = "@HD\tVN:1.4\tSO:queryname\n";
        resultWriter.writeAdd(header.c_str(), header.size(), 0);

        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            char *data = alnDbr.getData(i, 0);
            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                if (headerWritten[dbKey] == false) {
                    headerWritten[dbKey] = true;
                    unsigned int tId = tDbr->sequenceReader->getId(dbKey);
                    unsigned int seqLen = tDbr->sequenceReader->getSeqLen(tId);
                    unsigned int tHeaderId = tDbrHeader->sequenceReader->getId(dbKey);
                    const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, 0);
                    std::string targetId = Util::parseFastaHeader(tHeader);
                    int count = snprintf(buffer, sizeof(buffer), "@SQ\tSN:%s\tLN:%d\n", targetId.c_str(),
                                         (int32_t) seqLen);
                    if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                        Debug(Debug::WARNING) << "Truncated line in header " << i << "!\n";
                        continue;
                    }
                    resultWriter.writeAdd(buffer, count, 0);
                }
                resultWriter.writeEnd(0, 0, false, 0);
                data = Util::skipLine(data);
            }
        }
        delete[] headerWritten;
    } else if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
        // size_t dstSize = ZSTD_findDecompressedSize(result_viz_prelude_fs_html_zst, result_viz_prelude_fs_html_zst_len);
        // char* dst = (char*)malloc(sizeof(char) * dstSize);
        // size_t realSize = ZSTD_decompress(dst, dstSize, result_viz_prelude_fs_html_zst, result_viz_prelude_fs_html_zst_len);
        
        size_t dstSize = ZSTD_findDecompressedSize(vendor_js_zst, vendor_js_zst_len);
        char* dst = (char*)malloc(sizeof(char) * dstSize);
        size_t realSize = ZSTD_decompress(dst, dstSize, vendor_js_zst, vendor_js_zst_len);
        
        std::string mainJS(const_cast<char *>(reinterpret_cast<const char *>(main_js)), main_js_len);
        std::string htmlTemplate(
R"html(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta http-equiv="x-ua-compatible" content="ie=edge">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title><%= STRINGS.APP_NAME %> Search Server</title>
</head>
<div id="app"></div>
)html");
        
        resultWriter.writeData(htmlTemplate.c_str(), htmlTemplate.size(), 0, 0, false, false);

        std::string scriptStart = "<script>";
        std::string scriptEnd   = "</script>";

        // vendor.js
        resultWriter.writeData(scriptStart.c_str(), scriptStart.size(), 0, 0, false, false);
        resultWriter.writeData(dst, realSize, 0, 0, false, false);
        resultWriter.writeData(scriptEnd.c_str(), scriptEnd.size(), 0, 0, false, false);
        
        // main.js
        resultWriter.writeData(scriptStart.c_str(), scriptStart.size(), 0, 0, false, false);
        resultWriter.writeData(mainJS.c_str(), mainJS.size(), 0, 0, false, false);
        resultWriter.writeData(scriptEnd.c_str(), scriptEnd.size(), 0, 0, false, false);
        
        // Data <div>
        const char* dataStart = "<div id=\"data\" style=\"display: none;\">\n[";
        resultWriter.writeData(dataStart, strlen(dataStart), 0, 0, false, false);

        free(dst);
    } else if (addColumnHeaders == true && outcodes.empty() == false) {
        std::vector<std::string> outfmt = Util::split(par.outfmt, ",");
        std::string header(outfmt[0]);
        for(size_t i = 1; i < outfmt.size(); i++) {
            header.append(1, '\t');
            header.append(outfmt[i]);
        }
        header.append(1, '\n');
        resultWriter.writeData(header.c_str(), header.length(), 0, 0, false, false);
    }

    Debug::Progress progress(alnDbr.getSize());
    std::vector<std::string> qChains;
    std::vector<std::string> tChains;
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[1024];

        TMaligner *tmaligner = NULL;
        if(needTMaligner) {
            tmaligner = new TMaligner(
                    std::max(tDbr->sequenceReader->getMaxSeqLen() + 1, qDbr.sequenceReader->getMaxSeqLen() + 1), false, true);
        }
        LDDTCalculator *lddtcalculator = NULL;
        if(needLDDT) {
            lddtcalculator = new LDDTCalculator(qDbr.sequenceReader->getMaxSeqLen() + 1, tDbr->sequenceReader->getMaxSeqLen() + 1);
        }

        std::string result;
        result.reserve(1024*1024);

        std::string caStr;
        caStr.reserve(1024*1024);

        std::string queryProfData;
        queryProfData.reserve(1024);

        std::string queryBuffer;
        queryBuffer.reserve(1024);

        std::string queryHeaderBuffer;
        queryHeaderBuffer.reserve(1024);

        std::string targetProfData;
        targetProfData.reserve(1024);

        std::string newBacktrace;
        newBacktrace.reserve(1024);

        const TaxonNode * taxonNode = NULL;
        TMaligner::TMscoreResult tmres;

        Coordinate16 qcoords;
        Coordinate16 tcoords;

        std::string tmpBt;
        double rmsd = 0.0;
#pragma omp  for schedule(dynamic, 10)
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            progress.updateProgress();

            const unsigned int queryKey = alnDbr.getDbKey(i);
            char *querySeqData = NULL;
            size_t querySeqLen = 0;
            queryProfData.clear();
            if (needSequenceDB) {
                size_t qId = qDbr.sequenceReader->getId(queryKey);
                querySeqData = qDbr.sequenceReader->getData(qId, thread_idx);
                querySeqLen = qDbr.sequenceReader->getSeqLen(qId);
                if(sameDB && qDbr.sequenceReader->isCompressed()){
                    queryBuffer.assign(querySeqData, querySeqLen);
                    querySeqData = (char*) queryBuffer.c_str();
                }
                if (queryProfile) {
                    size_t queryEntryLen = qDbr.sequenceReader->getEntryLen(qId);
                    Sequence::extractProfileConsensus(querySeqData, queryEntryLen, *subMat, queryProfData);
                }
            }
            float *queryCaData = NULL;
            if (needQCA) {
                size_t qSeqId = qDbr.sequenceReader->getId(queryKey);
		        querySeqLen = qDbr.sequenceReader->getSeqLen(qSeqId);
                size_t qId = qcadbr->sequenceReader->getId(queryKey);
                char *qcadata = qcadbr->sequenceReader->getData(qId, thread_idx);
                size_t qCaLength = qcadbr->sequenceReader->getEntryLen(qId);
                queryCaData = qcoords.read(qcadata, querySeqLen, qCaLength);
            }
            size_t qHeaderId = qDbrHeader.sequenceReader->getId(queryKey);
            const char *qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
            size_t qHeaderLen = qDbrHeader.sequenceReader->getSeqLen(qHeaderId);
            std::string queryId = Util::parseFastaHeader(qHeader);
            if (sameDB && needFullHeaders) {
                queryHeaderBuffer.assign(qHeader, qHeaderLen);
                qHeader = (char*) queryHeaderBuffer.c_str();
            }
            if(needLDDT){
	            lddtcalculator->initQuery(querySeqLen, queryCaData, &queryCaData[querySeqLen], &queryCaData[querySeqLen+querySeqLen]);
            }

            if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
                const char* jsStart = "{\"query\": {\"header\": \"%s\",\"sequence\": \"";
                int count = snprintf(buffer, sizeof(buffer), jsStart, queryId.c_str(), querySeqData);
                if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                    Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                    continue;
                }
                result.append(buffer, count);
                if (queryProfile) {
                    result.append(queryProfData);
                } else {
                    result.append(querySeqData, querySeqLen);
                }
                result.append("\", \"qCa\": \"");
                caStr.clear();
                caToStr(queryCaData, querySeqLen, caStr);
                result.append(caStr, 0, caStr.size()-1);
                result.append("\"}, \"results\": [\n{\"db\": \"");
                result.append(par.db2);
                result.append("\", \"alignments\": [");
            }
            char *data = alnDbr.getData(i, thread_idx);
            Matcher::result_t res;
            while (*data != '\0') {
                const char *entry[255];
                Util::getWordsOfLine(data, entry, 255);
                ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                if (!retComplex.isValid) {
                    res = Matcher::parseAlignmentRecord(data, true);
                }
                data = Util::skipLine(data);
                if (res.backtrace.empty() && needBacktrace == true) {
                    Debug(Debug::ERROR) << "Backtrace cigar is missing in the alignment result. Please recompute the alignment with the -a flag.\n"
                                           "Command: mmseqs align " << par.db1 << " " << par.db2 << " " << par.db3 << " " << "alnNew -a\n";
                    EXIT(EXIT_FAILURE);
                }
                size_t tHeaderId = tDbrHeader->sequenceReader->getId(res.dbKey);
                const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, thread_idx);
                size_t tHeaderLen = tDbrHeader->sequenceReader->getSeqLen(tHeaderId);
                float *targetCaData = NULL;
                if (needTCA) {
                    size_t tId = tcadbr->sequenceReader->getId(res.dbKey);
                    char *tcadata = tcadbr->sequenceReader->getData(tId, thread_idx);
                    size_t tCaLength = tcadbr->sequenceReader->getEntryLen(tId);
                    targetCaData = tcoords.read(tcadata, res.dbLen, tCaLength);
                }

                std::string targetId = Util::parseFastaHeader(tHeader);

                unsigned int gapOpenCount = 0;
                unsigned int alnLen = res.alnLength;
                unsigned int missMatchCount = 0;
                unsigned int identical = 0;
                if (res.backtrace.empty() == false) {
                    size_t matchCount = 0;
                    alnLen = 0;
                    for (size_t pos = 0; pos < res.backtrace.size(); pos++) {
                        int cnt = 0;
                        if (isdigit(res.backtrace[pos])) {
                            cnt += Util::fast_atoi<int>(res.backtrace.c_str() + pos);
                            while (isdigit(res.backtrace[pos])) {
                                pos++;
                            }
                        }
                        alnLen += cnt;

                        switch (res.backtrace[pos]) {
                            case 'M':
                                matchCount += cnt;
                                break;
                            case 'D':
                            case 'I':
                                gapOpenCount += 1;
                                break;
                        }
                    }
//                res.seqId = X / alnLength;
                    identical = static_cast<unsigned int>(res.seqId * static_cast<float>(alnLen) + 0.5);
                    //res.alnLength = alnLength;
                    missMatchCount = static_cast<unsigned int>( matchCount - identical);
                } else {
                    const int adjustQstart = (res.qStartPos == -1) ? 0 : res.qStartPos;
                    const int adjustDBstart = (res.dbStartPos == -1) ? 0 : res.dbStartPos;
                    const float bestMatchEstimate = static_cast<float>(std::min(abs(res.qEndPos - adjustQstart), abs(res.dbEndPos - adjustDBstart)));
                    missMatchCount = static_cast<unsigned int>(bestMatchEstimate * (1.0f - res.seqId) + 0.5);
                }
                if(needTMaligner){
                    tmaligner->initQuery(queryCaData, &queryCaData[res.qLen], &queryCaData[res.qLen+res.qLen], NULL, res.qLen);
                    tmpBt = Matcher::uncompressAlignment(res.backtrace);
                    tmres = tmaligner->computeTMscore(targetCaData, &targetCaData[res.dbLen], &targetCaData[res.dbLen+res.dbLen], res.dbLen,
                                                      res.qStartPos, res.dbStartPos, tmpBt,
                                                      std::min(std::min(res.dbLen, res.qLen), static_cast<unsigned int>(tmpBt.size())));
                    rmsd = tmres.rmsd;

                }
                LDDTCalculator::LDDTScoreResult lddtres;
                if(needLDDT) {
                    lddtres = lddtcalculator->computeLDDTScore(res.dbLen, res.qStartPos, res.dbStartPos, Matcher::uncompressAlignment(res.backtrace), targetCaData, &targetCaData[res.dbLen], &targetCaData[res.dbLen+res.dbLen]);
                }
                switch (format) {
                    case Parameters::FORMAT_ALIGNMENT_BLAST_TAB: {
                        if (outcodes.empty()) {
                            int count = snprintf(buffer, sizeof(buffer),
                                                 "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                                                 queryId.c_str(), targetId.c_str(), res.seqId, alnLen,
                                                 missMatchCount, gapOpenCount,
                                                 res.qStartPos + 1, res.qEndPos + 1,
                                                 res.dbStartPos + 1, res.dbEndPos + 1,
                                                 res.eval, res.score);
                            if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                                Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                                continue;
                            }
                            result.append(buffer, count);
                        } else {
                            char *targetSeqData = NULL;
                            targetProfData.clear();
                            unsigned int taxon = 0;
                            if (needTaxonomy || needTaxonomyMapping) {
                                taxon = mapping->lookup(res.dbKey);
                                if (taxon == 0) {
                                    taxonNode = NULL;
                                } else if (needTaxonomy) {
                                    taxonNode = t->taxonNode(taxon, false);
                                }
                            }

                            if (needSequenceDB) {
                                size_t tId = tDbr->sequenceReader->getId(res.dbKey);
                                targetSeqData = tDbr->sequenceReader->getData(tId, thread_idx);
                                if (targetProfile) {
                                    size_t targetEntryLen = tDbr->sequenceReader->getEntryLen(tId);
                                    Sequence::extractProfileConsensus(targetSeqData, targetEntryLen, *subMat, targetProfData);
                                }
                            }
                            for(size_t i = 0; i < outcodes.size(); i++) {
                                switch (outcodes[i]) {
                                    case Parameters::OUTFMT_QUERY:
                                        result.append(queryId);
                                        break;
                                    case Parameters::OUTFMT_TARGET:
                                        result.append(targetId);
                                        break;
                                    case Parameters::OUTFMT_EVALUE:
                                        result.append(SSTR(res.eval));
                                        break;
                                    case Parameters::OUTFMT_GAPOPEN:
                                        result.append(SSTR(gapOpenCount));
                                        break;
                                    case Parameters::OUTFMT_FIDENT:
                                        result.append(SSTR(res.seqId));
                                        break;
                                    case Parameters::OUTFMT_PIDENT:
                                        result.append(SSTR(res.seqId*100));
                                        break;
                                    case Parameters::OUTFMT_NIDENT:
                                        result.append(SSTR(identical));
                                        break;
                                    case Parameters::OUTFMT_QSTART:
                                        result.append(SSTR(res.qStartPos + 1));
                                        break;
                                    case Parameters::OUTFMT_QEND:
                                        result.append(SSTR(res.qEndPos + 1));
                                        break;
                                    case Parameters::OUTFMT_QLEN:
                                        result.append(SSTR(res.qLen));
                                        break;
                                    case Parameters::OUTFMT_TSTART:
                                        result.append(SSTR(res.dbStartPos + 1));
                                        break;
                                    case Parameters::OUTFMT_TEND:
                                        result.append(SSTR(res.dbEndPos + 1));
                                        break;
                                    case Parameters::OUTFMT_TLEN:
                                        result.append(SSTR(res.dbLen));
                                        break;
                                    case Parameters::OUTFMT_ALNLEN:
                                        result.append(SSTR(alnLen));
                                        break;
                                    case Parameters::OUTFMT_RAW:
                                        result.append(SSTR(static_cast<int>(evaluer->computeRawScoreFromBitScore(res.score) + 0.5)));
                                        break;
                                    case Parameters::OUTFMT_BITS:
                                        result.append(SSTR(res.score));
                                        break;
                                    case Parameters::OUTFMT_CIGAR:
                                        if(isTranslatedSearch == true && targetNucs == true && queryNucs == true ){
                                            Matcher::result_t::protein2nucl(res.backtrace, newBacktrace);
                                            res.backtrace = newBacktrace;
                                        }
                                        result.append(SSTR(res.backtrace));
                                        newBacktrace.clear();
                                        break;
                                    case Parameters::OUTFMT_QSEQ:
                                        if (queryProfile) {
                                            result.append(queryProfData.c_str(), res.qLen);
                                        } else {
                                            result.append(querySeqData, res.qLen);
                                        }
                                        break;
                                    case Parameters::OUTFMT_TSEQ:
                                        if (targetProfile) {
                                            result.append(targetProfData.c_str(), res.dbLen);
                                        } else {
                                            result.append(targetSeqData, res.dbLen);
                                        }
                                        break;
                                    case Parameters::OUTFMT_QHEADER:
                                        result.append(qHeader, qHeaderLen);
                                        break;
                                    case Parameters::OUTFMT_THEADER:
                                        result.append(tHeader, tHeaderLen);
                                        break;
                                    case Parameters::OUTFMT_QALN:
                                        if (queryProfile) {
                                            structurePrintSeqBasedOnAln(result, queryProfData.c_str(), res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                        } else {
                                            structurePrintSeqBasedOnAln(result, querySeqData, res.qStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                                        }
                                        break;
                                    case Parameters::OUTFMT_TALN: {
                                        if (targetProfile) {
                                            structurePrintSeqBasedOnAln(result, targetProfData.c_str(), res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                        } else {
                                            structurePrintSeqBasedOnAln(result, targetSeqData, res.dbStartPos,
                                                               Matcher::uncompressAlignment(res.backtrace), true,
                                                               (res.dbStartPos > res.dbEndPos),
                                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                                        }
                                        break;
                                    }
                                    case Parameters::OUTFMT_MISMATCH:
                                        result.append(SSTR(missMatchCount));
                                        break;
                                    case Parameters::OUTFMT_QCOV:
                                        result.append(SSTR(res.qcov));
                                        break;
                                    case Parameters::OUTFMT_TCOV:
                                        result.append(SSTR(res.dbcov));
                                        break;
                                    case Parameters::OUTFMT_QSET:
                                        result.append(SSTR(qSetToSource[qKeyToSet[queryKey]]));
                                        break;
                                    case Parameters::OUTFMT_QSETID:
                                        result.append(SSTR(qKeyToSet[queryKey]));
                                        break;
                                    case Parameters::OUTFMT_TSET:
                                        result.append(SSTR(tSetToSource[tKeyToSet[res.dbKey]]));
                                        break;
                                    case Parameters::OUTFMT_TSETID:
                                        result.append(SSTR(tKeyToSet[res.dbKey]));
                                        break;
                                    case Parameters::OUTFMT_TAXID:
                                        result.append(SSTR(taxon));
                                        break;
                                    case Parameters::OUTFMT_TAXNAME:
                                        result.append((taxonNode != NULL) ? t->getString(taxonNode->nameIdx) : "unclassified");
                                        break;
                                    case Parameters::OUTFMT_TAXLIN:
                                        result.append((taxonNode != NULL) ? t->taxLineage(taxonNode, true) : "unclassified");
                                        break;
                                    case Parameters::OUTFMT_EMPTY:
                                        result.push_back('-');
                                        break;
                                    case Parameters::OUTFMT_QORFSTART:
                                        result.append(SSTR(res.queryOrfStartPos));
                                        break;
                                    case Parameters::OUTFMT_QORFEND:
                                        result.append(SSTR(res.queryOrfEndPos));
                                        break;
                                    case Parameters::OUTFMT_TORFSTART:
                                        result.append(SSTR(res.dbOrfStartPos));
                                        break;
                                    case Parameters::OUTFMT_TORFEND:
                                        result.append(SSTR(res.dbOrfEndPos));
                                        break;
                                    case LocalParameters::OUTFMT_QCA:
                                        caStr.clear();
                                        caToStr(queryCaData, res.qLen, caStr);
                                        result.append(caStr, 0, caStr.size()-1);
                                        break;
                                    case LocalParameters::OUTFMT_TCA:
                                        caStr.clear();
                                        caToStr(targetCaData, res.dbLen, caStr);
                                        result.append(caStr, 0, caStr.size()-1);
                                        break;
                                    case LocalParameters::OUTFMT_U:
                                        result.append(SSTR(tmres.u[0][0]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.u[0][1]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.u[0][2]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.u[1][0]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.u[1][1]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.u[1][2]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.u[2][0]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.u[2][1]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.u[2][2]));
                                        break;
                                    case LocalParameters::OUTFMT_T:
                                        result.append(SSTR(tmres.t[0]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.t[1]));
                                        result.push_back(',');
                                        result.append(SSTR(tmres.t[2]));
                                        break;
                                    case LocalParameters::OUTFMT_ALNTMSCORE:
                                        result.append(SSTR(tmres.tmscore));
                                        break;
                                    case LocalParameters::OUTFMT_QTMSCORE:
                                        tmres = tmaligner->computeTMscore(targetCaData, &targetCaData[res.dbLen], &targetCaData[res.dbLen+res.dbLen], res.dbLen,
                                                                          res.qStartPos, res.dbStartPos, Matcher::uncompressAlignment(res.backtrace),
                                                                          res.qLen);
                                        result.append(SSTR(tmres.tmscore));
                                        break;
                                    case LocalParameters::OUTFMT_TTMSCORE:
                                        tmres = tmaligner->computeTMscore(targetCaData, &targetCaData[res.dbLen], &targetCaData[res.dbLen+res.dbLen], res.dbLen,
                                                                          res.qStartPos, res.dbStartPos, Matcher::uncompressAlignment(res.backtrace),
                                                                          res.dbLen);
                                        result.append(SSTR(tmres.tmscore));
                                        break;
                                    case LocalParameters::OUTFMT_RMSD:
                                        result.append(SSTR(rmsd));
                                        break;
                                    case LocalParameters::OUTFMT_LDDT:
                                        // TODO: make SSTR_approx that outputs %2f, not %3f
                                        result.append(SSTR(lddtres.avgLddtScore));
                                        break;
                                    case LocalParameters::OUTFMT_LDDT_FULL:
                                        for(int i = 0; i < lddtres.scoreLength - 1; i++) {
                                            result.append(SSTR(lddtres.perCaLddtScore[i]));
                                            result.push_back(',');
                                        }
                                        result.append(SSTR(lddtres.perCaLddtScore[lddtres.scoreLength - 1]));
                                        break;
                                    case LocalParameters::OUTFMT_PROBTP:
                                        result.append(SSTR(CalcProbTP::calculate(res.score)));
                                        break;
                                    case LocalParameters::OUTFMT_Q_COMPLEX_TMSCORE:
                                        if (!retComplex.isValid) {
                                            Debug(Debug::ERROR) << "The column qcomplextmscore is only for scorecomplex result.\n";
                                            EXIT(EXIT_FAILURE);
                                        }
                                        result.append(SSTR(retComplex.qTmScore));
                                        break;
                                    case LocalParameters::OUTFMT_T_COMPLEX_TMSCORE:
                                        if (!retComplex.isValid) {
                                            Debug(Debug::ERROR) << "The column tcomplextmscore is only for scorecomplex result.\n";
                                            EXIT(EXIT_FAILURE);
                                        }
                                        result.append(SSTR(retComplex.tTmScore));
                                        break;
                                    case LocalParameters::OUTFMT_ASSIGN_ID:
                                        if (!retComplex.isValid) {
                                            Debug(Debug::ERROR) << "The column assignid is only for scorecomplex result.\n";
                                            EXIT(EXIT_FAILURE);
                                        }
                                        result.append(SSTR(retComplex.assId));
                                        break;
//                                    case LocalParameters::OUTFMT_COMPLEX_T:
//                                        if (!retComplex.isValid) {
//                                            Debug(Debug::ERROR) << "The column complext is only for scorecomplex result.\n";
//                                            EXIT(EXIT_FAILURE);
//                                        }
//                                        result.append(retComplex.tString);
//                                        break;
//                                    case LocalParameters::OUTFMT_COMPLEX_U:
//                                        if (!retComplex.isValid) {
//                                            Debug(Debug::ERROR) << "The column complexu is only for scorecomplex result.\n";
//                                            EXIT(EXIT_FAILURE);
//                                        }
//                                        result.append(retComplex.uString);
//                                        break;
                                }
                                if (i < outcodes.size() - 1) {
                                    result.push_back('\t');
                                }
                            }
                            result.push_back('\n');
                        }
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_BLAST_WITH_LEN: {
                        int count = snprintf(buffer, sizeof(buffer),
                                             "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%d\t%d\t%d\n",
                                             queryId.c_str(), targetId.c_str(), res.seqId, alnLen,
                                             missMatchCount, gapOpenCount,
                                             res.qStartPos + 1, res.qEndPos + 1,
                                             res.dbStartPos + 1, res.dbEndPos + 1,
                                             res.eval, res.score,
                                             res.qLen, res.dbLen);

                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                            continue;
                        }

                        result.append(buffer, count);
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_SAM: {
                        bool strand = res.qEndPos > res.qStartPos;
                        int rawScore = static_cast<int>(evaluer->computeRawScoreFromBitScore(res.score) + 0.5);
                        uint32_t mapq = -4.343 * log(exp(static_cast<double>(-rawScore)));
                        mapq = (uint32_t) (mapq + 4.99);
                        mapq = mapq < 254 ? mapq : 254;
                        int count = snprintf(buffer, sizeof(buffer), "%s\t%d\t%s\t%d\t%d\t",  queryId.c_str(), (strand) ? 16: 0, targetId.c_str(), res.dbStartPos + 1, mapq);
                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                            continue;
                        }
                        result.append(buffer, count);
                        if (isTranslatedSearch == true && targetNucs == true && queryNucs == true) {
                            Matcher::result_t::protein2nucl(res.backtrace, newBacktrace);
                            result.append(newBacktrace);
                            newBacktrace.clear();

                        } else {
                            result.append(res.backtrace);
                        }
                        result.append("\t*\t0\t0\t");
                        int start = std::min(res.qStartPos, res.qEndPos);
                        int end   = std::max(res.qStartPos, res.qEndPos);
                        if (queryProfile) {
                            result.append(queryProfData.c_str() + start, (end + 1) - start);
                        } else {
                            result.append(querySeqData + start, (end + 1) - start);
                        }
                        count = snprintf(buffer, sizeof(buffer), "\t*\tAS:i:%d\tNM:i:%d\n", rawScore, missMatchCount);
                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                            continue;
                        }
                        result.append(buffer, count);
                        break;
                    }
                    case LocalParameters::FORMAT_ALIGNMENT_PDB_SUPERPOSED:{
                        // rotate and translate the target Calpha
                        // and write the results as pdb file
                        std::string filename = resultWriter.getDataFileName() + queryId+"_"+targetId+".pdb";
                        FILE * fp = fopen(filename.c_str(), "w");
                        result.append("MODEL\n");
                        result.append("REMARK ");
                        result.append(queryId);
                        result.append(" ");
                        result.append(targetId);
                        result.append("\n");
                        for(unsigned int tpos = 0; tpos < res.dbLen; tpos++){
                            size_t tId = tDbr->sequenceReader->getId(res.dbKey);
                            char* targetSeqData  = (char*) tDbr->sequenceReader->getData(tId, thread_idx);
                            // printf for ATOM Calpha record pdb
                            int count = snprintf(buffer, sizeof(buffer),
                                   "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                   tpos+1, "CA", singleLetterToThree(targetSeqData[tpos]), "A", tpos+1,
                                   tmres.t[0] + targetCaData[tpos] * tmres.u[0][0] + targetCaData[res.dbLen+tpos] * tmres.u[0][1] + targetCaData[res.dbLen+res.dbLen+tpos] * tmres.u[0][2],
                                   tmres.t[1] + targetCaData[tpos] * tmres.u[1][0] + targetCaData[res.dbLen+tpos] * tmres.u[1][1] + targetCaData[res.dbLen+res.dbLen+tpos] * tmres.u[1][2],
                                   tmres.t[2] + targetCaData[tpos] * tmres.u[2][0] + targetCaData[res.dbLen+tpos] * tmres.u[2][1] + targetCaData[res.dbLen+res.dbLen+tpos] * tmres.u[2][2],
                                   1.0, 0.0);
                            if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                                Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                                continue;
                            }
                            result.append(buffer, count);
                        }
                        result.append("ENDMDL\n");
                        // use result size to write to fp
                        fwrite(result.c_str(), sizeof(char), result.size(), fp);
                        fclose(fp);
                        result.clear();
                        break;
                    }
                    case Parameters::FORMAT_ALIGNMENT_HTML: {
                        const char* jsAln = "{\"target\": \"%s\", \"prob\": %1.2f, \"seqId\": %1.3f, \"alnLength\": %d, \"mismatch\": %d, \"gapopen\": %d, \"qStartPos\": %d, \"qEndPos\": %d, \"dbStartPos\": %d, \"dbEndPos\": %d, \"eval\": %.2E, \"score\": %d, \"qLen\": %d, \"dbLen\": %d, \"qAln\": \"";
                        int count = snprintf(buffer, sizeof(buffer), jsAln,
                                             targetId.c_str(),
                                             CalcProbTP::calculate(res.score),
                                             res.seqId, alnLen,
                                             missMatchCount, gapOpenCount,
                                             res.qStartPos + 1, res.qEndPos + 1,
                                             res.dbStartPos + 1, res.dbEndPos + 1,
                                             res.eval, res.score,
                                             res.qLen, res.dbLen);
                        if (count < 0 || static_cast<size_t>(count) >= sizeof(buffer)) {
                            Debug(Debug::WARNING) << "Truncated line in entry" << i << "!\n";
                            continue;
                        }
                        result.append(buffer, count);
                        if (queryProfile) {
                            structurePrintSeqBasedOnAln(result, queryProfData.c_str(), res.qStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                        } else {
                            structurePrintSeqBasedOnAln(result, querySeqData, res.qStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), false, (res.qStartPos > res.qEndPos),
                                               (isTranslatedSearch == true && queryNucs == true), translateNucl);
                        }
                        result.append("\", \"dbAln\": \"");
                        size_t tId = tDbr->sequenceReader->getId(res.dbKey);
                        char* targetSeqData = tDbr->sequenceReader->getData(tId, thread_idx);
                        if (targetProfile) {
                            size_t targetEntryLen = tDbr->sequenceReader->getEntryLen(tId);
                            Sequence::extractProfileConsensus(targetSeqData, targetEntryLen, *subMat, targetProfData);
                            structurePrintSeqBasedOnAln(result, targetProfData.c_str(), res.dbStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), true,
                                               (res.dbStartPos > res.dbEndPos),
                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                        } else {
                            structurePrintSeqBasedOnAln(result, targetSeqData, res.dbStartPos,
                                               Matcher::uncompressAlignment(res.backtrace), true,
                                               (res.dbStartPos > res.dbEndPos),
                                               (isTranslatedSearch == true && targetNucs == true), translateNucl);
                        }
                        result.append("\", \"tCa\": \"");
                        caStr.clear();
                        caToStr(targetCaData, res.dbLen, caStr);
                        result.append(caStr, 0, caStr.size()-1);
                        
                        result.append("\", \"tSeq\": \"");
                        result.append(targetSeqData, 0, res.dbLen);

                        result.append("\" },\n");
                        break;
                    }
                    default:
                        Debug(Debug::ERROR) << "Not implemented yet";
                        EXIT(EXIT_FAILURE);
                }
            }

            if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
                result.erase(result.end() - 2); // remove ,\n of last entry
                result.append("]}]},\n");
            }
            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx, isDb);
            result.clear();
        }
        if(tmaligner != NULL){
            delete tmaligner;
        }
        if(lddtcalculator != NULL) {
            delete lddtcalculator;
        }
    }
    const char* htmlEndBlock = "]\n</div>";
    if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
        resultWriter.writeData(htmlEndBlock, strlen(htmlEndBlock), 0, localThreads - 1, false, false);
    }
    // tsv output
    resultWriter.close(true);

    if (format == Parameters::FORMAT_ALIGNMENT_HTML) {
        // replace last , to make this valid json
        FILE* handle = FileUtil::openFileOrDie(par.db4.c_str(), "r+b", true);
        // + newline + comma
        fseek(handle, -(strlen(htmlEndBlock) + 2), SEEK_END);
        const char space = ' ';
        fwrite(&space, 1, 1, handle);
        fclose(handle);
    }
    
    if (isDb == false) {
        FileUtil::remove(par.db4Index.c_str());
    }

    if (needQCA) {
        delete qcadbr;
    }
    if (needTCA && sameDB == false) {
        delete tcadbr;
    }

    if (needTaxonomy) {
        delete t;
    }
    if (mapping != NULL) {
        delete mapping;
    }
    alnDbr.close();
    if (sameDB == false) {
        delete tDbr;
        delete tDbrHeader;
    }
    if (needSequenceDB) {
        delete evaluer;
    }
    delete subMat;

    return EXIT_SUCCESS;
}



