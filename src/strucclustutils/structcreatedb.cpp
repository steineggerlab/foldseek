//
// Created by Martin Steinegger on 11/14/20.
//
#include "Command.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBWriter.h"
#include "FastSort.h"

#include "structureto3di.h"
#include "SubstitutionMatrix.h"
#include "GemmiWrapper.h"
#include "PulchraWrapper.h"
#include "microtar.h"
#include "PatternCompiler.h"

#include <iostream>
#include <dirent.h>

#ifdef HAVE_GCS
#include "google/cloud/storage/client.h"
#endif

#ifdef OPENMP
#include <omp.h>
#endif


#ifdef HAVE_ZLIB
#include <zlib.h>
static int structure_file_gzread(mtar_t *tar, void *data, size_t size) {
    size_t res = gzread((gzFile)tar->stream, data, size);
    return (res == size) ? MTAR_ESUCCESS : MTAR_EREADFAIL;
}

static int structure_file_gzseek(mtar_t *tar, long offset, int whence) {
    int res = gzseek((gzFile)tar->stream, offset, whence);
    return (res != -1) ? MTAR_ESUCCESS : MTAR_ESEEKFAIL;
}

static int structure_file_gzclose(mtar_t *tar) {
    gzclose((gzFile)tar->stream);
    return MTAR_ESUCCESS;
}

int structure_mtar_gzopen(mtar_t *tar, const char *filename) {
    // Init tar struct and functions
    memset(tar, 0, sizeof(*tar));
    tar->read = structure_file_gzread;
    tar->seek = structure_file_gzseek;
    tar->close = structure_file_gzclose;
    tar->isFinished = 0;
    // Open file
    tar->stream = gzopen(filename, "rb");
    if (!tar->stream) {
        return MTAR_EOPENFAIL;
    }

#if defined(ZLIB_VERNUM) && ZLIB_VERNUM >= 0x1240
    if (gzbuffer((gzFile)tar->stream, 1 * 1024 * 1024) != 0) {
        Debug(Debug::WARNING) << "Could not set gzbuffer size, performance might be bad\n";
    }
#endif

    return MTAR_ESUCCESS;
}
#endif

size_t
writeStructureEntry(SubstitutionMatrix & mat, GemmiWrapper & readStructure, StructureTo3Di & structureTo3Di,
                    PulchraWrapper & pulchra, std::vector<char> & alphabet3di, std::vector<char> & alphabetAA,
                    std::vector<float> & camol, std::string & header, std::string & name,
                    DBWriter & aadbw, DBWriter & hdbw, DBWriter & torsiondbw, DBWriter & cadbw, int chainNameMode,
                    float maskBfactorThreshold, size_t & tooShort, size_t &globalCnt, int thread_idx,
                    std::string & filename,  size_t  &fileidCnt,
                    std::map<std::string, size_t> & entrynameToFileId,
                    std::map<std::string, size_t> & filenameToFileId,
                    std::map<size_t, std::string> & fileIdToName) {
    size_t id = __sync_fetch_and_add(&globalCnt, readStructure.chain.size());
    size_t entriesAdded = 0;
    for(size_t ch = 0; ch < readStructure.chain.size(); ch++){
        size_t dbKey = id + ch;
        size_t chainStart = readStructure.chain[ch].first;
        size_t chainEnd = readStructure.chain[ch].second;
        size_t chainLen = chainEnd - chainStart;
        if(chainLen <= 3){
            tooShort++;
            continue;
        }

        // Detect if structure is Ca only
        if (std::isnan(readStructure.n[chainStart + 0].x) &&
            std::isnan(readStructure.n[chainStart + 1].x) &&
            std::isnan(readStructure.n[chainStart + 2].x) &&
            std::isnan(readStructure.n[chainStart + 3].x) &&
            std::isnan(readStructure.c[chainStart + 0].x) &&
            std::isnan(readStructure.c[chainStart + 1].x) &&
            std::isnan(readStructure.c[chainStart + 2].x) &&
            std::isnan(readStructure.c[chainStart + 3].x))
        {
            pulchra.rebuildBackbone(&readStructure.ca[chainStart],
                                    &readStructure.n[chainStart],
                                    &readStructure.c[chainStart],
                                    &readStructure.ami[chainStart],
                                    chainLen);
        }

        char * states = structureTo3Di.structure2states(&readStructure.ca[chainStart],
                                                        &readStructure.n[chainStart],
                                                        &readStructure.c[chainStart],
                                                        &readStructure.cb[chainStart],
                                                        chainLen);
        for(size_t pos = 0; pos < chainLen; pos++){
            if(readStructure.ca_bfactor[pos] < maskBfactorThreshold){
                alphabet3di.push_back(tolower(mat.num2aa[static_cast<int>(states[pos])]));
                alphabetAA.push_back(tolower(readStructure.ami[chainStart+pos]));
            }else{
                alphabet3di.push_back(mat.num2aa[static_cast<int>(states[pos])]);
                alphabetAA.push_back(readStructure.ami[chainStart+pos]);
            }
        }
        alphabet3di.push_back('\n');
        alphabetAA.push_back('\n');
        torsiondbw.writeData(alphabet3di.data(), alphabet3di.size(), dbKey, thread_idx);
        aadbw.writeData(alphabetAA.data(), alphabetAA.size(), dbKey, thread_idx);
        header.clear();
        header.append(readStructure.names[ch]);
        if(chainNameMode == LocalParameters::CHAIN_MODE_ADD ||
           (chainNameMode == LocalParameters::CHAIN_MODE_AUTO && readStructure.names.size() > 1)){
            header.push_back('_');
            header.append(readStructure.chainNames[ch]);
        }
        if(readStructure.title.size() > 0){
            header.push_back(' ');
            header.append(readStructure.title);
        }
        header.push_back('\n');
        std::string entryName = Util::parseFastaHeader(header.c_str());
#pragma omp critical
        {
            std::map<std::string, size_t>::iterator it = filenameToFileId.find(filename);
            size_t fileid;
            if(it != filenameToFileId.end()){
                fileid = it->second;
            }else{
                fileid = fileidCnt;
                filenameToFileId[filename] = fileid;
                fileIdToName[fileid] = filename;
                fileidCnt++;
            }
            entrynameToFileId[entryName] = fileid;
        }
        hdbw.writeData(header.c_str(), header.size(), dbKey, thread_idx);
        name.clear();
        for(size_t pos = 0; pos < chainLen; pos++){
            float val = (std::isnan(readStructure.ca[chainStart+pos].x))
                        ? 0.0 : readStructure.ca[chainStart+pos].x;
            camol.push_back(val);
        }
        for(size_t pos = 0; pos < chainLen; pos++) {
            float val = (std::isnan(readStructure.ca[chainStart+pos].y))
                        ? 0.0 : readStructure.ca[chainStart+pos].y;
            camol.push_back(val);
        }
        for(size_t pos = 0; pos < chainLen; pos++) {
            float val = (std::isnan(readStructure.ca[chainStart+pos].z))
                        ? 0.0 : readStructure.ca[chainStart+pos].z;
            camol.push_back(val);
        }
        cadbw.writeData((const char*)camol.data(), camol.size() * sizeof(float), dbKey, thread_idx);
        alphabet3di.clear();
        alphabetAA.clear();
        camol.clear();
        entriesAdded++;
    }
    return entriesAdded;
}

void sortDatafileByIdOrder(DBWriter & dbw,
                           DBReader<unsigned int> &dbr,  std::vector<std::pair<std::string, unsigned int>> & order) {
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

#pragma omp for schedule(static)
        for (size_t i = 0; i < dbr.getSize(); i++) {
            size_t id = order[i].second;
            char *data = dbr.getData(id, thread_idx);
            size_t length = dbr.getEntryLen(id);
            dbw.writeData(data, (length == 0 ? 0 : length - 1), dbr.getDbKey(id), thread_idx);
        }
    }
}



int createdb(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_COMMON);

    std::vector<std::string> filenames(par.filenames);
    std::string outputName = filenames.back();
    filenames.pop_back();
    if (filenames.size() == 1 && FileUtil::directoryExists(filenames.back().c_str())) {
        std::vector<std::string> dirs;
        dirs.push_back(filenames.back());
        filenames.pop_back();
        while(dirs.size() != 0){
            std::string dir = dirs.back();
            dirs.pop_back();
            DIR * dpdf = opendir(dir.c_str());
            if (dpdf != NULL) {
                while (dirent * epdf = readdir(dpdf)) {
                    std::string filename(epdf->d_name);
                    if(filename != "." && filename !=".."){
                        if (epdf->d_type == DT_DIR){
                            dirs.push_back(dir+"/"+filename);
                        }else{
                            filenames.push_back(dir+"/"+filename);
                        }
                    }
                }
                closedir(dpdf);
            }
        }
    }

    Debug(Debug::INFO) << "Output file: " << outputName << "\n";
    SORT_PARALLEL(filenames.begin(), filenames.end());

    DBWriter torsiondbw((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    torsiondbw.open();
    DBWriter hdbw((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
    hdbw.open();
    DBWriter cadbw((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
    cadbw.open();
    DBWriter aadbw((outputName).c_str(), (outputName+".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    aadbw.open();
    SubstitutionMatrix mat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    Debug::Progress progress(filenames.size());
    std::map<std::string, size_t> entrynameToFileId;
    std::map<size_t, std::string> fileIdToName;
    std::map<std::string, size_t> filenameToFileId;

    size_t globalCnt = 0;
    size_t globalFileidCnt = 0;
    size_t incorrectFiles = 0;
    size_t tooShort = 0;
    bool needsReorderingAtTheEnd = false;
    // Process tar files!
    for(size_t i = 0; i < filenames.size(); i++) {
        mtar_t tar;
        if (Util::endsWith(".tar.gz", filenames[i]) || Util::endsWith(".tgz", filenames[i])) {
#ifdef HAVE_ZLIB
            if (structure_mtar_gzopen(&tar, filenames[i].c_str()) != MTAR_ESUCCESS) {
                Debug(Debug::ERROR) << "Cannot open file " << filenames[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
#else
            Debug(Debug::ERROR) << "Foldseek was not compiled with zlib support. Cannot read compressed input.\n";
            EXIT(EXIT_FAILURE);
#endif
        } else if (Util::endsWith(".tar", filenames[i])) {
            if (mtar_open(&tar, filenames[i].c_str()) != MTAR_ESUCCESS) {
                Debug(Debug::ERROR) << "Cannot open file " << filenames[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
        } else {
            continue;
        }
        progress.updateProgress();
#ifdef OPENMP
        int localThreads = par.threads;
        if(localThreads > 1){
            needsReorderingAtTheEnd = true;
        }
#endif

#pragma omp parallel default(none) shared(tar, par, torsiondbw, hdbw, cadbw, aadbw, mat, filenames, progress, globalCnt, globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName) num_threads(localThreads) reduction(+:incorrectFiles, tooShort)
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#ifdef HAVE_ZLIB
            const unsigned int CHUNK = 128 * 1024;
            unsigned char in[CHUNK];
            unsigned char out[CHUNK];
            z_stream strm;
            memset(&strm, 0, sizeof(z_stream));
            strm.zalloc = Z_NULL;
            strm.zfree = Z_NULL;
            strm.opaque = Z_NULL;
            strm.next_in = in;
            strm.avail_in = 0;
            int status = inflateInit2(&strm, 15 | 32);
            if (status < 0) {
                Debug(Debug::ERROR) << "Cannot initialize zlib stream\n";
                //EXIT(EXIT_FAILURE);
            }
#endif
            //recon_related
            StructureTo3Di structureTo3Di;
            PulchraWrapper pulchra;
            GemmiWrapper readStructure;
            std::vector<char> alphabet3di;
            std::vector<char> alphabetAA;
            std::vector<float> camol;
            std::string header;
            std::string name;
            std::string pdbFile;
            mtar_header_t tarHeader;
            char buffer[4096];
            size_t bufferSize = 1024 * 1024;
            char *dataBuffer = (char *) malloc(bufferSize);
            size_t inflateSize = 1024 * 1024;
            char *inflateBuffer = (char *) malloc(inflateSize);
            bool proceed = true;
            PatternCompiler include(par.tarInclude.c_str());
            PatternCompiler exclude(par.tarExclude.c_str());
            while (proceed) {
                bool writeEntry = true;
#pragma omp critical
                {
                    if (tar.isFinished == 0 && (mtar_read_header(&tar, &tarHeader)) != MTAR_ENULLRECORD) {
                        // GNU tar has special blocks for long filenames
                        if (tarHeader.type == MTAR_TGNU_LONGNAME || tarHeader.type == MTAR_TGNU_LONGLINK) {
                            if (tarHeader.size > bufferSize) {
                                bufferSize = tarHeader.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, tarHeader.size) != MTAR_ESUCCESS) {
                                Debug(Debug::ERROR) << "Cannot read entry " << tarHeader.name << "\n";
                                //EXIT(EXIT_FAILURE);
                            }
                            name.assign(dataBuffer, tarHeader.size);
                            // skip to next record
                            if (mtar_read_header(&tar, &tarHeader) == MTAR_ENULLRECORD) {
                                Debug(Debug::ERROR) << "Tar truncated after entry " << name << "\n";
                                //EXIT(EXIT_FAILURE);
                            }
                        } else {
                            name = tarHeader.name;
                        }
                        if (tarHeader.type == MTAR_TREG || tarHeader.type == MTAR_TCONT ||
                            tarHeader.type == MTAR_TOLDREG) {
                            if (tarHeader.size > bufferSize) {
                                bufferSize = tarHeader.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, tarHeader.size) != MTAR_ESUCCESS) {
                                Debug(Debug::ERROR) << "Cannot read entry " << name << "\n";
                                //EXIT(EXIT_FAILURE);
                            }
                            if (include.isMatch(name.c_str()) == false || exclude.isMatch(name.c_str()) == true) {
                                proceed = true;
                                writeEntry = false;
                            } else {
                                proceed = true;
                                writeEntry = true;
                            }
                        } else {
                            proceed = true;
                            writeEntry = false;
                        }
                    } else {
                        tar.isFinished = 1;
                        proceed = false;
                        writeEntry = false;
                    }
                }

                if (writeEntry) {
                    pdbFile.clear();
                    if (Util::endsWith(".gz", name)) {
#ifdef HAVE_ZLIB
                        inflateReset(&strm);
                        strm.avail_in = tarHeader.size;
                        strm.next_in = (unsigned char *) dataBuffer;
                        do {
                            unsigned have;
                            strm.avail_out = CHUNK;
                            strm.next_out = out;
                            int err = inflate(&strm, Z_NO_FLUSH);
                            switch (err) {
                                case Z_OK:
                                case Z_STREAM_END:
                                case Z_BUF_ERROR:
                                    break;
                                default:
                                    inflateEnd(&strm);
                                    Debug(Debug::ERROR) << "Gzip error " << err << " entry " << name << "\n";
                                    //EXIT(EXIT_FAILURE);
                            }
                            have = CHUNK - strm.avail_out;
                            pdbFile.append((char *) out, have);
                        } while (strm.avail_out == 0);
#else
                        Debug(Debug::ERROR) << "MMseqs2 was not compiled with zlib support. Cannot read compressed input.\n";
                EXIT(EXIT_FAILURE);
#endif
                    } else {
                        pdbFile.append(buffer, bufferSize);
                    }
                    if (readStructure.loadFromBuffer(pdbFile.c_str(), pdbFile.size(), name) == false) {
                        incorrectFiles++;
                        continue;
                    }
                    writeStructureEntry(mat, readStructure, structureTo3Di, pulchra,
                                        alphabet3di, alphabetAA, camol, header, name, aadbw, hdbw, torsiondbw, cadbw,
                                        par.chainNameMode, par.maskBfactorThreshold, tooShort, globalCnt, thread_idx,
                                        name, globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName);
                }
            } // end while
            free(inflateBuffer);
            free(dataBuffer);
        } // end omp open
        mtar_close(&tar);
    } // end file for


    //===================== single_process ===================//__110710__//
#pragma omp parallel default(none) shared(par, torsiondbw, hdbw, cadbw, aadbw, mat, filenames, progress, globalCnt, globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName) reduction(+:incorrectFiles, tooShort)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        //recon_related
        StructureTo3Di structureTo3Di;
        PulchraWrapper pulchra;
        GemmiWrapper readStructure;
        std::vector<char> alphabet3di;
        std::vector<char> alphabetAA;
        std::vector<float> camol;
        std::string header;
        std::string name;


#pragma omp for schedule(static)
        for (size_t i = 0; i < filenames.size(); i++) {
            // skip tars since they were processed on top
            if (Util::endsWith(".tar.gz", filenames[i])
                || Util::endsWith(".tgz", filenames[i])
                || Util::endsWith(".tar", filenames[i])
#ifdef HAVE_GCS
                || Util::startWith("gcs://", filenames[i])
#endif
                ) {
                continue;
            }
            progress.updateProgress();

            if(readStructure.load(filenames[i]) == false){
                incorrectFiles++;
                continue;
            }
            // clear memory
            writeStructureEntry(mat, readStructure, structureTo3Di,  pulchra,
                                alphabet3di, alphabetAA, camol, header, name, aadbw, hdbw, torsiondbw, cadbw,
                                par.chainNameMode, par.maskBfactorThreshold, tooShort, globalCnt, thread_idx,
                                filenames[i], globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName);
        }

    }

#ifdef HAVE_GCS
    namespace gcs = ::google::cloud::storage;
    auto options = google::cloud::Options{}
        .set<gcs::ConnectionPoolSizeOption>(par.threads)
        .set<google::cloud::storage_experimental::HttpVersionOption>("2.0");
    auto client = gcs::Client(options);
    for (size_t i = 0; i < filenames.size(); i++) {
        if (Util::startWith("gcs://", filenames[i]) == false) {
            continue;
        }
        std::vector<std::string> parts = Util::split(filenames[i], "/");
        if (parts.size() == 1) {
            Debug(Debug::WARNING) << "Skipping invalid URI " << filenames[i] << "\n";
            continue;
        }
        std::string bucket_name = parts[1];

        char filter = '\0';
        if (parts.size >= 3) {
            filter = parts[2][0];
        }
        progress.reset(SIZE_MAX);
#pragma omp parallel default(none) shared(par, torsiondbw, hdbw, cadbw, aadbw, mat, filenames, progress, globalCnt, globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName, client, bucket_name, filter) reduction(+:incorrectFiles, tooShort)
        {
            StructureTo3Di structureTo3Di;
            PulchraWrapper pulchra;
            GemmiWrapper readStructure;
            std::vector<char> alphabet3di;
            std::vector<char> alphabetAA;
            std::vector<float> camol;
            std::string header;
            std::string name;

#pragma omp single
            for (auto&& object_metadata : client.ListObjects(bucket_name, gcs::Projection::NoAcl(), gcs::MaxResults(15000))) {
                std::string obj_name = object_metadata->name();
#pragma omp task firstprivate(obj_name, alphabet3di, alphabetAA, camol, header, name, filter) private(structureTo3Di, pulchra, readStructure)
                {
                    bool skipFilter = filter != '\0' && obj_name.length() >= 9 && obj_name[8] == filter;
                    bool allowedSuffix = Util::endsWith(".cif", obj_name) || Util::endsWith(".pdb", obj_name);
                    if (skipFilter && allowedSuffix) {
                        unsigned int thread_idx = 0;
#ifdef OPENMP
                        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
                        progress.updateProgress();

                        auto reader = client.ReadObject(bucket_name, obj_name);
                        if (!reader.status().ok()) {
                            Debug(Debug::ERROR) << reader.status().message() << "\n";
                        } else {
                            std::string contents{std::istreambuf_iterator<char>{reader}, {}};
                            if (readStructure.loadFromBuffer(contents.c_str(), contents.size(), obj_name) == false) {
                                incorrectFiles++;
                            } else {
                                writeStructureEntry(mat, readStructure, structureTo3Di,  pulchra,
                                        alphabet3di, alphabetAA, camol, header, name, aadbw, hdbw, torsiondbw, cadbw,
                                        par.chainNameMode, par.maskBfactorThreshold, tooShort, globalCnt, thread_idx,
                                        obj_name, globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName);
                            }
                        }
                    }
                }
            }
        }
    }
#endif

    torsiondbw.close(true);
    hdbw.close(true);
    cadbw.close(true);
    aadbw.close(true);


    if(needsReorderingAtTheEnd) {
        Debug(Debug::INFO) << "Reordering by identifier\n";
        DBReader<unsigned int> header_reorder((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        header_reorder.open(DBReader<unsigned int>::NOSORT);
        header_reorder.readMmapedDataInMemory();
        std::vector<std::pair<std::string, unsigned int>> mappingOrder(header_reorder.getSize());
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(static)
            for (size_t i = 0; i < header_reorder.getSize(); i++) {
                mappingOrder[i].first = Util::parseFastaHeader(header_reorder.getData(i,thread_idx));
                mappingOrder[i].second = i;
            }
        }
        SORT_PARALLEL(mappingOrder.begin(), mappingOrder.end());
        std::string lookupFile = outputName + ".lookup";
        FILE* file = FileUtil::openAndDelete(lookupFile.c_str(), "w");
        std::string buffer;
        buffer.reserve(2048);
        for (unsigned int id = 0; id < header_reorder.getSize(); id++) {
            DBReader<unsigned int>::LookupEntry entry;
            entry.id = id;
            entry.entryName = mappingOrder[id].first;
            entry.fileNumber = 0;
            header_reorder.lookupEntryToBuffer(buffer, entry);
            size_t written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), file);
            if (written != buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffer.clear();
        }
        fclose(file);

        DBWriter hdbw_reorder((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
        hdbw_reorder.open();
        sortDatafileByIdOrder(hdbw_reorder, header_reorder, mappingOrder);
        hdbw_reorder.close(true);
        header_reorder.close();

        DBReader<unsigned int> torsiondbr_reorder((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        torsiondbr_reorder.open(DBReader<unsigned int>::NOSORT);
        torsiondbr_reorder.readMmapedDataInMemory();
        DBWriter torsiondbw_reorder((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
        torsiondbw_reorder.open();
        sortDatafileByIdOrder(torsiondbw_reorder, torsiondbr_reorder, mappingOrder);
        torsiondbw_reorder.close(true);
        torsiondbr_reorder.close();

        DBReader<unsigned int> cadbr_reorder((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        cadbr_reorder.open(DBReader<unsigned int>::NOSORT);
        cadbr_reorder.readMmapedDataInMemory();
        DBWriter cadbw_reorder((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
        cadbw_reorder.open();
        sortDatafileByIdOrder(cadbw_reorder, cadbr_reorder, mappingOrder);
        cadbw_reorder.close(true);
        cadbr_reorder.close();

        DBReader<unsigned int> aadbr_reorder((outputName).c_str(), (outputName+".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        aadbr_reorder.open(DBReader<unsigned int>::NOSORT);
        aadbr_reorder.readMmapedDataInMemory();
        DBWriter aadbw_reorder((outputName).c_str(), (outputName+".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
        aadbw_reorder.open();
        sortDatafileByIdOrder(aadbw_reorder, aadbr_reorder, mappingOrder);
        aadbw_reorder.close(true);
        aadbr_reorder.close();
    }else{
        DBWriter::createRenumberedDB((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
        DBWriter::createRenumberedDB((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
        DBWriter::createRenumberedDB((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
        DBWriter::createRenumberedDB((outputName).c_str(), (outputName+".index").c_str(), "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
    }

    if (par.writeLookup == true) {
        DBReader<unsigned int> readerHeader((outputName + "_h").c_str(), (outputName + "_h.index").c_str(),
                                            1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        readerHeader.open(DBReader<unsigned int>::NOSORT);
        // create lookup file
        std::string sourceFilename = outputName + ".source";
        FILE *sourceFile = FileUtil::openAndDelete(sourceFilename.c_str(), "w");
        std::string lookupFile = outputName + ".lookup";
        FILE *file = FileUtil::openAndDelete(lookupFile.c_str(), "w");
        std::string buffer;
        size_t maxFileId = 0;
        buffer.reserve(2048);
        DBReader<unsigned int>::LookupEntry entry;
        for (unsigned int id = 0; id < readerHeader.getSize(); id++) {
            char *header = readerHeader.getData(id, 0);
            entry.id = readerHeader.getDbKey(id);
            entry.entryName = Util::parseFastaHeader(header);
            size_t fileId = entrynameToFileId[entry.entryName];
            entry.fileNumber = fileId;
            maxFileId = std::max(fileId, maxFileId);
            readerHeader.lookupEntryToBuffer(buffer, entry);
            size_t written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), file);
            if (written != buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffer.clear();
        }
        readerHeader.close();
        for(unsigned int fileId = 0; fileId <= maxFileId; fileId++) {
            std::string filename = fileIdToName[fileId];
            buffer.append(SSTR(fileId));
            buffer.push_back('\t');
            buffer.append(FileUtil::baseName(fileIdToName[fileId]));
            buffer.push_back('\n');
            size_t written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), sourceFile);
            if (written != buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffer.clear();
        }
    }


    Debug(Debug::INFO) << "Ignore " << tooShort+incorrectFiles << " out of " << globalCnt << ".\n";
    Debug(Debug::INFO) << "Too short: " << tooShort << ", incorrect  " << incorrectFiles << ".\n";
    return EXIT_SUCCESS;
}
