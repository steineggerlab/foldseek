/**
 * File: input_processor.h
 * Created: 2023-02-10 17:04:08
 * Author: Milot Mirdita (milot@mirdita.de)
 */

#pragma once

#include "microtar.h"
#include "utility.h"
#include "database_reader.h"

#include <utility>
#include <functional>
#include <vector>
#include <string>

#ifdef HAVE_GCS
#include "google/cloud/storage/client.h"
#endif

// #ifdef HAVE_AWS_S3
// #include <aws/core/Aws.h>
// #include <aws/s3/S3Client.h>
// #include <aws/s3/model/ListObjectsRequest.h>
// #include <aws/s3/model/GetObjectRequest.h>
// #endif

// OpenMP for parallelization
#ifdef OPENMP
#include <omp.h>
#endif

#include <zlib.h>
static int file_gzread(mtar_t *tar, void *data, size_t size) {
    size_t res = gzread((gzFile)tar->stream, data, size);
    return (res == size) ? MTAR_ESUCCESS : MTAR_EREADFAIL;
}

static int file_gzseek(mtar_t *tar, long offset, int whence) {
    int res = gzseek((gzFile)tar->stream, offset, whence);
    return (res != -1) ? MTAR_ESUCCESS : MTAR_ESEEKFAIL;
}

static int file_gzclose(mtar_t *tar) {
    gzclose((gzFile)tar->stream);
    return MTAR_ESUCCESS;
}

int mtar_gzopen(mtar_t *tar, const char *filename) {
    // Init tar struct and functions
    memset(tar, 0, sizeof(*tar));
    tar->read = file_gzread;
    tar->seek = file_gzseek;
    tar->close = file_gzclose;
    // Open file
    tar->stream = gzopen(filename, "rb");
    if (!tar->stream) {
        return MTAR_EOPENFAIL;
    }

#if defined(ZLIB_VERNUM) && ZLIB_VERNUM >= 0x1240
    gzbuffer((gzFile)tar->stream, 1 * 1024 * 1024);
#endif

    return MTAR_ESUCCESS;
}

using process_entry_func = std::function<bool(const char* name, const char* content, size_t size)>;

class Processor {
public:
    virtual ~Processor() {};
    virtual void run(process_entry_func, int) {};
};

class DirectoryProcessor : public Processor {
public:
    DirectoryProcessor(const std::string& input, bool recursive) {
        files = getFilesInDirectory(input, recursive);
    };
    DirectoryProcessor(std::vector<std::string> files) : files(std::move(files)) {};

    void run(process_entry_func func, int num_threads) override {
#pragma omp parallel shared(files) num_threads(num_threads)
        {
            char* dataBuffer;
            ssize_t bufferSize;
#pragma omp for
            for (size_t i = 0; i < files.size(); i++) {
                std::string name = files[i];
                FILE* file = fopen(name.c_str(), "r");
                dataBuffer = file_map(file, &bufferSize, 0);
                if (!func(name.c_str(), dataBuffer, bufferSize)) {
                    std::cerr << "[Error] processing dir entry " << name << " failed." << std::endl;
                    file_unmap(dataBuffer, bufferSize);
                    continue;
                }
                file_unmap(dataBuffer, bufferSize);
                fclose(file);
            }
        }
    }

private:
    std::vector<std::string> files;
};

class TarProcessor : public Processor {
public:
    TarProcessor(const std::string& input) {
        if (stringEndsWith(".gz", input) || stringEndsWith(".tgz", input)) {
            if (mtar_gzopen(&tar, input.c_str()) != MTAR_ESUCCESS) {
                std::cerr << "[Error] open tar " << input << " failed." << std::endl;
            }
        } else {
            if (mtar_open(&tar, input.c_str(), "r") != MTAR_ESUCCESS) {
                std::cerr << "[Error] open tar " << input << " failed." << std::endl;
            }
        }
    };

    ~TarProcessor() {
        mtar_close(&tar);
    };

    void run(process_entry_func func, int num_threads) override {
#pragma omp parallel shared(tar) num_threads(num_threads)
        {
            bool proceed = true;
            mtar_header_t header;
            size_t bufferSize = 1024 * 1024;
            char* dataBuffer = (char*)malloc(bufferSize);
            std::string name;
            while (proceed) {
                bool writeEntry = true;
#pragma omp critical
                {
                    if (tar.isFinished == 0 && (mtar_read_header(&tar, &header)) != MTAR_ENULLRECORD) {
                        // GNU tar has special blocks for long filenames
                        if (header.type == MTAR_TGNU_LONGNAME || header.type == MTAR_TGNU_LONGLINK) {
                            if (header.size > bufferSize) {
                                bufferSize = header.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                std::cerr << "[Error] cannot read entry " << header.name << std::endl;
                                goto done;
                            }
                            name.assign(dataBuffer, header.size);
                            // skip to next record
                            if (mtar_read_header(&tar, &header) == MTAR_ENULLRECORD) {
                                std::cerr << "[Error] tar truncated after entry " << name << std::endl;
                                goto done;
                            }
                        } else {
                            name = header.name;
                        }
                        if (header.type == MTAR_TREG || header.type == MTAR_TCONT || header.type == MTAR_TOLDREG) {
                            if (header.size > bufferSize) {
                                bufferSize = header.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                std::cerr << "[Error] cannot read entry " << name << std::endl;
                                goto done;
                            }
                            proceed = true;
                            writeEntry = true;
                        } else {
                            if (header.size > 0 && mtar_skip_data(&tar) != MTAR_ESUCCESS) {
                                std::cerr << "[Error] cannot skip entry " << name << std::endl;
                                goto done;
                            }
                            proceed = true;
                            writeEntry = false;
                        }
                    } else {
done:
                        tar.isFinished = 1;
                        proceed = false;
                        writeEntry = false;
                    }
                } // end read in
                if (proceed && writeEntry) {
                    if (!func(name.c_str(), dataBuffer, header.size)) {
                        std::cerr << "[Error] failed processing tar entry " << name << std::endl;
                        continue;
                    }
                }
            }
            free(dataBuffer);
        }
    }

private:
    mtar_t tar;
};

class DatabaseProcessor : public Processor {
public:
    DatabaseProcessor(const std::string& input) {
        std::string index = input + ".index";
        int mode = DB_READER_USE_DATA | DB_READER_USE_LOOKUP_REVERSE;
        handle = make_reader(input.c_str(), index.c_str(), mode);
    };

    DatabaseProcessor(const std::string& input, std::string& user_id_file) {
        std::string index = input + ".index";
        int mode = DB_READER_USE_DATA | DB_READER_USE_LOOKUP;
        handle = make_reader(input.c_str(), index.c_str(), mode);
        _read_id_list(user_id_file);
    };

    DatabaseProcessor(const std::string& input, const std::vector<std::string>& ids) {
        std::string index = input + ".index";
        int mode = DB_READER_USE_DATA | DB_READER_USE_LOOKUP;
        handle = make_reader(input.c_str(), index.c_str(), mode);
        user_ids = ids;
    };

    ~DatabaseProcessor() {
        free_reader(handle);
    };

    void run(process_entry_func func, int num_threads) override {
        size_t db_size = reader_get_size(handle);
#pragma omp parallel shared(handle) num_threads(num_threads)
        if (user_ids.size() == 0) { // process all entries in db
            {
#pragma omp for
                for (size_t i = 0; i < db_size; i++) {
                    uint32_t key = reader_get_key(handle, i);
                    const char* name = reader_lookup_name_alloc(handle, key);
                    // If name == "", throw warning
                    if (name && !name[0]) {
                        std::cerr << "[Warning] empty name for key " << key << std::endl;
                        free((void*)name);
                        continue;
                    }
                    if (!func(name, reader_get_data(handle, i), reader_get_length(handle, i))) {
                        std::cerr << "[Error] processing db entry " << name << " failed." << std::endl;
                        free((void*)name);
                        continue;
                    }
                    free((void*)name);
                }
            }
        } else { // process only entries in user_ids
            {
#pragma omp for
                for (size_t i = 0; i < user_ids.size(); i++) {
                    uint32_t key = reader_lookup_entry(handle, user_ids[i].c_str());
                    int64_t id = reader_get_id(handle, key);
                    if (id == -1 || key == UINT32_MAX) {
                        // NOT found
                        std::cerr << "[Warning] " << user_ids[i] << " not found in database." << std::endl;
                        continue;
                    }
                    if (!func(user_ids[i].c_str(), reader_get_data(handle, id), reader_get_length(handle, id))) {
                        std::cerr << "[Error] processing db entry " << user_ids[i] << " failed." << std::endl;
                        continue;
                    }
                }
            }
        }
    }

private:
    void* handle;
    std::vector<std::string> user_ids;

    void _read_id_list(std::string& file) {
        // Check if file exists
        if (!std::ifstream(file)) {
            std::cerr << "[Error] user id '" << file << "' does not exist." << std::endl;
            return;
        }
        std::ifstream infile(file);
        std::string line;
        while (std::getline(infile, line)) {
            user_ids.push_back(line);
        }
        infile.close();
    }
};

#ifdef HAVE_GCS
class GcsProcessor : public Processor {
public:
    namespace gcs = ::google::cloud::storage;
    GcsProcessor(const std::string& input) {
        auto options = google::cloud::Options{}
            .set<gcs::ConnectionPoolSizeOption>(num_threads)
            .set<google::cloud::storage_experimental::HttpVersionOption>("2.0");
        client = gcs::Client(options);
        bucket_name = input;
    };

    void run(process_entry_func func, int num_threads) override {
       #pragma omp parallel num_threads(num_threads)
        {
#pragma omp single
            // Get object list from gcs bucket
            for (auto&& object_metadata : client.ListObjects(bucket_name, gcs::Projection::NoAcl(), gcs::MaxResults(100000))) {
                std::string obj_name = object_metadata->name();
                // Set zero padding for ID with 4 digits
#pragma omp task firstprivate(obj_name)
                {
                    // Filter for splitting input into 10 different processes
                    // bool skipFilter = filter != '\0' && obj_name.length() >= 9 && obj_name[8] == filter;
                    bool skipFilter = true;
                    bool allowedSuffix = stringEndsWith(".cif", obj_name) || stringEndsWith(".pdb", obj_name);
                    if (skipFilter && allowedSuffix) {
                        auto reader = client.ReadObject(bucket_name, obj_name);
                        if (!reader.status().ok()) {
                            std::cerr << "Could not read object " << obj_name << std::endl;
                        } else {
                            std::string contents{ std::istreambuf_iterator<char>{reader}, {} };
                            func(obj_name.c_str(), contents.c_str(), contents.length());
                        }
                    }
                }
            }
        }
    }

private:
    google::cloud::storage::Client::Client client;
    std::string bucket_name;
};
#endif

// #ifdef HAVE_AWS_S3
// class S3Processor : public Processor {
// public:
//     S3Processor(const std::string& input) {
//         Aws::Client::ClientConfiguration clientConfig;
//         clientConfig.region = Aws::Region::US_WEST_2;  // Update the region if needed
//         client = std::make_shared<Aws::S3::S3Client>(clientConfig);
//         bucket_name = input;
//     };

//     void run(process_entry_func func, int num_threads) override {
// #pragma omp parallel num_threads(num_threads)
//         {
// #pragma omp single
//             // Get object list from S3 bucket
//             Aws::S3::Model::ListObjectsRequest objects_request;
//             objects_request.WithBucket(bucket_name);

//             auto list_objects_outcome = client->ListObjects(objects_request);

//             if (list_objects_outcome.IsSuccess()) {
//                 auto object_list = list_objects_outcome.GetResult().GetContents();
//                 for (auto const& s3_object : object_list) {
//                     std::string obj_name = s3_object.GetKey();
// #pragma omp task firstprivate(obj_name)
//                     {
//                         bool skipFilter = true;
//                         bool allowedSuffix = stringEndsWith(".cif", obj_name) || stringEndsWith(".pdb", obj_name);
//                         if (skipFilter && allowedSuffix) {
//                             Aws::S3::Model::GetObjectRequest object_request;
//                             object_request.WithBucket(bucket_name).WithKey(obj_name);

//                             auto get_object_outcome = client->GetObject(object_request);

//                             if (get_object_outcome.IsSuccess()) {
//                                 auto& retrieved_file = get_object_outcome.GetResultWithOwnership().GetBody();
//                                 std::string contents{ std::istreambuf_iterator<char>{retrieved_file}, {} };
//                                 func(obj_name.c_str(), contents.c_str(), contents.length());
//                             }
//                             else {
//                                 std::cerr << "Could not read object " << obj_name << std::endl;
//                             }
//                         }
//                     }
//                 }
//             }
//             else {
//                 std::cout << "ListObjects error: "
//                     << list_objects_outcome.GetError().GetExceptionName() << " - "
//                     << list_objects_outcome.GetError().GetMessage() << std::endl;
//             }
//         }
//     }

// private:
//     std::shared_ptr<Aws::S3::S3Client> client;
//     std::string bucket_name;
// };
// #endif
