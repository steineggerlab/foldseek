#ifndef PROSTT5_FORK_RUNNER_H
#define PROSTT5_FORK_RUNNER_H
#include "DBReader.h"
#include "DBWriter.h"
#include "ProstT5.h"

#include <unistd.h>
#include <signal.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/wait.h>

#ifdef OPENMP
#include <omp.h>
#endif

struct TaskMsg {
    long mtype;
    long idx;
};

void prostt5Forking(
    const std::string& modelWeights,
    unsigned int split_length,
    unsigned int minSplitLength,
    DBReader<unsigned int>& reader,
    const std::string& db,
    const std::string& index,
    int threads,
    int compressed
    ) {
#ifdef OPENMP
    // forking does not play well with OpenMP threads
    omp_set_num_threads(1);
#endif

    int procs = (threads + 3) / 4;
    int leftover = threads;
    int msgid = msgget(IPC_PRIVATE, IPC_CREAT | 0666);
    if (msgid == -1) {
        Debug(Debug::ERROR) << "Could not create SysV message queue!\n";
        EXIT(EXIT_FAILURE);
    }

    Debug::Progress progress(reader.getSize());
    std::vector<pid_t> children;
    int maxSplits = procs;
    for (int p = 0; p < procs; ++p) {
        int inner = std::min(4, leftover);
        leftover -= inner;
        switch (pid_t pid = fork()) {
            default:
                children.push_back(pid);
                break;
            case -1:
                Debug(Debug::ERROR) << "Could not fork worker process!\n";
                EXIT(EXIT_FAILURE);
            case 0: {
                std::string device = "none";
                ProstT5Model model(modelWeights, device);
                ProstT5 context(model, inner);
                std::string result;
                const char newline = '\n';

                std::pair<std::string, std::string> outDb = Util::createTmpFileNames(db, index, p);
                DBWriter writer(outDb.first.c_str(), outDb.second.c_str(), 1, compressed, reader.getDbtype());
                writer.open();
                while (true) {
                    TaskMsg msg;
                    if (msgrcv(msgid, &msg, sizeof(msg.idx), 0, 0) == -1) {
                        Debug(Debug::ERROR) << "msgrcv failed in child " << p << "\n";
                        _Exit(1);
                    }
                    if (msg.idx == -1) {
                        break;
                    }

                    size_t i = static_cast<size_t>(msg.idx);
                    unsigned int key = reader.getDbKey(i);
                    size_t length = reader.getSeqLen(i);
                    std::string seq = std::string(reader.getData(i, 0), length);
                    result.clear();
                    // splitting input sequences longer than ProstT5 attention (current cutoff 6000 AAs)
                    // split length of 0 will deactivate splitting
                    // Debug(Debug::INFO) << "split_length: " << split_length << " minSplitLength: " <<  minSplitLength << "\n";
                    // Debug(Debug::INFO) << "seq: " << seq << "\n";
                    if (split_length > 0 && length > split_length) {
                        unsigned int n_splits, overlap_length;
                        n_splits = int(length / split_length) + 1;
                        overlap_length = length % split_length;

                        // ensure minimum overlap length; adjustment length was not computed properly with ceil/ceilf now using simple int cast
                        if (overlap_length < minSplitLength) {
                            split_length -= int((minSplitLength - overlap_length) / (n_splits - 1)) + 1;
                        }

                        // loop over splits and predict
                        for (unsigned int i = 0; i < n_splits; i++){
                            unsigned int split_start = i * split_length;
                            result.append(context.predict(seq.substr(split_start, split_length)));
                        }
                    } else {
                        result.append(context.predict(seq));
                    }
                    // Debug(Debug::INFO) << "p: " << p << "pred: " << result << "\n";

                    writer.writeStart(0);
                    writer.writeAdd(result.c_str(), result.length(), 0);
                    writer.writeAdd(&newline, 1, 0);
                    writer.writeEnd(key, 0);
                    progress.updateProgress(i);
                }

                std::cout.setstate(std::ios_base::failbit);
                writer.close(true);
                fflush(NULL);
                sync();
                _Exit(0);
            }
        }
        if (leftover <= 0) {
            maxSplits = p + 1;
            break;
        }
    }

    for (size_t i = 0; i < reader.getSize(); ++i) {
        TaskMsg msg {1, static_cast<long>(i)};
        if (msgsnd(msgid, &msg, sizeof(msg.idx), 0) == -1) {
            Debug(Debug::ERROR) << "msgsnd failed for index " << i << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    for (int p = 0; p < procs; ++p) {
        TaskMsg quitMsg {1, -1};
        msgsnd(msgid, &quitMsg, sizeof(quitMsg.idx), 0);
    }

    for (const pid_t& child_pid : children) {
        int status = 0;
        while (waitpid(child_pid, &status, 0) == -1) {
            if (errno == EINTR) {
                continue;
            }
            perror("waitpid");
            break;
        }
    }
    msgctl(msgid, IPC_RMID, nullptr);
    fflush(NULL);
    sync();
    std::pair<std::string, std::string> outDb = std::make_pair(db, index);
    std::vector<std::pair<std::string, std::string>> splitFiles;
    for (int p = 0; p < maxSplits; ++p) {
        splitFiles.emplace_back(Util::createTmpFileNames(outDb.first, outDb.second, p));
    }
    DBWriter::mergeResults(outDb.first, outDb.second, splitFiles);
}
#endif