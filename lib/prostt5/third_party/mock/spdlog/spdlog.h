#ifndef MOCK_SPDLOG_H
#define MOCK_SPDLOG_H

#include <mutex>

namespace spdlog {

template <typename... Args>
void debug(Args&&...) {}

template <typename... Args>
void info(Args&&...) {}

template <typename... Args>
void warn(Args&&...) {}

template <typename... Args>
void error(Args&&...) {}

template <typename... Args>
void critical(Args&&...) {}

}

#endif
