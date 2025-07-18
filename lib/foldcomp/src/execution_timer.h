/**
 * File: execution_timer.h
 * Project: foldcomp
 * Created: 2022-09-29 16:45:16
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code is written as part of project "foldcomp".
 * ---
 * Last Modified: 2022-12-04 16:23:52
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 */
#include <chrono>
#include <type_traits>
#include <sstream>
/**
 * @brief General function to measure the running time of a function
 * https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
 * @param func
 * @return double
 */
auto static measureRunningTime = [](int n, auto && func, auto&&... params) {
    // get time before function invocation
    const auto& start = std::chrono::high_resolution_clock::now();
    // function invocation using perfect forwarding
    // n is the number of times to run the function; default is 100000
    for (auto i = 0; i < n; ++i) {
        std::forward<decltype(func)>(func)(std::forward<decltype(params)>(params)...);
    }
    // get time after function invocation
    const auto& stop = std::chrono::high_resolution_clock::now();
    return (stop - start) / n;
};

template<class Resolution = std::chrono::milliseconds>
class ExecutionTimer {
public:
    using Clock = std::conditional_t<std::chrono::high_resolution_clock::is_steady,
        std::chrono::high_resolution_clock,
        std::chrono::steady_clock>;
private:
    const Clock::time_point mStart = Clock::now();

public:
    ExecutionTimer() = default;
    ~ExecutionTimer() {}

    inline void stop() {
        const auto end = Clock::now();
        std::ostringstream strStream;
        strStream << "Stop Elapsed: "
            << std::chrono::duration<double>(end - mStart).count()
            << std::endl;
        std::cout << strStream.str() << std::endl;
    }

    inline double getElapsed() const {
        const auto end = Clock::now();
        return std::chrono::duration<double>(end - mStart).count();
    }

}; // ExecutionTimer

class TimerGuard {
private:
    ExecutionTimer<>* mTimer;
    const std::string& prefix;
    bool enabled;
public:
    TimerGuard(const std::string& prefix, bool enabled) : prefix(prefix), enabled(enabled) {
        if (enabled) {
            mTimer = new ExecutionTimer<>();
        }
    }
    ~TimerGuard() {
        if (!enabled) {
            return;
        }
        std::string elapsed = prefix;
        elapsed.append(1, '\t');
        elapsed.append(std::to_string(mTimer->getElapsed()));
        elapsed.append(1, '\n');
        std::cout << elapsed;
        delete mTimer;
    }
};