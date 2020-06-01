#ifndef STRUCCLUST_UTILS_H
#define STRUCCLUST_UTILS_H


// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#ifndef NDEBUG
#include <iostream>
#endif

#include <time.h>
#include <string>
#include <limits>
#include <algorithm>


template<typename T>
static inline T fast_atoi( const char * str )
{
    T val = 0;
    int sign=1;
    if(std::numeric_limits<T>::is_signed == true){
        if(*str == '-') {
            sign = -1;
            str++;
        }
    }
    while (*str >= '0' && *str <= '9') {
        val = val*10 + (*str++ - '0');
    }
    return sign * val;
}

namespace ba{



    static void to_upper(std::string &data){
        std::transform(data.begin(), data.end(), data.begin(), ::toupper);
    }

    static bool not_whitespace(int ch)
    {
        if (ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t') return false;
        return true;
    }

// trim from start
    static std::string & ltrim(std::string & s)
    {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_whitespace));
        return s;
    }

// trim from end
    static std::string & rtrim(std::string & s)
    {
        s.erase(std::find_if(s.rbegin(), s.rend(), not_whitespace).base(), s.end());
        return s;
    }

// trim from both ends
    static void trim(std::string & s)
    {
        ltrim(rtrim(s));
    }

    static std::string trim_copy(std::string s)
    {
        return ltrim(rtrim(s));
    }


    static bool starts_with(std::string &s, char * pattern){
        if (s.rfind(pattern, 0) == 0) {
            return true;
        }
        return false;
    }

}


// --------------------------------------------------------------------

// --------------------------------------------------------------------

#ifndef NDEBUG
struct stats
{
    stats() : m_max(0), m_count(0), m_cumm(0) {}
    ~stats();

    void operator()(uint32_t i)
    {
        if (m_max < i)
            m_max = i;
        ++m_count;
        m_cumm += i;
    }

    uint32_t m_max, m_count, m_cumm;
};
#endif

// --------------------------------------------------------------------

void WriteToFD(int inFD, const std::string& inText);



#endif //STRUCCLUST_UTILS_H
