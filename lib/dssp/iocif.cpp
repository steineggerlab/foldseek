//
// Created by Martin Steinegger on 29.09.18.
//

#include "iocif.h"

#include "mas.h"

#include <cassert>

#include <iostream>
#include <string>
#include <vector>

#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/algorithm/string.hpp>

#include "iocif.h"
#include "utils.h"

using namespace std;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

//	Our CIF implementation consists of flyweight classes.

namespace mmCIF
{

// skip routines to quickly position character pointer p at a next interesting location
    const char* skip_line(const char* p, const char* end);
    const char* skip_white(const char* p, const char* end);	// skip over white-space and comments
    const char* skip_value(const char* p, const char* end);

    string row::operator[](const char* inName) const
    {
        string result;

                foreach (const field& f, m_fields)
                    {
                        if (strncmp(inName, f.m_name, f.m_name_end - f.m_name) == 0)
                        {
                            result = f.value();
                            break;
                        }
                    }

        return result;
    }

    row record::front() const
    {
        row result;
        result.m_data = m_start;
        result.m_field = 0;

        const char* p = m_start;

        if (m_loop)
        {
            for (uint32 i = 0; i < m_field_count; ++i)
            {
                assert(*p == '_');
                assert(*(p + m_name.length()) == '.');

                field field = {};

                field.m_name = p = p + m_name.length() + 1;
                while (p != m_end and not isspace(*p))
                    ++p;

                field.m_name_end = p;

                p = skip_white(p, m_end);

                result.m_fields.push_back(field);
            }

                    foreach (field& fld, result.m_fields)
                        {
                            fld.m_data = skip_white(p, m_end);
                            fld.m_data_end = skip_value(fld.m_data, m_end);
                            p = skip_white(fld.m_data_end, m_end);
                        }
        }
        else
        {
            for (uint32 i = 0; i < m_field_count; ++i)
            {
                assert(*p == '_');
                assert(*(p + m_name.length()) == '.');

                field field;
                field.m_name = p = p + m_name.length() + 1;
                while (p != m_end and not isspace(*p))
                    ++p;

                field.m_name_end = p;

                p = skip_white(p, m_end);
                field.m_data = p;

                p = skip_value(p, m_end);
                field.m_data_end = p;
                p = skip_white(p, m_end);

                result.m_fields.push_back(field);
            }
        }

        return result;
    }

    record::iterator record::begin() const
    {
        return const_iterator(*this, front());
    }

    record::iterator record::end() const
    {
        row end = { m_start, -1 };
        return const_iterator(*this, end);
    }

    void record::advance(row& row) const
    {
        if (m_loop and not row.m_fields.empty())
        {
            const char* p = skip_white(row.m_fields.back().m_data_end, m_end);

            if (p >= m_end)
            {
                row.m_fields.clear();
                row.m_field = -1;
            }
            else
            {
                        foreach (field& fld, row.m_fields)
                            {
                                fld.m_data = skip_white(p, m_end);
                                fld.m_data_end = skip_value(fld.m_data, m_end);
                                p = skip_white(fld.m_data_end, m_end);
                            }

                row.m_field += 1;
            }
        }
        else
        {
            row.m_fields.clear();
            row.m_field = -1;
        }
    }

    string record::get_joined(const char* inName, const char* inDelimiter) const
    {
        string result;

        for (iterator i = begin(); i != end(); ++i)
        {
            string s = i->operator[](inName);
            ba::trim(s);
            result = (result.empty() ? result : result + inDelimiter) + s;
        }

        return result;
    }

    file::file(istream& is)
    {
        // first extract data into a buffer
        m_buffer.reserve(10 * 1024 * 1024);	// reserve 10 MB, should be sufficient for most

        io::copy(is, io::back_inserter(m_buffer));

        m_data = &m_buffer[0];
        m_end = m_data + m_buffer.size();

        m_buffer.push_back(0);				// end with a null character, makes coding a bit easier

        // CIF files are simple to parse

        const char* p = m_data;

        if (strncmp(p, "data_", 5) != 0)
            throw mas_exception("Is this an mmCIF file?");

        p = skip_line(p, m_end);

        record rec = { p };
        bool loop = false;

        while (p < m_end)
        {
            if (isspace(*p))	// skip over white space
            {
                ++p;
                continue;
            }

            if (*p == '#')	 // line starting with hash, this is a comment, skip
            {
                p = skip_line(p, m_end);
                continue;
            }

            if (strncmp(p, "loop_", 5) == 0)
            {
                if (not m_records.empty() and m_records.back().m_end == nullptr)
                    m_records.back().m_end = p;

                loop = true;
                rec.m_loop = false;
                p = skip_line(p + 5, m_end);

                continue;
            }

            const char* s = p;

            if (*p == '_')	// a label
            {
                // scan for first dot
                bool newName = loop;
                const char* n = rec.m_start;

                for (;;)
                {
                    if (not newName and *p != *n)
                        newName = true;

                    ++p;
                    ++n;

                    if (p == m_end or *p == '.' or isspace(*p))
                        break;
                }

                if (*p == '.')	// OK, found a record
                {
                    if (newName)
                    {
                        // store start as end for the previous record, if any
                        if (not m_records.empty() and m_records.back().m_end == nullptr)
                            m_records.back().m_end = s;

                        rec.m_start = s;
                        rec.m_end = nullptr;
                        rec.m_loop = loop;
                        rec.m_field_count = 1;
                        rec.m_name = string(s, p);

                        m_records.push_back(rec);
                    }
                    else
                        m_records.back().m_field_count += 1;

                    // skip over field name
                    while (p != m_end and not isspace(*p))
                        ++p;
                }
                else
                {
                    // store start as end for the previous record, if any
                    if (not m_records.empty() and m_records.back().m_end == nullptr)
                        m_records.back().m_end = s;

                    // a record without a field (is that possible in mmCIF?)
                    cerr << "record without field: " << string(s, p) << endl;

                    rec.m_start = s;
                    rec.m_end = nullptr;
                    rec.m_loop = loop;
                    rec.m_field_count = 0;
                    rec.m_name = string(s, p);

                    m_records.push_back(rec);
                }

                if (not rec.m_loop)
                    p = skip_value(p, m_end);

                loop = false;
                continue;
            }

            if (rec.m_loop == false)
            {
                // guess we should never reach this point
                throw mas_exception("invalid CIF file? (unexpected data, not in loop)");
            }

            p = skip_value(p, m_end);
            p = skip_white(p, m_end);

            // check for a new data_ block
            if (p != m_end and strncmp(p, "data_", 5) == 0)
                throw mas_exception("Multiple data blocks in CIF file");
        }

        if (not m_records.empty() and m_records.back().m_end == nullptr)
            m_records.back().m_end = p;

        sort(m_records.begin(), m_records.end());
    }

    record file::operator[](const char* inName) const
    {
        record result = {};
        result.m_name = inName;

        vector<record>::const_iterator i = lower_bound(m_records.begin(), m_records.end(), result);
        if (i != m_records.end() and i->m_name == inName)
            result = *i;

        return result;
    }

    string file::get(const char* inName) const
    {
        const char* p = strchr(inName, '.');
        assert(p != nullptr);
        if (p == nullptr)
            throw logic_error("incorrect name");

        record r = operator[](string(inName, p).c_str());
        return r.front()[string(p + 1).c_str()];
    }

    string file::get_joined(const char* inName, const char* inDelimiter) const
    {
        const char* p = strchr(inName, '.');
        assert(p != nullptr);
        if (p == nullptr)
            throw logic_error("incorrect name");

        record test;
        test.m_name.assign(inName, p);

        string result;

        vector<record>::const_iterator i = lower_bound(m_records.begin(), m_records.end(), test);
        if (i != m_records.end() and i->m_name == test.m_name)
            result = i->get_joined(p + 1, inDelimiter);

        return result;
    }

// skip to first character after the next NL character
    const char* skip_line(const char* p, const char* end)
    {
        while (p != end)
        {
            if (*p++ == '\n')
                break;
        }

        return p;
    }

// skip over white space and comments
    const char* skip_white(const char* p, const char* end)
    {
        while (p != end)
        {
            if (isspace(*p))
            {
                ++p;
                continue;
            }

            if (*p == '#')
            {
                do ++p; while (p < end and *p != '\n');
                continue;
            }

            break;
        }

        return p;
    }

// skip over values for a record
    const char* skip_value(const char* p, const char* end)
    {
        for (;;)
        {
            if (isspace(*p))
            {
                ++p;
                continue;
            }

            if (*p == ';' and *(p - 1) == '\n')
            {
                do p = skip_line(p, end); while (p < end and *p != ';');
                ++p;
                break;
            }

            if (*p == '\'')
            {
                do ++p; while (p != end and *p != '\'');
                ++p;
                break;
            }

            if (*p == '\"')
            {
                do ++p; while (p != end and *p != '\"');
                ++p;
                break;
            }

            if (*p == '#')
            {
                p = skip_line(p, end);
                continue;
            }

            while (p != end and not isspace(*p))
                ++p;

            break;
        }

        return p;
    }

}

//void ReadCIF(std::istream& in, MProtein& out)
//{
//	file cif(&buffer[0], buffer.size() - 1);
//	
//	
//
//	cout << "id: " << cif["_entry"].front()["id"].value() << endl;
//	
//	foreach (const row& row, cif["_atom_type"])
//	{
//		cout << row["symbol"].value() << endl;
//	}
//	
//	foreach (const row& row, cif["_atom_site"])
//	{
//		cout << "ATOM  " << row["Cartn_x"].value() << ' ' << row["Cartn_y"].value() << ' ' << row["Cartn_z"].value() << endl;
//	}
//}
//
