//
// Created by Martin Steinegger on 29.09.18.
//

#ifndef STRUCCLUST_IOCIF_H
#define STRUCCLUST_IOCIF_H


#include <iostream>

//	Our CIF implementation consists of flyweight classes.

namespace mmCIF
{

    struct field
    {
        std::string name() const
        {
            return std::string(m_name, m_name_end);
        }

        std::string value() const
        {
            std::string result;

            if (m_data[0] == '\'' and m_data_end > m_data and m_data_end[-1] == '\'')
                result = std::string(m_data + 1, m_data_end - 1);
            else if (m_data[0] == '"' and m_data_end > m_data and m_data_end[-1] == '"')
                result = std::string(m_data + 1, m_data_end - 1);
            else if (m_data[0] == ';' and m_data_end > m_data and m_data_end[-1] == ';')
                result = std::string(m_data + 1, m_data_end - 1);
            else
                result = std::string(m_data, m_data_end);

            return result;
        }

        const char*		m_name;
        const char*		m_name_end;
        const char*		m_data;
        const char*		m_data_end;
    };

    struct row
    {
        std::string operator[](const char* inName) const;

        bool operator==(const row& rhs) const
        {
            return m_data == rhs.m_data and m_field == rhs.m_field;
        }

        const char*			m_data;
        int32				m_field;
        std::vector<field>	m_fields;
    };

    struct record
    {
        std::string name() const
        {
            return m_name;
        }

        struct const_iterator : public std::iterator<std::forward_iterator_tag, const row>
        {
            typedef std::iterator<std::forward_iterator_tag, const row>	base_type;
            typedef base_type::reference								reference;
            typedef base_type::pointer									pointer;

            const_iterator(const record& rec, const row& row) : m_rec(rec), m_row(row) {}
            const_iterator(const const_iterator& iter) : m_rec(iter.m_rec), m_row(iter.m_row) {}
            const_iterator&	operator=(const const_iterator& iter)			{ m_row = iter.m_row; return *this; }

            reference		operator*() const								{ return m_row; }
            pointer			operator->() const								{ return &m_row; }

            const_iterator&	operator++()									{ m_rec.advance(m_row); return *this; }
            const_iterator	operator++(int)									{ const_iterator iter(*this); operator++(); return iter; }

            bool			operator==(const const_iterator& iter) const	{ return m_row == iter.m_row; }
            bool			operator!=(const const_iterator& iter) const	{ return not operator==(iter); }

        private:
            const record&	m_rec;
            row				m_row;
        };

        typedef const_iterator iterator;

        row front() const;
        row back() const;

        const_iterator begin() const;
        const_iterator end() const;

        void advance(row& row) const;	// update pointers to next data row, if any

        bool operator<(const record& rhs) const
        {
            return m_name < rhs.m_name;
        }

        std::string		get_joined(const char* inFieldName, const char* inDelimiter) const;

        const char*		m_start;
        const char*		m_end;
        bool			m_loop;
        uint32			m_field_count;
        std::string		m_name;
    };

    class file
    {
    public:
        file(std::istream& is);

        record operator[](const char* inName) const;

        std::string get(const char* inName) const;
        std::string get_joined(const char* inName, const char* inDelimiter) const;

    private:
        std::vector<char>	m_buffer;
        std::vector<record>	m_records;
        const char*			m_data;
        const char*			m_end;
    };

}


#endif //STRUCCLUST_IOCIF_H
