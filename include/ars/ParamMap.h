#ifndef PARAM_MAP_H
#define PARAM_MAP_H

#include <iostream>
#include <fstream>
#include <sstream>
//#include <unordered_map>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "definitions.h"

namespace ars {

    /** Reads and stores parameters from a string, a file, etc.
     */
    class ParamMap {
    public:
        //typedef std::unordered_map<std::string, std::string> table_type;
        typedef std::map<std::string, std::string> table_type; // a map stores lexically ordered parameters (nicer to view!)
        typedef table_type::iterator iterator;
        typedef table_type::const_iterator const_iterator;

        const static char COMMENT = '#';
        const static unsigned int MAX_LEN = 2000;

        /** Default constructor.
         */
        ParamMap() : table_() {
        }

        /** Destructor.
         */
        ~ParamMap() {
        }

        /** Clears all the content of param table.
         */
        void clear() {
            table_.clear();
        }

        /** Reads params from an input stream in the format:
         *   key1 value1
         *   key2 value2
         *   ...
         */
        bool read(std::istream& in);

        /** Reads params from an input file (format as above).
         */
        bool read(std::string& filename);

        /** Reads from a command line. Required format
         *   ...
         *   argv[i] = "-key1"      // starting with "-"
         *   argv[i+1] = "value1"
         */
        bool read(int argc, char** argv);

        /** Writes the pairs (key,value).
         */
        bool write(std::ostream& out) const;

        /** Writes the parameters to an output file (format as above).
         */
        bool write(std::string& filename) const;

        /** Sets the param (as a set of string).
         */
        void setParamString(std::string paramName, std::string paramValue);

        /** Sets the param (as a set of string).
         */
        template <typename Value>
        void setParam(std::string paramName, const Value& paramValue) {
            std::stringstream sstr;
            sstr << paramValue;
            table_.erase(paramName);
            table_.insert(std::make_pair(paramName, sstr.str()));
        }

        /** Casts the string value of a given parameters to the desired value.
         */
        template <typename Value>
        bool getParam(std::string paramName, Value& value, const Value& defaultValue) {
            const_iterator v = table_.find(paramName);
            if (v != table_.end()) {
                try {
                    value = boost::lexical_cast<Value>(v->second);
                } catch (boost::bad_lexical_cast const&) {
                    std::cerr << __FILE__ << "," << __LINE__ << ": Error: cannot cast string \""
                            << v->second << "\" to type \"" << typeid (Value).name() << "\" for variable \"" << v->first << "\"" << std::endl;
                }
            } else {
                //            std::cerr << "Parameter " << paramName << " not found." << std::endl;
                value = defaultValue;
                setParam(paramName, defaultValue);
                return false;
            }
            return true;
        }

        template <typename Value, typename Iterator>
        bool getParamContainer(std::string paramName, Iterator beg, Iterator end, std::string defaultString, const Value &defaultValue, std::string delim = "[],") {
            Iterator cit;
            // Initializes the vector values with default value
            for (cit = beg; cit != end; ++cit) {
                *cit = defaultValue;
            }
            const_iterator v = table_.find(paramName);
            if (v != table_.end()) {
                fillWithTokens(v->second, beg, end, defaultValue, delim);
                return true;
            } else {
                fillWithTokens(defaultString, beg, end, defaultValue, delim);
                setParam(paramName, defaultString);
                return false;
            }
            return true;
        }

    protected:
        table_type table_;

        static bool isOption(std::string str);

        template <typename Value, typename Iterator>
        static void fillWithTokens(const std::string& valueString, Iterator beg, Iterator end, const Value& valueDefault, const std::string& delim) {
            try {
                // Splits the value into tokens, e.g. "[1,2,3]" with delim "[,]" should become tokens "1", "2" and "3"
                boost::char_separator<char> sep(delim.c_str());
                boost::tokenizer<boost::char_separator<char> > tokens(valueString, sep);
                // Casts each token into a value
                auto sit = tokens.begin();
                for (Iterator vit = beg; vit != end; ++vit) {
                    if (sit != tokens.end()) {
                        *vit = boost::lexical_cast<Value>(*sit);
                        ++sit;
                    } else {
                        *vit = valueDefault;
                    }
                }
            } catch (boost::bad_lexical_cast const&) {
                ARS_ERROR("Error: cannot cast \"" << valueString << "\" to container values");
                return;
            }
        }
    };

} // end of namespace 

#endif

