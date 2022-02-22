#include "ars/mpeg7_io.h"

namespace mpeg7io {
    // ----------------------------------------------
    // I/O OPERATIONS
    // ----------------------------------------------

    void glob(const std::string globPath, std::vector<std::string>& matchingFiles) {
        glob_t glob_result;
        matchingFiles.clear();

        // glob struct resides on the stack
        memset(&glob_result, 0, sizeof (glob_result));

        ::glob(globPath.c_str(), GLOB_TILDE, NULL, &glob_result);
        for (unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
            matchingFiles.push_back(std::string(glob_result.gl_pathv[i]));
        }
        globfree(&glob_result);
    }

    void getDirectoryFiles(const std::string& dirPath, std::vector<std::string>& matchingFiles) {
        std::cout << "Looking in directory \"" << dirPath << "\"" << std::endl;
        std::experimental::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
        for (std::experimental::filesystem::directory_iterator i(dirPath); i != end_itr; ++i) {
            // Skip if not a file
            if (std::experimental::filesystem::is_regular_file(i->status())) {
                matchingFiles.push_back(i->path().string());
            }
        }
        std::sort(matchingFiles.begin(), matchingFiles.end());
    }

    std::string generateStampedString(const std::string prefix, const std::string postfix) {
        boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
        std::ostringstream formatter;
        std::string formatstring = prefix + "%Y%m%d_%H%M_%S" + postfix;
        formatter.imbue(std::locale(std::cout.getloc(), new boost::posix_time::time_facet(formatstring.c_str())));
        formatter << now;
        return formatter.str();
    }

}