/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cFiles/file.h to edit this template
 */

/* 
 * File:   mpeg7_io.h
 * Author: rimlab
 *
 * Created on February 22, 2022, 11:56 AM
 */

#ifndef MPEG7_IO_H
#define MPEG7_IO_H

#include <glob.h>
#include <experimental/filesystem>
#include "boost/date_time/posix_time/posix_time.hpp"



namespace mpeg7io {

    // ----------------------------------------------
    // I/O OPERATIONS
    // ----------------------------------------------

    /** Returns a list of files based on Unix-like GLOB. 
     */
    void glob(const std::string globPath, std::vector<std::string>& matchingFiles);

    /** Returns the list of files in the given directory. 
     */
    void getDirectoryFiles(const std::string& dirPath, std::vector<std::string>& matchingFiles);

    /** Generates a filename dependent on date and time.
     */
    std::string generateStampedString(const std::string prefix = "", const std::string postfix = "");

}



#endif /* MPEG7_IO_H */

