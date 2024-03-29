#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "data_structures.hpp"
#include <fstream>
#include <Eigen/Dense>

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

template<typename T>
void writeToCSVfile(std::string name, Eigen::MatrixBase<T>& matrix, std::string const header = "")
{
    std::ofstream file(name.c_str());
    if (header != "") {
        file << "#" << header << "\n";
    }
    file << matrix.format(CSVFormat);
    file.close();
}

#endif