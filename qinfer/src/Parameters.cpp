#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <fstream>
#include <regex>
#include <vector>

#include "Parameters.h"

void
Parameters::parseParameters(int argc, char* argv[]) {

  // We need three input files (+ program name as first arg)
  if (argc < 4) {
    std::cerr << "At least 3 arguments required. "
              << argc - 1
              << " given. Abort!"
              << std::endl;
    throw std::runtime_error("Not enough parameters!");
  }

  for(size_t i = 3; i < (size_t)argc; ++i) {

    std::string arg(argv[i]);
    std::regex rgx1("-([A-Za-z0-9-]*)=([A-Za-z0-9-_\\.]*)");
    std::regex rgx2("-([A-Za-z0-9-]*)");
    std::smatch matches;
    std::string parameterString = "";
    bool setValue = false;

    if(std::regex_search(arg, matches, rgx1)) {
      parameterString = matches[1].str();
      setValue = true;
    } else if(std::regex_search(arg, matches, rgx2)) {
      parameterString = matches[1].str();
    } else if(i > 3) {
      throw std::runtime_error("Illegal parameter format");
    }

    if(parameterString != ""){
      if(parameterString == "e" || parameterString == "-epsilon"){
        m_restrictedY = true;
        if(setValue){
          m_epsilon = std::stod(matches[2].str());
        }

      } else if(parameterString == "d" || parameterString == "-disable-quartet"){
        m_disableQuartet =true;

      } else if(parameterString == "w" || parameterString == "-weighted"){
        m_weightedMode = true;

      } else if(parameterString == "s" || parameterString == "-subtree-files"){
        m_subtreeFiles = true;

      } else if(parameterString == "a" || parameterString == "-all-outgroups"){
        m_relativeOutgroups = true;
        if(setValue){
          m_incongruenceThreshold = std::stod(matches[2].str());
        }

      } else if(parameterString == "l" || parameterString == "-outgroup-limit"){
        if(setValue){
          m_outgroupLimit = std::stoi(matches[2].str());
        }

      } else if(parameterString == "b" || parameterString == "-benchmark"){
        m_benchmark = true;
        if(setValue){
          m_benchmarkFile = matches[2].str();
        }
      } else {
        throw std::runtime_error("Illegal parameter: " + parameterString);
      }
    }
  }

  m_matrixFile = std::string(argv[1]);
  m_speciesFile = std::string(argv[2]);
  if(!m_disableQuartet){
    m_treeFile = std::string(argv[3]);
  }
}

bool
Parameters::checkIntegrity() {
  // some logic checks
  if(m_subtreeFiles && m_relativeOutgroups){
    std::cerr << "All outgroups option is not available in combination "
              << "with subtree files! Use a Newick-tree file instead."
              << std::endl;
    return false;
  }

  if(m_incongruenceThreshold < 0.0){
    m_incongruenceThreshold = 0.0;
  } else if(m_incongruenceThreshold > 1.0){
    m_incongruenceThreshold = 1.0;
  }

  if(m_outgroupLimit < 1){
    std::cerr << "Outgroup limit too low: "
              << m_outgroupLimit
              << std::endl;
    return false;
  }

  // check if files exist
  std::vector<std::string> fileNames = {m_matrixFile, m_speciesFile, m_treeFile};
  const size_t fileCount = (m_disableQuartet) ? 2 : 3;
  auto fileError = false;
  for(size_t i = 0; i < fileCount; ++i) {
    const std::filesystem::path filepath(fileNames[i]);

    if(!std::filesystem::exists(filepath) || !std::filesystem::is_regular_file(filepath)) {
      std::cerr << "File "
                << filepath
                << " does not exist or is no file!"
                << std::endl;
      fileError = true;
    }
  }

  return !fileError;
}
