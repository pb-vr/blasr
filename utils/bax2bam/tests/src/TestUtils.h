// Author: Derek Barnett

#include "SMRTSequence.hpp"
#include <pbbam/BamRecord.h>
#include <string>
#include <vector>

void RemoveFile(const std::string& filename);
void RemoveFiles(const std::vector<std::string>& filenames);

int RunBax2Bam(const std::vector<std::string>& baxFilenames,
               const std::string& outputType,
               const std::string& additionalArgs = std::string());
