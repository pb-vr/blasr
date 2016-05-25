#include "MetadataWriter.h"

std::string internal::Replace(const std::string & in_str,
                              const std::string & to_find,
                              const std::string & to_replace) {
    // Replace the first occurrence of to_find by to_replace.
    std::string ret = in_str;
    std::size_t pos = ret.find(to_find);
    if (pos != std::string::npos) {
        ret.replace(pos, to_find.size(), to_replace);
    }
    return ret;
}

MetadataWriter::MetadataWriter(const std::string & filename, 
                               const PacBio::BAM::ReadGroupInfo & rg,
                               const std::string & analysisDir) {
    MetadataWriter(filename, 
                   rg.BasecallerVersion(),
                   rg.SequencingKit(),
                   rg.BindingKit(),
                   analysisDir);
}

MetadataWriter::MetadataWriter(const std::string & filename, 
                               const std::string & basecallerVersion,
                               const std::string & sequencingKit,
                               const std::string & bindingKit,
                               const std::string & analysisDir) {
    assert(analysisDir.find('/') == std::string::npos);
    std::ofstream ofile; 
    ofile.open(filename, std::ofstream::out);

    std::string to_print = internal::META_CONTENT;
    to_print = internal::Replace(to_print, "__BASECALLERVERSION__", basecallerVersion);
    to_print = internal::Replace(to_print, "__SEQUENCINGKIT__", sequencingKit);
    to_print = internal::Replace(to_print, "__BINDINGKIT__", bindingKit);
    to_print = internal::Replace(to_print, "__ANALYSISDIR__", analysisDir);

    ofile << to_print << std::endl;
    ofile.close();
}
