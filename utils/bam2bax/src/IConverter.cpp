#include "IConverter.h"

IConverter::IConverter(Settings & settings)
:settings_(settings) {}

IConverter::~IConverter(void) {}

std::vector<std::string> IConverter::Errors(void) const {
    return errors_;
}

