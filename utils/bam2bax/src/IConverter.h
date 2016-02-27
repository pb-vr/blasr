// Author: Yuan Li

#ifndef BAM2BAX_ICONVERTER_H_
#define BAM2BAX_ICONVERTER_H_

#include "Settings.h"

class IConverter {
public:
    virtual ~IConverter(void);

public:
    virtual std::vector<std::string> Errors(void) const final;
    virtual bool Run(void) = 0;

protected:
    IConverter(Settings & settings);

    void AddErrorMessage(const std::string & errmsg) {
        errors_.push_back(errmsg);
    }

protected:
    // protected variables
    Settings & settings_;
    std::vector<std::string> errors_;
};
#endif
