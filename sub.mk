
SRCDIR:=$(dir $(realpath $(firstword ${MAKEFILE_LIST})))
-include ${CURDIR}/defines.mk
-include ${SRCDIR}/rules.mk

export

configure-submodule:
	cd libcpp && ./configure.py
