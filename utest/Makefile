.PHONY: all gtest clean

all:
	@(if [ -e common/Makefile ]; \
	then $(MAKE) -C common -f Makefile; \
	fi)

gtest:
	@(if [ -e common/Makefile ]; \
	then $(MAKE) -C common -f Makefile gtest; \
	fi)

clean:
	@(if [ -e common/Makefile ]; \
	then $(MAKE) -C common -f Makefile clean; \
	fi)

