include ../../common.mk
SUB_DIRS = alignment alignmentset anchoring bwt \
           mapping matrix metagenome reads repmask \
		   saf sequence suffixarray tuplelists

all:
	@(for subdir in $(SUB_DIRS); do \
		if [ -e $$subdir/Makefile ];\
		then $(MAKE) -C $$subdir -f Makefile; \
		fi; \
	done)

gtest:
	@(for subdir in $(SUB_DIRS); do \
		if [ -e $$subdir/Makefile ];\
		then $(MAKE) -C $$subdir -f Makefile gtest; \
		fi; \
	done)

clean:
	@(for subdir in $(SUB_DIRS); do \
		if [ -e $$subdir/Makefile ];\
		then $(MAKE) -C $$subdir -f Makefile clean; \
		fi; \
	done)
