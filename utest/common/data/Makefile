include ../../common.mk

SUB_DIRS = hdf amos

all:
	@(for subdir in $(SUB_DIRS); do \
		if [ -e $$subdir/Makefile ];\
		then $(MAKE) -C $$subdir -f Makefile || true ; \
		fi; \
	done)

gtest:
	@(for subdir in $(SUB_DIRS); do \
		if [ -e $$subdir/Makefile ];\
		then $(MAKE) -C $$subdir -f Makefile gtest || true; \
		fi; \
	done)

clean:
	@(for subdir in $(SUB_DIRS); do \
		if [ -e $$subdir/Makefile ];\
		then $(MAKE) -C $$subdir -f Makefile clean || true; \
		fi; \
	done)
