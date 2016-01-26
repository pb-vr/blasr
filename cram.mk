FAST_CTESTS := \
ctest/affineAlign.t            ctest/bamOut.t           ctest/ccsH5.t           ctest/filtercriteria.t  ctest/m0-5.t             ctest/samNM.t \
ctest/aggressiveIntervalCut.t  ctest/fofn.t             ctest/multipart.t \
ctest/alignScore.t             ctest/bug25741.t         ctest/ecoli.t           ctest/hitpolicy.t       ctest/noSplitSubreads.t \
ctest/bamIn.t                  ctest/fastMaxInterval.t  ctest/open_fail.t       ctest/verbose.t

MILD_CTESTS := \
	ctest/useccsallBestN1.t ctest/concordant.t ctest/bug25766.t ctest/holeNumbers.t

SLOW_CTESTS := ctest/bug25328.t ctest/useccsallLargeGenome.t

#BLASR_PATH=/mnt/secondary/builds/full/3.0.0/prod/current-build_smrtanalysis/private/otherbins/internalall/bin/
#export BLASR_PATH


cramfast:
	cram -v --shell=/bin/bash ${FAST_CTESTS}

crammild:
	cram -v --shell=/bin/bash ${MILD_CTESTS}

cramslow:
	cram -v --shell=/bin/bash ${SLOW_CTESTS}

cramtests:
	cram -v --shell=/bin/bash ${FAST_CTESTS} ${MILD_CTESTS} ${SLOW_CTESTS}

cramqu:
	for test in ${FAST_CTESTS}; do \
		qsub -pe smp 15 -V -cwd -b y -N cramqu $@cram -v --shell=bin/bash $$test;\
	done

clean:
	rm -f cramqu.* ctest/*.err
