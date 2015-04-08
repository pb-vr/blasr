#!/bin/bash

# ---- required subroutines
set_globals() {
    # required globals:
    g_name="blasr"

    # other globals:
    g_make_rootdir=$(readlink -f "${g_progdir}/../../../../prebuilt.tmpout/make/make_3.81-8.2ubuntu3/usr")
#    g_gcc_rootdir=$(readlink -f "${g_progdir}/../../../../prebuilt.tmpout/gcc-toolchain/gcc-toolchain-4.8.2/usr")
    g_gcc_rootdir=$(readlink -f "${g_progdir}/../../../../prebuilt.out/3.x/gcc/gcc-4.8.4/libc-2.5")

    g_make_cmd='
            ${g_make_rootdir}/bin/make -C "${g_progdir}/.."
                CXX="${g_gcc_rootdir}/bin/g++"
                AR="${g_gcc_rootdir}/bin/ar"
    '
    g_make_cmd=$(echo $g_make_cmd)

    g_outdir_name="_output";
    g_outdir_abs="$(readlink -m "$g_progdir/../$g_outdir_name")"
    g_outdir=$g_outdir_abs;
    g_installbuild_dir_abs="${g_outdir}/install-build"
    g_installbuild_dir="${g_installbuild_dir_abs}"
}

# ---- build targets

clean() {
    echo "Running $g_name 'clean' target..."
    # Clean the build artifacts
    # NOTE: don't need HDF5_INC/LIB for clean, supply /dev/null to inhibit
    #       errors
    eval "$g_make_cmd" clean \
	HDF5_INC="/dev/null" \
	HDF5_LIB="/dev/null" \
	${1+"$@"}
    # Remove the _output directory
    rm -rf "${g_outdir}"
}
cleanall() {
    echo "Running $g_name 'cleanall' target..."
    clean;
}
build() {
    echo "Running $g_name 'build' target..."

    # Create dependency links
    rm -rf "${g_outdir}/deplinks"
    mkdir -p "${g_outdir}/deplinks"
    ln -s "../../../../../prebuilt.out/zlib/zlib-1.2.5/centos-5" "${g_outdir}/deplinks/zlib"
    ln -s "../../../../../prebuilt.out/hdf5/hdf5-1.8.12/centos-5" "${g_outdir}/deplinks/hdf5"
    ln -s "../../../../lib/cpp/pbdata/${g_outdir_name}/install-build" "${g_outdir}/deplinks/libpbdata"
    ln -s "../../../../lib/cpp/hdf/${g_outdir_name}/install-build" "${g_outdir}/deplinks/libpbihdf"
    ln -s "../../../../lib/cpp/alignment/${g_outdir_name}/install-build" "${g_outdir}/deplinks/libblasr"
    ln -s "../../../../staging/PostPrimary/pbbam/${g_outdir_name}/install-build" "${g_outdir}/deplinks/pbbam"
    ln -s "../../../../staging/PostPrimary/pbbam/third-party/htslib/${g_outdir_name}/install-build" "${g_outdir}/deplinks/htslib"
    ln -s "../../../../../prebuilt.out/boost/boost_1_55_0" "${g_outdir}/deplinks/boost"

    # build
    eval "$g_make_cmd" \
	NO_SUBMAKES=true \
	ZLIB_ROOT="${g_outdir_name}/deplinks/zlib" \
	HDF5_INC="${g_outdir_name}/deplinks/hdf5/include" \
        LIBPBDATA_INCLUDE="${g_outdir_name}/deplinks/libpbdata/include" \
	LIBPBIHDF_INCLUDE="${g_outdir_name}/deplinks/libpbihdf/include/hdf" \
	LIBBLASR_INCLUDE="${g_outdir_name}/deplinks/libblasr/include/alignment" \
	PBBAM_INCLUDE="${g_outdir_name}/deplinks/pbbam/include" \
	HTSLIB_INCLUDE="${g_outdir_name}/deplinks/htslib/include" \
	BOOST_INCLUDE="${g_outdir_name}/deplinks/boost" \
	HDF5_LIB="${g_outdir_name}/deplinks/hdf5/lib" \
        LIBPBDATA_LIB="${g_outdir_name}/deplinks/libpbdata/lib" \
	LIBPBIHDF_LIB="${g_outdir_name}/deplinks/libpbihdf/lib" \
	LIBBLASR_LIB="${g_outdir_name}/deplinks/libblasr/lib" \
	PBBAM_LIB="${g_outdir_name}/deplinks/pbbam/lib" \
	HTSLIB_LIB="${g_outdir_name}/deplinks/htslib/lib" \
	${1+"$@"}
}
install_build() {
    if ! $opt_no_sub_targets; then
	build;
    fi

    echo "Running $g_name 'install-build' target..."

    # clean install dir
    rm -rf "$g_installbuild_dir";
    mkdir -p "$g_installbuild_dir";

    # install bin executables
    mkdir "$g_installbuild_dir/bin"
    cp -a "${g_progdir}/../${g_name}"  "$g_installbuild_dir/bin"

}
install_prod() {
    echo "Running $g_name 'install-prod' target..."
}
publish_build() {
    if ! $opt_no_sub_targets; then
	install_build;
    fi

    echo "Running $g_name 'publish-build' target..."

}
publish_prod() {
    if ! $opt_no_sub_targets; then
	install_prod;
    fi
    echo "Running $g_name 'cleanall' target..."

}


# ---- End of Module-specific code
# Common code from here on out, do not modify...

# ---- error handling
set -o errexit;
set -o posix;
set -o pipefail;
set -o errtrace;
unexpected_error() {
    local errstat=$?
    echo "${g_prog:-$(basename $0)}: Error! Encountered unexpected error at 'line $(caller)', bailing out..." 1>&2
    exit $errstat;
}
trap unexpected_error ERR;


g_prog=$(basename $0);
g_progdir=$(dirname $0);

# ---- usage

usage() {
  local exitstat=2;
  if [[ ! -z "$2" ]] ; then
      exitstat=$2;
  fi

  # Only redirect to stderr on non-zero exit status
  if [[ $exitstat -ne 0 ]] ; then
      exec 1>&2;
  fi

  if [[ ! -z "$1" ]] ; then
      echo "$g_prog: Error! $1" 1>&2;
  fi

  echo "Usage: $g_prog [--help] \\"
#  echo "              -t|--target buildtarget";
#  echo "         -t|--target     -- chef target to build (e.g. 'cookbookname::build')";
  echo "         --help          -- print this usage";
  echo "";

  # bash only:
  if [[ $exitstat -ne 0 ]] ; then
      echo "  at: $(caller)";
  fi
  exit $exitstat;
}

# ---- argument parsing

# Save off the original args, use as "${g_origargs[@]}" (with double quotes)
declare -a g_origargs;
g_origargs=( ${1+"$@"} )

opt_target_exist_check=false;
opt_no_sub_targets=false;
opt_process_all_deps=false;
declare -a opt_additional_options;
declare -a opt_targets;
while [[ $# != 0 ]]; do
    opt="$1"; shift;
    case "$opt" in
	# Flag with no argument example:
	#   --flag|--fla|--fl|--f)
	#     opt_flag=true;;
	# Option with argument example:
	#   --arg|--ar|--a)
	#     [[ $# -eq 0 ]] && usage;
	#     opt_somearg=$1; shift;;
	-e|--exists|--exist-check|--target-exist-check) opt_target_exist_check=true;;
	-s|--no-sub|--no-subs|--no-sub-targets|--single) opt_no_sub_targets=true;;
	-d|--deps|--process-all-deps|--all-deps|-all) opt_process_all_deps=true;;
	-o) 
	    [[ $# -eq 0 ]] && usage;
	    opt_additional_options=( "${opt_additional_options[@]}" "$1" );
	    shift;;
	-h|-help|--help|--hel|--he|--h) usage "" 0;;
	--*) opt_targets=( "${opt_targets[@]}" "$opt" );;
	-*) usage "Unrecognized option: $opt";;
	*)  usage "Extra trailing arguments: $opt $@";;
    esac
done

# ---- error functions
merror() {
    echo "$g_prog: Error! ""$@" 1>&2;
    exit 1;
}
minterror() {
    echo "$g_prog: Internal Error! ""$@" 1>&2;
    exit 1;
}
mwarn() {
    echo "$g_prog: Warning! ""$@" 1>&2;
}

# ---- globals

# ---- subroutines

munge_target() {
    local target=$1; shift;
    local mtarget=$target;
    
    mtarget=${mtarget#--}
    mtarget=${mtarget//-/_}
    echo "$mtarget"
}

# ---- main

set_globals;

warnings=false;
for target in "${opt_targets[@]}"; do
    mtarget=$(munge_target "$target");
    if ! declare -f -F "$mtarget" > /dev/null; then
	if $opt_strict; then
	    mwarn "target '$target' does not exist"
	    warnings=true;
	else
	    echo "$g_prog: target '$target' does not exist"
	fi
    fi
done
if $warnings; then
    merror "Detected warnings, bailing out..."
fi	

if ! $opt_target_exist_check; then
    for target in "${opt_targets[@]}"; do
	mtarget=$(munge_target "$target");
	eval "$mtarget" "${opt_additional_options[@]}"
    done
fi

exit 0;

