#!/usr/bin/env python
"""Configure the build.

- Create defines.mk
"""
import commands
import contextlib
import optparse
import os
import sys

#DEFAULTCXXFLAG := -O3
#DEBUGCXXFLAG := -g -ggdb -fno-inline
#PROFILECXXFLAG := -Os -pg 
#GCXXFLAG := -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 

ROOT = os.path.abspath(os.path.dirname(__file__))

def log(msg):
    sys.stderr.write(msg)
    sys.stderr.write('\n')

def shell(cmd):
    log(cmd)
    status, output = commands.getstatusoutput(cmd)
    if status:
        raise Exception('%d <- %r' %(status, cmd))
    return output

def system(cmd):
    log(cmd)
    status = os.system(cmd)
    if status:
        raise Exception('%d <- %r' %(status, cmd))
    return

def mkdirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

@contextlib.contextmanager
def cd(nwd):
    cwd = os.getcwd()
    log('cd %r -> %r' %(cwd, nwd))
    os.chdir(nwd)
    yield
    os.chdir(cwd)
    log('cd %r <- %r' %(cwd, nwd))

def update_content(fn, content):
    current_content = open(fn).read() if os.path.exists(fn) else None
    if content != current_content:
        log('writing to %r:' %fn)
        log('"""\n' + content + '"""\n----')
        open(fn, 'w').write(content)

def get_OS_STRING():
    G_BUILDOS_CMD = """bash -c 'set -e; set -o pipefail; id=$(lsb_release -si | tr "[:upper:]" "[:lower:]"); rel=$(lsb_release -sr); case $id in ubuntu) printf "$id-%04d\n" ${rel/./};; centos) echo "$id-${rel%%.*}";; *) echo "$id-$rel";; esac' 2>/dev/null"""
    return shell(G_BUILDOS_CMD)
def get_PREBUILT():
    cmd = 'cd ../../../../prebuilt.out 2>/dev/null && pwd || echo -n notfound'
    return shell(cmd)
def ifenvf(env, key, func):
    if key in env:
        return env[key]
    else:
        return func()
def setifenvf(envout, envin, key, func):
    envout[key] = ifenvf(envin, key, func)
def setifenv(envout, envin, key, val):
    envout[key] = envin.get(key, val)
def setenv(envout, key, val):
    envout[key] = val
def update_env_if(envout, envin, keys):
    for key in keys:
        if key in envin:
            envout[key] = envin[key]
def compose_defs_env(env):
    # We disallow env overrides for anything with a default from GNU make.
    nons = ['CXX', 'CC', 'AR'] # 'SHELL'?
    ovr    = ['%-20s ?= %s' %(k, v) for k,v in env.items() if k not in nons]
    nonovr = ['%-20s := %s' %(k, v) for k,v in env.items() if k in nons]
    return '\n'.join(ovr + nonovr + [''])
def compose_defines_pacbio(envin):
    """
    This is used by mobs via buildcntl.sh.
    """
    env = dict()
    setenv(env, 'SHELL', 'bash')
    #setifenvf(env, envin, 'OS_STRING', get_OS_STRING)
    #setifenvf(env, envin, 'PREBUILT', get_PREBUILT)
    nondefaults = set([
            'CXX',
            'LIBPBDATA_INCLUDE', 'LIBPBDATA_LIB', 'LIBPBDATA_LIBFLAGS',
            'LIBPBIHDF_INCLUDE', 'LIBPBIHDF_LIB', 'LIBPBIHDF_LIBFLAGS',
            'LIBBLASR_INCLUDE', 'LIBBLASR_LIB', 'LIBBLASR_LIBFLAGS',
            'HDF5_INCLUDE', 'HDF5_LIB', 'HDF5_LIBFLAGS',
            'PBBAM_INCLUDE', 'PBBAM_LIB', 'PBBAM_LIBFLAGS',
            'HTSLIB_INCLUDE', 'HTSLIB_LIB', 'HTSLIB_LIBFLAGS',
            'BOOST_INCLUDE',
            'ZLIB_LIB', 'ZLIB_LIBFLAGS',
            'PTHREAD_LIBFLAGS',
            'DL_LIBFLAGS',
    ])
    update_env_if(env, envin, nondefaults)
    return compose_defs_env(env)


def configure_pacbio(envin, shared, build_dir):
    content1 = compose_defines_pacbio(envin)
    if shared:
        content1 += 'LDLIBS+=-lrt\n'
    else:
        content1 += 'LDFLAGS+=-static\n'
    content1 += 'SUB_CONF_FLAGS+=--shared\n'
    update_content(os.path.join(build_dir, 'defines.mk'), content1)

def set_defs_submodule_defaults(env, nopbbam):
    subdir = os.path.join(ROOT, 'libcpp')
    defaults = {
        'LIBPBDATA_INCLUDE': os.path.join(subdir, 'pbdata'),
        'LIBBLASR_INCLUDE':  os.path.join(subdir, 'alignment'),
        #'LIBPBIHDF_INCLUDE': '' if nopbbam else os.path.join(subdir, 'hdf'),
        'LIBPBDATA_LIB': os.path.join(subdir, 'pbdata'),
        'LIBBLASR_LIB':  os.path.join(subdir, 'alignment'),
        #'LIBPBIHDF_LIB': '' if nopbbam else os.path.join(subdir, 'hdf'),
    }
    for k in defaults:
        if k not in env:
            env[k] = defaults[k]

def set_defs_defaults(env, nopbbam):
    # OS := $(shell uname)
    # if Darwin, -lsz (for static builds?)
    defaults = {
        'LIBBLASR_INCLUDE':  os.path.join(ROOT, 'libcpp', 'alignment'),
        'LIBPBDATA_INCLUDE':  os.path.join(ROOT, 'libcpp', 'pbdata'),
        'LIBPBIHDF_INCLUDE':  os.path.join(ROOT, 'libcpp', 'hdf'),
        'LIBBLASR_LIB':  os.path.join(ROOT, 'libcpp', 'alignment'),
        'LIBPBDATA_LIB':  os.path.join(ROOT, 'libcpp', 'pbdata'),
        'LIBPBIHDF_LIB':  os.path.join(ROOT, 'libcpp', 'hdf'),
        'LIBBLASR_LIBFLAGS':  '-lblasr',
        'LIBPBDATA_LIBFLAGS': '-lpbdata',
        'LIBPBIHDF_LIBFLAGS': '-lpbihdf',
        'HDF5_LIBFLAGS': '-lhdf5_cpp -lhdf5',
        'ZLIB_LIBFLAGS': '-lz',
        'PTHREAD_LIBFLAGS': '-lpthread',
        'DL_LIBFLAGS': '-ldl', # neeeded by HDF5 always
        'SHELL': 'bash -xe',
    }
    #setifenvf(defaults, env, 'OS_STRING', get_OS_STRING)
    #setifenvf(defaults, env, 'PREBUILT', get_PREBUILT)
    pbbam_defaults = {
        'PBBAM_LIBFLAGS': '-lpbbam',
        'HTSLIB_LIBFLAGS': '-lhts',
        'ZLIB_LIBFLAGS': '-lz',
        'PTHREAD_LIBFLAGS': '-lpthread',
        'DL_LIBFLAGS': '-ldl', # neeeded by HDF5 always
    }
    if not nopbbam:
        defaults.update(pbbam_defaults)
    for k in defaults:
        if k not in env:
            env[k] = defaults[k]

def get_make_style_env(envin, args):
    envout = dict()
    for arg in args:
        if '=' in arg:
            k, v = arg.split('=')
            envout[k] = v
    envout.update(envin)
    return envout

def parse_args(args):
    parser = optparse.OptionParser()
    parser.add_option('--no-pbbam', action='store_true',
            help='Avoid compiling anything which would need pbbam.')
    parser.add_option('--submodules', action='store_true',
            help='Set variables to use our git-submodules, which must be pulled and built first. (Implies --no-pbbam.)')
    parser.add_option('--shared', action='store_true',
            help='Build for dynamic linking.')
    parser.add_option('--build-dir',
            help='Can be different from source directory, but only when *not* also building submodule.')
    return parser.parse_args(list(args))

def write_makefile(build_dir_root, src_dir_root, makefilename, relpath):
    src_dir = os.path.join(src_dir_root, relpath)
    build_dir = os.path.join(build_dir_root, relpath)
    content = """\
include %(src_dir)s/%(makefilename)s
VPATH:=%(src_dir)s
""" %dict(makefilename=makefilename, src_dir=src_dir)
    mkdirs(build_dir)
    fn = os.path.join(build_dir, makefilename)
    update_content(fn, content)

def write_makefiles(build_dir):
    write_makefile(build_dir, ROOT, 'makefile', '.')
    write_makefile(build_dir, ROOT, 'makefile', 'utils')
    write_makefile(build_dir, ROOT, 'makefile', 'extrautils')

def main(prog, *args):
    """We are still deciding what env-vars to use, if any.
    """
    # Set up an alias, until everything uses one consistently.
    if 'HDF5_INC' in os.environ and 'HDF5_INCLUDE' not in os.environ:
        os.environ['HDF5_INCLUDE'] = os.environ['HDF5_INC']
    conf, makevars = parse_args(args)
    if conf.build_dir is not None:
        write_makefiles(conf.build_dir)
    else:
        conf.build_dir = '.'
    conf.build_dir = os.path.abspath(conf.build_dir)
    envin = get_make_style_env(os.environ, makevars)
    if conf.submodules:
        set_defs_submodule_defaults(envin, conf.no_pbbam)
        conf.no_pbbam = True
    set_defs_defaults(envin, conf.no_pbbam)
    configure_pacbio(envin, conf.shared, conf.build_dir)


if __name__=="__main__":
    main(*sys.argv)
