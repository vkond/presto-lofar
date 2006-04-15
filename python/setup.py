from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc, get_python_lib
import os
import sys
import numpy

###################################################################
# build the extension
#

define_macros = []
undef_macros = []
extra_compile_args = []
include_dirs = []

ppgplot_libraries = ["cpgplot", "pgplot", "X11", "png", "m", "g2c"]
ppgplot_library_dirs = ["/usr/X11R6/lib"]

presto_libraries = ["presto", "fftw3f", "m"]
presto_library_dirs = []

# Use NumPy instead of Numeric or numarray
make_extension = Extension
include_dirs.append(numpy.get_numpy_include())
ppgplot_include_dirs = include_dirs
presto_include_dirs = include_dirs
undef_macros.append('USE_NUMARRAY')
print >>sys.stderr, "using Scipy_core..."

if os.name == "posix":
    if os.environ.has_key("PGPLOT_DIR"):
        ppgplot_library_dirs.append(os.environ["PGPLOT_DIR"])
        ppgplot_include_dirs.append(os.environ["PGPLOT_DIR"])
    else:
        print >>sys.stderr, "PGPLOT_DIR env var not defined!"
    if os.environ.has_key("PRESTO"):
        presto_library_dirs.append(os.path.join(os.environ["PRESTO"], "lib"))
        presto_include_dirs.append(os.path.join(os.environ["PRESTO"], "include"))
    else:
        print >>sys.stderr, "PRESTO env var not defined!"
else:
    raise Exception, "os not supported"

ext_ppgplot = make_extension('_ppgplot',
			[os.path.join('ppgplot_src','_ppgplot.c')],
			include_dirs=ppgplot_include_dirs,
			libraries=ppgplot_libraries,
			library_dirs=ppgplot_library_dirs,
			define_macros=define_macros,
			extra_compile_args=extra_compile_args)

ext_presto = make_extension('_presto',
			[os.path.join('presto_src', 'presto_wrap.c')],
			include_dirs=presto_include_dirs,
			libraries=presto_libraries,
			library_dirs=presto_library_dirs,
			define_macros=define_macros,
			extra_compile_args=extra_compile_args)



###################################################################
# the package
#

setup(name="PRESTO",
      version="0.99",
      description="Python interfaces to PGPLOT and PRESTO",
      author="Scott Ransom (ppgplot from Nick Patavlis)",
      author_email="ransom@physics.mcgill.ca",
      packages=['ppgplot', "presto"],
      package_dir={'ppgplot':'ppgplot_src','presto':'presto_src'},
      ext_modules=[ext_ppgplot,ext_presto])
