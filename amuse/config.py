#
# src/amuse/config.py.  Generated from configpy.in by configure.
#

class interpreters(object):
    python = r'/home/tgarrow/local/amuse/prerequisites//bin/python3.6' 

class compilers(object):
    cxx = 'g++' 
    cc  = 'gcc'
    fc = 'gfortran'
    
    cxx_flags = '-g -O2 -fPIC'
    cc_flags  = '-I/usr/local/include -fPIC'
    fc_flags = '-g -O2 -fPIC'
    ld_flags = '-L/usr/local/lib'
    
    found_fftw = 'yes'
    fftw_flags = '-I/home/tgarrow/local/amuse/prerequisites//include'
    fftw_libs = '-L/home/tgarrow/local/amuse/prerequisites//lib -lfftw3  -lfftw3_threads'
    
    found_gsl = 'yes'
    gsl_flags = ' '
    gsl_libs = '-lgsl -lgslcblas -lm  '
    
    gfortran_version = '4.8.5'
    ifort_version = ''

    fc_iso_c_bindings = 'yes'=='yes'
    
    cython = '/home/bovy/local/bin/cython'
    pythondev_cflags = '-I/usr/include/python2.7 -I/usr/include/python2.7 -fno-strict-aliasing -O2 -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv -DNDEBUG -O2 -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv'
    pythondev_ldflags = '-lpthread -ldl -lutil -lm -lpython2.7 -Xlinker -export-dynamic'

class mpi(object):
    is_enabled = 'yes'=='yes'
    mpicxx = 'mpicxx' 
    mpicc  = 'mpicc'
    mpif95 = 'mpif90'
    mpiexec = '/home/tgarrow/local/amuse/prerequisites//bin/mpiexec'

class java(object):
    is_enabled = 'yes'=='yes'
    java = '/usr/bin/java'
    javac = ''
    jar = ''

class cuda(object):
    is_enabled   = 'yes'=='yes'
    compiler     = '/usr/local/cuda/bin/nvcc'
    toolkit_path = '/usr/local/cuda-10.0'
    sdk_path     = '@CUDA_SDK@'
    cuda_libs = '-L/usr/lib -lcuda  -L/usr/local/cuda-10.0/lib64 -lcudart'
    sapporo_version = 'light'
    
class openmp(object):
    is_enabled   = 'yes'=='yes'
    fcflags = '-fopenmp'
    cflags = '-fopenmp' 
    
class modules(object):
    have_matplotlib = 0==1
    have_h5py       = 0==1
