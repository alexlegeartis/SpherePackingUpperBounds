"""Exact compiled integer subdivision; no sign or coverage decisions.

The C routine implements the same recurrence as certificate_core.split.
Compilation is local from the included source. Inputs and intermediates
are bounded below2^59; only signed integer arithmetic is performed.
"""
from pathlib import Path
import ctypes, hashlib, subprocess, tempfile
import numpy as np

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'fast_split.cpp'
SOURCE_SHA256=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
# Build anew for each independent parent process. A cached binary is never
# accepted based only on its filename or the source digest in that name.
BUILD_DIRECTORY=Path(tempfile.mkdtemp(prefix='fast_split_build_',dir=ROOT))
LIBRARY=BUILD_DIRECTORY/'fast_split.so'
subprocess.run(['g++','-std=c++20','-O3','-march=native','-shared','-fPIC',str(SOURCE),'-o',str(LIBRARY)],check=True)
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA256
LIBRARY_SHA256=hashlib.sha256(LIBRARY.read_bytes()).hexdigest()
_lib=ctypes.CDLL(str(LIBRARY))
_f=_lib.split_bernstein_i64
_ptr=ctypes.POINTER(ctypes.c_int64)
_f.argtypes=[_ptr,_ptr,_ptr,ctypes.c_size_t,ctypes.c_size_t,ctypes.c_int]
_f.restype=ctypes.c_int

def split(a,axis,upper):
    assert isinstance(a,np.ndarray) and a.dtype==np.dtype(np.int64)
    assert type(axis) is int and 0<=axis<a.ndim and type(upper) is bool
    assert a.size>0 and all(x>0 for x in a.shape)
    assert int(a.min())>-(1<<59) and int(a.max())<(1<<59)
    w=np.array(np.moveaxis(a,axis,0),dtype=np.int64,order='C',copy=True)
    rows=w.shape[0];columns=w.size//rows
    left=np.empty_like(w);right=np.empty_like(w)
    result=_f(w.ctypes.data_as(_ptr),left.ctypes.data_as(_ptr),right.ctypes.data_as(_ptr),rows,columns,int(upper))
    assert result==0
    return np.moveaxis(left,0,axis),np.moveaxis(right,0,axis)

def source_contract():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA256
    assert hashlib.sha256(LIBRARY.read_bytes()).hexdigest()==LIBRARY_SHA256
    return {'source':str(SOURCE),'source_sha256':SOURCE_SHA256,
            'library':str(LIBRARY),'library_sha256':LIBRARY_SHA256,
            'arithmetic':'signed int64, exact floor/ceil midpoint averages',
            'role':'subdivision only; no sign or coverage decisions'}
