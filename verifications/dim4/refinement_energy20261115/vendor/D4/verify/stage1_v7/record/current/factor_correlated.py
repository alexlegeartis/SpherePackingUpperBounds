"""Optional exact correlated-denominator Bernstein acceptance.

For a branch numerator P and positive factor H, compare each Bernstein
coefficient of P with correction_lower times the corresponding coefficient
of H.  It is a sufficient test only.  An inconclusive comparison returns
False and the caller retains every existing subdivision and failure gate.
"""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
import ctypes,hashlib,math,subprocess,tempfile
import numpy as np

ROOT=Path(__file__).resolve().parent
WRAPPER=Path(__file__).resolve()
WRAPPER_SHA256=hashlib.sha256(WRAPPER.read_bytes()).hexdigest()
SOURCE=ROOT/'factor_correlated.cpp'
SOURCE_SHA256=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
BUILD_DIRECTORY=Path(tempfile.mkdtemp(prefix='factor_correlated_build_',dir=ROOT))
LIBRARY=BUILD_DIRECTORY/'factor_correlated.so'
subprocess.run(['g++','-std=c++20','-O3','-march=native','-shared','-fPIC',str(SOURCE),'-o',str(LIBRARY)],check=True)
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA256
LIBRARY_SHA256=hashlib.sha256(LIBRARY.read_bytes()).hexdigest()
_lib=ctypes.CDLL(str(LIBRARY));_f=_lib.factor_correlated_i64
_iptr=ctypes.POINTER(ctypes.c_int64);_uptr=ctypes.POINTER(ctypes.c_uint64)
_f.argtypes=[_iptr,_uptr,_uptr,_uptr,ctypes.c_size_t,ctypes.c_size_t,ctypes.c_size_t,ctypes.c_uint64]
_f.restype=ctypes.c_int
BITS=20;UNIT=1<<BITS

def exact_q(value):
    assert isinstance(value,(int,Q)) and not isinstance(value,bool)
    return Q(int(value.numerator),int(value.denominator)) if isinstance(value,Q) else Q(value)

@lru_cache(maxsize=256)
def _axis_integer_weights(n,L):
    assert type(n) is int and type(L) is int and 0<=L<=n
    H=math.lcm(*(math.comb(n,j) for j in range(L+1)))
    weights=tuple(tuple(math.comb(L,j)*math.comb(i,j)*(H//math.comb(n,j))
                        for j in range(min(i,L)+1)) for i in range(n+1))
    return weights,H

@lru_cache(maxsize=65536)
def normalized_axis_coefficients(n,L,left,right):
    """floor(2^20 times degree-n Bernstein coeffs of (x/right)^L).

    x=left+(right-left)t.  With left/right=A/D, the exact numerator is
    sum_j binom(L,j)binom(i,j)H/binom(n,j)*A^(L-j)*(D-A)^j,
    over the positive denominator D^L H.  Every term is an integer.
    """
    left,right=exact_q(left),exact_q(right)
    assert 0<left<=right
    weights,H=_axis_integer_weights(n,L)
    ratio=left/right;A,D=ratio.numerator,ratio.denominator
    monomial=[A**(L-j)*(D-A)**j for j in range(L+1)]
    denominator=D**L*H
    result=tuple(UNIT*sum(w*monomial[j] for j,w in enumerate(row))//denominator
                 for row in weights)
    assert all(type(x) is int and 0<=x<=UNIT for x in result)
    return result

def prepare(high,scale,bounds,kind,L,correction_lower):
    assert isinstance(high,np.ndarray) and high.dtype==np.dtype(np.int64) and high.ndim==3
    assert all(x>0 for x in high.shape)
    assert int(high.min())>-(1<<59) and int(high.max())<(1<<59)
    scale=exact_q(scale);correction_lower=exact_q(correction_lower)
    assert scale>0 and correction_lower>=0
    assert kind in ['full','outside'] and type(L) is int and L>=0
    assert len(bounds)==3
    bounds=tuple((exact_q(a),exact_q(b)) for a,b in bounds)
    assert all(a<=b for a,b in bounds)
    active=3 if kind=='full' else 2
    assert all(bounds[i][0]>0 for i in range(active))
    factor=Q(2)
    for i in range(active):factor*=bounds[i][1]
    Hmax=factor**L
    available=scale*correction_lower*Hmax
    T=min((1<<59)-1,available.numerator//available.denominator)
    q=[]
    for i in range(3):
        values=normalized_axis_coefficients(high.shape[i]-1,L,*bounds[i]) if i<active else (UNIT,)*high.shape[i]
        q.append(np.array(values,dtype=np.uint64))
    return np.ascontiguousarray(high),T,q

def accept(high,scale,bounds,kind,L,correction_lower):
    """True proves P/H<=correction_lower throughout the entire closed box."""
    high,T,q=prepare(high,scale,bounds,kind,L,correction_lower)
    answer=_f(high.ctypes.data_as(_iptr),*(x.ctypes.data_as(_uptr) for x in q),
              *high.shape,T)
    assert answer in [0,1]
    return answer==1

def reference_accept(high,scale,bounds,kind,L,correction_lower):
    """Independent exact Python-integer comparison for controls."""
    high,T,q=prepare(high,scale,bounds,kind,L,correction_lower)
    for index in np.ndindex(high.shape):
        h=int(high[index])
        if h>0 and h*(1<<60)>T*math.prod(int(q[j][index[j]]) for j in range(3)):
            return False
    return True

def clear_cache():
    normalized_axis_coefficients.cache_clear()
    _axis_integer_weights.cache_clear()

def source_contract():
    assert hashlib.sha256(WRAPPER.read_bytes()).hexdigest()==WRAPPER_SHA256
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA256
    assert hashlib.sha256(LIBRARY.read_bytes()).hexdigest()==LIBRARY_SHA256
    return dict(wrapper=str(WRAPPER),wrapper_sha256=WRAPPER_SHA256,
                source=str(SOURCE),source_sha256=SOURCE_SHA256,
                library=str(LIBRARY),library_sha256=LIBRARY_SHA256,
                normalized_axis_bits=BITS,
                comparison='hi*2^60 <= T*q0*q1*q2 for every positive numerator coefficient',
                role='Additional sufficient sign test only; all old coverage and failure gates remain required.')
