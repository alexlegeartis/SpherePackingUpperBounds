"""Complete-domain Bernstein proof with nonnegative local Gram multiplier.
Acceptance H+lambda*G<0 and G>=0 implies H<0. Lambda is a nonnegative
binary64 value, interpreted exactly as a dyadic rational; its selection is
heuristic only, while acceptance is entirely outward-rounded arithmetic.
"""
if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

from bernstein import *

def gram_init(d):
    n=d+1;delta=CAP+1
    P=[-1+delta*fmpq(i,d) for i in range(n)]
    S=[1-2*delta*fmpq(i,d)+delta**2*fmpq(i*(i-1),d*(d-1)) for i in range(n)]
    gl=np.empty((n,n,n));gh=np.empty_like(gl)
    for i,j,k in np.ndindex(gl.shape):
        val=arb(1+2*P[i]*P[j]*P[k]-S[i]-S[j]-S[k]);gl[i,j,k]=np.nextafter(float(val.lower()),-np.inf);gh[i,j,k]=np.nextafter(float(val.upper()),np.inf)
    return gl,gh

def gram_combination(h,gh):
    # Candidate lambda interval from upper coefficients. This computation does
    # not establish proof; final multiplication and sum are rounded upward.
    neg=gh<0;pos=gh>0;zero=~(neg|pos)
    if np.any(h[zero]>=0):return None
    low=max(0.,float(np.max(h[neg]/(-gh[neg])))) if np.any(neg) else 0.
    high=float(np.min(-h[pos]/gh[pos])) if np.any(pos) else np.inf
    if not (low<high and np.isfinite(low)):return None
    lam=(low+high)*.5 if np.isfinite(high) else np.nextafter(low,np.inf)
    if not (lam>=0 and np.isfinite(lam)):return None
    hu=np.nextafter(h+np.nextafter(lam*gh,np.inf),np.inf)
    if np.max(hu)<0:return lam,float(np.max(hu))
    return None

def run_gram(limit=200000,seconds=600):
    lo,hi=init();gl,gh=gram_init(lo.shape[0]-1)
    stack=[(lo,hi,gl,gh,((Fraction(0),Fraction(1)),)*3,0)]
    checked=passed=passed_multiplier=excluded_sym=excluded_gram=0;maxdepth=0;fail=None;peakstack=1;minmargin=np.inf;maxlam=0.;accepted_volume=Fraction(0);excluded_volume=Fraction(0);transcript=hashlib.sha256();events=[]
    while stack:
        l,h,gl,gh,box,depth=stack.pop();checked+=1;maxdepth=max(maxdepth,depth)
        volume=math.prod(b-a for a,b in box)
        if np.max(h)<0:
            passed+=1;accepted_volume+=volume;minmargin=min(minmargin,-float(np.max(h)));transcript.update(('P'+str(box)).encode());events.append({'kind':'H'});continue
        if box[0][0]>box[1][1] or box[1][0]>box[2][1]:
            excluded_sym+=1;excluded_volume+=volume;transcript.update(('S'+str(box)).encode());events.append({'kind':'S'});continue
        if np.max(gh)<0:
            excluded_gram+=1;excluded_volume+=volume;transcript.update(('G'+str(box)).encode());events.append({'kind':'G'});continue
        proof=gram_combination(h,gh)
        if proof is not None:
            lam,upper=proof;passed_multiplier+=1;accepted_volume+=volume;minmargin=min(minmargin,-upper);maxlam=max(maxlam,lam);transcript.update(('M'+str(box)+repr(lam)).encode());events.append({'kind':'M','lambda_rational':list(lam.as_integer_ratio())});continue
        if checked>=limit or depth>=100 or time.monotonic()-start>seconds:
            fail={'box':[[str(a),str(b)]for a,b in box],'depth':depth,'upper':float(np.max(h)),'remaining':len(stack),'reason':'bounded work limit'};break
        mid=(l+h)*.5;scores=[float(np.max(np.abs(np.diff(mid,axis=i)))) for i in range(3)];axis=int(np.argmax(scores));ll,lh,rl,rh=split(l,h,axis);gll,glh,grl,grh=split(gl,gh,axis)
        events.append({'kind':'split','axis':axis});a,b=box[axis];m=(a+b)/2;bx1=list(box);bx2=list(box);bx1[axis]=(a,m);bx2[axis]=(m,b)
        stack.extend([(rl,rh,grl,grh,tuple(bx2),depth+1),(ll,lh,gll,glh,tuple(bx1),depth+1)]);peakstack=max(peakstack,len(stack))
        if checked%1000==0:print('PROGRESS',checked,'passed',passed,'multiplier',passed_multiplier,'gram',excluded_gram,'sym',excluded_sym,'stack',len(stack),'depth',maxdepth,'seconds',time.monotonic()-start,flush=True)
    verified=fail is None and not stack
    if verified:assert accepted_volume+excluded_volume==1
    out={'verified':verified,'source_sha256':SHAS,'checked':checked,'passed_direct':passed,'passed_gram_multiplier':passed_multiplier,'excluded_symmetry':excluded_sym,'excluded_gram':excluded_gram,'max_depth':maxdepth,'peak_stack':peakstack,'minimum_acceptance_margin':minmargin,'largest_nonnegative_multiplier':maxlam,'covered_unit_cube_volume_exact':str(accepted_volume+excluded_volume),'accepted_volume_exact':str(accepted_volume),'excluded_volume_exact':str(excluded_volume),'partition_transcript_sha256':transcript.hexdigest(),'seconds':time.monotonic()-start,'failure':fail,'method':'Exact rational cyclic polynomial and tensor-Bernstein initialization, coordinate degree16. Exact rational degree-elevated Gram polynomial G. Outward binary64 dyadic de Casteljau subdivision. Accept H<0 or H+lambda*G<0 for exactly nonnegative dyadic lambda, which implies H<0 wherever G>=0. Exclude only strictly Gram-negative boxes or boxes strictly outside sorted chamber. Exact permutation symmetry extends to full cube. No samples used for acceptance.'}
    (ROOT/'proof_events.json').write_text(json.dumps({'input_sha256':SHAS,'polynomial_sha256':hashlib.sha256((ROOT/'polynomial.json').read_bytes()).hexdigest(),'events':events},separators=(',',':')))
    (ROOT/'bernstein_gram_result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
if __name__=='__main__':run_gram(int(sys.argv[1]) if len(sys.argv)>1 else 200000,600)
