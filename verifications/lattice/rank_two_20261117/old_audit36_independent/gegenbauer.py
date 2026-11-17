from fractions import Fraction as F

def gegen(n,k):
    a=F(n-2,2); prev=[F(1)]; cur=[F(0),2*a]
    if k==0:return prev
    for j in range(2,k+1):
        out=[F(0)]*(j+1)
        for i,x in enumerate(cur):out[i+1]+=2*(j+a-1)*x/j
        for i,x in enumerate(prev):out[i]-=(j+2*a-2)*x/j
        prev,cur=cur,out
    den=sum(cur)
    return [x/den for x in cur]
