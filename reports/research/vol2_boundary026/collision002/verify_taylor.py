from pathlib import Path
import hashlib,json
prior=Path('reports/research/vol2_boundary026/verify_quantum_contraction.py')
exec(prior.read_text().split('counts = []')[0])
covariance={(0,1):F(1),(0,2):F(1,9),(1,2):F(1,4)}
K,D,expn,projection,H=operators(covariance)
normal=lambda a:expn(a,-1)
def parity(a):
 degrees={key[2].bit_count()%2 for key in a}
 assert len(degrees)<=1
 return next(iter(degrees),0)
def f2(a,b):return add(normal(multiply(a,b)),scale(multiply(normal(a),normal(b)),-1))
def bracket(a,b):return add(D(multiply(a,b)),scale(multiply(D(a),b),-1),scale(multiply(a,D(b)),-(-1)**parity(a)))
def f3(a,b,c):
 return add(normal(multiply(multiply(a,b),c)),scale(multiply(normal(multiply(a,b)),normal(c)),-1),scale(multiply(normal(multiply(a,c)),normal(b)),-(-1)**(parity(b)*parity(c))),scale(multiply(normal(a),normal(multiply(b,c))),-1),scale(multiply(multiply(normal(a),normal(b)),normal(c)),2))
x=[add(variable(i,'delta'),variable(i,'y')) for i in range(3)]
theta=[basis(mask=1<<i) for i in range(3)]
samples=x+theta+[multiply(x[0],x[0]),multiply(theta[0],x[1]),multiply(theta[0],theta[1])]
c2=c3=0
for a in samples:
 for b in samples:
  left=add(d0(f2(a,b)),scale(f2(D(a),b),-1),scale(f2(a,D(b)),-(-1)**parity(a)))
  assert left==normal(bracket(a,b)),('binary',a,b)
  c2+=1
  for c in samples:
   left=add(d0(f3(a,b,c)),scale(f3(D(a),b,c),-1),scale(f3(a,D(b),c),-(-1)**parity(a)),scale(f3(a,b,D(c)),-(-1)**(parity(a)+parity(b))))
   right=add(f2(bracket(a,b),c),scale(f2(bracket(a,c),b),(-1)**(parity(b)*parity(c))),scale(f2(bracket(b,c),a),(-1)**(parity(a)*(parity(b)+parity(c)))))
   assert left==right,('ternary',a,b,c)
   c3+=1
assert not f3(*x)
assert f3(multiply(x[0],x[0]),x[1],x[2])==scale(basis(hbar=2),2*covariance[0,1]*covariance[0,2])
assert f2(x[0],x[1])==scale(basis(hbar=1),-covariance[0,1])
assert f2(D(theta[0]),x[1])==scale(bracket(theta[0],x[1]),-1)
print(json.dumps({'result':'PASS','binary_identities':c2,'ternary_identities':c3,'arithmetic':'exact rational','python':platform.python_version(),'new_examples':['three linear currents vanish','quadratic-linear-linear cumulant equals 2 hbar^2 p12 p13','binary contact cancellation'],'helper_source_sha256':hashlib.sha256(prior.read_bytes()).hexdigest(),'scope':'Finite checks of the new Taylor identities, not repeated contraction checks.'},indent=2))
