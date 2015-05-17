import scipy as sp
from mycwt import cwt

pi = sp.pi
mu = [100.0, 500.0, 900.0]
sigma = [5.0, 10.0, 20.0]
a = [3.0, 1.0, 0.5]
t = sp.arange(0, 1000)*1.0
x = sp.zeros(t.shape)
for i in range(0, len(mu)):
    x += 1/sp.sqrt(2*pi)/sigma[i]*sp.exp(-0.5*((t-mu[i])/sigma[i])**2)

smax = 128
wname = 'bior2.6'
scales = sp.arange(1, smax+1)*1.0
coefs = cwt(x, scales, wname, bplot=True)
