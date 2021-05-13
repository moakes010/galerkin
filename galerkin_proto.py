'''
Galerkin Method Example
'''
import numpy as np
import scipy.integrate as integrate
from functools import partial
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
 
mpl.style.use('seaborn-paper')
'''
class GalerkinSolver:
   def init(self, gfunc, bfunc, operator):
       self.initialzed = false
       self.bfunc = bfunc
       self.gfunc = gfunc
       self.iterations = 3
 
   def solve(self):
 
       gm = np.zeros((N,1))
       Imn = np.zeros((N,N))
       for idx in range(N):
           integrand = self.gfunc
           gm[idx,:] = sp.integrate((x - x**((idx + 1)+1))*G, (x, 0, 1))
'''
      
def exact_solution(x):
   return (5/6)*x -x**2/2 - x**4/3
def approx_solution(x , cs):
   ret=0
   for idx, cval in enumerate(cs):
       N=idx+1
       ret+=cval*(x-x**(N+1))
 
   return ret
if __name__ == "__main__":
   x = sp.Symbol('x')
   n = sp.Symbol('n')
   m = sp.Symbol('m')
   G = 1 + 4*x**2
   '''
   w_n = z - z**(n+1)
   w_m = z - z**(m+1)
   dz_dx = sp.simplify(w_n.diff('z'))
  
   gm = sp.simplify(sp.integrate(F*w_n, ('z', 0, 1)))
   inn = w_m*d2Fn_dx2
   print(sp.simplify(inn))
   Imn = sp.simplify(sp.integrate(inn, ('z', 0, 1)))
   print(gm)
   print(Imn)
   '''
 
   #w_m = lambda w: w - w**(N+1)
   #w_n = lambda y: M*(M+1)*y**(M-1)
   w_m = lambda w: w**N
   w_n = lambda y: n(1-n)*y**(M-2)
   G = lambda x: 1 + 4*x**2
  
   T = 4
   gm = np.zeros((T,1))
   Imn = np.zeros((T,T))
   for idx in range(T):
       #pidx = idx +1
       #w_n = lambda w: w - w**(pidx+1)
       w_n = lambda w: w**idx
       F = lambda x: w_n(x)*G(x)
       gm[idx,:] = integrate.quad(F, 0, 1)[0]
   w_n = x - x**(n+1)
   w_m = x - x**(m+1)
 
   T=4
   #print(d2Fn_dx2)
   #d2Fn_dx2 = n*x**n + x**n - 1
   for j in range(T):
       for i in range(T):
           #M = i+1
           #N = j+1
           #w_m = lambda w: w - w**(N+1)
           #w_n = lambda y: M*(M+1)*y**(M-1)
           w_m = lambda w: w**(j)
           w_n = lambda y: i*(1-i)*y**(i-2)
           F = lambda x: w_m(x)*w_n(x)
           out = integrate.quad(F, 0, 1)
           print(out)
           Imn[i,j] = out[0]
   #Cn*Imn = gm
   print(Imn)
   Cn = np.dot(np.linalg.inv(Imn),gm)
   print(Cn)
 
   xs = np.linspace(0,1,100)
   plt.suptitle("Galerkin Method")
   plt.title("N=2")
   plt.plot(xs, exact_solution(xs))
   plt.plot(xs, approx_solution(xs, Cn.flatten()))
 
   err = np.linalg.norm(approx_solution(xs, Cn.flatten()) - exact_solution(xs))/np.linalg.norm(exact_solution(xs))
 
   print(err)
   plt.grid()
   plt.show()
