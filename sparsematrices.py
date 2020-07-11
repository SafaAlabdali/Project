# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 09:17:03 2019

@author: E1260826
"""

import numpy as np
from scipy.sparse import csr_matrix
import scipy.sparse.linalg as sp
import matplotlib.pyplot as plt

plt.close('all')

#Open file to read as binary ('rb')
f=open('dumpMatrix000.ymat','rb')
#Construct an array of data from binary file, interpretting as 32 bit integers (np.uint32)
a=np.fromfile(f,dtype=np.uint32)

#Use the given information about what is contained within the file to read the data
#Note that the data is not given in standard CSR (compressed sparse row) form - it must be converted later
Ns=a[0]         #System size
na=a[1]         #Number of diagonal entries
nija=a[2]       #Size of topology data (size of matrix values val and row pointers/column indices ija)
nnz=nija-na-1   #Number of non-zero off-diagonal entries

ija=a[3:3+nija]                                 #First na+1 are row pointers, last nnz are column indices for off-diagonal entries
val=np.uint64(a[3+nija:3+nija+nija])            #Matrix values: first na are diagonals, next is a dummy, last nnz are offdiagonals corresponding to column indices
RHS=np.uint64(a[3+nija+nija:3+nija+nija+na])    #RHS of the matrix equation

#If all this has been read correctly, then the following two numbers should be equal...
print('Does',len(RHS)+len(val)+len(ija)+3,'=',3+nija+nija+na,'?')
#...which they are :)

#Converting into CSR form
#We need 3 arrays: value, index, pointer
#These are calculated from ija and val as follows:

#Matrix values, length=total number of non-zero entries=na+nnz
value=[]
for k in range(na):
    value=np.concatenate((value,[val[k]]),axis=None)
    value=np.concatenate((value,val[ija[k]:ija[k+1]]),axis=None)

#Column indices, length=total number of non-zero entries=na+nnz
index=[]
for k in range(na):
    index=np.concatenate((index,[k]),axis=None)
    index=np.concatenate((index,ija[ija[k]:ija[k+1]]),axis=None)

#Row pointers, length=number of rows (ie. number of diagonal entries) + 1
pointer=[]
for k in range(na+1):
    pointer=np.concatenate((pointer,[k+ija[k]-ija[0]]),axis=None)

#Putting the matrix in CSR form    
csr=csr_matrix((value,index,pointer))

#Some tests on the matrix:

#Norm
norm=sp.norm(csr)
print('Frobenius norm =',norm)

#Lower bound of 1-norm
onenorm=sp.onenormest(csr)
print('Lower bound of 1-norm =',onenorm)

#Solving the system
#Note that sp.spsolve doesn't seem to work. It says that 'colind and rwoptr must be of type cint', so use lsmr instead
soln=sp.lsmr(csr,RHS)   #Least squares approx
print('Least squares solution =',soln[0])

#Condition number
print('Condition number =',soln[6])

#Plotting the sparsity pattern
fig=plt.figure()
ax=fig.add_subplot(111)
plot=ax.spy(csr,markersize=1)
plt.title('Sparsity pattern')