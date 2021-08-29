import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt


def gauss_seidel(A,b,x0,max_iter,tol):
    D = np.diag(np.diag(A))
    L = np.tril(A,k=-1)
    U = np.triu(A,k=1)
    M = D + L
    N = U
    B = -np.dot(npl.inv(M),N)
    C = np.dot(npl.inv(M),b)
    xN = x0
    r = npl.norm(np.dot(A,xN)-b)
    iter = 0
    while r > tol and iter < max_iter:
        xN = np.dot(B,xN) + C
        iter = iter + 1
        r = npl.norm(np.dot(A,xN)-b)    
    if iter==max_iter:
       print('max_iter reached')
    return([xN,iter])


def jacobi(A,b,x0,max_iter,tol):
    D = np.diag(np.diag(A))
    L = np.tril(A,k=-1)
    U = np.triu(A,k=1)
    M = D 
    N = U + L
    B = -np.dot(npl.inv(M),N)
    C = np.dot(npl.inv(M),b)
    xN = x0
    r = npl.norm(np.dot(A,xN)-b) 
    iter = 0
    while r > tol and iter < max_iter:
        xN = np.dot(B,xN) + C
        iter = iter + 1
        r = npl.norm(np.dot(A,xN)-b)    
    if iter==max_iter:
       print('max_iter reached')
    return([xN,iter])
 

b = np.array([[1],[1]])
x0 = np.array([[0],[0]])
max_iter = 10
tol = 10**-10

A1 = np.array([[1,0],[0,2]])
A2 = np.array([[1,1],[0,2]])

#test of the jacobi and gauss seidel functions for matrix A1 and matrix A2.
A = A1
solJ = jacobi(A,b,x0,max_iter,tol)
solGS = gauss_seidel(A,b,x0,max_iter,tol)
print ('solution by Jacobi method for matrix A1:'+str(solJ[0]))
print ('solution by Gauss Seidel method for matrix A1:'+str(solGS[0]))


A = A2
solJ = jacobi(A,b,x0,max_iter,tol)
solGS = gauss_seidel(A,b,x0,max_iter,tol)
print ('solution by Jacobi method for matrix A2:'+str(solJ[0]))
print ('solution by Gauss Seidel method for matrix A2:'+str(solGS[0]))


#In the while loop the condition is det (A) == 0, to compute
#matrix A again if it turns out to be singular.
def matrix_and_vector(n):
    A = np.random.rand(n, n)
    while npl.det(A)==0:   
      A = np.random.rand(n, n)
    b = np.random.rand(n, 1)
    return([A,b])

trials_J = []
trials_GS = []
Ntrials = 10**3
trials = 0
x0 = np.transpose(np.array([[0,0,0]]))
max_iter = 10**6
tol = 10**-5
while trials<Ntrials:
    Ab = matrix_and_vector(3)
    A = Ab [0]
    b = Ab [1]
    rGS= gauss_seidel(A,b,x0,max_iter,tol)
    rJ = jacobi(A,b,x0,max_iter,tol)
    if rJ[1]!=max_iter and rGS[1]!=max_iter:
       trials_J.append(rJ[1])
       trials_GS.append(rGS[1])
       trials = trials + 1
       #print(trials)


plt.scatter(trials_J,trials_GS,c=trials_J,cmap=plt.cm.get_cmap("cool",7))
plt.ylabel('Gauss-Seidel')
plt.xlabel('Jacobi')
plt.yscale('log')
plt.xscale('log')
plt.title('number of steps')
plt.savefig('JvsGS')
plt.show()


