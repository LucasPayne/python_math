import numpy as np

def make_rotation_matrix(u, v):
    u = u/np.linalg.norm(u)
    v = v - u*np.dot(v, u)
    v = v/np.linalg.norm(v)
    w = np.cross(u, v)
    w = w/np.linalg.norm(w) #-maybe to prevent error
    return np.vstack((u,v,w))

Q = make_rotation_matrix(np.array([1,0,0]),np.array([2,3,4]))
R = make_rotation_matrix(np.array([-1,1,0]),np.array([4,1,0]))
print("Q =", Q)
print("R =", R)
A = np.linalg.inv(Q).dot(R)
evals, evecs = np.linalg.eig(A)
print(A)
print(np.diag(evals))
print(evecs.dot(np.diag(evals)).dot(np.linalg.inv(evecs)))

for t in np.linspace(0,1,4):
    print("t = {}".format(t))
    evalpows = np.exp(np.log(evals)*t)
    Apow = evecs.dot(np.diag(evalpows)).dot(np.linalg.inv(evecs))
    print(Q.dot(Apow))
