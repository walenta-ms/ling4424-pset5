import sys
import math
import numpy as np

modelfile = open(sys.argv[1])
testfile = open(sys.argv[2])

Q = []
O = []
Qindices = {}
Oindices = {}
transes = {}
emits = {}
seq = []

def dcadd(dic, key, amt=1):
    if key in dic:
        if type(dic[key]) == type(0) or type(dic[key]) == type(0.0):
            dic[key] += amt
    else:
        dic[key] = amt
def lshash(ls, elmt):
    if elmt not in ls:
        ls += [elmt]
def initialize(mfile, tfile):
    global Q
    global O
    global Qindices
    global Oindices
    global transes
    global emits
    global seq
    for line in tfile:
        if line.strip() != '':
            seq = line.strip().split()
    # for line in seq:
    #     line[0:0] = ['<s>']
    #     line += ['<f>']
    for line in mfile:
        l = line.strip().split()
        if l[0] == 'T':
            lshash(Q, l[1])
            transes[(l[1], l[2])] = l[-1]
        else:
            lshash(O, l[2])
            emits[(l[1], l[2])] = l[-1]
    for i in range(len(Q)-2):
        Qindices[Q[i+1]] = i
    for i in range(len(O)):
        Oindices[O[i]] = i
initialize(modelfile, testfile)


A = np.zeros((len(Q), len(Q)))
for i in range(len(Q)):
    for j in range(len(Q)):
        A[i, j] = transes[(Q[i], Q[j])]

B = np.zeros((len(Q)-2, len(O)))
for i in range(1, len(Q)-1):
    for j in range(len(O)):
        B[i-1, j] = emits[(Q[i], O[j])]


def forward(seq):
    N = len(Q) - 2
    T = len(seq)
    alpha = np.zeros((N, T))
    for i in range(N):
        alpha[i, 0] = alph(i, 0)
    for j in range(1, T):
        for i in range(N):
            alpha[i, j] = alph(i, j, alpha, True)
    return alpha

def alph(n, t, alpha=None, recursing=False):
    N = len(Qindices)
    if not recursing:
        return A[t, n+1]*B[n, Oindices[seq[t]]]
    else:
        summ = 0
        Bprob = B[n, Oindices[seq[t]]]
        for i in range(N):
            s = alpha[i, t-1] * A[i+1, n+1] * Bprob
            summ += s
        return summ

def terminatef(alpha):
    N = len(alpha)
    summ = 0
    for i in range(N):
        summ += alpha[i, -1]*A[i+1, -1]
    return math.log(summ)

alpha = forward(seq)
fdprob = terminatef(alpha)
out = ''
out += 'Forward:\n'
out += str(fdprob) + '\n'
out += str(alpha) + '\n'

def backward(seq):
    N = len(Q) - 2
    T = len(seq)
    beta = np.zeros((N, T))
    for i in range(N):
        beta[i, -1] = bet(i, -1)
    for j in range(T-2, -1, -1):
    # for j in range(1, T):
        for i in range(N):
            beta[i, j] = bet(i, j, beta, True)
            # beta[i, -(j+1)] = bet(i, j, beta, True)
    return beta

def bet(n, t, beta=None, recursing=False):
    N = len(Qindices)
    if not recursing:
        return A[n+1, -1]
    else:
        summ = 0
        for i in range(N):
            # summ += A[i+1, n+1] * B[i, Oindices[seq[-(t+1)]]] * beta[i, -(t)]
            a = A[n+1, i+1]
            b1 = B[i, Oindices[seq[t+1]]]
            bt1 = beta[i, t+1]
            summ += a*b1*bt1
        return summ

def terminateb(beta):
    N = len(beta)
    summ = 0
    for i in range(N):
        summ += A[0, i+1]*B[i, Oindices[seq[0]]]*beta[i, 0]
    return math.log(summ)

beta = backward(seq)
bdprob = terminateb(beta)
out += str("Backward:") + '\n'
out += str(bdprob) + '\n'
out += str(beta)

f = open('hmm.likelihood', 'w')
f.write(out)
f.close
# print(out)
