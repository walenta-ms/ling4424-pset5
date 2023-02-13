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
seqs = []

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
    global seqs
    for line in tfile:
        if line.strip() != '':
            seqs += [line.strip().split()]
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

A0 = np.zeros((len(Q), len(Q)))
for i in range(len(Q)):
    for j in range(len(Q)):
        A0[i, j] = transes[(Q[i], Q[j])]

B0 = np.zeros((len(Q)-2, len(O)))
for i in range(1, len(Q)-1):
    for j in range(len(O)):
        B0[i-1, j] = emits[(Q[i], O[j])]

# AFTER THIS POINT, DO NOT USE: transes, emits

def forward(seq, A, B, longAlph=False):
    N = len(Q) - 2
    T = len(seq)
    if not longAlph:
        alpha = np.zeros((N, T))
        for i in range(N):
            alpha[i, 0] = alph(seq, i, 0, A, B)
        for j in range(1, T):
            for i in range(N):
                alpha[i, j] = alph(seq, i, j, A, B, alpha, True)
    else:
        alpha = np.zeros((N, T+1))
        for i in range(N):
            alpha[i, 0] = 1
            alpha[i, 1] = alph(seq, i, 0, A, B)
        for j in range(2, T+1):
            for i in range(N):
                alpha[i, j] = alph(seq, i, j-1, A, B, alpha, True)
    return alpha

def alph(seq, n, t, A, B, alpha=None, recursing=False):
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

def terminatef(seq, alpha, A):
    N = len(alpha)
    summ = 0
    for i in range(N):
        summ += alpha[i, -1]*A[i+1, -1]
    if summ != 0:
        return math.log(summ)
    else:
        return summ

def backward(seq, A, B, longBet=False):
    N = len(Q) - 2
    T = len(seq)
    if not longBet:
        beta = np.zeros((N, T))
        for i in range(N):
            beta[i, -1] = bet(seq, i, -1, A, B)
        for j in range(T-2, -1, -1):
            for i in range(N):
                beta[i, j] = bet(seq, i, j, A, B, beta, True)
    else:
        beta = np.zeros((N, T+1))
        for i in range(N):
            beta[i, -1] = 1
        for i in range(N):
            beta[i, -2] = bet(seq, i, -1, A, B)
        for j in range(T-2, -1, -1):
            for i in range(N):
                beta[i, j] = bet(seq, i, j, A, B, beta, True)
    return beta

def bet(seq, n, t, A, B, beta=None, recursing=False):
    N = len(Qindices)
    if not recursing:
        return A[n+1, -1]
    else:
        summ = 0
        for i in range(N):
            a = A[n+1, i+1]
            b1 = B[i, Oindices[seq[t+1]]]
            bt1 = beta[i, t+1]
            summ += a*b1*bt1
        return summ

def terminateb(seq, beta, A, B):
    N = len(beta)
    summ = 0
    for i in range(N):
        summ += A[0, i+1]*B[i, Oindices[seq[0]]]*beta[i, 0]
    return math.log(summ)

# sets j and t inputs for B to standard
def getB(seq, j, t, B):
    if Q[j] == '<s>' or t == 0:
        if Q[j] == '<s>' and t == 0:
            return 1
        else:
            return 0
    elif Q[j] == '<f>' or t == len(seq)+1:
        if Q[j] == '<f>' and t == len(seq)+1:
            return 1
        else:
            return 0
    else:
        return B[j-1, Oindices[seq[t-1]]]

def getAlpha(seq, i, t, alpha):
    if Q[i] == '<s>' or t == 0:
        if Q[i] == '<s>' and t == 0:
            return 1
        else:
            return 0
    elif Q[i] == '<f>' or t == len(seq)+1:
        if Q[i] == '<f>' and t == len(seq)+1:
            return 1
        else:
            return 0
    else:
        return alpha[i-1, t-1]

def getBeta(seq, j, t, beta):
    if Q[j] == '<s>' or t == 0:
        if Q[j] == '<s>' and t == 0:
            return 1
        else:
            return 0
    elif Q[j] == '<f>' or t == len(seq)+1:
        if Q[j] == '<f>' and t == len(seq)+1:
            return 1
        else:
            return 0
    else:
        return beta[j-1, t-1]

def buildXi(seq, A, B, alpha, beta):
    N = len(Q) - 2
    T = len(seq)
    xi = np.zeros((T+1, N+2, N+2))
    prob = math.e**terminatef(seq, alpha, A)
    for t in range(0, T+1):
        for i in range(0, N+2):
            for j in range(0, N+2):
                ap = getAlpha(seq, i, t, alpha)
                bb = getB(seq, j, t+1, B)
                bt = getBeta(seq, j, t+1, beta)
                xi[t, i, j] = ap*A[i, j]*bb*bt/prob
    return xi

def xiList(seqs, A, B):
    xis = []
    for seq in seqs:
        alpha = forward(seq, A, B)
        beta = backward(seq, A, B)
        xis += [buildXi(seq, A, B, alpha, beta)]
    return xis

def buildGamma(seq, A, B, alpha, beta):
    N = len(Q) - 2
    T = len(seq)
    gamma = np.zeros((N, T))
    prob = math.e**terminatef(seq, alpha, A)
    for t in range(0, T):
        for j in range(0, N):
            ap = alpha[j, t]
            bt = beta[j, t]
            gamma[j, t] = ap*bt/prob
    return gamma

def gammaList(seqs, A, B):
    gammas = []
    for seq in seqs:
        alpha = forward(seq, A, B)
        beta = backward(seq, A, B)
        gammas += [buildGamma(seq, A, B, alpha, beta)]
    return gammas

def getGamma(seq, j, t, gamma):
    if Q[i] == '<s>' or t == 0:
        if Q[i] == '<s>' and t == 0:
            return 1 # OR RETURN 1/PROB!!!!!!!!!!!!!!!!
        else:
            return 0
    elif Q[i] == '<f>' or t == len(seq)+1:
        if Q[i] == '<f>' and t == len(seq)+1:
            return 1 # OR RETURN 1/PROB!!!!!!!!!!!!!!!!
        else:
            return 0
    else:
        return gamma[i-1, t-1]

def betterA(xis):
    N = len(Q) - 2
    A = np.zeros((N+2, N+2))
    for xi in xis:
        Aprime = np.sum(xi, 0)
        A += Aprime
    divisors = np.zeros(N+2)
    for xi in xis:
        divisors += np.sum(xi, (0, 2))
    for i in range(len(A)):
        for j in range(len(A[i])):
            if A[i, j] != 0:
                A[i, j] = A[i, j]/divisors[i]
    return A

def betterB(seqs, gammas):
    N = len(Q) - 2
    B = np.zeros((N, len(Oindices)))
    divisors = np.sum(gammas[0], 1)
    for g in range(1, len(seqs)):
        divisors += np.sum(gammas[g], 1)
    for g in range(len(seqs)):
        for i in range(len(gammas[g])):
            for o in range(len(gammas[g][i])):
                for t in range(len(seqs[g])):
                    if Oindices[seqs[g][t]] == o:
                        B[i, o] += gammas[g][i, t]
    for i in range(len(B)):
        for j in range(len(B[i])):
            B[i, j] /= divisors[i]
    return B


def bestModel(seqs, Ao, Bo):
    A = Ao
    B = Bo
    current = 0
    for seq in seqs:
        current += terminatef(seq, forward(seq, A, B), A)
    ago1 = current
    ago2 = current - 20
    ago3 = current - 40
    while ago3 + .1 < current:
        temp = current
        xl = xiList(seqs, A, B)
        gl = gammaList(seqs, A, B)
        A = betterA(xl)
        B = betterB(seqs, gl)
        current = 0
        for seq in seqs:
            current += terminatef(seq, forward(seq, A, B), A)
        ago3 = ago2
        ago2 = ago1
        ago1 = temp
    return (A, B)


xl = xiList(seqs, A0, B0)
gl = gammaList(seqs, A0, B0)

BESTpair = bestModel(seqs, A0, B0)
A = BESTpair[0]
B = BESTpair[1]
out = ''
model = open('hmm.learning', 'w')
for i in range(len(A)):
    for j in range(len(A[i])):
        out += 'T ' + Q[i] + ' ' + Q[j] + ' : '+ str(A[i][j]) + '\n'
for i in range(len(B)):
    for j in range(len(B[i])):
        out += 'E ' + Q[i+1] + ' ' + O[j] + ' : ' + str(B[i][j]) + '\n'
model.write(out)
model.close()


