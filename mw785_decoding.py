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
    # seq += ['<f>']
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

def viterbi(seq):
    N = len(Q) - 2
    T = len(seq)
    V = np.zeros((N, T+1))
    vb = np.zeros((N, T), dtype='object')
    # vbx = []
    for i in range(N):
        V[i, 0] = alph(i, 0)
    for j in range(1, T+1):
        # ls1 = []
        # ls2 = []
        for i in range(N):
            pr = alph(i, j, V, True, vb)
            # maxx = alph(i, j, V, True)
            # print(pr, Q[pr[1]+1])
            # V[i, j] = pr[0]
            V[i, j] = pr[0]
            # ls1 += [pr[0]]
            # ls2 += [pr[1]]
            vb[i, j-1] = str(Q[pr[1]+1])
            # print(vb, end='\n'+str(j)+' '+str(i)+'\n')
        # print(ls1, ls2)
        # vbx += [Q[1+ls1.index(max(ls1))]]
    return (V, vb)
    # return (V, vbx)



def alph(n, t, V=None, recursing=False, vb=None):
    N = len(Qindices)
    if not recursing:
        return A[0, n+1]*B[n, Oindices[seq[t]]]
    else:
        summ = []
        if t == len(seq):
            Bprob = 1
            state = -1
        else:
            Bprob = B[n, Oindices[seq[t]]]
            state = n+1
        # if t == 1:
        #     ls = []
        #     for i in range(len(V)):
        #         ls += [V[i, 0] * A[i+1, state] * Bprob]
        #         maxx = max(ls)
        #         prevQ = ls.index(maxx)
        # else:
        #     # vb[]
        #     ls = []
        #     for i in range(N):
        #         ls += [V[i, t-1] * A[i+1, state] * Bprob]
        #     maxx = max(ls)
        #     prevQ = ls.index(maxx)
        #     print(t, ls)
        for i in range(N):
            s = V[i, t-1] * A[i+1, state] * Bprob
            # print(V[i, t-1], ' '*(12-len(str(V[i, t-1]))), A[i+1, state], ' '*(12-len(str(A[i+1, state]))), Bprob)
            summ += [s]
        maxx = max(summ)
        prevQ = summ.index(maxx)
        # print(summ)
        return (maxx, prevQ)

def assemblevb(vbx, V):
    # return (vbx, V)
    vbf = []
    prob = 0
    # print()
    vbf = vbDive(vbx, V)
    # for j in range(0, len(vbx[0])):
    #     ls = []
    #     for i in range(len(V)):
    #         ls += [V[i, j]]
    #     # print(ls)
    #     big = ls.index(max(ls))
    #     vbf += [vbx[big, j]]
    ls = []
    for i in range(len(V)):
        ls += [V[i, -1]]
    # print(ls)
    big = ls.index(max(ls))
    if V[big, -1] != 0:
        prob = math.log(V[big, -1])
    else:
        prob = 0
    return (vbf, prob)

def vbDive(vb, V, ind=-1):
    if len(vb[0]) == 1:
        return [vb[ind, 0]]
    elif ind == -1:
        fstate = vb[ind, -1]
        newInd = Qindices[fstate]
        return vbDive(vb[:, 0:-1], V[:, 0:-1], newInd) + [fstate]
    else:
        fstate = vb[ind, -1]
        # slicer1 = vb[:, -1]
        # slicer2 = list(V[:, -1])
        # newInd = slicer2.index(max(slicer2))
        newInd = Qindices[fstate]
        return vbDive(vb[:, 0:-1], V[:, 0:-1], newInd) + [fstate]

pair = viterbi(seq)
V = pair[0]
vb = pair[1]
pair2 = assemblevb(vb, V)
vbf = pair2[0]
prob = pair2[1]

out = 'Best Sequence : '
for v in vbf:
    out += str(v)+ ' '
out += '\nBest LL : '
out += str(prob)
# print(out)
f = open('hmm.decoding', 'w')
f.write(out)
f.close()
