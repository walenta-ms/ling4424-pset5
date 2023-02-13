import sys
import numpy as np

trainingfile = open(sys.argv[1])

oannotations = {}
oemissions = {}
opairs = {}
otrans = {}
ostarts = {}
oends = {}
otranscount = {}
pairs = []
# potentially implement a function that increments everything 

# Helper functions
def dcadd(dic, key, amt=1):
    if key in dic:
        if type(dic[key]) == type(0) or type(dic[key]) == type(0.0):
            dic[key] += amt
    else:
        dic[key] = amt

# def printmat(matrix, axis, end=''):
#     print('  ', end='')
#     for e in axis:
#         print(e, end=' '*(11-len(e)))
#     print()
#     print(matrix)
#     print(end, end='')

# Processing the input
def initialize(tfile):
    global oannotations
    global oemissions
    global opairs
    global otrans
    global ostarts
    global oends
    global otranscount
    global pairs
    index = -1
    for line in tfile:
        if line.strip() != '':
            l = line.strip().split(';')
            pairs += [[]]
            index += 1
            for i in range(len(l)):
                pair = l[i]
                p = pair.split()
                p = tuple(p)
                pairs[index] += [p]
                dcadd(oannotations, p[0])
                dcadd(oemissions, p[1])
                dcadd(opairs, p)
                if i < len(l) - 1:
                    p2 = tuple(l[i+1].split())
                    trans = (p[0], p2[0])
                    if len(trans) == 2:
                        dcadd(otranscount, trans[0])
                    dcadd(otrans, trans)
                if i == 0:
                    t = ('<s>', p[0])
                    dcadd(otranscount, t[0])
                    dcadd(otrans, t)
                    dcadd(ostarts, p[0])
                elif i + 1 == len(l):
                    t = (p[0], '<f>')
                    dcadd(otranscount, t[0])
                    dcadd(otrans, t)
                    dcadd(oends, p[0])

initialize(trainingfile)

Q = list(oannotations)
Q[0:0] = ['<s>']
Q += ['<f>']
for state0 in Q:
    for state1 in Q:
        if (state0, state1) not in otrans:
            otrans[(state0, state1)] = 0
    if state0 not in otranscount:
        otranscount[state0] = 0
O = list(oemissions)

A = np.zeros((len(Q), len(Q)))
for i in range(len(Q)):
    for j in range(len(Q)):
        A[i, j] = 0 if otranscount[Q[i]] == 0 else otrans[(Q[i], Q[j])]/otranscount[Q[i]]

B = np.zeros((len(Q)-2, len(O)))
for i in range(1, len(Q)-1):
    for j in range(len(O)):
        p = (Q[i], O[j])
        B[i-1, j] = 0 if p not in opairs else opairs[p]/oannotations[Q[i]]

# printmat(A, Q, '\n')
# printmat(B, O)

out = ''
model = open('model.hmm', 'w')
for i in range(len(A)):
    for j in range(len(A[i])):
        out += 'T ' + Q[i] + ' ' + Q[j] + ' : '+ str(A[i][j]) + '\n'
for i in range(len(B)):
    for j in range(len(B[i])):
        out += 'E ' + Q[i+1] + ' ' + O[j] + ' : ' + str(B[i][j]) + '\n'
model.write(out)
model.close()
# print(out)
