#function HammingG
#input: a number r
#output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code
def hammingGeneratorMatrix(r):
    n = 2**r-1
    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    return G
#print(hammingGeneratorMatrix(3))



#Helper Function
def decimalToVector(i,l):
    vec = []
    while i != 0:
        vec += [i % 2]
        i = i // 2
    x = l - len(vec)
    vec += [0]*x
    return vec[::-1]

#non-examinable function
def vectorToDecimal(j):
    j = int("".join(map(str, j[0:len(j)])))
    l = 0
    r = len(str(j))+1
    for i in range(0, r):
        x = j%10
        j = j//10
        l = l + (x*(2**i))
    return l


#Functions for repetition codes
def repetitionEncoder(m,n):
    return m*n


def repetitionDecoder(v):
    z = []
    o = []
    for bit in v:
        if bit == 0:
            z += [0]
        elif bit == 1:
            o += [1]
    if len(z) > len(o):
        return [0]
    elif len(o) > len(z):
        return [1]
    else:
        return []


#Functions for hamming codes
def message(a):
    l = len(a)
    r = 2
    while 2**r - 2*r - 1 < l:
        r += 1
    k = 2**r - r - 1
    m = []
    m += decimalToVector(l,r) + a + [0]*(k-(r+l))
    return m


def hammingEncoder(m):
    length = len(m)
    r = 2
    k = 2**r - r - 1
    while length != k:
        if length < k:
            return []
        r += 1
        k = 2**r - r - 1
    n = 2**r - 1
    G = hammingGeneratorMatrix(r)
    v = []
    for i in range(0, n):
        temp = 0
        for j in range(0, k):
            temp += m[j]*G[j][i]
        if temp % 2 == 0:     
            v += [0]  
        else: 
            v += [1]
    return v


def hammingDecoder(v):
    l = len(v)
    power = 0
    positions = [1]
    while positions[power] < l:
        power += 1
        positions += [2**power]
    positions.pop()
    length = l - len(positions)
    r = 2
    k = 2**r - r - 1
    while length != k:
        if length < k:
            return []
        r += 1
        k = 2**r - r - 1
    HT = []
    vHT = []
    for i in range(1, l+1):
        HT += [decimalToVector(i,r)]
    for i in range(0, r):
        temp = 0
        for j in range(0, l):
            temp += v[j]*HT[j][i]
        if temp % 2 == 0:     
            vHT += [0]  
        else: 
            vHT += [1]
    i = vectorToDecimal(vHT)
    c = v
    if i == 0:
        return c
    else:
        i -= 1
        if c[i] == 0:
            c[i] = 1
        else:
            c[i] = 0
        return c


def messageFromCodeword(c):
    l = len(c)
    power = 0
    positions = [1]
    while positions[power] < l:
        power += 1
        positions += [2**power]
    positions.pop()
    length = l - len(positions)
    r = 2
    k = 2**r - r - 1
    while length != k:
        if length < k:
            return []
        r += 1
        k = 2**r - r - 1
    m = []
    for i in range(0,l):
        if i+1 not in positions:
            m += [c[i]]
    return m


def dataFromMessage(m):
    length = len(m)
    r = 2
    k = 2**r - r - 1
    while length != k:
        if length < k:
            return []
        r += 1
        k = 2**r - r - 1
    l = vectorToDecimal(m[0:r])
    if l > len(m[r:k]):
        return []
    a = m[r:r+l]
    return a

