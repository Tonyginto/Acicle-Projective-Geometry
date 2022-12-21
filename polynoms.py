import Uni


# Generating placements
def NextSet(a, n, m):
    j = m - 1
    while (j >= 0) & (a[j] == n):
        j=j-1
    if (j < 0): 
        return False
    if (a[j] >= n):
        j=j-1
    a[j]=a[j]+1
    if (j == m - 1):
        return True
    for k in range(j+1,m):
        a[k] = 0
    return True
    
# Function for checking the simplicity of the number p
def is_prime(p):
    if (p == 2):
        return True
    if (p % 2 == 0):
        return False
    for i in range(3, int(p ** (1 / 2)), 2):
        if (p % i == 0):
            return False
    return True

# Decomposition of the number n into Prime factors
def PrimMultipliers(n):
    rez = []
    ch = n
    i = 2
    while i <= n:
        if (ch%i==0) & is_prime(i):
            while (ch%i==0):
                ch = ch/i
            rez.append(i)
        i=i+1
    return rez   


# Test for irreducibility of a polynomial in the field Fp^l
def LidlNiderreiterF(f):
    n = len(f)-1
    if (n==1):
        return 1
    s = PrimMultipliers(n)
    k = len(s)
    x = [0,1]
    ft = x
    for t in range(1,n):
        ft = Uni.delzeroend(Uni.powerFpl(ft, Uni.p))
        for i in range(0,k):
            if ((t==n/s[i]) & (ft == x)) | (len(Uni.EuclidFpX(f,Uni.subFpX(ft,x))[2])!=1):
                return 0
    fn = Uni.powerFpl(ft, Uni.p)
    if (fn == x):
        return 1
    else:
        return 0


# Test for primitivity of a polynomial in the field Fp^l
def isPrimitiveF(f):
    Uni.f = f
    n = len(f)-1
    ordf = Uni.p**n-1
    s = PrimMultipliers(ordf)
    x = [0,1]
    for si in s:
        if Uni.powerFpl(x,int(ordf/si))==[1]:
            return 0
    return 1


# Generating all irreducible and primitive polynomials in the field Fp^l in increasing weight order
def genF():
    a = []
    for i in range(Uni.l):
        a.append(0)
    a.append(1)
    rez = [[],[]]
    Uni.f = a
    if LidlNiderreiterF(a)==1:
        rez[0].append(Uni.copyelem(a))
        if isPrimitiveF(a)==1:
            rez[1].append(Uni.copyelem(a))
    while(NextSet(a, Uni.p-1, Uni.l)):
        Uni.f = a
        if LidlNiderreiterF(a)==1:
            rez[0].append(Uni.copyelem(a))
            if isPrimitiveF(a)==1:
                rez[1].append(Uni.copyelem(a))
    return rez


# Test for irreducibility of a polynomial in the field Fp^ld
def LidlNiderreiterG(g):
    n = len(g)-1
    s = PrimMultipliers(n)
    k = len(s)
    x = [[0],[1]]
    gt = x
    for t in range(1,n):
        gt = Uni.clearNulls(Uni.powerFpld(gt, Uni.p**Uni.l))
        for i in range(0,k):
            if ((t==n/s[i]) & (gt == x)) | (len(Uni.EuclidFplY(g,Uni.subFplY(gt,x))[2])!=1):
                return 0
    gn = Uni.powerFpld(gt, Uni.p**Uni.l)
    if (gn == x):
        return 1
    else:
        return 0


# Test for primitivity of a polynomial in the field Fp^l
def isPrimitiveG(g):
    Uni.mod = g
    n = len(g)-1
    ordg = (Uni.p**Uni.l)**n-1
    s = PrimMultipliers(ordg)
    x = [[0],[1]]
    for si in s:
        if Uni.clearNulls(Uni.powerFpld(x,int(ordg/si)))==[[1]]:
            return 0
    return 1


# Decomposition of numbers by powers of p
def razlP(num):
    rez=[]
    for i in range(Uni.l):
        rez.append(0)
    s = Uni.p
    i = Uni.l-1
    while num !=0:
        if (s**i<=num):
            rez[i] = int(num/(s**i))
            num = num - s**i*rez[i]
        i=i-1
    return rez


# Generating all irreducible and primitive polynomials in the Fp^ld field in increasing weight order
def genG():
    rez = [[],[]]
    a = []
    for i in range(Uni.d):
        a.append(0)
    a.append(1)
    el = []
    for i in a:
        el.append(Uni.delzeroend(razlP(i)))
    Uni.mod = el
    if LidlNiderreiterG(el)==1:
        rez[0].append(Uni.copyelem(el))
        if isPrimitiveG(el)==1:
            rez[1].append(Uni.copyelem(el))
    while(NextSet(a, Uni.p**Uni.l-1, Uni.d)):
        el = []
        for i in a:
            el.append(Uni.delzeroend(razlP(i)))
        Uni.mod = el
        if LidlNiderreiterG(el)==1:
            rez[0].append(Uni.copyelem(el))
            if isPrimitiveG(el)==1:
                rez[1].append(Uni.copyelem(el))
    return rez


# Generating first primitive polynomial in the field Fp^l
def genFirstF():
    a = []
    for i in range(Uni.l):
        a.append(0)
    a.append(1)
    Uni.f = a
    if LidlNiderreiterF(a)==1:
        if isPrimitiveF(a)==1:
            return Uni.copyelem(a)
    while(NextSet(a, Uni.p-1, Uni.l)):
        Uni.f = a
        if LidlNiderreiterF(a)==1:
            if isPrimitiveF(a)==1:
                return Uni.copyelem(a)


# Generating first primitive polynomial in the field Fp^ld
def genFirstG():
    a = []
    for i in range(Uni.d):
        a.append(0)
    a.append(1)
    el = []
    for i in a:
        el.append(Uni.delzeroend(razlP(i)))
    Uni.mod = el
    if LidlNiderreiterG(el)==1:
        if isPrimitiveG(el)==1:
            return Uni.copyelem(el)
    while(NextSet(a, Uni.p**Uni.l-1, Uni.d)):
        el = []
        for i in a:
            el.append(Uni.delzeroend(razlP(i)))
        Uni.mod = el
        if LidlNiderreiterG(el)==1:
            if isPrimitiveG(el)==1:
                return Uni.copyelem(el)


# Output of the polynomial in the field Fp^l
def outF(el):
    st = ""
    if el[0]!=0:
        st = str(el[0])
    for i in range(1,len(el)):
        if el[i]!=0:
            if st!="":
                st = st + " + " 
            if el[i]!=1:
                st = st + str(el[i]) + "*"
            st = st + "z"
            if i!=1:
                st = st + "^" + str(i)
    return(st)


# Output of the polynomial in the field Fp^ld
def outG(el):
    st = outF(el[0])
    for i in range(1,len(el)):
        if el[i]!=[0]:
            st = st + " + " 
            if el[i]!=[1]:
                if len(outF(el[i]))>1:
                    st = st + "("+outF(el[i]) + ")*"
                else:
                    st = st + outF(el[i]) + "*"
            st = st + "x"
            if i!=1:
                st = st + "^" + str(i)
    return(st)
