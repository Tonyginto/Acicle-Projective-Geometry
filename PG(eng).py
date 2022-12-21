import Uni
import math
import polynoms

global n,k,v,d,s

globalpar= [2, 1, 4]
Uni.p=globalpar[0]
Uni.l=globalpar[1]
Uni.d=globalpar[2]
Uni.n=Uni.p**Uni.l
npol = 0
if Uni.l == 1:
    npol = 1
Uni.f=polynoms.genF()[1][npol]
x = [0,1,0]
n=Uni.n
d=Uni.d+1
# kolichestvo blokov
v=int((Uni.powerZ(Uni.n, Uni.d+1)-1)/(Uni.n-1))
# kolichestco elementov v bloke
k=int((v-1)/Uni.n)


# funktsiya vozvedeniya v stepen' elementa polya
def xPow(i):
    global x
    if i==0:
        return [0]
    else:
        return Uni.powerFpl(x,i)

# funktsiya opredeleniya indeksa elementa v pole
def ind(ell):
    if ell==[0]:
        return 0
    else:
        i=1
        while not xPow(i)==ell:
            i=i+1
    return i

# polucheniye nomera elementa psi dlya proyektivnoy geometrii (dlya lyuboy geometrii, d natural'noye chislo)
def Psi(el):
    global n
    flag = 0
    ret = 0
    for i in range(0, Uni.d+1):
        if flag == 0 and i == Uni.d:
            ret = 0
        if flag == 0 and el[i] == [1]:
            flag = 1
            for j in range(0, Uni.d-i):
                ret+=Uni.powerZ(Uni.n,j)
        else:
            if el[i] != [0]:
                ret+=Uni.powerZ(Uni.n,Uni.d-i)*ind(el[i])
    return ret

# funktsiya razlozheniya chisla po stepenyam n, s koeffitsiyentami ot 1 do n vklyuchitel'no
def razlN(num):
    global n
    ret=[]
    co=-1
    while num != 0:
        pr = Uni.divZ(num,n)[1]
        if pr == 0:
            pr = n
        ret.append(pr)
        num = Uni.divZ(num-pr,n)[0]
        co+=1
    return ret, co


# funktsiya preobrazovaniya tselogo neotritsatel'nogo chisla v yelement psi
def PsiInverse(N):
    B=[]
    if N==0:
        for i in range(0, Uni.d):
            B.append([0])
            if i==Uni.d-1:
                B.append([1])
    else:
        smth=razlN(N)
        for i in range(0, Uni.d-smth[1]-1):
            B.append([0])
        B.append([1])
        for i in range(0,smth[1]+1):
            B.append(xPow(smth[0][smth[1]-i]-1))
    return B


# funktsiya summy vektorov
def vectorsum(v,w):
    s=[]
    for i in range(0, Uni.d+1):
        s.append((Uni.addFpl(v[i],w[i])))
    return s

# funktsiya proizvedeniya vektorov
def vectorproduct(v,w):
    s=[]
    for i in range(0, Uni.d+1):
        s.append((Uni.mulFpl(v[i],w[i])))
    return s

# funktsiya raznosti vektorov
def vectorsub(v,w):
    a=[]
    for i in range(0, Uni.d+1):
        a.append((Uni.subFpl(v[i],w[i])))
    return a

# funktsiya umnozheniya vectora na chislo
def mulbyscalar(v,w):
    a=[]
    for i in range(0, Uni.d+1):
        a.append((Uni.mulFpl(v[i],w)))
    return a

# postroyeniye atsiklicheskogo bloka
def AciclBlock():
    global s
    B=[]
    D=[]
    # Postroyeniye bloka
    razl=razlN(s)
    q=1
    for i in range(0, Uni.d):
        pr = [[0]] * (Uni.d + 1)
        if i<Uni.d-razl[1]-1:
            pr[Uni.d-i]=[1]
        else:
            pr[razl[1]+1]=xPow(razl[0][razl[1]-q+1]-1)
            pr[razl[1]+1-q]=[1]
            q+=1
        B.append(pr)
        D.append(Psi(pr))
        if i > 0:
            cnt = len(B)
            for ii in range(cnt - 1):
                for iii in range(1, Uni.n):
                    smthq = vectorsum(B[cnt-1], mulbyscalar(B[ii], xPow(iii)))
                    B.append(smthq)
                    D.append(Psi(smthq))
    return B, D

print('PG( d =',Uni.d,', n =',Uni.n,')')
print()
print('elementy blokov v forme mnogochlenov')
for i in range(0,v):
    s=i
    print(i, AciclBlock()[0])

print()
print('elementy blokov v forme chisel')
for i in range(0,v):
    s=i
    print(i, AciclBlock()[1])

