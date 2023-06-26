import random
import math


def isPrime(q):
    if (q > 1):
        for i in range(2, int(math.sqrt(q)) + 1):
            if (q % i == 0):
                return False
        return True
    else:
        return False


# To generate random prime less than N
def randPrime(N):
    primes = []
    for q in range(2, N + 1):
        if (isPrime(q)):
            primes.append(q)
    return primes[random.randint(0, len(primes) - 1)]

#value assign to ? and for aplhabets
def asci(s):
    if s=="?":
        return 0
    else:
        return (ord(s) - 64)

# now if we have to reduce the false positive cases in the computation we have formed a findN function which gives the most appropraite prime N where eps is the upper bound on your error PR.
#we have use the formula n/log2(n) >= (2*m*log2(26))//eps
# setting the (n/log2(n))>=squaret(n) that means finding any k3 <=(squartroot(n) will satsifying our requirement which we can do using binary search in time log(a) a<=(m//eps)
def findN(eps, m):
     a = (9.5*(math.log2(m/eps)))**2
     return math.ceil(a)


# main function which
def dd(d,m,p,q):
    ans = 0
    d=26
    for i in range(0,m):
        if i==0:
            ans = (asci(p[i]))%q
        else:
            ans = (ans*d + asci(p[i]))%q
    return ans
# a type of function which calcualte the hash function for it using only O(log(Q)) bits of space and time complexity will be O(mlog(q))

# return 26**(m-1) whenever required for the computation of the hash value during text scanning
def power(d,m,q):
    h=1
    for i in range(m - 1):
        h = (h * d) % q
    return h
# this function efficiently reduce the sapce taken for compuatation by doing modular arithmetic operation every time

def modPatternMatch(q, p, x):
    n = len(x)
    m = len(p)
    hash_p = 0
    hash_x = 0
    d = 26
    ans = []
    h = power(26, m, q)       # takes O(mlog(q))time
    hash_p = dd(d,m,p,q)     # takes O(mlog(q))time
    hash_x = dd(d,m,x,q)     # takes O(mlog(q))time
    for i in range(n-m+1):   # time taken is O(nlogq) as n>=m for this loop
        if hash_x == hash_p:
            ans.append(i)
        if i < n - m:
            hash_x = (d * (hash_x - asci(x[i])*h) + asci(x[i + m])) % q  # use of this trick takes O(logq) time as we adding and subtarcting some values from hash and this also prevents overflow
            if hash_x<0:                                     # handle when hash value become negative
                hash_x = hash_x+q
    return(ans)
#ime complexity is O(mlogq+nlogq)  and space complexity O(k+logn+logq)  logn is the space taken by the saccning index and k is the taken by output list and log(q) is the sapce taken by we taking the modulus of the hash value


def randPatternMatch(eps, p, x):
    N = findN(eps, len(p))
    q = randPrime(N)
    return (modPatternMatch(q, p, x))
# as we have to ignore space and time taken by randprime algorithm
# now as q<=N and N is the minimum number satsfying error probality bound for the equality so hence q<(m//eps)
# now the major time is taken by modPatternmatch which is (m+n)logq <= (m+n)log(m//eps) hence time complexity is O((m+n)log(m//eps))
# the major space taken is by modPatternmatch which is k+lohn+logq < k+logn+log(m//eps) ,hence space complexity is O(k+logn+log(m//eps))


def modPatternMatchWildcard(q, p, x):
    n = len(x)
    m = len(p)
    hash_p= 0
    d = 26
    ans = []
    hash_x =0
    h= power(d,m,q)           # O(m) time taken
    a = 0
    r = 1
    for i in range(m):                    # finding the index where ? occurs in the pattern
         hash_p= (d * hash_p + asci(p[i])) % q
         if p[i]=="?":                    # we take hash_value of all the alphabets except "?"
             a = i
    if a==(m-1):                     # r gives  the how many time 26 is raised some power at index a
        r = 1                        # and further multiply it with asci(x[a])
    else:
        for i in range(0,m-a-1):
            r = (r*d)%q
    hash_x = dd(d,m,x,q)
    for i in range(n-m+1):          # takes O(nlogq) time
        hash_x2 = (hash_x - r*asci(x[a+i]))%q          # removing the value of the substring  at index (a) from hash_x and storing it in hash_x2,then comapring it with hash_p ,if matches append to the list
        if hash_x2==hash_p:
            ans.append(i)
        if i<(n-m):
            hash_x = (d * (hash_x - asci(x[i]) * h) + asci(x[i + m])) % q     # use of this trick takes O(logq) time as we adding and subtarcting some values from hash and this alsa prevents overflow
            if hash_x < 0:                       # handle when hash value becomes negative
                hash_x = hash_x + q
    return (ans)
# time complexity will be same as modpatternmatch as we are just doing one extra operation which takes O(logq) time in the second loop
# space compexity is same as modpatternmatch and some explaination as we have created one variable which takes O(logq) sapce
def randPatternMatchWildcard(eps,p,x):
	N = findN(eps,len(p))
	q = randPrime(N)
	return modPatternMatchWildcard(q,p,x)
# same explanation as randPatternmatch above
