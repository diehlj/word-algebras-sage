from word_algebras import *
import itertools




words = []
for w in itertools.product( range(1,3+1), repeat=3 ):
    s = ''.join(map(str,w))
    words.append( s )



R = PolynomialRing(QQ, [ 'a'+s for s in words ] )
NSA = NilpotentShuffleAlgebra(R,3,'y',nilpotency_order=6)

P = NSA.zero()
for i in range(len(words)):
    s = words[i]
    a = R.gens()[i]
    P += a * from_str( s, NSA )

det = from_str('11',NSA) * from_str('22',NSA) * from_str('33',NSA) \
    + from_str('21',NSA) * from_str('32',NSA) * from_str('13',NSA) \
    + from_str('12',NSA) * from_str('23',NSA) * from_str('31',NSA) \
    - from_str('31',NSA) * from_str('22',NSA) * from_str('13',NSA) \
    - from_str('32',NSA) * from_str('23',NSA) * from_str('11',NSA) \
    - from_str('21',NSA) * from_str('12',NSA) * from_str('33',NSA)

vol = from_str('123', NSA) - from_str('132', NSA)\
    - from_str('213', NSA) + from_str('231', NSA)\
    + from_str('312', NSA) - from_str('321', NSA)

assert 2**3 * det == vol**2 # Sanity check.

target = as_vector(det)

assert len(as_vector(det)) == 1 + 3 + 3**2 + 3**3 + 3**4 + 3**5 + 3**6

ansatz = as_vector( P**2 )

n = 0
for i in range(len(target)):
    if target[i] == 0 and ansatz[i] == 0:
        pass
    else:
        #print 'target=', target[i], '    ansatz=', ansatz[i]
        print target[i] - ansatz[i]
        n += 1

print 'total nr of equations=', n
