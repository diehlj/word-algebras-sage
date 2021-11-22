from word_algebras import *
import itertools


def shuffle_conjecture():
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
            print( target[i] - ansatz[i] )
            n += 1

    print( 'total nr of equations=', n )

import numpy as np
def volume_identities():
    #a, b, c = map(from_str('{}'.format(i),NSA) from 
    vol = vol_bracket

    NSA = NilpotentShuffleAlgebra(QQ,3,'y',nilpotency_order=3)
    a = from_str('1', NSA)
    b = from_str('2', NSA)
    c = from_str('3', NSA)

    vs = []
    for i,j,k in itertools.permutations( (a,b,c) ):
        vs.append( as_vector_hom( vol(i,j,k), 3 ) )

    #M = matrix( vs )
    M = np.array( vs )
    np.save('M-3', M)
    #assert 1 == M.rank()
    assert 1 == np.linalg.matrix_rank(M)

    #print( as_vector_hom(vol(a,b,c), 3 ) )
    #print( len(as_vector_hom(vol(a,b,c), 3 )) )

    NSA = NilpotentShuffleAlgebra(QQ,5,'y',nilpotency_order=5)
    a = from_str('1', NSA)
    b = from_str('2', NSA)
    c = from_str('3', NSA)
    d = from_str('4', NSA)
    e = from_str('5', NSA)
    letters = [a,b,c,d,e]
    S = Set( {0,1,2,3,4} )
    vs = []
    for A in Subsets(S, 3):
        i,j,k = map( lambda idx: letters[idx], sorted(list(A)) )
        ell,m = map( lambda idx: letters[idx],sorted(list(S-A)) )
        vs.append( as_vector_hom( vol(vol(i,j,k),ell,m), 5 ) )
        #vs.append( as_vector( vol(vol(i,j,k),ell,m) ) )
    # M = matrix( vs )
    M = np.array( vs )
    np.save('M-5', M)
    #assert binomial(5,3) == M.rank()
    assert binomial(5,3) == np.linalg.matrix_rank(M)

    # TOO SLOW:
    NSA = NilpotentShuffleAlgebra(QQ,7,'y',nilpotency_order=7)
    a = from_str('1', NSA)
    b = from_str('2', NSA)
    c = from_str('3', NSA)
    d = from_str('4', NSA)
    e = from_str('5', NSA)
    f = from_str('6', NSA)
    g = from_str('7', NSA)
    letters = [a,b,c,d,e,f,g]
    S = Set( {0,1,2,3,4,5,6} )
    vs = []
    for A in Subsets(S, 3):
        print(A)
        for B in Subsets(S-A,2):
            C = S-A-B
            i,j,k = map( lambda idx: letters[idx], sorted(list(A)) )
            ell,m = map( lambda idx: letters[idx],sorted(list(B)) )
            n,o   = map( lambda idx: letters[idx],sorted(list(C)) )
            vs.append( as_vector_hom( vol(vol(vol(i,j,k),ell,m),n,o), 7 ) )
        for B in Subsets(S-A,3):
            C = S-A-B
            i,j,k   = map( lambda idx: letters[idx], sorted(list(A)) )
            ell,m,n = map( lambda idx: letters[idx],sorted(list(B)) )
            o       = list(map( lambda idx: letters[idx],sorted(list(C)) ))[0]
            vs.append( as_vector_hom( vol(vol(i,j,k),vol(ell,m,n),o), 7 ) )
    # M = matrix( vs )
    M = np.array( vs )
    np.save('M-7', M)
    print( npp.linalg.matrix_rank(M) )

    # Not sure how to count the expected kernel here:
    #vs = []
    #for i,j,k,ell,m in itertools.permutations( (a,b,c,d,e) ):
    #    vs.append( as_vector( vol( vol(i,j,k), ell, m) ) )
    #    vs.append( as_vector( vol( ell, vol(i,j,k), m) ) )
    #    vs.append( as_vector( vol( ell, m, vol(i,j,k) ) ) )
    #M = matrix( vs )
    #print(M.rank())

if __name__ == '__main__':
    #shuffle_conjecture()
    volume_identities()

