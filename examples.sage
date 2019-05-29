from word_algebras import *

NSA = NilpotentShuffleAlgebra(QQ,4,'y',nilpotency_order=4)

el_1 = from_str('12', NSA)
el_2 = from_str('34', NSA)

assert el_1 * el_2 == from_str('1234 + 1324 + 1342 + 3124 + 3142 + 3412', NSA)

# Collecting coefficients.
R = PolynomialRing(QQ,'a1,a2,a3')
a1, a2, a3 = R.gens()

NSA = NilpotentShuffleAlgebra(R,2,'y',nilpotency_order=4)

el_1 = a1 * from_str('12', NSA)
el_2 = a2 * from_str('2', NSA)
el_3 = a3 * from_str('22', NSA)

print 'the element=', el_1 * el_2 + el_3
print 'as vector=', as_vector(el_1 * el_2 + el_3)
