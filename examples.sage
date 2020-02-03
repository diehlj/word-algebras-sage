from word_algebras import *

NSA = NilpotentShuffleAlgebra(QQ,4,'y',nilpotency_order=4)

el_1 = from_str('12', NSA)
el_2 = from_str('34', NSA)

assert el_1 * el_2 == from_str('1234 + 1324 + 1342 + 3124 + 3142 + 3412', NSA)
print('')
print('EXAMPLE WITH QQ-COEFFICIENTS')
print('the shuffle of ', el_1, ' and ', el_2, ' is =', el_1 * el_2)
print('as vector (of dimension 1 + 4 + 4^2 + 4^3 + 4^4)=', as_vector(el_1 * el_2 ))
assert 1 + 4 + 4**2 + 4**3 + 4**4 == len(as_vector(el_1 * el_2 ))



R = PolynomialRing(QQ,'a1,a2,a3')
a1, a2, a3 = R.gens()

NSA = NilpotentShuffleAlgebra(R,2,'y',nilpotency_order=3)

el_1 = a1 * from_str('12', NSA)
el_2 = a2 * from_str('2', NSA)
el_3 = a3 * from_str('22', NSA)

print('')
print('EXAMPLE WITH POLYNOMIAL COEFFICIENTS')
print('the element=', el_1 * el_2 + el_3)
print('as vector (of dimension 1 + 2 + 2^2 + 2^3)=', as_vector(el_1 * el_2 + el_3))
assert 1 + 2 + 2**2 + 2**3 == len(as_vector(el_1 * el_2 + el_3))
