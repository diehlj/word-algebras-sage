from __future__ import print_function

from sage.algebras.free_algebra import FreeAlgebra_generic, FreeAlgebraElement
from sage.modules.with_basis.indexed_element import IndexedFreeModuleElement
#from sage.structure.element import AlgebraElement
from sage.combinat.words.word import Word
from sage.all_cmdline import factorial, vector, Integer, flatten
import operator
import pyparsing as pp
import six
import functools

#def test_misc():
#    # TODO CONTINUE HERE
#    NFA = NilpotentFreeAlgebra(QQ,2,'y',nilpotency_order=3)
#    print(NFA.dimension())
#    #print(NFA._dense_free_module())
#    y0, y1 = NFA.gens()
#    el = y0**2 + y1
#    print(el)
#    print(vector(el))


#####################
# CONVENIENCE METHODS
#####################
def from_str(s, ALGEBRA, sympy_coefficients=False):

    element_parser = pp.Word(pp.nums)
    ALGEBRA_GENS = [ list(g)[0][0] for g in ALGEBRA.gens() ]
    def monomoial_from_str(s):
        return ALGEBRA.monomial( functools.reduce( operator.mul, [ ALGEBRA_GENS[int(i)-1] for i in s ] ) )


    if sympy_coefficients:
        from sympy.parsing.sympy_parser import parse_expr
        coeff_s = pp.QuotedString("[",endQuoteChar="]")
        coeff_s.setParseAction(lambda t: [parse_expr(t[0])])
        coeff = pp.Optional(coeff_s,1)
    else:        
        coeff_i=pp.Suppress("[")+pp.Word(pp.nums)+pp.Suppress("]")
        coeff_i.setParseAction(lambda t: [int(t[0])])
        coeff_f=pp.Suppress("[")+pp.Combine(pp.Optional(pp.Word(pp.nums))+
                                            "."+
                                            pp.Optional(pp.Word(pp.nums)))+pp.Suppress("]")
        coeff_f.setParseAction(lambda t: [float(t[0])])
        coeff=pp.Optional(coeff_i|coeff_f,1)
    if six.PY2:
        minus = pp.Literal("-")
    else:
        #In python 3, where str is unicode, it is easy to allow the minus sign character.
        #This means you can copy from a formula in a pdf
        minus = pp.Literal("-")|pp.Literal(chr(0x2212))
        minus.setParseAction(lambda t:["-"])
    firstTerm=pp.Optional(minus,"+")+coeff+pp.Optional(element_parser,"")
    otherTerm=(pp.Literal("+")|minus)+coeff+pp.Optional(element_parser,"")
    expn = pp.Group(firstTerm)+pp.ZeroOrMore(pp.Group(otherTerm))
    exp = expn.parseString(s,True)
    x = [(b if a=="+" else -b)*monomoial_from_str(c) for a,b,c in exp]
    out = functools.reduce(operator.add,x)
    return out  

def test_from_str():
    DIM = ORDER = 3
    CONCAT = NilpotentFreeAlgebra(QQ,DIM,'y',ORDER)
    SHUFFLE = NilpotentShuffleAlgebra(QQ,DIM,'y',ORDER)

    #print concatenation_word_from_str('112')
    print( from_str('112', CONCAT) )
    print( from_str('112 + 23', SHUFFLE) )
    #from_str('112')
    # TODO Put asserts.
#test_from_str()

def word_to_shuffle_monomial(SHUFF, w): # XXX name
    SHUFFLE_LETTERS = [ list(g)[0][0] for g in SHUFF.gens() ] # XXX
    if len(w) == 0:
        return SHUFF.one()
    return SHUFF.monomial( functools.reduce(operator.mul, map(lambda i: SHUFFLE_LETTERS[i], w)) ) # XXX starting at 0 vs starting at 1 ..

def shuffle_monomial_to_word(SHUFF, m):
    SHUFFLE_LETTERS = [ list(g)[0][0] for g in SHUFF.gens() ] # XXX
    return flatten( [ [SHUFFLE_LETTERS.index(ell)] * power for ell, power in m] )

def shuffle_to_concat(x): # TODO
    pass

def concat_to_shuffle(x): # TODO
    pass

def safe_add(d,key,value):
    if key in d:
        d[key] += value
    else:
        d[key] = value

def to_word( m ): # XXX naming
    ell = []
    for v, power in m:
        ell += [v] * power
    return Word( ell )

def to_monomial( w ): # XXX misnomer, this is still FreeMonoid
    if len(w) == 0:
        raise "Not implemented yet." # XXX
    else:
        return functools.reduce(operator.mul, w)

class NilpotentShuffleElement(FreeAlgebraElement): #IndexedFreeModuleElement, AlgebraElement):

    def __mul__(self,y):
        A = self.parent()
        z_elt = {}
        for mx, cx in self:
            for my, cy in y:
                if not A.nilpotency_order or len(my) + len(mx) <= A.nilpotency_order:
                    if len(mx) == 0:
                        safe_add(z_elt, my, cx * cy)
                    elif len(my) == 0:
                        safe_add(z_elt, mx, cx * cy)
                    else:
                        shuffle = to_word(mx).shuffle( to_word(my) )
                        for w in shuffle:
                            key = to_monomial(w)
                            safe_add(z_elt, key, cx * cy)
                            if not z_elt[key]:
                                del z_elt[key]
        return A._from_dict(z_elt)

    def succ(self,y):
        A = self.parent()
        z_elt = {}
        for mx, cx in self:
            for my, cy in y:
                left = to_word(mx)
                right = to_word(my)
                inner_shuffle = left.shuffle( right[:-1 ] )
                for w in inner_shuffle:
                    key = to_monomial(w) * right[-1 ]
                    safe_add(z_elt, key, cx * cy)
                    if not z_elt[key]:
                        del z_elt[key]
        return A._from_dict(z_elt)

class NilpotentShuffleAlgebra(FreeAlgebra_generic):

    Element = NilpotentShuffleElement # This is respected somewhere in a parent.

    def __init__(self, R, n, names, nilpotency_order=None):
        """
        The free shuffle algebra on `n` generators over a base ring.

        TESTS:
            sage: SHUFF = NilpotentShuffleAlgebra(QQ,2,'x',3)
            sage: x0,x1 = SHUFF.gens()
            sage: x0 * x1
            x0*x1 + x1*x0
            sage: x0 * x1 * x1
            2*x0*x1^2 + 2*x1*x0*x1 + 2*x1^2*x0
            sage: x0 * x1 * x1 * x1 # Nilpotency.
            0
        """
        
        super(NilpotentShuffleAlgebra, self).__init__(R, n, names)
        self.nilpotency_order = nilpotency_order


class NilpotentFreeAlgebraElement(FreeAlgebraElement): #IndexedFreeModuleElement, AlgebraElement):

    #def __mul__(self,y):
    #    A = self.parent()
    #    z_elt = {}
    #    for mx, cx in self:
    #        for my, cy in y:
    #            if not A.nilpotency_order or len(my) + len(mx) <= A.nilpotency_order:
    #                shuffle = to_word(mx).shuffle( to_word(my) )
    #                for w in shuffle:
    #                    key = to_monomial(w)
    #                    safe_add(z_elt, key, cx * cy)
    #                    if not z_elt[key]:
    #                        del z_elt[key]
    #    return A._from_dict(z_elt)

    def _mul_(self, y):
        # Copied from /Applications/SageMath/local/lib/python2.7/site-packages/sage/algebras/free_algebra_element.py
        A = self.parent()
        z_elt = {}
        for mx, cx in self:
            for my, cy in y:
                key = mx*my
                if not A.nilpotency_order or len(key) <= A.nilpotency_order:
                    if key in z_elt:
                        z_elt[key] += cx*cy
                    else:
                        z_elt[key] = cx*cy
                    if not z_elt[key]:
                        del z_elt[key]
        return A._from_dict(z_elt)

class NilpotentFreeAlgebra(FreeAlgebra_generic):

    Element = NilpotentFreeAlgebraElement # This is respected somewhere in a parent.

    def __init__(self, R, n, names, nilpotency_order=None):
        """
        The free algebra on `n` generators over a base ring.

        TESTS:
            sage: CONCAT = NilpotentFreeAlgebra(QQ,2,'x',3)
            sage: x0,x1 = CONCAT.gens()
            sage: x0 * x1
            x0*x1
            sage: x0 * x1 * x1
            x0*x1^2
            sage: x0 * x1 * x1 * x1 # Nilpotency.
            0
        """
        
        super(NilpotentFreeAlgebra, self).__init__(R, n, names)
        self.nilpotency_order = nilpotency_order

#from sage.rings import QQ
from sage.rings.rational_field import RationalField
from sage.rings.complex_field import ComplexField


def as_vector(el):
    """Convert an element for NilpotentShuffleAlgebra or NilpotentFreeAlgebra into a vector.

        TESTS:
        sage: NFA = NilpotentFreeAlgebra(QQ,2,'x',nilpotency_order=3)
        sage: x0, x1 = NFA.gens()
        sage: el = 1 + 2 * x0 + 3 * x1 + 4 * x0**2 + 5 * x0 * x1 + 6 * x1 * x0 + 7 * x1**2\
            + 8 * x0**3 + 9 * x0*x0*x1 + 10 * x0 * x1 * x0 + 11 * x0 * x1 * x1 \
            + 12 * x1*x0*x0 + 13 * x1*x0*x1 + 14 * x1 * x1 * x0 + 15 * x1 * x1 * x1
        sage: as_vector(el)
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
    """
    gens = el.parent().gens()
    DIM = len(gens)
    UPTO_LEVEL = el.parent().nilpotency_order
    assert UPTO_LEVEL
    SIZE = sum( [DIM**i for i in range(UPTO_LEVEL+1) ] )
    #v = vector( [0] * SIZE, RationalField() )
    #v = vector( [0] * SIZE, ComplexField() ) # XXX what?
    #v = vector( [0] * SIZE, el.parent() )
    #v = vector( el.parent(), SIZE )
    v = [0] * SIZE
    for m, c in list(el):
        flat = list(map( lambda x: gens.index(x), to_word( m ) ))
        i = 0
        current_level = 1
        for f in reversed(flat):
            i += DIM**(current_level-1) * (f+1)
            current_level += 1
        v[i] = c
        #try:
        #    v[i] = c.complex_embedding()
        #except:
        #    v[i] = c
    return v

def left_action_permutation_on_list(sigma, ell):
    """
    TESTS:
    sage: assert [77,99,8] == left_action_permutation_on_list( Permutation([3,2,1]), [8,99,77] )
    """
    #assert len(sigma) == len(ell)
    sigma_inverse = sigma.inverse()
    #print( 'ell=', ell )
    #print( 'sigma_inverse=', sigma_inverse, type(sigma_inverse) )
    return [ ell[ sigma_inverse(i) - 1 ] for i in range(1,len(ell)+1) ]

def left_action_permutation_on_words(tau, el):
    """
    TESTS:
    sage: perm = Permutation( [2,1] )
    sage: A = GroupAlgebra(SymmetricGroup(2), QQ)
    sage: tau = A(perm)
    sage: CONCAT = NilpotentFreeAlgebra(QQ,2,'x',2)
    sage: x0, x1 = CONCAT.gens()
    sage: left_action_permutation_on_words( tau, x0 * x1 )
    x1*x0
    """
    P = el.parent()
    def _left_action(tau):
        def f(x):
            w = to_word( x )
            #print( 'w=', w, type(w), type(w[0]) , w[0] * w[1] )
            for sigma, c in tau:
                yield ( to_monomial(left_action_permutation_on_list(sigma, w)), c )
        return f
    return apply_linear_map_gen(_left_action(tau), el )


def apply_linear_map_gen(fn, el): # XXX name
    """
    TESTS:
    sage: NFA = NilpotentFreeAlgebra(QQ,2,'y',nilpotency_order=3)
    sage: y0, y1 = NFA.gens()
    sage: def f(x):
    ...      yield (x*x, 3)
    sage: res = apply_linear_map_gen(f, 5*y0)
    sage: res
    15*y0^2
    sage: F = NFA.module_morphism( on_basis=lambda b: 3 * NFA.monomial(b*b), codomain=NFA )
    sage: F(5*y0)
    15*y0^2
    """
    di = {}
    A = el.parent()
    for x1, c1 in list(el):
        for x2, c2 in fn(x1):
            safe_add( di, x2, c1 * c2 )
    return A._from_dict( di )


def _r(m):
    def _r_helper(w):
        if len(w) == 0:
            pass
        elif len(w) == 1:
            yield (w,1)
        else:
            for w1, c1 in _r_helper( w[1:] ):
                yield ( w[:1] + w1, c1 )
                yield ( w1 + w[:1], -c1 )
    for m1, c1 in _r_helper(to_word(m)):
        yield (to_monomial(m1), c1)

def r(el):
    """Dynkin map, p.XXX in Reutenauer - Free Lie algebras.
    
    TESTS:
    sage: NFA = NilpotentFreeAlgebra(QQ,2,'y',nilpotency_order=3)
    sage: y0, y1 = NFA.gens()
    sage: assert y0 * y1 - y1 * y0 == r( y0 * y1 )
    sage: assert 0 == r( NFA.one() )
    """
    return apply_linear_map_gen(_r, el)

def _rho(m):
    def _rho_helper(w):
        if len(w) == 0:
            pass
        elif len(w) == 1:
            yield (w,1)
        else:
            a = w[:1]
            b = w[-1:]
            u = w[1:-1]
            for prev in _rho_helper( u + b ):
                yield (a + prev[0], prev[1])
            for prev in _rho_helper( a + u ):
                yield (b + prev[0], -prev[1])
    for m1, c1 in _rho_helper(to_word(m)):
        yield (to_monomial(m1), c1)

def rho(el):
    """The adjoint map to r."""
    return apply_linear_map_gen(_rho, el)

def id_otimes_r(el):
    def _id_otimes_r(t):
        left, right = t
        for w1, c1 in _r(right):
            yield ( (left, w1), c1 )
    return apply_linear_map_gen(_id_otimes_r, el)

def id_otimes_D(el):
    def _id_otimes_D(t):
        left, right = t
        assert len(left) == len(right)
        yield ( (left, right), len(right) )
    return apply_linear_map_gen(_id_otimes_D, el)

def lie_bracket(a,b):
    return a*b - b*a

def area_bracket(a,b):
    return a.succ(b) - b.succ(a)

def exp(el):
    ret = el.base_ring().one()
    ORDER = el.parent().nilpotency_order
    for i in range(1,ORDER+1):
        add = el**i
        add *= el.base_ring().one() / factorial(i)
        ret += add
    return ret

def inverse(el):
    ret = el.parent().zero()
    ORDER = el.parent().nilpotency_order
    for i in range(0,ORDER+1):
        add = (el-1)**i
        add *= Integer( (-1)**(i) )
        ret += add
    return ret

def log(el):
    """
    TESTS:
    sage: NFA = NilpotentFreeAlgebra(QQ,2,'x',3)
    sage: x0,x1 = NFA.gens()
    sage: assert x0 == log(exp(x0))
    """
    ret = 0
    ORDER = el.parent().nilpotency_order
    for i in range(1,ORDER+1):
        add = (el-1)**i
        add *= Integer( (-1)**(i+1) ) / i
        ret += add
    return ret

if __name__ == '__main__':
    test_misc()
    #examples()
