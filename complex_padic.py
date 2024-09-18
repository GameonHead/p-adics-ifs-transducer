from numbers import Complex

import p_adic as pa
from math import lcm, gcd

class Complex_pAdic:
    def __init__(self, p: int, real: pa.pAdic, imaginary: pa.pAdic):
        assert p > 1
        assert p == real.p and p == imaginary.p
        self.p = p
        self.re = real
        self.im = imaginary

    def to_rational(self) -> tuple[complex, int]:
        ra, rb = self.re.to_rational()
        ia, ib = self.im.to_rational()
        g = lcm(rb, ib)
        return (ra * g//rb) + 1j * (ia * g//ib), g

    @staticmethod
    def zero(p):
        return Complex_pAdic(p, pa.pAdic.zero(p), pa.pAdic.zero(p))

    @staticmethod
    def to_p_adic(p, num: complex, denom: complex = 1):
        assert denom != 0
        if num == 0:
            return Complex_pAdic.zero(p)
        # Make real denominator
        a = num * denom.conjugate()
        b = denom * denom.conjugate()
        return Complex_pAdic(p, pa.pAdic.to_p_adic(p, int(a.real), int(b.real)), pa.pAdic.to_p_adic(p, int(a.imag), int(b.real)))

    @staticmethod
    def _extend_to_complex(a: pa.pAdic):
        return Complex_pAdic(a.p, a, pa.pAdic.zero(a.p))

    @staticmethod
    def from_digit_sequence(p, sequence):
        reals = tuple(int(a.real) for a in sequence)
        imags = tuple(int(a.imag) for a in sequence)
        return Complex_pAdic(p, pa.pAdic(p, whole_part=reals), pa.pAdic(p, whole_part=imags))

    def __add__(self, other):
        if type(other) == pa.pAdic:
            other = Complex_pAdic._extend_to_complex(other)
        return Complex_pAdic(self.p, self.re + other.re, self.im + other.im)

    def __radd__(self, other):
        if type(other) == pa.pAdic:
            other = Complex_pAdic._extend_to_complex(other)
        return Complex_pAdic(self.p, self.re + other.re, self.im + other.im)

    def __neg__(self):
        return Complex_pAdic(self.p, -self.re, -self.im)

    def condense(self):
        return Complex_pAdic(self.p, self.re.condense(), self.im.condense())

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __lshift__(self, other):
        assert type(other) is int
        assert other >= 0
        if other == 0:
            return self
        return Complex_pAdic(self.p, self.re << other, self.im << other)

    def __rshift__(self, other):
        assert type(other) is int
        assert other >= 0
        if other == 0:
            return self
        return Complex_pAdic(self.p, self.re >> other, self.im >> other)

    def __getitem__(self, index):
        return self.re[index] + self.im[index] * 1j

    def __eq__(self, other):
        return self.re == other.re and self.im == other.im

    def __repr__(self):
        return repr(self.re) + " + i" + repr(self.im)

    def __str__(self):
        return str(self.re) + " + i" + str(self.im)

    def __xor__(self, other):
        assert type(other) is int
        final = other % 4
        match final:
            case 0:
                return self
            case 1:
                return Complex_pAdic(self.p, -self.im, self.re)
            case 2:
                return -self
            case 3:
                return Complex_pAdic(self.p, self.im, -self.re)
            case _:
                raise Exception('Impossible state reached')

    def __hash__(self):
        return hash((self.re, self.im))