from math import lcm, gcd
from typing import Tuple
from itertools import batched


class pAdic:
    digits: tuple[tuple, tuple, tuple]
    precision_tolerance = 1e-6

    def __init__(self, p, repeating_digits: tuple = (0,), whole_part: tuple = tuple(), frac_part: tuple = tuple()):
        """
Represents a p-adic number, for a given prime number p.
        :param p: Must be a prime number
        :param repeating_digits: The repeating digits in the p-adic.
        :param whole_part: The non-repeating digits in the p-adic.
        :param frac_part: The fractional part of the p-adic.
        """
        self.p = p
        self.digits = (repeating_digits, whole_part, frac_part)

    @staticmethod
    def zero(p):
        return pAdic(p)

    def __repr__(self):
        return f"[{"_".join(map(str, self.digits[0]))}]{"_".join(map(str, self.digits[1]))}.{"_".join(map(str, self.digits[2]))}"

    def __str__(self):
        """
Returns a readable version of the p-adic number after simplification.
        :return: str representation of pAdic
        """
        return repr(self.condense())

    def __hash__(self):
        reduced = self.condense()
        return hash((self.p, reduced.digits[0], reduced.digits[1], reduced.digits[2]))

    def __eq__(self, other):
        return self.p == other.p and self.digits == other.digits

    def __abs__(self):
        reduced = self.condense()
        if reduced == pAdic.zero(self.p):
            return 0
        if len(reduced.digits[2]) > 0:
            return self.p ** (len(reduced.digits[2]))
        try:
            aperiodic_part = tuple(map(lambda l: l > 0, reversed(reduced.digits[1])))
            k = aperiodic_part.index(True)
            return self.p ** (-k)
        except ValueError:
            k = len(reduced.digits[1])
            periodic_part = tuple(map(lambda l: l > 0, reversed(reduced.digits[0])))
            k += periodic_part.index(True)
            return self.p ** (-k)

    @staticmethod
    def p_adic_valuation(num, denom, p):
        def highest_power_of_p(x):
            assert x % 1 == 0
            i = 0
            while x % p == 0:
                i += 1
                x = x / p
            return i

        return highest_power_of_p(num) - highest_power_of_p(denom)

    @staticmethod
    def to_p_adic(p, num: int, denom: int = 1):
        """
Converts a given rational number to a p-adic.
        :param p: Must be a prime number
        :param num: Must be an integer
        :param denom: Must be a non-zero integer
        :return: pAdic number
        """
        assert denom != 0
        if num == 0:
            return pAdic.zero(p)
        available_digits = list(range(p))
        n_0 = pAdic.p_adic_valuation(num, denom, p)

        # shift over if not a p-adic integer
        if n_0 < 0:
            denom *= p ** n_0

        # this allows for the repeating structure to be captured
        visited: dict[int, int] = {}
        place_value = min(0, n_0)
        digits: dict[int, int] = {}
        a = num
        highest_index = 0
        while a not in visited:
            for digit in available_digits:
                if (a - digit * denom) % p != 0:
                    continue
                digits[place_value] = digit
                if place_value >= 0:
                    visited[a] = place_value
                    highest_index += 1
                a = (a - digit * denom) // p
                break
            place_value += 1

        # separate it into the component parts
        fractional_part = tuple([digits[j] for j in range(-1, n_0 - 1, -1)])
        aperiodic_part = tuple([digits[j] for j in range(visited[a] - 1, -1, -1)])
        periodic_part = tuple([digits[j] for j in range(highest_index - 1, visited[a] - 1, -1)])
        return pAdic(p, periodic_part, aperiodic_part, fractional_part)

    @staticmethod
    def _loop_over(x, amount=1):
        assert len(x) > 0
        length = len(x)
        if amount > length:
            amount %= length
        return tuple(list(x[-amount:]) + list(x[:-amount]))

    @staticmethod
    def _spew(cycle, remainder, amount=1):
        """
        Adds more to remainder from cycle
        :param cycle: cyclic part (tuple)
        :param remainder: remaining part (tuple)
        :param amount: how much to add to remainder
        :return: cycle in new state, extended remainder
        """
        assert amount > 0
        loops, cycling_amount = divmod(amount, len(cycle))
        remainder = tuple(list(cycle[-cycling_amount:]) + list(cycle) * loops + list(remainder))
        cycle = pAdic._loop_over(cycle, cycling_amount)
        return cycle, remainder

    def to_rational(self) -> tuple[int, int]:
        """
        Returns a tuple representing a rational number equivalent to the p-adic number
        :return: Numerator, denominator.
        """
        p = self.p
        a = 0
        b = 1
        for i in self.digits[0]:
            a *= p
            b *= p
            a += i
        a = -a
        b -= 1
        for i in self.digits[1]:
            a *= p
            a += i * b
        b_0 = b
        for i in self.digits[2]:
            a *= p
            b *= p
            a += i * b_0
        g = gcd(a, b)
        a //= g
        b //= g
        return a, b

    def __add__(self, other):
        """
Adds together two p-adic numbers, for the same p.
        :param other: p-adic number
        :return: p-adic number
        """
        if isinstance(other, int):
            other = pAdic.to_p_adic(self.p, other)
        assert isinstance(other, pAdic)
        # Cannot add incompatible p-adic numbers
        if self.p != other.p:
            raise ValueError
        p = self.p

        lhs_periodic_part, lhs_aperiodic_part, lhs_fractional_part = self.digits
        if lhs_periodic_part == tuple():
            lhs_periodic_part = (0,)

        rhs_periodic_part, rhs_aperiodic_part, rhs_fractional_part = other.digits
        if rhs_periodic_part == tuple():
            rhs_periodic_part = (0,)

        # Need periodic parts to be of common size
        common_period = lcm(len(lhs_periodic_part), len(rhs_periodic_part))
        lhs_periodic_part *= common_period // len(lhs_periodic_part)
        rhs_periodic_part *= common_period // len(rhs_periodic_part)

        # Need aperiodic parts to be of common size
        delay_of_period = max(len(lhs_aperiodic_part), len(rhs_aperiodic_part))
        l_difference, r_difference = delay_of_period - len(lhs_aperiodic_part), delay_of_period - len(
            rhs_aperiodic_part)
        if l_difference > 0:
            lhs_periodic_part, lhs_aperiodic_part = pAdic._spew(lhs_periodic_part, lhs_aperiodic_part,
                                                                l_difference + common_period)
        if r_difference > 0:
            rhs_periodic_part, rhs_aperiodic_part = pAdic._spew(rhs_periodic_part, rhs_aperiodic_part,
                                                                r_difference + common_period)
        lhs_aperiodic_part = tuple(list(lhs_periodic_part) + list(lhs_aperiodic_part))
        rhs_aperiodic_part = tuple(list(rhs_periodic_part) + list(rhs_aperiodic_part))

        # Need fractional parts to be of common size
        if len(lhs_fractional_part) < len(rhs_fractional_part):
            difference = len(rhs_fractional_part) - len(lhs_fractional_part)
            lhs_fractional_part = tuple(list(lhs_fractional_part) + [0] * difference)
        elif len(rhs_fractional_part) < len(lhs_fractional_part):
            difference = len(lhs_fractional_part) - len(rhs_fractional_part)
            rhs_fractional_part = tuple(list(rhs_fractional_part) + [0] * difference)

        # Add with carrying
        frac_list = []
        carry = 0
        for a, b in zip(reversed(lhs_fractional_part), reversed(rhs_fractional_part)):
            carry, c = divmod(a + b + carry, p)
            frac_list.append(c)
        aper_list = []
        for a, b in zip(reversed(lhs_aperiodic_part), reversed(rhs_aperiodic_part)):
            carry, c = divmod(a + b + carry, p)
            aper_list.append(c)
        periodic_list = []
        for a, b in zip(reversed(lhs_periodic_part), reversed(rhs_periodic_part)):
            carry, c = divmod(a + b + carry, p)
            periodic_list.append(c)
        return pAdic(p, tuple(reversed(periodic_list)), tuple(reversed(aper_list)),
                     tuple(reversed(frac_list))).condense()

    def condense(self):
        """
Simplifies a p-adic number to remove trailing zeroes, repetition in the periodic part, and parts of the periodic part in the aperiodic part.
        :return: equivalent p-adic number
        """
        p, a, f = self.digits
        # Collapse anything into periodic part
        while len(a) > 0 and p[0] == a[0]:
            p = pAdic._loop_over(p, -1)
            a = a[1:]

        # If the periodic part can be condensed, do so.
        p_len = len(p)
        for i in range(1, p_len):
            if p_len % i != 0:
                continue
            if all(p[:i] == batch for batch in batched(p, i)):
                p = p[:i]
                break

        # Remove trailing 0s in the fraction
        while len(f) > 0 and f[-1] == 0:
            f = f[:-1]

        return pAdic(self.p, p, a, f)

    def __neg__(self):
        p = self.p
        r, a, f = self.digits
        r = tuple(map(lambda i: p - 1 - i, r))
        a = tuple(map(lambda i: p - 1 - i, a))
        f = tuple(map(lambda i: p - 1 - i, f))
        result: pAdic
        if len(f) > 0:
            result = pAdic(self.p, r, a, f) + pAdic(self.p, frac_part=tuple([0] * (len(f) - 2) + [1]))
        else:
            one = pAdic(p, whole_part=(1,))
            result = pAdic(self.p, r, a, f) + one
        return result

    def __sub__(self, other):
        return self + (- other)

    def __mul__(self, other):
        """
        Multiplies if rational. If this is extended into the irrationals this will need a total rewrite.
        :param other:
        :return: self * other
        """
        assert isinstance(other, pAdic)
        assert other.p == self.p
        a, b = self.to_rational()
        c, d = other.to_rational()
        return pAdic.to_p_adic(self.p, a * c, b * d)

    def __lshift__(self, other):
        """
Equivalent to multiplying by p^other
        :param other: non-negative integer, represents how much to shift digits
        :return: p-adic number = self * p^other
        """
        assert type(other) is int
        assert other >= 0
        if other == 0:
            return self
        p, a, f = self.digits
        f = tuple(list(f) + [0] * other)
        a = tuple(list(a) + list(f[:other]))
        f = f[other:]
        return pAdic(self.p, p, a, f).condense()

    def __rshift__(self, other):
        """
        Equivalent to dividing by p^other
        :param other: non-negative integer, represents how much to shift digits
        :return: p-adic number = self * p^-other
        """
        assert type(other) is int
        assert other >= 0
        if other == 0:
            return self
        p, a, f = self.digits
        p, a = pAdic._spew(p, a, other)
        f = tuple(list(a[-other:]) + list(f))
        a = a[:-other]
        return pAdic(self.p, p, a, f).condense()

    def __getitem__(self, index):
        if isinstance(index, int):
            if index < -len(self.digits[2]):
                raise IndexError('Out of bounds')
            if index < 0:
                return self.digits[2][-1 - index]
            if index < len(self.digits[1]):
                return self.digits[1][len(self.digits[1]) - 1 - index]
            else:
                index -= len(self.digits[1])
                index %= len(self.digits[0])
                return self.digits[0][len(self.digits[0]) - 1 - index]
        else:
            raise NotImplementedError
