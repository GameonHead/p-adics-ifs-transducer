from typing import Callable, Any
from math import gcd

def prod(*args, default=1):
    if len(args) == 0:
        return default
    if len(args) == 1:
        return args[0]
    cumulative = args[0]
    for i in args[1:]:
        cumulative *= i
    return cumulative

class Polynomial[T]:
    def __init__(self, *coefficients: T, default = 0):
        self.coefficients = coefficients
        self.discard_zero()
        self.degree = len(self.coefficients) - 1
        self.default = default

    def discard_zero(self):
        j = len(self.coefficients) - 1
        while j > 0:
            if self.coefficients[j]:
                return
            self.coefficients = self.coefficients[:j]
            j -= 1

    def __repr__(self):
        return " + ".join(map(lambda x: repr(self.coefficients[x]) + (f" * x^{x}" if x > 1 else f" * x" if x == 1 else ""), range(len(self.coefficients))))

    def __str__(self):
        return repr(self)

    def __getitem__(self, item):
        if item < len(self.coefficients):
            return self.coefficients[item]
        return self.default

    def __add__(self, other):
        if not isinstance(other, Polynomial):
            return Polynomial(*([self.coefficients[0] + other] + list(self.coefficients[1:])))
        if self.degree < other.degree:
            return other + self
        return Polynomial(*[self[i] + other[i] for i in range(len(self.coefficients))], default=self.default)

    def __eq__(self, other):
        return all(map(lambda x, y: x == y, self.coefficients, other.coefficients))

    def __neg__(self):
        return Polynomial(*map(lambda x: -x, self.coefficients), default=self.default)

    def __sub__(self, other):
        return self + (-other)

    def __lshift__(self, other):
        if other == 0:
            return self
        return Polynomial(*([self.default] * other + list(self.coefficients)), default=self.default)

    def __rshift__(self, other):
        if other == 0:
            return self
        return Polynomial(*(self.coefficients[other:]), default=self.default)

    def __mul__(self, other):
        if not isinstance(other, Polynomial):
            return Polynomial(*map(lambda x: other * x, self.coefficients), default=self.default)
        return sum([(self * other[i]) << i for i in range(len(other.coefficients))], start=Polynomial(self.default, default=self.default))

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        assert not isinstance(other, type(self))
        return Polynomial(*map(lambda x: x / other, self.coefficients), default=self.default)

    def __hash__(self):
        return hash(self.coefficients)


def make_algebraic_class[T](labels: list[tuple[str, str]], rule: Callable[[int, tuple[T, ...]], tuple[T, ...]], length: int | None = None, default=None):
    if length == None:
        length = len(labels)
    class Algebraic:
        def __init__(self, *parts: T):
            self.parts = ([i for i in parts] + [default]*(length - len(parts)))
        def __add__(self, other):
            return Algebraic(*map(lambda x,y: x + y, self.parts, other.parts))

        def __sub__(self, other):
            return Algebraic(*map(lambda x,y: x - y, self.parts, other.parts))

        def __lshift__(self, other):
            return Algebraic(*map(lambda x,y: x << y, self.parts, other.parts))

        def __rshift__(self, other):
            return Algebraic(*map(lambda x,y: x >> y, self.parts, other.parts))

        def apply_rule(self, other: int):
            return Algebraic(*rule(other, self.parts))

        def __repr__(self):
            return "[" + ", ".join(map(repr, self.parts)) + "]"

        def __str__(self):
            return " + ".join([str(self.parts[i]) + labels[i][1] for i in range(len(self.parts)) if self.parts[i]])

        def multiply_simple(self, other: T):
            return Algebraic(*tuple(p * other for p in self.parts))

        def __mul__(self, other):
            if type(self) != type(other):
                raise TypeError
            a = self.apply_rule(0)
            a = a.multiply_simple(other.parts[0])
            for i in range(1, len(other.parts)):
                a += self.apply_rule(i).multiply_simple(other.parts[i])
            return a

        def __rmul__(self, other):
            return self * other

        def __getattr__(self, item):
            return self.parts[list(map(lambda x: x[0] == item, labels)).index(True)]

        def signed_permutation(self, *arrangement: int):
            return Algebraic(signed_permute(tuple(self.parts), arrangement))

        def __neg__(self):
            return self.signed_permutation(*(-(i + 1) for i in range(len(self.parts))))

    return Algebraic


def signed_permute(to_move: tuple[Any,...], arrangement: tuple[int,...]):
    assert len(arrangement) == len(to_move)
    return tuple(map(lambda i: to_move[i - 1] if i > 0 else -to_move[1 + i], range(len(arrangement))))

def nth_root_rule[T](x: T):
    def rule(k: int, i: tuple[T, ...]) -> tuple[T, ...]:
        if k == 0:
            return i
        return *(map(lambda u: x * u, i[-k:])), *i[:k+1]
    return rule

def nth_root_class[T](n: int, x: T):
    def rule(k: int, i: tuple[T, ...]) -> tuple[T, ...]:
        if k == 0:
            return i
        return *(map(lambda u: x * u, i[-k:])), *i[:k+1]
    names = [("e_0", "")] + [(f"e_{i}", "*" + str(x) + f"^{i}/{n}") for i in range(1,n)]

    class NthRootClass(make_algebraic_class(names, rule)):
        def __init__(self, *parts: T, default=None):
            super().__init__(*parts, default=default)

        def __truediv__(self, other):
            assert isinstance(other, type(self))

def pseudo_division(A: Polynomial[Any], B: Polynomial[Any], iteration_limit=10):
    r = A
    q = Polynomial(A.default, default=A.default)
    d = B[B.degree]
    e = A.degree - B.degree + 1
    i = 0
    while not r.degree < B.degree:
        if i > iteration_limit:
            raise RecursionError
        s = Polynomial(*([A.default] * (r.degree - B.degree) + [r[r.degree]]), default=A.default)
        q = (q * d) + s
        r = r * d - s * B
        e -= 1
        i += 1
    small_q = prod(*[d for i in range(e)])
    return q*small_q, r*small_q

def sub_resultant_gcd(A: Polynomial[Any], B: Polynomial[Any], g=1, h=1, One=Polynomial(1), gcd_function=gcd, iteration_limit=10):
    if B.degree > A.degree:
        return sub_resultant_gcd(B, A, g, h, One, gcd_function, iteration_limit)
    if B == Polynomial(B.default, default=B.default):
        return A
    a = gcd_function(*A.coefficients)
    b = gcd_function(*B.coefficients)
    print(a, b)
    d = gcd_function(a,b)
    A = A/a
    B = B/b
    Q, R = pseudo_division(A, B)
    print(f"{A} = {a} * ({Q})({B}) + {a*b} * ({R})")
    delta = A.degree - B.degree
    i = 0
    while i < iteration_limit:
        if R.degree == 0:
            if R.coefficients[0]:
                B = One
            gcd_in = gcd_function(*B.coefficients)
            return B * d / gcd_in
        A = B
        B = R/(g * prod(*(h for i in range(delta))))
        g = A[A.degree]
        h = prod(*(g for i in range(delta)), default=One[0])/prod(*(h for i in range(delta - 1)), default=One[0])
        Q, R = pseudo_division(A, B)
        print(f"{A} = ({Q})({B}) + {R}")
        delta = A.degree - B.degree
        i += 1
    raise RecursionError


def mix_rules(*rule_tuple: tuple[int, Callable[[int, tuple[Any, ...]], tuple[Any, ...]]]) -> Callable[[int, tuple[Any, ...]], tuple[Any, ...]]:
    states = tuple(i[0] for i in rule_tuple)
    rules = tuple(i[1] for i in rule_tuple)
    def rule(k: int, i: tuple[Any, ...]) -> tuple[Any, ...]:
        cur = i
        gap_size = 1
        length = len(i)
        block_count = length
        for l in range(len(rule_tuple)):
            k, u = divmod(k, states[l])
            block_count //= states[l]
            v = [cur[x: x + gap_size * states[l]: gap_size] for y in range(block_count) for x in range(gap_size * states[l] * y, gap_size * states[l] * y + gap_size)]
            w = [rules[l](u, a) for a in v]
            w = [list(zip(*w[i*gap_size:i*gap_size+gap_size])) for i in range(block_count)]
            cur = tuple(a for c in w for b in c for a in b)
            gap_size *= states[l]
        return cur
    return rule

def multiple_roots_class[T](*roots: tuple[int, T], default: None):
    indices = tuple(u for u, _ in roots)
    def multi_index(n) -> tuple[int, ...]:
        result = []
        for i in indices:
            n, k = divmod(n, i)
            result.append(k)
        return tuple(result)
    names = [("a" + "_".join(["0" for i in roots]), "")] + [("a" + "_".join(map(str, multi_index(i))), " * " + " * ".join(map(lambda x, y: f"{y[1]}^{x}/{y[0]}" if x != 0 else "\b\b\b", multi_index(i), roots))) for i in range(1, prod(*indices))]
    return make_algebraic_class(names, mix_rules(*((n, nth_root_rule(root)) for n, root in roots)), default=default)
