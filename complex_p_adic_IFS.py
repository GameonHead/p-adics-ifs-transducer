import complex_padic as cpa
import typing
from math import lcm, gcd

type State = tuple[cpa.Complex_pAdic, int]
type Result = tuple[State, tuple]


class Complex_pAdicFunction:
    def __init__(self, p, name='', d: None | cpa.Complex_pAdic = None, k=1, rotations=0):
        assert k > 0
        self.p = p
        self.d = d if d is not None else cpa.Complex_pAdic.zero(p)
        self.k = k
        self.name = name
        self.epsilon = rotations

    def as_func(self):
        def f(x: cpa.Complex_pAdic):
            return ((x << self.k) ^ self.epsilon) + self.d
        return f

    def __str__(self):
        name_part = self.name + ': ' if self.name != '' else ''
        sign = ('-' if (self.epsilon // 2) % 2 == 1 else '') + ('i * ' if self.epsilon % 2 == 1 else '')
        exponent_part = str(self.p) + (f'^{self.k} * ' if self.k > 1 else ' * ')
        added_part = ' + ' + str(self.d) if self.d != cpa.Complex_pAdic.zero(self.p) else ''
        return f"{name_part}{sign}{exponent_part}x{added_part}"

    def __repr__(self):
        return str(self)

class Transducer:
    def __init__(self, p, i: State, *functions: Complex_pAdicFunction):
        self.p = p
        self.functions = functions
        self.i = i
        self.states: dict[State, dict[Complex_pAdicFunction, Result]] = {i: {}}
        self.nodes: set[State] = {i}

    def simplify(self):
        f_tuples = [(f.d.to_rational(), f.epsilon, f.k, f.name) for f in self.functions]
        z, d = f_tuples[0][0]
        denom = lcm(*(rational[1] for rational, _, _, _ in f_tuples))


    def shift(self, x: cpa.Complex_pAdic, shift_count: int = 1) -> tuple[cpa.Complex_pAdic, tuple]:
        d = tuple(x[i] for i in range(shift_count))
        d_pAdic = cpa.Complex_pAdic.from_digit_sequence(self.p, tuple(reversed(d)))
        s = (x - d_pAdic) >> shift_count
        return s, d

    def apply_function(self, function: Complex_pAdicFunction, state: State) -> Result:
        output: tuple
        next_state: State
        next_state_part, output = self.shift((state[0] + function.d)^state[1], function.k)
        next_state = (next_state_part, (state[1] + function.epsilon) % 4)
        return next_state, output

    def find_all_from(self, state: State):
        assert state in self.states
        for f in self.functions:
            new_state, output = self.apply_function(f, state)
            if new_state not in self.nodes:
                self.states[new_state] = {}
                self.nodes.add(new_state)
            self.states[state][f] = (new_state, output)

    def create_transducer(self):
        visited = set()
        while visited != self.nodes:
            m = self.nodes.difference(visited)
            this_state = m.pop()
            self.find_all_from(this_state)
            visited.add(this_state)