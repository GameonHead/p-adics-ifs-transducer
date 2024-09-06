import p_adic as pa
import typing as t

type State = tuple[pa.pAdic, int]
type Result = tuple[State, tuple]


class pAdicFunction:
    def __init__(self, p, name='', d=None, k=1, sign='+'):
        assert k > 0
        assert sign == '-' or sign == '+'
        self.p = p
        self.d = d if d is not None else pa.pAdic(p)
        self.k = k
        self.name = name
        self.sign = 1 if sign == '+' else -1

    def as_func(self):
        def f(x: pa.pAdic):
            return (x << self.k) + self.d if self.sign == 1 else -(x << self.k) + self.d

        return f

    def __repr__(self):
        return f"{self.name}{': ' if self.name != '' else ''}{'' if self.sign > 0 else '-'}{self.p}{('^' + str(self.k)) if self.k != 1 else ''}x + {self.d}"


class Transducer:
    def __init__(self, p, i: State, *functions: pAdicFunction):
        self.p = p
        self.functions = functions
        self.i = i
        self.states: dict[State, dict[pAdicFunction, Result]] = {i: {}}
        self.nodes: set[State] = {i}

    def shift(self, x: pa.pAdic, shift_count: int = 1) -> tuple[pa.pAdic, tuple]:
        d = tuple(x[i] for i in range(shift_count))
        pAdic_d = -pa.pAdic(self.p, whole_part=tuple(reversed(d)))
        s = (x + pAdic_d) >> shift_count
        return s, d

    def apply_function(self, function: pAdicFunction, state: State) -> Result:
        output: tuple
        next_state: State
        if state[1] == 1:
            next_state_part, output = self.shift(state[0] + function.d, function.k)
            next_state = (next_state_part, state[1] * function.sign)
        else:
            next_state_part, output = self.shift(-(state[0] + function.d), function.k)
            next_state = (next_state_part, state[1] * function.sign)
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
