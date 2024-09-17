import logging
from graphviz import Digraph
import p_adic_IFS as pIFS
import p_adic as pa
import numpy as np


class MyGraph:
    def __init__(self):
        self.graph: dict[str, list[tuple[str, str]]] = {}

    def add_node(self, node: str):
        assert node not in self.graph.keys()
        self.graph[node] = []

    def add_edge(self, tail: str, head: str, label: str = ""):
        assert tail in self.graph.keys()
        assert head in self.graph.keys()
        self.graph[tail].append((label, head))

    def to_graphviz(self):
        g = Digraph()
        for n in self.graph:
            g.node(n, label=n)
        for n in self.graph:
            for l, m in self.graph[n]:
                g.edge(n, m, l)
        return g

    def remove_node(self, node):
        assert node in self.graph.keys()
        del self.graph[node]
        for u in self.graph:
            self.graph[u] = [(a, b) for a, b in self.graph[u] if b != node]

    def edge_count(self, tail, head):
        if tail in self.graph:
            return [a[1] == head for a in self.graph[tail]].count(True)
        return 0

    def adjacency_matrix(self):
        ordering = list(self.graph.keys())
        matrix = [[self.edge_count(t,h) for h in ordering] for t in ordering]
        return np.array(matrix)


def node_name(node: pIFS.State) -> str:
    num, denom = node[0].to_rational()
    return f"({num}{f'/{denom}' if denom != 1 else ''}, {node[1]})"


def arc_label(f: pIFS.pAdicFunction, o: pIFS.Result) -> str:
    display = f.name if f.name != '' else f.__repr__()
    output = ','.join(map(str, o[1]))
    return f'{display}/{output}'

def ndfa_arcs(g: Digraph, tail: str, head: str, o: pIFS.Result, counter: int, q: MyGraph | None = None) -> int:
    output = o[1]
    if len(output) == 1:
        g.edge(tail, head, label=str(output[0]))
        if q is not None:
            q.add_edge(tail, head, str(output[0]))
        return counter
    prev = tail
    for i in output[:-1]:
        logging.info(f"NDFA: creating node {i}")
        g.node(str(counter), label=str(counter))
        g.edge(prev, str(counter), label=str(i))
        if q is not None:
            q.add_node(str(counter))
            q.add_edge(prev, str(counter), str(i))
        prev = str(counter)
        counter += 1
    g.edge(prev, head, label=str(output[-1]))
    if q is not None:
        q.add_edge(prev, head, str(i))
    return counter


def view_transducer(transducer: pIFS.Transducer):
    transducer.create_transducer()
    graph = Digraph()
    for i in transducer.nodes:
        graph.node(node_name(i), label=node_name(i))
    for tail, transitions in transducer.states.items():
        for f, result in transitions.items():
            graph.edge(node_name(tail), node_name(result[0]), label=arc_label(f, result))
    return graph

def make_ndfa(transducer: pIFS.Transducer, suppress_output=False) -> MyGraph:
    counter = 0
    g = MyGraph()
    transducer.create_transducer()
    graph = Digraph()
    for i in transducer.nodes:
        logging.info(f"NDFA: creating node {i}")
        graph.node(node_name(i), label=node_name(i))
        g.add_node(node_name(i))
    for tail, transitions in transducer.states.items():
        for f, result in transitions.items():
            counter = ndfa_arcs(graph, node_name(tail), node_name(result[0]), result, counter, g)
    if not suppress_output:
        print(graph.source)
    return g

def make_dfa(transducer: pIFS.Transducer):
    inputs = list(map(str, range(transducer.p)))
    unexplored_states: list[set[str]] = []
    traversed_states: list[set[str]] = []
    unexplored_states.append({node_name(transducer.i)})
    r = make_ndfa(transducer, True)
    state_nodes: dict[str, dict[str, set[str]]] = {}
    for node in r.graph:
        state_nodes[node] = {i: set() for i, _ in r.graph[node]}
        for i in state_nodes[node]:
            state_nodes[node][i] = {j for L, j in r.graph[node] if L == i}

    def combine_states(*args):
        assert all(map(lambda x: x in state_nodes, args))
        state_transitions: dict[str, set[str]] = {}
        for arg in args:
            for u in state_nodes[arg]:
                if u not in state_transitions:
                    state_transitions[u] = set()
                state_transitions[u] = state_transitions[u].union(state_nodes[arg][u])
        return state_transitions

    dfa = MyGraph()
    dfa.add_node(str(tuple(sorted(unexplored_states[0]))))
    while len(unexplored_states) > 0:
        cur_state = unexplored_states[0]
        states_from_here = combine_states(*cur_state)
        for i in states_from_here:
            if states_from_here[i] not in unexplored_states and states_from_here[i] not in traversed_states:
                dfa.add_node(str(tuple(sorted(states_from_here[i]))))
                unexplored_states.append(states_from_here[i])
            dfa.add_edge(str(tuple(sorted(cur_state))), str(tuple(sorted(states_from_here[i]))), i)
        traversed_states.append(cur_state)
        unexplored_states = unexplored_states[1:]
    return dfa

def hausdorff_dimension(transducer: pIFS.Transducer):
    adjacency_matrix = make_dfa(transducer).adjacency_matrix()
    spectral_radius = max(map(np.abs, np.linalg.eig(adjacency_matrix)[0]))
    return np.log(spectral_radius)/np.log(transducer.p)
