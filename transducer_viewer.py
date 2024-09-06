from graphviz import Digraph
import p_adic_IFS as pIFS
import p_adic as pa


def node_name(node: pIFS.State) -> str:
    num, denom = node[0].to_rational()
    return f"({num}{f'/{denom}' if denom != 1 else ''}, {node[1]})"


def arc_label(f: pIFS.pAdicFunction, o: pIFS.Result) -> str:
    display = f.name if f.name != '' else f.__repr__()
    output = ','.join(map(str, o[1]))
    return f'{display}/{output}'


def view_transducer(transducer: pIFS.Transducer):
    transducer.create_transducer()
    graph = Digraph()
    for i in transducer.nodes:
        graph.node(node_name(i), label=node_name(i))
    for tail, transitions in transducer.states.items():
        for f, result in transitions.items():
            graph.edge(node_name(tail), node_name(result[0]), label=arc_label(f, result))
    print(graph.source)
    return graph


a = pIFS.pAdicFunction(5, 'A', pa.pAdic.to_p_adic(5,1,2))
b = pIFS.pAdicFunction(5, 'B', pa.pAdic.to_p_adic(5,1,3))
T = pIFS.Transducer(5,(pa.pAdic.zero(5),1), a, b)
view_transducer(T)