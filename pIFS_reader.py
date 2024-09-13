import p_adic as pa
import p_adic_IFS as pIFS
import transducer_viewer as tv
from p_adic_IFS import Transducer

with open('functions.txt', 'r+t') as f:
    lines = f.readlines()

    # Line 1 defines the p
    p_definition_line = lines[0]
    if not p_definition_line.startswith("p:"):
        raise LookupError("Unable to determine p")
    p = int(p_definition_line.removeprefix("p:").strip())

    # All other lines determine the functions
    function_list: list[pIFS.pAdicFunction] = []
    for line in lines[2:]:
        line = line.removesuffix('\n')
        if len(line) == 0:
            continue
        function_name, _, expression = line.partition(':')
        function_name = function_name.strip()
        expression = expression.strip()
        coefficient, _, constant = expression.partition('x')
        sign, _, exponent = coefficient.partition(str(p))
        _, _, exponent = exponent.partition('^')
        exponent = exponent.replace('*', '')
        exponent = exponent.strip()
        k_i = int(exponent) if exponent != '' else 1
        sign = sign.strip()
        e_i = sign if sign != '' else '+'
        numerator, _, denominator = constant.partition('/')
        numerator = numerator.replace(' ','')
        a_i = int(numerator) if numerator != '' else 0
        denominator = denominator.strip()
        b_i = int(denominator) if denominator != '' else 1
        this_function = pIFS.pAdicFunction(p, function_name, pa.pAdic.to_p_adic(p, a_i, b_i), k_i, e_i)
        function_list.append(this_function)
    transducer: pIFS.Transducer = pIFS.Transducer(p, (pa.pAdic.zero(p), 1), *function_list)

    # Line 2 decides what to do:
    if lines[1].strip() == 'DFS':
        print(tv.make_dfa(transducer).to_graphviz())
    elif lines[1].strip() == 'NDFS':
        tv.make_ndfa(transducer)
    else:
        print(tv.view_transducer(transducer))