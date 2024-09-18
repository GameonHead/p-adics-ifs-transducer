import complex_padic as cpa
import complex_p_adic_IFS as cpIFS
import transducer_viewer as tv
import re

cplx_finder_regex = '''^((?:\+|-)\d+(?:\/\d+)?)((?:(?:\+|-)(?:i|j)\*?\d+(?:\/\d+)?)|(?:(?:\+|-)(?:\d+)?\*?(?:i|j)(?:\/\d+)?)|(?:(?:\+|-)(?:\d+(?:\/\d+)?)?\*?(?:i|j)))?$|^((?:(?:\+|-)(?:i|j)\*?\d+(?:\/\d+)?)|(?:(?:\+|-)(?:\d+)?\*?(?:i|j)(?:\/\d+)?)|(?:(?:\+|-)(?:\d+(?:\/\d+)?)?\*?(?:i|j)))((?:\+|-)\d+(?:\/\d+)?)?$'''

with open('functions.txt', 'r+t') as f:
    lines = f.readlines()

    # Line 1 defines the p
    p_definition_line = lines[0]
    if not p_definition_line.startswith("p:"):
        raise LookupError("Unable to determine p")
    p = int(p_definition_line.removeprefix("p:").strip())

    # All other lines determine the functions
    function_list: list[cpIFS.Complex_pAdicFunction] = []
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
        sign = sign.replace('*', '')
        sign = sign.strip()
        e_i = 0
        if '-' in sign:
            e_i += 2
        if 'i' in sign or 'j' in sign:
            e_i += 1
        print(e_i)
        constant = constant.replace(' ', '')
        a_i = 0j
        b_i = 1+0j
        if constant != '':
            print(constant)
            m = re.match(cplx_finder_regex, constant, re.M)
            capture_groups = m.groups()
            real_str = m.group(1) if m.group(1) is not None else m.group(4) if m.group(4) is not None else '0'
            imag_str = m.group(2) if m.group(2) is not None else m.group(3) if m.group(3) is not None else '0i'
            imag_str = imag_str.replace('i', '').replace('j', '').replace('*', '')
            re_num_str, _, re_denom_str = real_str.partition('/')
            if re_denom_str == '':
                re_denom_str = '1'
            if re_num_str == '':
                re_num_str = '0'
            im_num_str, _, im_denom_str = imag_str.partition('/')
            if im_denom_str == '':
                im_denom_str = '1'
            if im_num_str == '':
                im_num_str = '0'
            a_i = int(re_num_str) * int(im_denom_str) + int(im_num_str) * int(re_denom_str) * 1j
            b_i = int(re_denom_str) * int(im_denom_str)
        this_function = cpIFS.Complex_pAdicFunction(p, function_name, cpa.Complex_pAdic.to_p_adic(p, a_i, b_i), k_i, e_i)
        function_list.append(this_function)
    print(function_list)
    transducer: cpIFS.Transducer = cpIFS.Transducer(p, (cpa.Complex_pAdic.zero(p), 1), *function_list)

    # Line 2 decides what to do:
    if lines[1].strip().upper() == 'DFA':
        print(tv.make_dfa(transducer).to_graphviz())
    elif lines[1].strip().upper() == 'NDFA':
        tv.make_ndfa(transducer)
    elif lines[1].strip().upper() == 'A':
        print(tv.make_dfa(transducer).adjacency_matrix())
    elif lines[1].strip().upper() == 'DIMENSION':
        print(tv.hausdorff_dimension(transducer))
    elif lines[1].strip().upper() == 'SIMPLIFY':
        simple_t = transducer.simplify()
        print(simple_t.functions)
        print(f'Hausdorff Dimension: {tv.hausdorff_dimension(simple_t)}')
        print(f'DFA:')
        print(tv.make_dfa(simple_t).to_graphviz())
    else:
        print(tv.view_transducer(transducer))