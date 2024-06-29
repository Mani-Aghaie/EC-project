import sympy as sp
import signal

#-----------------------------------------------FUNCTIONS------------------------------------------------------

def simplify_D(expr):
    """
    Simplifies the expression by eliminating D-operator in the denominator.
    
    Parameters:
    expr (sympy Eq): The equation to be simplified.
    
    Returns:
    sympy Eq: The simplified equation with no D-operator in the denominator.
    """
    # Continue simplifying while D-operator is in the denominator
    while is_D_in_denominator(expr):
        # Expand both sides of the equation after multiplying by D
        expr = sp.Eq(sp.expand(D * expr.lhs), sp.expand(D * expr.rhs))

    # Expand for those with no Ds
    return sp.Eq(sp.expand(expr.lhs), sp.expand(expr.rhs))

def is_D_in_denominator(expr):
    """
    Checks if the D-operator is present in the denominator of the equation.
    
    Parameters:
    expr (sympy Eq): The equation to be checked.
    
    Returns:
    bool: True if D-operator is found in the denominator, False otherwise.
    """
    # Get the denominator of the left-hand side and right-hand side of the equation
    _, denominator_lhs = sp.fraction(sp.simplify(sp.expand(expr.lhs)))
    #_, denominator_rhs = sp.fraction(sp.simplify(sp.expand(expr.rhs)))

    # Check if D-operator is in the denominator of either side
    return denominator_lhs.has(D) #or denominator_rhs.has(D)

def D_form_to_differential_form(D_form_equation, f_y_t, f_w_t):
    """
    Converts an equation in terms of the D operator to its corresponding differential equation form.
    
    Parameters:
    D_form_equation (sympy Equation): An equation which includes powers of D.
    f_y_t (sympy Function): The desired output signal as a function of time.
    f_w_t (sympy Function): The input waveform as a function of time.
    
    Returns:
    sympy Eq: Differential equation form of the same equation.
    """
    # Split the equation into left-hand side (lhs) and right-hand side (rhs)
    lhs = D_form_equation.lhs
    rhs = D_form_equation.rhs
    
    # Break the lhs and rhs into their additive components
    lhs_terms = lhs.as_ordered_terms()
    rhs_terms = rhs.as_ordered_terms()
    
    # Replacing Ds with derivatives on the left-hand side
    for i, term in enumerate(lhs_terms):
        D_counter = 0
        # While term contains D, count the power of D and simplify the term
        while term.has(D):
            term = term / D
            term = sp.expand(term)
            D_counter += 1
        # Replace D^n with nth derivative of y_t
        term = sp.expand((term / f_y_t) * sp.Derivative(f_y_t, (t, D_counter)))
        lhs_terms[i] = term
    
    # Combine the modified terms back into a single expression
    lhs = sum(lhs_terms)
    
    # Replacing Ds with derivatives/integrals on the right-hand side
    for i, term in enumerate(rhs_terms):
        D_counter = 0
        # While term contains D or D^-1, count the power of D and simplify the term
        while term.has(D) or term.has(D**(-1)):
            if term.has(D**(-1)):
                term = term * D
                term = sp.expand(term)
                D_counter -= 1
            else:
                term = term / D
                term = sp.expand(term)
                D_counter += 1
        # Replace D^n with nth derivative or integral based on D_counter
        if D_counter > 0:
            term = sp.expand((term / f_w_t) * sp.Derivative(f_w_t, (t, D_counter)))
        elif D_counter < 0:
            term = sp.expand(sp.integrate(term, (t, -D_counter)))
        rhs_terms[i] = term
    
    # Combine the modified terms back into a single expression
    rhs = sum(rhs_terms)
    
    # Return the resulting differential equation
    return sp.Eq(lhs, rhs)

def analyzer(A, Y, W, N, M, w_t, y_t, return_D_form_equation = False):
    """
    Develops the differential equation relating y(t) to w(t) for a given LTI circuit.
    
    Parameters:
    A (list of lists): n x n matrix whose elements are multiples of {D, D^-1, 1}
    Y (list): n x 1 vector of circuit variables
    W (list): n x 1 vector, where the Mth element is w(t) and others are 0
    N (int): Index of the desired signal y(t) in Y
    M (tuple): Indices of the input signal w(t) in W
    w_t (sympy function): Input waveform w(t)
    y_t (sympy function): output response y(t)
    
    Returns:
    sympy Eq: Differential equation relating y(t) to w(t)
    """

    # Convert A into sympy matrices
    A = sp.Matrix(A)

    # Calculating the minor matrix
    minor_matrix = A.minorMatrix(M, N)

    # Calculating the required determinants
    det_A = A.det()
    det_minor = minor_matrix.det()

    # Extract the equation and Turn it to sympy eq
    D_form_equation = simplify_D(sp.Eq(det_A * y_t, det_minor * w_t))

    if return_D_form_equation :
        return D_form_equation

    diff_equation = D_form_to_differential_form(D_form_equation, Y[N], W[M])

    return diff_equation 

def D_form_to_H_s(D_form_equation, f_y_t, f_w_t):
    """
    Converts an equation in terms of the D operator to its corresponding transfer function H(s).

    Parameters:
    D_form_equation (sympy Equation): An equation which includes powers of D.
    f_y_t (sympy Function): The desired output signal as a function of time.
    f_w_t (sympy Function): The input waveform as a function of time.

    Returns:
    sympy Function: The transfer function H(s) in terms of s.
    """
    # Isolate y(t) and w(t) terms
    y_terms = sp.collect(D_form_equation.lhs, f_y_t)
    w_terms = sp.collect(D_form_equation.rhs, f_w_t)
    
    # Extract coefficients of y(t) and w(t)
    coeff_y = y_terms.coeff(f_y_t)
    coeff_w = w_terms.coeff(f_w_t)
    
    # Solve for w(t)/y(t) to find the transfer function H(D)
    H_D = coeff_w / coeff_y

    # Substitute s for D to convert H(D) to H(s)
    s = sp.symbols('s')
    H_s = sp.Function('H')(s)
    H_s = H_D.replace(D, s)
    
    return H_s

def find_H_s(A, Y, W, N, M, w_t, y_t):
    """
    Finds the transfer function H(s) of a given LTI system.

    Parameters:
    A (list of lists): n x n matrix whose elements are multiples of {D, D^-1, 1}.
    Y (list): n x 1 vector of circuit variables.
    W (list): n x 1 vector, where the Mth element is w(t) and others are 0.
    N (int): Index of the desired signal y(t) in Y.
    M (tuple): Indices of the input signal w(t) in W.
    w_t (sympy Function): Input waveform w(t).
    y_t (sympy Function): Output response y(t).

    Returns:
    sympy Function: The transfer function H(s) in terms of s.
    """
    # Get the D-form equation using the analyzer function
    D_form = analyzer(A, Y, W, N, M, w_t, y_t, True)
    
    # Convert the D-form equation to H(s) using D_form_to_H_s
    return D_form_to_H_s(D_form, y_t, w_t)

def find_impulse_response_idea_2(A, Y, W, N, M, w_t, y_t):
    """
    Finds the impulse response h(t) of an LTI circuit given its transfer function.

    Parameters:
    A (list of lists): n x n matrix whose elements are multiples of {D, D^-1, 1}.
    Y (list): n x 1 vector of circuit variables.
    W (list): n x 1 vector, where the Mth element is w(t) and others are 0.
    N (int): Index of the desired signal y(t) in Y.
    M (tuple): Indices of the input signal w(t) in W.
    w_t (sympy Function): Input waveform w(t).
    y_t (sympy Function): Output response y(t).

    Returns:
    sympy Function: Impulse response h(t) as a function of time t.
    """
    # subsitute the delta function to w_t
    w_t = sp.DiracDelta(t)
    
    # Find the coresponding differential equation
    diff_eq = analyzer(A, Y, W, N, M, w_t, y_t)

    # Generate initial conditions, assuming they are all zero
    ics = {y_t.subs(t, 0): 0}
    
    # Determine the order of the differential equation by finding the highest derivative
    # Iterate over terms and find the highest derivative
    highest_order = 0
    for term in sp.preorder_traversal(diff_eq.lhs):
        if isinstance(term, sp.Derivative):
            highest_order = max(highest_order, term.derivative_count)
    
    # Add initial conditions for each derivative up to the order of the equation
    for i in range(1, highest_order + 1):
        ics[sp.Derivative(y_t, t, i).subs(t, 0)] = 0


    # solve 
    try:
        solution = sp.dsolve(diff_eq, y_t, inc = ics)
    except ValueError():
        print("Couldn't find the proper initial condition")
        solution = sp.dsolve(diff_eq, y_t)

    return sp.simplify(solution)

def find_impulse_response(A, Y, W, N, M, w_t, y_t):
    """
    Finds the impulse response h(t) of an LTI circuit given its transfer function.

    Parameters:
    A (list of lists): n x n matrix whose elements are multiples of {D, D^-1, 1}.
    Y (list): n x 1 vector of circuit variables.
    W (list): n x 1 vector, where the Mth element is w(t) and others are 0.
    N (int): Index of the desired signal y(t) in Y.
    M (tuple): Indices of the input signal w(t) in W.
    w_t (sympy Function): Input waveform w(t).
    y_t (sympy Function): Output response y(t).

    Returns:
    sympy Function: Impulse response h(t) as a function of time t.
    """
    
    # Find the coresponding differential equation
    diff_eq_homogenous = analyzer(A, Y, W, N, M, 0, y_t)
    diff_eq = analyzer(A, Y, W, N, M, w_t, y_t)

    # # Generate initial conditions, assuming they are all zero
    # ics = {y_t.subs(t, 0): 0}
    
    # Determine the order of the differential equation by finding the highest derivative
    # Iterate over terms and find the highest derivative
    n = 0
    for term in sp.preorder_traversal(diff_eq.lhs):
        if isinstance(term, sp.Derivative):
            n = max(n, term.derivative_count)
    
    m = 0
    for term in sp.preorder_traversal(diff_eq.rhs):
        if isinstance(term, sp.Derivative):
            m = max(m, term.derivative_count)

    # # Add initial conditions for each derivative up to the order of the equation
    # for i in range(1, n + 1):
    #     ics[sp.Derivative(y_t, t, i).subs(t, 0)] = 0
    
    y_homogenous = sp.dsolve(diff_eq_homogenous, y_t)

    # Create a list of symbolic constants C1, C2, ..., Cm+1 or C1, C2, ..., Cn
    if m - n > -1:
        h_t = y_homogenous.rhs * sp.Heaviside(t) 
        constants , constants_values = [], []
        for i in range(1, m + 2):
            constants.append(sp.symbols(f'C{i}'))
            constants_values.append(0)
            if i > n:
                h_t +=  h_t + constants[i-1] * sp.DiracDelta(t).diff((t, n-1))
    else:
        constants = [sp.symbols(f'C{i}') for i in range(1, n + 1)]
        constants_values = [0 for i in range(n)]
        h_t = y_homogenous.rhs * sp.Heaviside(t)

    
    # solve 
    

    lhs_terms = diff_eq.lhs.as_ordered_terms()
    rhs_terms = diff_eq.rhs.as_ordered_terms()

    eq_system = []
    rhs_terms_counter = 0


    for i in range(len(constants)):
        #rhs of teh new equations
        delta_coeff = 0
        
        # to be made in the loop
        lhs_equation_expression = 0
            
        derivative_count = 0
        # you can be sure they are ordered properly this is to check if one of them is zero
        try:
            if isinstance(rhs_terms[rhs_terms_counter].as_coeff_Mul()[1], sp.Derivative):
                derivative_count = rhs_terms[rhs_terms_counter].as_coeff_Mul()[1].derivative_count
        except IndexError:
            print("index error happened line 314")
            derivative_count = 0
        if derivative_count == i:
            delta_coeff = rhs_terms[rhs_terms_counter].as_coeff_Mul()[0]
            rhs_terms_counter += 1
            
        for j in range(n + 1):
            derivative_count = 0
            # you can be sure they are ordered properly this is to check if one of them is zero
            if isinstance(lhs_terms[j].as_coeff_Mul()[1], sp.Derivative):
                derivative_count = lhs_terms[j].as_coeff_Mul()[1].derivative_count
            
            if derivative_count != j:
                continue
        
            if m - n >= i - j and j - (i + 1)< 0:
                lhs_equation_expression += lhs_terms[j].as_coeff_Mul()[0] * constants[n + i -j]
            elif  j - (i + 1) >= 0:
                lhs_equation_expression += lhs_terms[j].as_coeff_Mul()[0] * sp.expand(y_homogenous.rhs.diff((t, j - (i + 1))).replace(t, 0))
        eq_system.append(sp.Eq(lhs_equation_expression, delta_coeff))
    
    constants_values = sp.solve(eq_system, constants)
    print(constants_values)
    print(sp.expand(sp.simplify(h_t)))
    for i in range(len(constants)):
        h_t = h_t.replace(constants[i], constants_values[constants[i]])


    return sp.expand(sp.simplify(h_t))

    
    

    


#------------------------------------------------TEST THE APP HERE----------------------------------------------------------


# Define the symbolic variables
t = sp.symbols('t')
D = sp.symbols('D')
s = sp.symbols('s')
h_t = sp.Function('h')(t)
w_t = sp.Function('w')(t)
x_t = sp.Function('x')(t)
y_t = sp.Function('y')(t)
z_t = sp.Function('z')(t)

# Define the system matrices for the given problem
A = [
    [D**(-1) , D, 3],
    [D , 2 *D, 3*D - 3],
    [D**(-1) , D, 2*D - 3]
]
Y = [y_t, z_t]
W = [0, w_t]
M = 1
N = 0


diff_eq = analyzer(A, Y, W, N, M, W[M], Y[N])
H_s = find_H_s(A, Y, W, N, M, W[M], Y[N])
h_t = find_impulse_response(A, Y, W, N, M, W[M], Y[N])

sp.pprint(diff_eq)
sp.pprint(H_s)
sp.pprint(h_t)
# Calculate the minimal differential equation
# minimal_diff_eq = minimal_differential_equation(A, Y, W, N, M, w_t)
# print("Minimal Differential Equation:")
# sp.pprint(minimal_diff_eq)
