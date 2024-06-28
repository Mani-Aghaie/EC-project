import sympy as sp

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
    _, denominator_rhs = sp.fraction(sp.simplify(sp.expand(expr.rhs)))

    # Check if D-operator is in the denominator of either side
    return denominator_lhs.has(D) or denominator_rhs.has(D)

def D_form_to_differential_form(D_form_equation):
    """
    Applies the D operator and turns all Ds to derivative    
     
    prameteters:
    D_form_equation (sympy Equation): An equation which includes powers of D

    Returns:
    sympy Eq: Differential equation form of the same equation
    """
    #splititng the function 
    lhs = D_form_equation.lhs
    rhs = D_form_equation.rhs
    lhs_terms = lhs.as_ordered_terms()
    rhs_terms = rhs.as_ordered_terms()
    
    # Replacing Ds with derivatives
    #left hand side
    i = 0
    for term in lhs_terms:
        D_counter = 0
        while(term.has(D)):
            term = term / D
            term = sp.expand(term)
            D_counter += 1
        term = sp.expand((term / y_t) * sp.Derivative(y_t, (t, D_counter)))
        lhs_terms[i] = term
        i += 1
    
    lhs = sum(lhs_terms)

    #right hand side
    i = 0
    for term in rhs_terms:
        D_counter = 0
        while(term.has(D)):
            term = term / D
            term = sp.expand(term)
            D_counter += 1
        term = sp.expand((term / w_t) * sp.Derivative(w_t, (t, D_counter)))
        rhs_terms[i] = term
        i += 1
    
    rhs = sum(rhs_terms)

    #return the differential equation
    return  sp.Eq(lhs, rhs)

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

    diff_equation = D_form_to_differential_form(D_form_equation)

    return diff_equation 



# Define the symbolic variables
t = sp.symbols('t')
D = sp.symbols('D')
w_t = sp.Function('w')(t)
y_t = sp.Function('y')(t)
z_t = sp.Function('z')(t)

# Define the system matrices for the given problem
A = [
    [D + 2 , 1 - D],
    [3 , D**(-1)]
]
Y = [y_t, z_t]
W = [0, w_t]
# W = [3 * w_t, D1(w_t, t)]
M = 1
N = 0


diff_eq = analyzer(A, Y, W, N, M, w_t, y_t)
sp.pprint(diff_eq)

# Calculate the minimal differential equation
# minimal_diff_eq = minimal_differential_equation(A, Y, W, N, M, w_t)
# print("Minimal Differential Equation:")
# sp.pprint(minimal_diff_eq)
