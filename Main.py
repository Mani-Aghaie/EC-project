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
    #_, denominator_rhs = sp.fraction(sp.simplify(sp.expand(expr.rhs)))

    # Check if D-operator is in the denominator of either side
    return denominator_lhs.has(D) #or denominator_rhs.has(D)

def D_form_to_differential_form(D_form_equation, f_y_t, f_w_t):
    """
    Converts an equation in terms of the D operator to its corresponding differential equation form.
    
    Parameters:
    D_form_equation (sympy Equation): An equation which includes powers of D.
    y_t (sympy Function): The desired output signal as a function of time.
    w_t (sympy Function): The input waveform as a function of time.
    
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



# Define the symbolic variables
t = sp.symbols('t')
D = sp.symbols('D')
w_t = sp.Function('w')(t)
x_t = sp.Function('x')(t)
y_t = sp.Function('y')(t)
z_t = sp.Function('z')(t)

# Define the system matrices for the given problem
A = [
    [D , 1 - D**(-1), 1 + D],
    [3 , 2 * D, 1],
    [D, D**(-1) + 1, D**(-1)]
]
Y = [y_t, z_t, x_t]
W = [0, w_t, 0]
# W = [3 * w_t, D1(w_t, t)]
M = 1
N = 0


diff_eq = analyzer(A, Y, W, N, M, W[M], Y[N])
D_eq = analyzer(A, Y, W, N, M, W[M], Y[N], True)
sp.pprint(D_eq)
sp.pprint(diff_eq)

# Calculate the minimal differential equation
# minimal_diff_eq = minimal_differential_equation(A, Y, W, N, M, w_t)
# print("Minimal Differential Equation:")
# sp.pprint(minimal_diff_eq)
