import sympy as sp
import signal
from sympy import I 

#----------------------------------------------- FUNCTIONS ------------------------------------------------------

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
    _ , denominator_lhs = sp.fraction(sp.simplify(sp.expand(expr.lhs)))
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
    det_minor = ((-1)**(N + M)) * minor_matrix.det()

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

def find_impulse_response_idea_1(A, Y, W, N, M, w_t, y_t):
    """
    Attempts to find the impulse response h(t) of an LTI circuit using the Laplace transform.

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

    Note:
    This attempt was abandoned due to the slow computation of the inverse Laplace transform.
    """

    # Find the transfer function H(s) of the system
    H_s = find_H_s(A, Y, W, N, M, w_t, y_t)
    
    # Compute the inverse Laplace transform to find the impulse response h(t)
    h_t = sp.inverse_laplace_transform(H_s, s, t)
    
    return h_t

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

    Note:
    This attempt was abandoned because SymPy could not be forced to calculate the constants correctly.
    """

    # Substitute the Dirac delta function for the input w(t)
    w_t = sp.DiracDelta(t)
    
    # Find the corresponding differential equation for the system
    diff_eq = analyzer(A, Y, W, N, M, w_t, y_t)

    # Generate initial conditions, assuming they are all zero
    ics = {y_t.subs(t, 0): 0}
    
    # Determine the order of the differential equation by finding the highest derivative
    highest_order = 0
    for term in sp.preorder_traversal(diff_eq.lhs):
        if isinstance(term, sp.Derivative):
            highest_order = max(highest_order, term.derivative_count)
    
    # Add initial conditions for each derivative up to the order of the equation
    for i in range(1, highest_order + 1):
        ics[sp.Derivative(y_t, t, i).subs(t, 0)] = 0

    # Solve the differential equation with initial conditions
    try:
        solution = sp.dsolve(diff_eq, y_t, ics=ics)
    except ValueError:
        print("Couldn't find the proper initial conditions")
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
    
    # Find the corresponding differential equations for the system
    diff_eq_homogeneous = analyzer(A, Y, W, N, M, 0, y_t)
    diff_eq = analyzer(A, Y, W, N, M, w_t, y_t)

    # Determine the order of the differential equation by finding the highest derivative on both sides
    n = 0
    for term in sp.preorder_traversal(diff_eq.lhs):
        if isinstance(term, sp.Derivative):
            n = max(n, term.derivative_count)
    
    m = 0
    for term in sp.preorder_traversal(diff_eq.rhs):
        if isinstance(term, sp.Derivative):
            m = max(m, term.derivative_count)

    # Solve the homogeneous differential equation
    y_homogeneous = sp.dsolve(diff_eq_homogeneous, y_t)

    # Create a list of symbolic constants C1, C2, ..., Cm+1 or C1, C2, ..., Cn
    if m - n > -1:
        h_t = y_homogeneous.rhs * sp.Heaviside(t)
        constants, constants_values = [], []
        for i in range(1, m + 2):
            constants.append(sp.symbols(f'C{i}'))
            constants_values.append(0)
            if i > n:
                h_t += constants[i-1] * sp.DiracDelta(t).diff((t, n-1))
    else:
        constants = [sp.symbols(f'C{i}') for i in range(1, n + 1)]
        constants_values = [0 for i in range(n)]
        h_t = y_homogeneous.rhs * sp.Heaviside(t)
    
    # Extract the left-hand side and right-hand side terms of the differential equation
    lhs_terms = diff_eq.lhs.as_ordered_terms()
    rhs_terms = diff_eq.rhs.as_ordered_terms()

    # Initialize an empty list to store the system of equations
    eq_system = []
    rhs_terms_counter = 0

    # Construct the system of equations to solve for the constants
    for i in range(len(constants)):
        delta_coeff = 0
        lhs_equation_expression = 0
        derivative_count = 0
        
        # Determine the derivative count for the right-hand side terms
        try:
            if isinstance(rhs_terms[rhs_terms_counter].as_coeff_Mul()[1], sp.Derivative):
                derivative_count = rhs_terms[rhs_terms_counter].as_coeff_Mul()[1].derivative_count
        except IndexError:
            print("Index error occurred")
            derivative_count = 0
        if derivative_count == i:
            delta_coeff = rhs_terms[rhs_terms_counter].as_coeff_Mul()[0]
            rhs_terms_counter += 1
            
        # Construct the left-hand side expression for the current equation
        for j in range(n + 1):
            derivative_count = 0
            if isinstance(lhs_terms[j].as_coeff_Mul()[1], sp.Derivative):
                derivative_count = lhs_terms[j].as_coeff_Mul()[1].derivative_count
            
            if derivative_count != j:
                continue
        
            if m - n >= i - j and j - (i + 1) < 0:
                lhs_equation_expression += lhs_terms[j].as_coeff_Mul()[0] * constants[n + i - j]
            elif j - (i + 1) >= 0:
                lhs_equation_expression += lhs_terms[j].as_coeff_Mul()[0] * sp.expand(y_homogeneous.rhs.diff((t, j - (i + 1))).replace(t, 0))
        eq_system.append(sp.Eq(lhs_equation_expression, delta_coeff))
    
    # Solve the system of equations for the constants
    constants_values = sp.solve(eq_system, constants)
    
    # Substitute the constants back into the homogeneous solution to find h(t)
    for i in range(len(constants)):
        h_t = h_t.replace(constants[i], constants_values[constants[i]])

    return sp.expand(sp.simplify(h_t))

def find_frequency_response(A, Y, W, N, M, w_t, y_t):
    """
    Finds the frequency response H(jω) of an LTI system given its transfer function.

    Parameters:
    A (list of lists): n x n matrix whose elements are multiples of {D, D^-1, 1}.
    Y (list): n x 1 vector of circuit variables.
    W (list): n x 1 vector, where the Mth element is w(t) and others are 0.
    N (int): Index of the desired signal y(t) in Y.
    M (tuple): Indices of the input signal w(t) in W.
    w_t (sympy Function): Input waveform w(t).
    y_t (sympy Function): Output response y(t).

    Returns:
    tuple: (H_jomega, amplitude(H_jomega), phase(H_jomega))
        - H_jomega: Frequency response H(jω).
        - amplitude(H_jomega): Amplitude response of H(jω).
        - phase(H_jomega): Phase response of H(jω).
    """
    # Find the transfer function H(s) of the system
    H_s = find_H_s(A, Y, W, N, M, w_t, y_t)
    
    # Replace s with jω to compute the frequency response H(jω)
    H_jomega = H_s.replace(s, sp.I * omega)

    return H_jomega, find_amplitude(H_jomega), find_phase(H_jomega)
    
def find_amplitude(F_jw):
    """
    Computes the amplitude of a complex function F(jω).

    Parameters:
    F_jw (sympy expression): A complex function of ω.

    Returns:
    sympy expression: The amplitude of the function F(jω).
    """
    # Get the numerator and denominator of the function
    numerator, denominator = sp.fraction(F_jw)
    
    # Compute the amplitude of the numerator and denominator
    numerator_amplitude = sp.Abs(numerator)
    denominator_amplitude = sp.Abs(denominator)
    
    # Return the amplitude response
    return sp.expand(numerator_amplitude) / sp.expand(denominator_amplitude)

def find_phase(F_jw):
    """
    Computes the phase of a complex function F(jω).

    Parameters:
    F_jw (sympy expression): A complex function of ω.

    Returns:
    sympy expression: The phase of the function F(jω).
    """
    # Get the numerator and denominator of the function
    numerator, denominator = sp.fraction(F_jw)
    
    # Compute the phase of the numerator and denominator
    numerator_phase = sp.arg(numerator)
    denominator_phase = sp.arg(denominator)
    
    # Return the phase response
    return sp.expand(numerator_phase) - sp.expand(denominator_phase)

    


#-------------------------------------------------------------------- GOLOBAL VARIABLES --------------------------------------------------------------------------


# Define the symbolic variables
t = sp.symbols('t')
D = sp.symbols('D')
s = sp.symbols('s')
j = I
omega = sp.symbols('\\omega') # For better looks replace with w
h_t = sp.Function('h')(t)
w_t = sp.Function('w')(t)
x_t = sp.Function('x')(t)
y_t = sp.Function('y')(t)
z_t = sp.Function('z')(t)


#-------------------------------------------------------------------- INPUTS OF THE APP ----------------------------------------------------------------------------

# Define the system matrices for the given problem
A = [
    [D**(-1), D, 3],
    [D, 2 * D, 3 * D - 3],
    [D**(-1), D, 2 * D - 3]
]
Y = [y_t, z_t, x_t]
W = [w_t, 0, 0]
M = 0
N = 0



#------------------------------------------------------------------ OUTPUTS AND CALCULATIONS --------------------------------------------------------------------------




# Analyze the system to get the differential equation
diff_eq = analyzer(A, Y, W, N, M, W[M], Y[N])

# Find the transfer function H(s)
H_s = find_H_s(A, Y, W, N, M, W[M], Y[N])

# Find the impulse response h(t)
h_t = find_impulse_response(A, Y, W, N, M, W[M], Y[N])

# Find the frequency response H(jω), its amplitude, and phase
H_jomega, H_jomega_amplitude, H_jomega_phase = find_frequency_response(A, Y, W, N, M, W[M], Y[N])

# Print the differential equation in LaTeX and pretty-print formats
print("Differential equation is: (it is minimized)")
# print(sp.latex(diff_eq))
# print("or")
sp.pprint(diff_eq)
print(" ")

# Print the impulse response h(t) in LaTeX and pretty-print formats
print("h(t) = ")
# print(sp.latex(h_t))
# print("or")
sp.pprint(h_t)
print(" ")

# Print the frequency response H(jω) in LaTeX and pretty-print formats
print("H(j\\omega) = ")
# print(sp.latex(H_jomega))
# print("or")
sp.pprint(H_jomega)
print(" ")

# Print the amplitude of H(jω) in LaTeX and pretty-print formats
print("amplitude(H(j\\omega)) = ")
# print(sp.latex(H_jomega_amplitude))
# print("or")
sp.pprint(H_jomega_amplitude)
print(" ")

# Print the phase of H(jω) in LaTeX and pretty-print formats
print("phase(H(j\\omega)) = ")
# print(sp.latex(H_jomega_phase))
# print("or")
sp.pprint(H_jomega_phase)
print(" ")