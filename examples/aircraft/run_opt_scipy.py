"""A scipy based optimization to trim an aircraft for an aicraft using optvl"""
import openmdao.api as om
import numpy as np
from scipy.optimize import minimize
from optvl import AVLSolver

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False)

# setup OptVL 
avl_solver.execute_run()

# Define your custom objective function
def custom_function(x):
    return x[0]**2 + x[1]**2 + 3*x[0]*x[1]

# Define the gradient of the custom objective function
def custom_gradient(x):
    df_dx0 = 2*x[0] + 3*x[1]
    df_dx1 = 2*x[1] + 3*x[0]
    return np.array([df_dx0, df_dx1])

# Define equality constraint: h(x) = 0
def eq_constraint(x):
    return x[0] + x[1] - 1

# Define the gradient of the equality constraint
def eq_constraint_jac(x):
    # Gradient of x[0] + x[1] - 1 is [1, 1]
    return np.array([1, 1])

# Define inequality constraint: g(x) >= 0
def ineq_constraint(x):
    return x[0] - 0.5

# Define the gradient of the inequality constraint
def ineq_constraint_jac(x):
    # Gradient of x[0] - 0.5 is [1, 0]
    return np.array([1, 0])

# Initial guess for the variables
x0 = np.array([0.6, 0.4])

# Define constraints with their gradients
constraints = [
    {'type': 'eq', 'fun': eq_constraint, 'jac': eq_constraint_jac},   # Equality constraint
    {'type': 'ineq', 'fun': ineq_constraint, 'jac': ineq_constraint_jac}  # Inequality constraint
]

# Call the minimize function with constraints and their gradients
result = minimize(
    fun=custom_function,       # The objective function to minimize
    jac=custom_gradient,       # The gradient function of the objective
    x0=x0,                     # Initial guess
    constraints=constraints,   # Constraints with gradients
    method='SLSQP',            # Optimization method that supports constraints
    options={'disp': True}     # Display convergence messages
)

# Print the result
print("Optimization result:")
print("x:", result.x)              # Optimized variables
print("fun:", result.fun)          # Final value of the objective function
print("success:", result.success)  # Whether the optimizer exited successfully
print("message:", result.message)  # Description of the cause of termination