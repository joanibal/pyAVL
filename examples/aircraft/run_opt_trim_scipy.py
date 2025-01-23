"""A scipy based optimization to trim an aircraft for an aicraft using optvl"""
import openmdao.api as om
import numpy as np
from scipy.optimize import minimize
from optvl import AVLSolver

avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False)

# setup OptVL 
avl_solver.set_case_parameter("Mach", 0.0)
avl_solver.add_constraint("alpha", 5.0)

# Define your custom objective function with outputs from OptVL
def custom_function(x):
    avl_solver.add_constraint("Elevator", x[0])
    
    avl_solver.set_surface_params({"Wing":{"aincs":x[1:]}})
    
    avl_solver.execute_run()
    cd = avl_solver.get_case_total_data()['CD']
    return cd

# Define the gradient (Jacobian) of the objective function
def custom_gradient(x):
    # Partial derivatives of the custom_function
    sens = avl_solver.execute_run_sensitivies(['CD'])
    dcd_dele = sens['CD']['Elevator']
    dcd_daincs = sens['CD']['Wing']['aincs']
    # concatinate the two and return the derivs
    return 

# Initial guess for the variables
x0 = np.concatenate(([0], [0]*5))
# Call the minimize function
result = minimize(
    fun=custom_function,       # The objective function to minimize
    x0=x0,                     # Initial guess
    jac=custom_gradient,       # The gradient function
    method='BFGS',             # Optimization method that supports gradient
    options={'disp': True}     # Display convergence messages
)

# Print the result
print("Optimization result:")
print("x:", result.x)              # Optimized variables
print("fun:", result.fun)          # Final value of the objective function
print("success:", result.success)  # Whether the optimizer exited successfully
print("message:", result.message)  # Description of the cause of termination