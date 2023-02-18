from pyavl import AVLSolver


avl_solver = AVLSolver(geo_file="aircraft.avl", debug=False)

avl_solver.addConstraint("alpha", 10)
avl_solver.executeRun()


