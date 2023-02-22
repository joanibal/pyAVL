"""A openmdao based optimization for sigma aicraft using pyavl"""
import openmao.api as om
from pyavl import AVLSolver

class AVLSolverComp(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('geom_file', types=str)
    
    def setup(self):
        geom_file = self.options['geom_file']
        self.avl_solver = AVLSolver(geo_file=geom_file, debug=False)
        
        # add the control surfaces as inputs 
        self.control_names = self.get_control_names()
        for c_name in self.control_names:
            self.add_input(c_name, val=0.0, units='deg', tag='geom')
        
        # add the geometric parameters as inputs 
        surf_data = self.avl_solver.get_surface_params()
        
        for surf in surf_data:
            for key in surf_data[surf]:
                self.add_input(f"{surf}:{key}", val=surf_data[surf][key])
        
        # add the outputs
        # forces
        self.add_output('CL', val=0.0)
        self.add_output('CD', val=0.0)
        
        # self.add_output('Cl', val=0.0, desc="about stability axis")
        self.add_output('Cm', val=0.0, desc="about stability axis")
        # self.add_output('Cn', val=0.0, desc="about stability axis")
    
    def compute(self, inputs, outputs):
        
        for c_name in self.control_names:
            self.avl_solver.addConstraint(c_name, inputs[c_name])
        
        # update the surface parameters
        surf_data = self.avl_solver.get_surface_params()
        
        for input_var in self.list_inputs(tags='geom'):
            # split the input name into surface name and parameter name
            surf, param = input_var.split(':')
            # update the corresponding parameter in the surface data
            surf_data[surf][param] = inputs[input_var]
        
        self.avl_solver.resetData()    
        self.avl_solver.executeRun()
        
        # set the outputs
        
        outputs['CL'] = self.avl_solver.CL
        outputs['CD'] = self.avl_solver.CD
        outputs['CM'] = self.avl_solver.CM

        
model = om.Group()
model.add_subsystem('avlsolver', AVLSolverComp(geom_file="aircraft.avl"))

prob = om.Problem(model)
prob.setup()