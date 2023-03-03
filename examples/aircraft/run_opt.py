"""A openmdao based optimization for sigma aicraft using pyavl"""
import openmdao.api as om
from pyavl import AVLSolver


class AVLSolverComp(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("geom_file", types=str)

    def setup(self):
        geom_file = self.options["geom_file"]
        self.avl_solver = AVLSolver(geo_file=geom_file, debug=False)

        # add the control surfaces as inputs
        self.control_names = self.avl_solver.get_control_names()
        for c_name in self.control_names:
            self.add_input(c_name, val=0.0, units="deg")

        # add the geometric parameters as inputs
        surf_data = self.avl_solver.get_surface_params(geom_only=True)

        # HACK: for now, just add the aincs as an input. need to talk with opendao team about not computing full jacobian
        # for surf in surf_data:
        #     for key in surf_data[surf]:
        # self.add_input(f"{surf}:{key}", val=surf_data[surf][key], tags="geom")

        surf = "Wing"
        key = "aincs"
        self.add_input(f"{surf}:{key}", val=surf_data[surf][key], tags="geom")

        # add the outputs
        # forces
        self.add_output("CL", val=0.0)
        self.add_output("CD", val=0.0)

        # self.add_output('Cl', val=0.0, desc="about stability axis")
        self.add_output("Cm", val=0.0, desc="about stability axis")
        # self.add_output('Cn', val=0.0, desc="about stability axis")

        # # the cowards way out
        self.declare_partials("*", "*", method="fd")

    def compute(self, inputs, outputs):
        for c_name in self.control_names:
            self.avl_solver.add_constraint(c_name, inputs[c_name][0])

        # update the surface parameters
        surf_data = self.avl_solver.get_surface_params()

        for input_var in self.list_inputs(tags="geom", val=False, out_stream=None):
            # split the input name into surface name and parameter name
            surf, param = input_var[0].split(":")
            # update the corresponding parameter in the surface data
            surf_data[surf][param] = inputs[input_var[0]]


        # self.avl_solver.execute_run()
        # print('aincs', inputs["Wing:aincs"], f"CL:{self.avl_solver.CL} CD:{self.avl_solver.CD} Cm:{self.avl_solver.CM}")
        self.avl_solver.set_surface_params(surf_data)


        self.avl_solver.execute_run()

        # set the outputs
        run_data = self.avl_solver.get_case_total_data()
        outputs["CL"] = run_data["CL"]
        outputs["CD"] = run_data["CD"]
        outputs["Cm"] = run_data["Cm"]
        print(
            "aincs",
            inputs["Wing:aincs"],
            inputs["Elevator"],
            f"CL:{self.avl_solver.CL} CD:{self.avl_solver.CD} Cm:{self.avl_solver.CM}",
        )


model = om.Group()
model.add_subsystem("des_vars", om.IndepVarComp())
model.des_vars.add_output(
    "Wing:aincs",
    shape_by_conn=True,
    val=0.2,
)
model.des_vars.add_output(
    "Elevator",
    val=0.0,
)
model.connect("des_vars.Elevator", "avlsolver.Elevator")
model.add_subsystem("avlsolver", AVLSolverComp(geom_file="aircraft.avl"))
model.connect("des_vars.Wing:aincs", "avlsolver.Wing:aincs")

model.add_design_var("des_vars.Wing:aincs", lower=-10, upper=10)
model.add_design_var("des_vars.Elevator", lower=-10, upper=10)
model.add_constraint("avlsolver.CL", equals=0.5)
model.add_constraint("avlsolver.Cm", equals=0.0)
model.add_objective("avlsolver.CD")
prob = om.Problem(model)
prob.driver = om.pyOptSparseDriver()
prob.driver.options["optimizer"] = "SLSQP"
prob.driver.opt_settings = {
    "IFILE": "SLSQP_print.out",
}
prob.setup()
om.n2(prob, show_browser=False, outfile="vlm_opt.html")
prob.run_driver()
