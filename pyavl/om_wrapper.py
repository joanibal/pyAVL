import os
import openmdao.api as om
from pyavl import AVLSolver
import numpy as np
import copy
import time

class AVLGroup(om.Group):
    def initialize(self):
        self.options.declare("geom_file", types=str)
        self.options.declare("mass_file", default=None)
        self.options.declare("write_grid", types=bool, default=False)
        self.options.declare("output_dir", types=str, recordable=False, default=".")

        self.options.declare("input_param_vals", types=bool, default=False)
        self.options.declare("input_ref_vals", types=bool, default=False)
        
        self.options.declare("output_stabililty_derivs", types=bool, default=False)
        self.options.declare("output_con_surf_derivs", types=bool, default=False)

    def setup(self):
        geom_file = self.options["geom_file"]
        mass_file = self.options["mass_file"]
        
        input_param_vals = self.options["input_param_vals"]
        input_ref_vals = self.options["input_ref_vals"]
        
        output_stabililty_derivs = self.options["output_stabililty_derivs"]
        output_con_surf_derivs = self.options["output_con_surf_derivs"]

        avl = AVLSolver(geo_file=geom_file, mass_file=mass_file, debug=False)

        self.add_subsystem("solver", AVLSolverComp(avl=avl, 
                                    input_param_vals=input_param_vals,
                                    input_ref_vals=input_ref_vals),
                                    promotes=["*"])
        self.add_subsystem("funcs", AVLFuncsComp(avl=avl,
                                    input_param_vals=input_param_vals,
                                    input_ref_vals=input_ref_vals,
                                    output_stabililty_derivs=output_stabililty_derivs,
                                    output_con_surf_derivs=output_con_surf_derivs),
                                    promotes=["*"])
        if self.options["write_grid"]:
            self.add_subsystem(
                "postprocess", AVLPostProcessComp(avl=avl, output_dir=self.options["output_dir"]),
                                                    input_param_vals=input_param_vals,
                                                    input_ref_vals=input_ref_vals,
                                                    promotes=["*"]
                                                )

# helper functions used by the AVL components
def add_avl_controls_as_inputs(self, avl):
    # add the control surfaces as inputs
    self.control_names = avl.get_control_names()
    for c_name in self.control_names:
        self.add_input(c_name, val=0.0, units="deg", tags="con_surf")
    return self.control_names

def add_avl_geom_vars(self, avl, add_as="inputs"):
    # add the geometric parameters as inputs
    surf_data = avl.get_surface_params()

    for surf in surf_data:
        for key in surf_data[surf]:
            geom_key = f"{surf}:{key}"
            if add_as == "inputs":
                self.add_input(geom_key, val=surf_data[surf][key], tags="geom")
            elif add_as == "outputs":
                self.add_output(geom_key, val=surf_data[surf][key], tags="geom")
                
def add_avl_conditions_as_inputs(sys, avl):
    # TODO: add all the condition constraints
        
    sys.add_input("alpha", val=0.0, units="deg", tags="flt_cond")
    sys.add_input("beta", val=0.0, units="deg" , tags="flt_cond")


def add_avl_params_as_inputs(sys, avl):
    # TODO: add all par vals with the analysis is supported
    
    # only adding the ones people would use for now
    for param in ["CD0", "Mach", "X cg", "Y cg", "Z cg"]:
        val = avl.get_case_parameter(param)
        sys.add_input(param, val=val, tags="param")

def add_avl_refs_as_inputs(sys, avl):
    ref_data = avl.get_reference_data()
    
    for key, val in ref_data.items():
        sys.add_input(key, val=val, tags="ref_val")


def om_input_to_surf_dict(sys, inputs):
    geom_inputs = sys.list_inputs(tags="geom", val=False, out_stream=None)
    # convert to a list witout tuples
    geom_inputs = [x[0] for x in geom_inputs]

    surf_data = {}
    for input_var in inputs:
        if input_var in geom_inputs:
            # split the input name into surface name and parameter name
            surf, param = input_var.split(":")

            # update the corresponding parameter in the surface data
            if surf not in surf_data:
                surf_data[surf] = {}
            if inputs[input_var].size == 1:
                # if the input is a scalar, convert to a 1D array for numpy
                surf_data[surf][param] = inputs[input_var][0]
            else:
                surf_data[surf][param] = inputs[input_var]
    return surf_data


def om_surf_dict_to_input(surf_dict):
    input_data = {}
    for surf_key in surf_dict:
        for geom_key in surf_dict[surf_key]:
            input_var = ":".join([surf_key, geom_key])

            input_data[input_var] = surf_dict[surf_key][geom_key]

    return input_data

def om_set_avl_inputs(sys, inputs):
    for c_name in sys.control_names:
        sys.avl.add_constraint(c_name, inputs[c_name][0])

    sys.avl.add_constraint("alpha", inputs["alpha"][0])
    sys.avl.add_constraint("beta", inputs["beta"][0])

    # add the parameters to the run
    for param in sys.avl.param_idx_dict:
        if param in inputs:
            val = inputs[param][0]
            sys.avl.set_case_parameter(param, val)
    
    # add the parameters to the run
    for ref in sys.avl.ref_var_to_fort_var:
        if ref in inputs:
            val = inputs[ref][0]
            sys.avl.set_reference_data({ref: val})
    


class AVLSolverComp(om.ImplicitComponent):
    """
    OpenMDAO component that wraps pyAVL
    """

    def initialize(self):
        self.options.declare("avl", types=AVLSolver, recordable=False)
        self.options.declare("input_param_vals", types=bool, default=False)
        self.options.declare("input_ref_vals", types=bool, default=False)
        

    def setup(self):
        self.avl = self.options["avl"]
        input_param_vals = self.options["input_param_vals"]
        input_ref_vals = self.options["input_ref_vals"]
        
        self.num_states = self.avl.get_mesh_size()
        self.num_cs = self.avl.get_num_control_surfs()

        self.add_output("gamma", val=np.zeros(self.num_states))
        self.add_output("gamma_d", val=np.zeros((self.num_cs, self.num_states)))
        
        add_avl_conditions_as_inputs(self, self.avl)
        
        if input_param_vals:
            add_avl_params_as_inputs(self, self.avl)
        
        if input_ref_vals:
            add_avl_refs_as_inputs(self, self.avl)
        
        self.control_names = add_avl_controls_as_inputs(self, self.avl)
        add_avl_geom_vars(self, self.avl, add_as="inputs")

    def apply_nonlinear(self, inputs, outputs, residuals):
        om_set_avl_inputs(self, inputs)
        
        surf_data = om_input_to_surf_dict(self, inputs)
        self.avl.set_surface_params(surf_data)

        gam_arr = outputs["gamma"]
        gam_d_arr = outputs["gamma_d"]

        self.avl.set_avl_fort_arr("VRTX_R", "GAM", gam_arr, slicer=(slice(0, gam_arr.size),))
        self.avl.set_avl_fort_arr(
            "VRTX_R", "GAM_D", gam_d_arr, slicer=(slice(0, self.num_cs), slice(0, self.num_states))
        )

        # propogate the seeds through without resolving
        self.avl.avl.update_surfaces()
        self.avl.avl.get_res()

        res_slice = (slice(0, self.num_states),)
        res_d_slice = (slice(0, self.num_cs), slice(0, self.num_states))

        res = copy.deepcopy(self.avl.get_avl_fort_arr("VRTX_R", "RES", slicer=res_slice))
        residuals["gamma"] = res
        res_d = copy.deepcopy(self.avl.get_avl_fort_arr("VRTX_R", "RES_D", slicer=res_d_slice))
        residuals["gamma_d"] = res_d
        # this shouldn't happen

    def solve_nonlinear(self, inputs, outputs):
        start_time = time.time()
        om_set_avl_inputs(self, inputs)

        # update the surface parameters
        surf_data = om_input_to_surf_dict(self, inputs)
        self.avl.set_surface_params(surf_data)

        # def_dict = self.avl.get_control_deflections()
        print('executing avl run')
        self.avl.execute_run()

        gam_arr = self.avl.get_avl_fort_arr("VRTX_R", "GAM", slicer=(slice(0, self.num_states),))

        outputs["gamma"] = copy.deepcopy(gam_arr)

        gam_d_arr = self.avl.get_avl_fort_arr(
            "VRTX_R", "GAM_D", slicer=(slice(0, self.num_cs), slice(0, self.num_states))
        )
        outputs["gamma_d"] = copy.deepcopy(gam_d_arr)

        # run_data = self.avl.get_case_total_data()
        # for func_key in run_data:
        #     print(func_key, run_data[func_key])
        # func_key = "CL"
        # print(func_key, run_data[func_key])
        # print("AVL solve time: ", time.time() - start_time)

    def apply_linear(self, inputs, outputs, d_inputs, d_outputs, d_residuals, mode):
        if mode == "fwd":
            con_seeds = {}
            for con_key in ["alpha", "beta"]:
                if con_key in d_inputs:
                    con_seeds[con_key] = d_inputs[con_key]
            for con_key in self.control_names:
                if con_key in d_inputs:
                    con_seeds[con_key] = d_inputs[con_key]
                    
            geom_seeds = om_input_to_surf_dict(self, d_inputs)
            
            
            param_seeds = {}
            for param in self.avl.param_idx_dict:
                if param in d_inputs:
                    param_seeds[param] = d_inputs[param]
                
            ref_seeds = {}
            for ref in self.avl.ref_var_to_fort_var:
                if ref in d_inputs:
                    ref_seeds[ref] = d_inputs[ref]

            _, res_seeds, _, res_d_seeds = self.avl.execute_jac_vec_prod_fwd(con_seeds=con_seeds, geom_seeds=geom_seeds,
                                                                           param_seeds=param_seeds, ref_seeds=ref_seeds)

            d_residuals["gamma"] += res_seeds
            d_residuals["gamma_d"] += res_seeds

        if mode == "rev":
            if "gamma" in d_residuals:
                self.avl.clear_ad_seeds_fast()
                res_seeds = d_residuals["gamma"]
                res_d_seeds = d_residuals["gamma_d"]

                con_seeds, geom_seeds, gamma_seeds, gamma_d_seeds, param_seeds, ref_seeds = self.avl.execute_jac_vec_prod_rev(
                    res_seeds=res_seeds, res_d_seeds=res_d_seeds
                )

                if "gamma" in d_outputs:
                    d_outputs["gamma"] += gamma_seeds

                if "gamma_d" in d_outputs:
                    d_outputs["gamma_d"] += gamma_d_seeds

                d_input_geom = om_surf_dict_to_input(geom_seeds)

                for d_input in d_inputs:
                    if d_input in d_input_geom:
                        d_inputs[d_input] += d_input_geom[d_input]
                    elif d_input in ["alpha", "beta"]:
                        d_inputs[d_input] += con_seeds[d_input]
                    elif d_input in self.control_names:
                        d_inputs[d_input] += con_seeds[d_input]
                    elif d_input in param_seeds:
                        d_inputs[d_input] += param_seeds[d_input]
                    elif d_input in ref_seeds:
                        d_inputs[d_input] += ref_seeds[d_input]

                  

    def solve_linear(self, d_outputs, d_residuals, mode):
        if mode == "rev":
            self.avl.set_gamma_ad_seeds(d_outputs["gamma"])
            self.avl.set_gamma_d_ad_seeds(d_outputs["gamma_d"])
            # start_time = time.time()
            self.avl.avl.solve_adjoint()
            # print("OM Solve adjoint time: ", time.time() - start_time)
            d_residuals["gamma"] = self.avl.get_residual_ad_seeds()
            d_residuals["gamma_d"] = self.avl.get_residual_d_ad_seeds()
        elif mode == "fwd":
            raise NotImplementedError("only reverse mode derivaties implemented. Use prob.setup(mode='rev')")


class AVLFuncsComp(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("avl", types=AVLSolver, recordable=False)
        self.options.declare("output_stabililty_derivs", types=bool, default=False)
        self.options.declare("output_con_surf_derivs", types=bool, default=False)
        self.options.declare("input_param_vals", types=bool, default=False)
        self.options.declare("input_ref_vals", types=bool, default=False)
        

    def setup(self):
        self.avl = self.options["avl"]
        self.num_states = self.avl.get_mesh_size()
        self.num_cs = self.avl.get_num_control_surfs()
        input_param_vals = self.options["input_param_vals"]
        input_ref_vals = self.options["input_ref_vals"]
        
        self.add_input("gamma", val=np.zeros(self.num_states))
        self.add_input("gamma_d", val=np.zeros((self.num_cs, self.num_states)))

        add_avl_conditions_as_inputs(self, self.avl)

        if input_param_vals:
            add_avl_params_as_inputs(self, self.avl)
        
        if input_ref_vals:
            add_avl_refs_as_inputs(self, self.avl)
            
        self.control_names = add_avl_controls_as_inputs(self, self.avl)
        add_avl_geom_vars(self, self.avl, add_as="inputs")

        # add the outputs
        for func_key in self.avl.case_var_to_fort_var:
            self.add_output(func_key)

        self.output_con_surf_derivs = self.options["output_con_surf_derivs"]
        if self.output_con_surf_derivs:
            for func_key in self.avl.case_derivs_to_fort_var:
                for con_name in self.control_names:
                    var_name = f"d{func_key}_d{con_name}"
                    self.add_output(var_name)
        
        self.output_stabililty_derivs = self.options["output_stabililty_derivs"]
        if self.output_stabililty_derivs:
            deriv_dict = self.avl.case_stab_derivs_to_fort_var
            for func_key in deriv_dict:
                for var in deriv_dict[func_key]:
                    var_name = f"d{func_key}_d{var}"
                    self.add_output(var_name)        

    def compute(self, inputs, outputs):
        # self.avl.set_gamma(inputs['gamma'])
        # for c_name in self.control_names:
        #     self.avl.add_constraint(c_name, inputs[c_name][0])
        # def_dict = self.avl.get_control_deflections()

        # TODO: add_constraint does not correctly do derives yet
        start_time = time.time()
        om_set_avl_inputs(self, inputs)

        # update the surface parameters
        surf_data = om_input_to_surf_dict(self, inputs)
        self.avl.set_surface_params(surf_data)

        gam_arr = inputs["gamma"]
        gam_d_arr = inputs["gamma_d"]

        self.avl.set_avl_fort_arr("VRTX_R", "GAM", gam_arr, slicer=(slice(0, gam_arr.size),))
        self.avl.set_avl_fort_arr(
            "VRTX_R", "GAM_D", gam_d_arr, slicer=(slice(0, self.num_cs), slice(0, self.num_states))
        )
        
        # TODO: only update what you need to. 
        # residuals (and AIC?) do not need to be calculated 
        # but in get_res, alpha and beta are set
        self.avl.avl.update_surfaces()
        self.avl.avl.get_res()
        self.avl.avl.velsum()
        self.avl.avl.aero()

        run_data = self.avl.get_case_total_data()

        for func_key in run_data:
            # print(f' {func_key} {run_data[func_key]}')
            outputs[func_key] = run_data[func_key]
        # print(f" CD {run_data['CD']} CL {run_data['CL']}")

        if self.output_con_surf_derivs:
            consurf_derivs_seeds = self.avl.get_case_coef_derivs()
            for func_key in consurf_derivs_seeds:
                for con_name in consurf_derivs_seeds[func_key]:
                    var_name = f"d{func_key}_d{con_name}"
                    outputs[var_name] = consurf_derivs_seeds[func_key][con_name]

        if self.output_stabililty_derivs:
            stab_derivs = self.avl.get_case_stab_derivs()
            for func_key in stab_derivs:
                for var in stab_derivs[func_key]:
                    var_name = f"d{func_key}_d{var}"
                    outputs[var_name] = stab_derivs[func_key][var]
                
        # print("Funcs Compute time: ", time.time() - start_time)

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        if mode == "fwd":
            con_seeds = {}
            for con_key in ["alpha", "beta"]:
                if con_key in d_inputs:
                    con_seeds[con_key] = d_inputs[con_key]
            
            for con_key in self.control_names:
                if con_key in d_inputs:
                    con_seeds[con_key] = d_inputs[con_key]
                    
            if "gamma" in d_inputs:
                gamma_seeds = d_inputs["gamma"]
            else:
                gamma_seeds = None

            if "gamma_d" in d_inputs:
                gamma_d_seeds = d_inputs["gamma_d"]
            else:
                gamma_d_seeds = None
                
            param_seeds = {}
            for param in self.avl.param_idx_dict:
                if param in d_inputs:
                    param_seeds[param] = d_inputs[param]
                
            ref_seeds = {}
            for ref in self.avl.ref_var_to_fort_var:
                if ref in d_inputs:
                    ref_seeds[ref] = d_inputs[ref]

            geom_seeds = self.om_input_to_surf_dict(self, d_inputs)

            func_seeds, _, csd_seeds, _ = self.avl.execute_jac_vec_prod_fwd(
                con_seeds=con_seeds, geom_seeds=geom_seeds, gamma_seeds=gamma_seeds, gamma_d_seeds=gamma_d_seeds,
                param_seeds=param_seeds, ref_seeds=ref_seeds
            )

            for func_key in func_seeds:
                d_outputs[func_key] += func_seeds[func_key]

            for func_key in csd_seeds:
                for con_name in csd_seeds[func_key]:
                    var_name = f"d{func_key}_d{con_name}"
                    d_outputs[var_name] = csd_seeds[func_key][con_name]

        if mode == "rev":
            self.avl.clear_ad_seeds_fast()

            # add the outputs
            func_seeds = {}
            for func_key in self.avl.case_var_to_fort_var:
                if func_key in d_outputs:
                    func_seeds[func_key] = d_outputs[func_key]
                    if np.abs(func_seeds[func_key]) > 0.0:
                        print(f'  running rev mode derivs for {func_key}')

            csd_seeds = {}
            con_names = self.avl.get_control_names()
            for func_key in self.avl.case_derivs_to_fort_var:
                csd_seeds[func_key] = {}
                for con_name in con_names:
                    var_name = f"d{func_key}_d{con_name}"

                    if var_name in d_outputs:
                        csd_seeds[func_key][con_name] = d_outputs[var_name]
                        
                        if np.abs(csd_seeds[func_key][con_name]) > 0.0:
                            # print(var_name, csd_seeds[func_key][con_name])
                            print(f'  running rev mode derivs for {func_key}:{con_name}')
                        
            con_seeds, geom_seeds, gamma_seeds, gamma_d_seeds, param_seeds, ref_seeds = self.avl.execute_jac_vec_prod_rev(
                func_seeds=func_seeds, consurf_derivs_seeds=csd_seeds
            )

            if "gamma" in d_inputs:
                d_inputs["gamma"] += gamma_seeds

            if "gamma_d" in d_inputs:
                d_inputs["gamma_d"] += gamma_d_seeds

            d_input_geom = om_surf_dict_to_input(geom_seeds)

            for d_input in d_inputs:
                if d_input in d_input_geom:
                    d_inputs[d_input] += d_input_geom[d_input]
                elif d_input in ["alpha", "beta"] or d_input in self.control_names:
                    d_inputs[d_input] += con_seeds[d_input]
                elif d_input in param_seeds:
                    d_inputs[d_input] += param_seeds[d_input]
                elif d_input in ref_seeds:
                    d_inputs[d_input] += ref_seeds[d_input]


# Optional components
class AVLPostProcessComp(om.ExplicitComponent):
    def initialize(self):
        self.options.declare("avl", types=AVLSolver, recordable=False)
        self.options.declare("output_dir", types=str, recordable=False, default=".")
        self.options.declare("input_param_vals", types=bool, default=False)
        self.options.declare("input_ref_vals", types=bool, default=False)
        
    def setup(self):
        self.avl = self.options["avl"]
        self.num_states = self.avl.get_mesh_size()
        self.num_cs = self.avl.get_num_control_surfs()
        input_param_vals = self.options["input_param_vals"]
        input_ref_vals = self.options["input_ref_vals"]
        

        self.add_input("gamma", val=np.zeros(self.num_states))
        self.add_input("gamma_d", val=np.zeros((self.num_cs, self.num_states)))

        add_avl_conditions_as_inputs(self, self.avl)
        
        if input_param_vals:
            add_avl_params_as_inputs(self, self.avl)
        
        if input_ref_vals:
            add_avl_refs_as_inputs(self, self.avl)
            
        add_avl_controls_as_inputs(self, self.avl)
        add_avl_geom_vars(self, self.avl, add_as="inputs")

        self.iter_count = 0

    def compute(self, inputs, outputs):
        # self.avl.set_gamma(inputs['gamma'])
        # for c_name in self.control_names:
        #     self.avl.add_constraint(c_name, inputs[c_name][0])
        # def_dict = self.avl.get_control_deflections()

        om_set_avl_inputs(self, inputs)

        # update the surface parameters
        surf_data = om_input_to_surf_dict(self, inputs)
        self.avl.set_surface_params(surf_data)

        gam_arr = inputs["gamma"]

        self.avl.set_avl_fort_arr("VRTX_R", "GAM", gam_arr, slicer=(slice(0, gam_arr.size),))

        gam_d_arr = inputs["gamma_d"]
        self.avl.set_avl_fort_arr(
            "VRTX_R", "GAM_D", gam_d_arr, slicer=(slice(0, gam_d_arr.shape[0]), slice(0, gam_d_arr.shape[1]))
        )

        file_name = f"vlm_{self.iter_count:03d}.avl"
        output_dir = self.options["output_dir"]
        self.avl.write_geom_file(os.path.join(output_dir, file_name))

        self.iter_count += 1

class AVLMeshReader(om.ExplicitComponent):
    """
        This class is moslty used to provide an initial set of coordinates for custom paramerization components
    """
    
    def initialize(self):
        self.options.declare("geom_file", types=str)
        self.options.declare("mass_file", default=None)

    def setup(self):
        geom_file = self.options["geom_file"]
        mass_file = self.options["mass_file"]

        avl = AVLSolver(geo_file=geom_file, mass_file=mass_file, debug=False)
        add_avl_geom_vars(self, avl, add_as="outputs")
