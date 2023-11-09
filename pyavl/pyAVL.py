"""
pyAVL

pyAVL is a wrapper for Mark Drela's AVL code. The purpose of this
class is to provide an easy to use wrapper for avl for intergration
into other projects.

Developers:
-----------
- Josh Anibal (JLA)

History
-------
    v. 0.0 - Initial Class Creation (JLA, 08 2016)
    v. 0.1 - Added CDFF and sectional properties for all surfaces (JLA, 12 2017)
    v. 1.0 - major refactor (I learned a lot more programming in 7 years) (JLA, 02 2023)
"""

# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
import copy
from pprint import pprint
from typing import Dict, List, Tuple, Union, Any
import warnings
import glob
from typing import Optional

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np


# =============================================================================
# Extension modules
# =============================================================================
from . import MExt


class AVLSolver(object):
    con_var_to_fort_var = {
        "alpha": ["CASE_R", "ALFA"],
        "beta": ["CASE_R", "BETA"],
    }

    param_idx_dict = {
        "alpha": 0,
        "beta": 1,
        "roll rate": 2,
        "pitch rate": 3,
        "yaw rate": 4,
        "CL": 5,
        "CD0": 6,
        "bank": 7,
        "elevation": 8,
        "heading": 9,
        "Mach": 10,
        "velocity": 11,
        "density": 12,
        "grav.acc.": 13,
        "turn rad.": 14,
        "load fac.": 15,
        "X cg": 16,
        "Y cg": 17,
        "Z cg": 18,
        "mass": 19,
        "Ixx": 20,
        "Iyy": 21,
        "Izz": 22,
        "Ixy": 23,
        "Iyz": 24,
        "Izx": 25,
        "visc CL_a": 26,
        "visc CL_u": 27,
        "visc CM_a": 28,
        "visc CM_u": 29,
    }
    conval_idx_dict = {
        "alpha": 0,
        "beta": 1,
        "roll rate": 2,
        "pitch rate": 3,
        "yaw rate": 4,
        "CL": 5,
        "CY": 6,
        "CR BA": 7,
        "CM": 8,
        "CR": 9,
    }

    # fmt: off
    # This dict has the following structure:
    # python key: [common block name, fortran varaiable name]
    case_var_to_fort_var = {
        # lift and drag from surface integration (wind frame)
        "CL": ["CASE_R", "CLTOT"],
        "CD": ["CASE_R", "CDTOT"],
        "CDv": ["CASE_R", "CDVTOT"], # viscous drag
        
        # lift and drag calculated from farfield integration
        "CLff": ["CASE_R", "CLFF"],
        "CYff": ["CASE_R", "CYFF"],
        "CDi": ["CASE_R", "CDFF"], # induced drag
        
        # non-dimensionalized forces 
        "CX": ["CASE_R", "CXTOT"], 
        "CY": ["CASE_R", "CYTOT"], 
        "CZ": ["CASE_R", "CZTOT"],   
        
        # non-dimensionalized moments (body frame)
        "CR BA": ["CASE_R", "CRTOT"],
        "CM": ["CASE_R", "CMTOT"],
        "CN BA": ["CASE_R", "CNTOT"],

        # non-dimensionalized moments (stablity frame)
        "CR SA": ["CASE_R", "CRSAX"],
        # "CM SA": ["CASE_R", "CMSAX"], # This is the same in both frames
        "CN SA": ["CASE_R", "CNSAX"],
        
        # spanwise efficiency
        "e": ["CASE_R", "SPANEF"],
    }
    
    case_derivs_to_fort_var = {
        # derivative of coefficents wrt control surface deflections
        "CL": ["CASE_R", "CLTOT_D"],
        "CD": ["CASE_R", "CDTOT_D"],
        "CX": ["CASE_R", "CXTOT_D"],
        "CY": ["CASE_R", "CYTOT_D"],
        "CZ": ["CASE_R", "CZTOT_D"],
        "CR": ["CASE_R", "CRTOT_D"],
        "CM": ["CASE_R", "CMTOT_D"],
        "CN": ["CASE_R", "CNTOT_D"],
    }
    
    case_stab_derivs_to_fort_var = {
        # derivative of coefficents wrt alpha
        "CL": {
            "alpha": ["CASE_R", "CLTOT_AL"],
            "beta": ["CASE_R", "CLTOT_BE"],
            "roll rate": ["CASE_R", "CLTOT_RX"],
            "pitch rate": ["CASE_R", "CLTOT_RY"],
            "yaw rate": ["CASE_R", "CLTOT_RZ"],
            },
        "CD": {
            "alpha": ["CASE_R", "CDTOT_AL"],
            "beta": ["CASE_R", "CDTOT_BE"],
            "roll rate": ["CASE_R", "CDTOT_RX"],
            "pitch rate": ["CASE_R", "CDTOT_RY"],
            "yaw rate": ["CASE_R", "CDTOT_RZ"],
            },
        "CY": {
            "alpha": ["CASE_R", "CYTOT_AL"],
            "beta": ["CASE_R", "CYTOT_BE"],
            "roll rate": ["CASE_R", "CYTOT_RX"],
            "pitch rate": ["CASE_R", "CYTOT_RY"],
            "yaw rate": ["CASE_R", "CYTOT_RZ"],
            },
        "CR SA": {
            "alpha": ["CASE_R", "CRTOT_AL"],
            "beta": ["CASE_R", "CRTOT_BE"],
            "roll rate": ["CASE_R", "CRTOT_RX"],
            "pitch rate": ["CASE_R", "CRTOT_RY"],
            "yaw rate": ["CASE_R", "CRTOT_RZ"],
            },
        "CM": {
            "alpha": ["CASE_R", "CMTOT_AL"],
            "beta": ["CASE_R", "CMTOT_BE"],
            "roll rate": ["CASE_R", "CMTOT_RX"],
            "pitch rate": ["CASE_R", "CMTOT_RY"],
            "yaw rate": ["CASE_R", "CMTOT_RZ"],
            },
        "CN SA": {
            "alpha": ["CASE_R", "CNTOT_AL"],
            "beta": ["CASE_R", "CNTOT_BE"],
            "roll rate": ["CASE_R", "CNTOT_RX"],
            "pitch rate": ["CASE_R", "CNTOT_RY"],
            "yaw rate": ["CASE_R", "CNTOT_RZ"],
            },
    }
    
    # This dict has the following structure:
    # python key: [common block name, fortran varaiable name]
    case_surf_var_to_fort_var = {
        # surface contributions to total lift and drag from surface integration (wind frame)
        "CL": ["SURF_R", "CLSURF"],
        "CD": ["SURF_R", "CDSURF"],
        "CDv": ["SURF_R", "CDVSURF"],  # viscous drag
        
        # non-dimensionalized forces 
        "CX": ["SURF_R", "CXSURF"], 
        "CY": ["SURF_R", "CYSURF"], 
        "CZ": ["SURF_R", "CZSURF"],   
        
        # non-dimensionalized moments (body frame)
        "CR": ["SURF_R", "CRSURF"],
        "CM": ["SURF_R", "CMSURF"],
        "CN": ["SURF_R", "CNSURF"],
        
        # forces non-dimentionalized by surface quantities
        # uses surface area instead of sref and takes moments about leading edge
        "CL surf" : ["SURF_R", "CL_SRF"],
        "CD surf" : ["SURF_R", "CD_SRF"],
        "CMLE surf" : ["SURF_R", "CMLE_SRF"],
        
        #TODO: add CF_SRF(3,NFMAX), CM_SRF(3,NFMAX)
    }
    # fmt: on

    ad_suffix = "_DIFF"

    def __init__(self, geo_file=None, mass_file=None, debug=False, timing=False):

        if timing:
            start_time = time.time()

        # MExt is important for creating multiple instances of the AVL solver that do not share memory
        # It is very gross, but I cannot figure out a better way (maybe use install_name_tool to change the dynamic library path to absolute).
        # increment this counter for the hours you wasted on trying find a better way
        # 7 hours

        module_dir = os.path.dirname(os.path.realpath(__file__))
        module_name = os.path.basename(module_dir)
        avl_lib_so_file = glob.glob(os.path.join(module_dir, "libavl*.so"))[0]
        # # get just the file name
        avl_lib_so_file = os.path.basename(avl_lib_so_file)
        self.avl = MExt.MExt("libavl", module_name, "pyavl_wrapper", lib_so_file=avl_lib_so_file, debug=debug)._module

        # this way doesn't work with mulitple isntances fo AVLSolver
        # from . import libavl
        # self.avl = libavl

        if not (geo_file is None):
            try:
                # check to make sure files exist
                file = geo_file
                f = open(geo_file, "r")
                f.close()

                if not (mass_file is None):
                    file = mass_file
                    f = open(mass_file, "r")
                    f.close()
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Could not open the file '{file}' from python. This is usually an issue with the specified file path"
                )

            self.avl.avl()
            if debug:
                self.set_avl_fort_arr("CASE_L", "LVERBOSE", True)

            if timing:
                self.set_avl_fort_arr("CASE_L", "LTIMING", True)

            self.avl.loadgeo(geo_file)

            if mass_file is not None:
                self.avl.loadmass(mass_file)

        else:
            raise ValueError("neither a geometry file or aircraft object was given")

        self.bref = self.get_avl_fort_arr("CASE_R", "BREF")
        self.cref = self.get_avl_fort_arr("CASE_R", "CREF")
        self.sref = self.get_avl_fort_arr("CASE_R", "SREF")

        control_names = self.get_control_names()
        self.control_variables = {}
        for idx_c_var, c_name in enumerate(control_names):
            self.control_variables[c_name] = f"D{idx_c_var+1}"

        #  the case parameters are stored in a 1d array,
        # these indices correspond to the position of each parameter in that arra
        self._init_surf_data()

        if timing:
            print(f"AVL init took {time.time() - start_time} seconds")

    def _init_surf_data(self):
        self.surf_geom_to_fort_var = {}
        surf_names = self.get_surface_names()

        self.surf_pannel_to_fort_var = {}
        self.con_surf_to_fort_var = {}

        for idx_surf, surf_name in enumerate(surf_names):
            idx_surf = surf_names.index(surf_name)
            slice_idx_surf = (idx_surf,)
            slice_surf_all = (idx_surf, slice(None))

            # only set unduplicated sufaces
            if self.get_avl_fort_arr("SURF_I", "IMAGS", slicer=slice_idx_surf) < 0:
                # this is a duplicated surface, skip it
                continue

            num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC", slicer=slice_idx_surf)

            slice_surf_secs = (idx_surf, slice(None, num_sec))
            slice_surf_secs_all = slice_surf_secs + (slice(None),)

            nasec = self.get_avl_fort_arr("SURF_GEOM_I", "NASEC", slicer=slice_surf_secs)
            # I don't see a case in the code where nasec is not the same for all sections
            # but just in case, we'll do a test and throw an error if not
            if np.unique(nasec).size != 1:
                raise RuntimeError("nasec is not the same for all sections")

            nasec = nasec[0]

            slice_surf_secs_nasec = slice_surf_secs + (slice(None, nasec),)

            self.surf_geom_to_fort_var[surf_name] = {
                "scale": ["SURF_GEOM_R", "XYZSCAL", slice_surf_all],
                "translate": ["SURF_GEOM_R", "XYZTRAN", slice_surf_all],
                "angle": ["SURF_GEOM_R", "ADDINC", slice_idx_surf],
                "xyzles": ["SURF_GEOM_R", "XYZLES", slice_surf_secs_all],
                "chords": ["SURF_GEOM_R", "CHORDS", slice_surf_secs],
                "aincs": ["SURF_GEOM_R", "AINCS", slice_surf_secs],
                "xasec": ["SURF_GEOM_R", "XASEC", slice_surf_secs_nasec],
                "sasec": ["SURF_GEOM_R", "SASEC", slice_surf_secs_nasec],
                "tasec": ["SURF_GEOM_R", "TASEC", slice_surf_secs_nasec],
                "clcdsec": ["SURF_GEOM_R", "CLCDSEC", slice_surf_secs_all],
                "claf": ["SURF_GEOM_R", "CLAF", slice_surf_secs],
            }
            self.surf_pannel_to_fort_var[surf_name] = {
                "nchordwise": ["SURF_GEOM_I", "NVC", slice_idx_surf],
                "cspace": ["SURF_GEOM_R", "CSPACE", slice_idx_surf],
                "nspan": ["SURF_GEOM_I", "NVS", slice_idx_surf],
                "sspace": ["SURF_GEOM_R", "SSPACE", slice_idx_surf],
                "sspaces": ["SURF_GEOM_R", "SSPACES", slice_surf_secs],
                "nspans": ["SURF_GEOM_I", "NSPANS", slice_surf_secs],
                "yduplicate": ["SURF_GEOM_R", "YDUPL", slice_idx_surf],
            }

            icontd_slices = []
            idestd_slices = []
            xhinged_slices = []
            vhinged_slices = []
            gaind_slices = []
            refld_slices = []
            gaing_slices = []

            for idx_sec in range(num_sec):
                slice_surf = (idx_surf, idx_sec)

                num_con_surf = self.get_avl_fort_arr("SURF_GEOM_I", "NSCON", slicer=slice_surf)
                slice_surf_con = slice_surf + (slice(None, num_con_surf),)
                slice_surf_con_hinge = slice_surf_con + (slice(None),)

                num_des_var = self.get_avl_fort_arr("SURF_GEOM_I", "NSDES", slicer=slice_surf)
                slice_surf_des_var = slice_surf + (slice(None, num_des_var),)

                icontd_slices.append(slice_surf_con)
                xhinged_slices.append(slice_surf_con)
                vhinged_slices.append(slice_surf_con_hinge)
                gaind_slices.append(slice_surf_con)
                refld_slices.append(slice_surf_con)
                idestd_slices.append(slice_surf_des_var)
                gaing_slices.append(slice_surf_des_var)

            self.con_surf_to_fort_var[surf_name] = {
                "icontd": ["SURF_GEOM_I", "ICONTD", icontd_slices],
                "xhinged": ["SURF_GEOM_R", "XHINGED", xhinged_slices],
                "vhinged": ["SURF_GEOM_R", "VHINGED", vhinged_slices],
                "gaind": ["SURF_GEOM_R", "GAIND", gaind_slices],
                "refld": ["SURF_GEOM_R", "REFLD", refld_slices],
                "idestd": ["SURF_GEOM_I", "IDESTD", idestd_slices],
                "gaing": ["SURF_GEOM_R", "GAING", gaing_slices],
            }

    def add_constraint(self, var, val, con_var=None):
        avl_variables = {
            "alpha": "A",
            "beta": "B",
            "roll rate": "R",
            "pitch rate": "P",
            "yaw rate": "Y",
        }

        avl_con_variables = copy.deepcopy(avl_variables)
        avl_con_variables.update(
            {
                "CL": "C",
                "CY": "S",
                "Cl roll moment": "RM",
                "Cm pitch moment": "PM",
                "Cn yaw moment": "YM",
            }
        )

        if var in avl_variables:
            # save the name of the avl_var
            avl_var = avl_variables[var]
        elif var in self.control_variables.keys():
            avl_var = self.control_variables[var]
        elif var in self.control_variables.values():
            avl_var = var
        else:
            raise ValueError(
                f"specified variable `{var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.control_variables.items()]}. Constraints that must be implicitly satisfied (such as `CL`) are set with `add_trim_constraint`."
            )

        if con_var is None:
            avl_con_var = avl_var
        elif con_var in avl_con_variables:
            avl_con_var = avl_con_variables[con_var]
        elif con_var in self.control_variables.keys():
            avl_con_var = self.control_variables[con_var]
        elif con_var in self.control_variables.values():
            avl_con_var = con_var
        else:
            raise ValueError(
                f"specified contraint variable `{con_var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.control_variables.items()]}."
            )

        # check that the type of val is correct
        if not isinstance(val, (int, float, np.floating, np.integer)):
            raise TypeError(f"contraint value must be a int or float for contraint {var}. Got {type(val)}")

        self.avl.conset(avl_var, f"{avl_con_var} {val} \n")

    def add_trim_condition(self, variable, val):
        options = {
            "bankAng": ["B"],
            "CL": ["C"],
            "velocity": ["V"],
            "mass": ["M"],
            "dens": ["D"],
            "G": ["G"],
            "X cg": ["X"],
            "Y cg": ["Y"],
            "Z cg": ["Z"],
        }

        if not (variable in options):
            raise ValueError(
                f"constraint variable `{variable}` not a valid option. Must be one of the following {[key for key in options]} "
            )

        self.avl.trmset("C1", "1 ", options[variable][0], (str(val) + "  \n"))

    def get_case_total_data(self) -> Dict[str, float]:
        """Get the aerodynamic data for the last run case and return it as a dictionary.

        Returns:
            Dict[str, float]: Dictionary of aerodynamic data. The keys the aerodyanmic coefficients.
        """

        total_data = {}

        for key, avl_key in self.case_var_to_fort_var.items():
            val = self.get_avl_fort_arr(*avl_key)
            # [()] because all the data is stored as a ndarray.
            # for scalars this results in a 0-d array.
            # It is easier to work with floats so we extract the value with [()]
            total_data[key] = val[()]

        return total_data

    def get_case_coef_derivs(self) -> Dict[str, Dict[str, float]]:
        deriv_data = {}

        control_names = self.get_control_names()

        for key, avl_key in self.case_derivs_to_fort_var.items():
            slicer = slice(0, len(control_names))
            val_arr = self.get_avl_fort_arr(*avl_key, slicer=slicer)
            deriv_data[key] = {}
            for idx_control, val in enumerate(val_arr):
                control = control_names[idx_control]
                deriv_data[key][control] = val[()]

        return deriv_data

    def get_case_stab_derivs(self) -> Dict[str, Dict[str, float]]:
        deriv_data = {}

        for func_key, var_dict in self.case_stab_derivs_to_fort_var.items():
            deriv_data[func_key] = {}

            for var_key, avl_key in var_dict.items():
                val_arr = self.get_avl_fort_arr(*avl_key)
                deriv_data[func_key][var_key] = val_arr[()]

        return deriv_data

    def get_avl_fort_arr(self, common_block, variable, slicer=None):
        # this had to be split up into two steps to work

        # get the corresponding common block object.
        # it must be lowercase because of f2py
        common_block = getattr(self.avl, common_block.upper())

        # get the value of the variable from the common block
        val = getattr(common_block, variable.upper())

        # convert from fortran ordering to c ordering
        val = val.ravel(order="F").reshape(val.shape[::-1], order="C")

        # Apply slicer if provided
        if slicer is not None:
            val = val[slicer]

        return val

    def set_avl_fort_arr(self, common_block, variable, val, slicer=None) -> None:
        # convert from fortran ordering to c ordering
        if isinstance(val, np.ndarray):
            val = val.ravel(order="C").reshape(val.shape[::-1], order="F")

        # this had to be split up into two steps to work
        # get the corresponding common block object.
        # it must be lowercase because of f2py
        common_block_obj = getattr(self.avl, common_block.upper())

        # get the value of the variable from the common block
        if slicer is None:
            setattr(common_block_obj, variable.upper(), val)
        else:
            # flip the order of the slicer to match the cordinates of the val
            new_slicer = slicer[::-1]

            original_val = getattr(common_block_obj, variable.upper())
            original_val[new_slicer] = val
            setattr(common_block_obj, variable.upper(), original_val)

        return

    def get_case_surface_data(self) -> Dict[str, Dict[str, float]]:
        surf_names = self.get_surface_names()

        # add a dictionary for each surface that will be filled later
        surf_data = {}
        for surf in surf_names:
            surf_data[surf] = {}

        for key, avl_key in self.case_surf_var_to_fort_var.items():
            vals = self.get_avl_fort_arr(*avl_key)

            # add the values to corresponding surface dict
            for idx_surf, surf_name in enumerate(surf_names):
                surf_data[surf_name][key] = vals[idx_surf]

        return surf_data

    def get_case_parameter(self, param_key: str) -> float:
        """
        analogous to ruinont Modify parameters for the oper menu to view parameters.
        """
        parvals = self.get_avl_fort_arr("CASE_R", "PARVAL")

        # [0] because pyavl only supports 1 run case
        param_val = parvals[0][self.param_idx_dict[param_key]]

        return param_val

    def get_case_constraint(self, con_key: str) -> float:
        """ """
        convals = self.get_avl_fort_arr("CASE_R", "CONVAL")

        # [0] because pyavl only supports 1 run case
        con_val = convals[0][self.conval_idx_dict[con_key]]

        return con_val

    def set_case_parameter(self, param_key: str, param_val: float) -> None:
        """
        analogous to ruinont Modify parameters for the oper menu to view parameters.
        """
        # warn the user that alpha, beta,
        if param_key is ["alpha", "beta", "pb/2V", "qc/2V", "rb/2V", "CL"]:
            raise ValueError(
                "alpha, beta, pb/2V, qc/2V, rb/2V, and CL are not allowed to be set,\n\
                             they are calculated during each run based on the contraints. to specify\n\
                             one of these values use the add_contraint method."
            )

        parvals = self.get_avl_fort_arr("CASE_R", "PARVAL")
        # [0] because pyavl only supports 1 run case
        parvals[0][self.param_idx_dict[param_key]] = param_val

        self.set_avl_fort_arr("CASE_R", "PARVAL", parvals)

        # (1) here because we want to set the first runcase with fortran indexing (the only one)
        self.avl.set_params(1)

    def get_control_deflections(self) -> Dict[str, float]:
        control_surfaces = self.get_control_names()

        def_arr = copy.deepcopy(self.get_avl_fort_arr("CASE_R", "DELCON"))

        def_dict = {}
        for idx_con, con_surf in enumerate(control_surfaces):
            def_dict[con_surf] = def_arr[idx_con]

        return def_dict

    def get_hinge_moments(self) -> Dict[str, float]:
        """
        get the hinge moments from the fortran layer and return them as a dictionary
        """
        hinge_moments = {}

        control_surfaces = self.get_control_names()
        mom_array = self.get_avl_fort_arr("CASE_R", "CHINGE")

        for idx_con, con_surf in enumerate(control_surfaces):
            hinge_moments[con_surf] = mom_array[idx_con]

        return hinge_moments

    def get_hinge_moments(self) -> Dict[str, float]:
        """
        get the hinge moments from the fortran layer and return them as a dictionary
        """
        hinge_moments = {}

        control_surfaces = self.get_control_names()
        mom_array = self.get_avl_fort_arr("CASE_R", "CHINGE")

        for idx_con, con_surf in enumerate(control_surfaces):
            hinge_moments[con_surf] = mom_array[idx_con]

        return hinge_moments

    def get_strip_data(self) -> Dict[str, Dict[str, np.ndarray]]:
        # fmt: off
        var_to_fort_var = {
            # geometric quantities
            "chord": ["STRP_R", "CHORD"],
            "width": ["STRP_R", "WSTRIP"],
            "XYZ LE": ["STRP_R", "RLE"],  # control point leading edge coordinates
            "twist": ["STRP_R", "AINC"],
            
            # strip contributions to total lift and drag from strip integration
            "CL": ["STRP_R", "CLSTRP"],
            "CD": ["STRP_R", "CDSTRP"],
            
            # strip contributions to non-dimensionalized forces 
            "CX": ["STRP_R", "CXSTRP"], 
            "CY": ["STRP_R", "CYSTRP"], 
            "CZ": ["STRP_R", "CZSTRP"],   
            
            
            # strip contributions to total moments (body frame)
            "CM": ["STRP_R", "CMSTRP"],
            "CN": ["STRP_R", "CNSTRP"],
            "CR": ["STRP_R", "CRSTRP"],
            
            
            # forces non-dimentionalized by strip quantities
            "CL strip" : ["STRP_R", "CL_LSTRP"],
            "CD strip" : ["STRP_R", "CD_LSTRP"],
            "CF strip" : ["STRP_R", "CF_STRP"], # forces in 3 directions
            "CM strip" : ["STRP_R", "CF_STRP"], # moments in 3 directions
        }    
        # fmt: on
        surf_names = self.get_surface_names()

        # add a dictionary for each surface that will be filled later
        strip_data = {}
        for surf in surf_names:
            strip_data[surf] = {}

        for key, avl_key in var_to_fort_var.items():
            vals = self.get_avl_fort_arr(*avl_key)

            # add the values to corresponding surface dict
            for idx_surf, surf_name in enumerate(surf_names):
                idx_srp_beg, idx_srp_end = self.get_surface_strip_indices(idx_surf)
                strip_data[surf_name][key] = vals[idx_srp_beg:idx_srp_end]

        cref = self.get_avl_fort_arr("CASE_R", "CREF")
        for surf_key in strip_data:
            # add sectional lift and drag
            strip_data[surf_key]["lift dist"] = strip_data[surf_key]["CL"] * strip_data[surf_key]["chord"] / cref
            strip_data[surf_key]["drag dist"] = strip_data[surf_key]["CD"] * strip_data[surf_key]["chord"] / cref

        return strip_data

    def get_surface_strip_indices(self, idx_surf):
        num_strips = np.trim_zeros(self.get_avl_fort_arr("SURF_I", "NJ"))
        idx_srp_beg = np.sum(num_strips[:idx_surf])
        idx_srp_end = np.sum(num_strips[: idx_surf + 1])

        return idx_srp_beg, idx_srp_end

    def executeRun(self):
        warnings.warn("executeRun is deprecated, use execute_run instead")
        self.execute_run()

    def execute_run(self):
        # run the analysis (equivalent to the avl command `x` in the oper menu)
        self.avl.oper()

    def CLSweep(self, start_CL, end_CL, increment=0.1):
        CLs = np.arange(start_CL, end_CL + increment, increment)

        for cl in CLs:
            self.add_trim_condition("CL", cl)
            self.execute_run()

    def get_control_names(self) -> List[str]:
        fort_names = self.get_avl_fort_arr("CASE_C", "DNAME")
        control_names = self._convertFortranStringArrayToList(fort_names)
        return control_names

    def get_surface_names(self, remove_dublicated=False) -> List[str]:
        """get the surface names from the geometry"""
        fort_names = self.get_avl_fort_arr("CASE_C", "STITLE")
        surf_names = self._convertFortranStringArrayToList(fort_names)

        if remove_dublicated:
            imags = self.get_avl_fort_arr("SURF_I", "IMAGS")
            unique_surf_names = []
            for idx_surf, surf_name in enumerate(surf_names):
                # get surfaces that have not been duplicated
                if imags[idx_surf] > 0:
                    unique_surf_names.append(surf_names[idx_surf])

            return unique_surf_names
        else:
            return surf_names

    def get_con_surf_param(self, surf_name, idx_slice, param):
        # the control surface and design variables need to be handeled differently because the number at each section is variable
        if param in self.con_surf_to_fort_var[surf_name].keys():
            fort_var = self.con_surf_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has control surface data {self.con_surf_to_fort_var[surf_name].keys()}"
            )
        param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2][idx_slice])
        return param

    def set_con_surf_param(self, surf_name, idx_slice, param, val):
        # the control surface and design variables need to be handeled differently because the number at each section is variable
        if param in self.con_surf_to_fort_var[surf_name].keys():
            fort_var = self.con_surf_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has control surface data {self.con_surf_to_fort_var[surf_name].keys()}"
            )
        # param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2][idx_slice])
        self.set_avl_fort_arr(fort_var[0], fort_var[1], val, slicer=fort_var[2][idx_slice])

    def get_surface_param(self, surf_name, param):
        # check that param is in self.surf_geom_to_fort_var
        if param in self.surf_geom_to_fort_var[surf_name].keys():
            fort_var = self.surf_geom_to_fort_var[surf_name][param]
        elif param in self.surf_pannel_to_fort_var[surf_name].keys():
            fort_var = self.surf_pannel_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has geom data {list(self.surf_geom_to_fort_var[surf_name].keys()) + list(self.surf_pannel_to_fort_var[surf_name].keys())}"
            )

        param = self.get_avl_fort_arr(fort_var[0], fort_var[1], slicer=fort_var[2])
        return param

    def set_surface_param(self, surf_name, param, val):
        if param in self.surf_geom_to_fort_var[surf_name].keys():
            fort_var = self.surf_geom_to_fort_var[surf_name][param]
        elif param in self.surf_pannel_to_fort_var[surf_name].keys():
            fort_var = self.surf_pannel_to_fort_var[surf_name][param]
        elif param in self.con_surf_to_fort_var[surf_name].keys():
            # the control surface and design variables need to be handeled differently because the number at each section is variable
            pass

        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has geom data {list(self.surf_geom_to_fort_var[surf_name].keys()) + list(self.surf_pannel_to_fort_var[surf_name].keys())}"
            )

        self.set_avl_fort_arr(fort_var[0], fort_var[1], val, slicer=fort_var[2])

    def get_surface_params(
        self,
        include_geom: bool = True,
        include_panneling: bool = False,
        include_con_surf: bool = False,
        include_airfoils: bool = False,
    ) -> Dict[str, Dict[str, Any]]:
        """get the surface level parameters from the geometry

        geom_only: only return the geometry parameters of the surface
            - xle, yle, zle, chord, twist
        """

        surf_names = self.get_surface_names()
        unique_surf_names = self.get_surface_names(remove_dublicated=True)
        surf_data = {}

        for surf_name in unique_surf_names:
            surf_data[surf_name] = {}
            if include_geom:
                for var in self.surf_geom_to_fort_var[surf_name]:
                    surf_data[surf_name][var] = self.get_surface_param(surf_name, var)

            idx_surf = surf_names.index(surf_name)
            if include_panneling:
                # add panneling parameters if requested
                for var in self.surf_pannel_to_fort_var[surf_name]:
                    surf_data[surf_name][var] = self.get_surface_param(surf_name, var)

                if not self.get_avl_fort_arr("SURF_GEOM_L", "LDUPL")[idx_surf]:
                    surf_data[surf_name].pop("yduplicate")

            if include_con_surf:
                # add control surface parameters if requested
                for var in self.con_surf_to_fort_var[surf_name]:
                    num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]

                    slice_data = []
                    for idx_sec in range(num_sec):
                        tmp = self.get_con_surf_param(surf_name, idx_sec, var)
                        slice_data.append(tmp)

                    surf_data[surf_name][var] = slice_data

            if include_airfoils:
                afiles = []
                num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]

                for idx_sec in range(num_sec):
                    afile = self.__decodeFortranString(self.avl.CASE_C.AFILES[idx_sec, idx_surf])
                    afiles.append(afile)

                surf_data[surf_name]["afiles"] = afiles

        return surf_data

    def set_surface_params(self, surf_data: Dict[str, Dict[str, any]]) -> None:
        """set the give surface data into the current avl object.
        ASSUMES THE CONTROL SURFACE DATA STAYS AT THE SAME LOCATION"""
        surf_names = self.get_surface_names()
        unique_surf_names = self.get_surface_names(remove_dublicated=True)

        for surf_name in surf_data:
            if surf_name not in unique_surf_names:
                raise ValueError(
                    f"surface name, {surf_name}, not found in the current avl object. Note duplicated surfaces can not be set directly"
                )

            for var in surf_data[surf_name]:
                # do not set the data this way if it is a control surface
                if var not in self.con_surf_to_fort_var[surf_name]:
                    self.set_surface_param(surf_name, var, surf_data[surf_name][var])
                else:
                    idx_surf = surf_names.index(surf_name)
                    num_sec = self.get_avl_fort_arr("SURF_GEOM_I", "NSEC")[idx_surf]
                    slice_data = []
                    for idx_sec in range(num_sec):
                        self.set_con_surf_param(surf_name, idx_sec, var, surf_data[surf_name][var][idx_sec])

        self.avl.update_surfaces()

    def write_geom_file(self, filename: str):
        """write the current avl geometry to a file"""
        with open(filename, "w") as fid:
            # write the header
            fid.write("# generated using pyAVL\n")
            self.__write_header(fid)

            surf_data = self.get_surface_params(
                include_geom=True, include_panneling=True, include_con_surf=True, include_airfoils=True
            )
            for surf_name in surf_data:
                self.__write_surface(fid, surf_name, surf_data[surf_name])

    def __write_fort_vars(self, fid, common_block: str, fort_var: str, newline: bool = True) -> None:
        var = self.get_avl_fort_arr(common_block, fort_var)

        out_str = ""
        # loop over the variables list and recursively convert the variables to a string and add to the output string
        if isinstance(var, np.ndarray):
            if var.size == 1:
                out_str += str(var[()])
            else:
                out_str += " ".join([str(item) for item in var])
        else:
            out_str += str(var)
        out_str += " "

        if newline:
            out_str += "\n"

        fid.write(out_str)

    def __write_header(self, fid):
        """write the header to a file"""
        # write the name of the aircraft

        self.__write_banner(fid, "Header")
        title_array = self.get_avl_fort_arr("case_c", "title")
        title = self._convertFortranStringArrayToList(title_array)
        fid.write(f"{title}\n")

        fid.write("#Mach\n")
        self.__write_fort_vars(fid, "case_r", "mach0")

        fid.write("#IYsym   IZsym   Zsym\n")
        self.__write_fort_vars(fid, "case_i", "iysym", newline=False)
        self.__write_fort_vars(fid, "case_i", "izsym", newline=False)
        self.__write_fort_vars(fid, "case_r", "zsym")

        fid.write("#Sref    Cref    Bref\n")
        self.__write_fort_vars(fid, "case_r", "sref", newline=False)
        self.__write_fort_vars(fid, "case_r", "cref", newline=False)
        self.__write_fort_vars(fid, "case_r", "bref")

        fid.write("#Xref    Yref    Zref\n")
        self.__write_fort_vars(fid, "case_r", "XYZREF")

        fid.write("#CD0\n")
        self.__write_fort_vars(fid, "case_r", "CDREF0")

        # fid.write(f" {self.get_avl_fort_arr('case_r', 'sref')}")

    def __write_banner(self, fid, header, line_width: int = 80):
        header = " " + header + " "  # pad with spaces

        width = line_width - 1

        banner = f"#{'='*(width)}\n" f"#{header.center(width,'-')}\n" f"#{'='*(width)}\n"
        fid.write(banner)

        # ======================================================
        # ------------------- Geometry File --------------------
        # ======================================================

    def __write_surface(self, fid, surf_name, data):
        """write a surface to a file"""
        # TODO add NACA and CLAF keyword support

        def __write_data(key_list, newline: bool = True):
            out_str = ""
            for key in key_list:
                val = data[key]

                if isinstance(val, np.ndarray):
                    if val.size == 1:
                        out_str += str(val[()])
                    else:
                        out_str += " ".join([str(item) for item in val])
                else:
                    out_str += str(val)
                out_str += " "

            if newline:
                out_str += "\n"

            fid.write(out_str)

        # start with the banner
        self.__write_banner(fid, surf_name)
        fid.write(f"SURFACE\n")
        fid.write(f"{surf_name}\n")

        fid.write(f"#Nchordwise  Cspace  [Nspanwise  Sspace]\n")
        __write_data(["nchordwise", "cspace", "nspan", "sspace"])

        if "yduplicate" in data:
            fid.write("YDUPLICATE\n")
            __write_data(["yduplicate"])

        fid.write("SCALE\n")
        __write_data(["scale"])

        fid.write("TRANSLATE\n")
        __write_data(["translate"])

        fid.write("ANGLE\n")
        __write_data(["angle"])

        fid.write("#---------------------------------------\n")

        num_sec = data["chords"].size
        control_names = self.get_control_names()

        for idx_sec in range(num_sec):
            fid.write("SECTION\n")
            fid.write("#Xle      Yle      Zle      | Chord    Ainc     Nspan  Sspace\n")
            fid.write(
                f" {data['xyzles'][idx_sec, 0]:.6f} "
                f"{data['xyzles'][idx_sec, 1]:.6f} "
                f"{data['xyzles'][idx_sec, 2]:.6f}   "
                f"{data['chords'][idx_sec]:.6f} "
                f"{data['aincs'][idx_sec]:.6f} "
                f"{data['nspans'][idx_sec]}      "
                f"{data['sspaces'][idx_sec]}\n"
            )

            afile = data["afiles"][idx_sec]

            if afile:
                fid.write(" AFILE\n")
                fid.write(f" {afile}\n")

            # check for control surfaces

            for idx_local_cont_surf, idx_cont_surf in enumerate(data["icontd"][idx_sec]):
                fid.write(" CONTROL\n")
                fid.write("#surface   gain xhinge       hvec       SgnDup\n")
                fid.write(f" {control_names[idx_cont_surf-1]} ")
                fid.write(f" {data['gaind'][idx_sec][idx_local_cont_surf]}")
                fid.write(f" {data['xhinged'][idx_sec][idx_local_cont_surf]}")
                vhinge = data["vhinged"][idx_sec][idx_local_cont_surf]
                fid.write(f" {vhinge[0]:.6f} {vhinge[1]:.6f} {vhinge[2]:.6f}")
                fid.write(f" {data['refld'][idx_sec][idx_local_cont_surf]}\n")

    def __decodeFortranString(self, fort_string) -> str:
        # TODO: need a more general solution for |S<variable> type

        if fort_string.dtype == np.dtype("|S0"):
            # there are no characters in the sting to add
            return ""
        if fort_string.dtype == np.dtype("|S1"):
            py_string = b"".join(fort_string).decode().strip()
        elif fort_string.dtype == np.dtype("<U1"):
            py_string = "".join(fort_string).strip()
        elif fort_string.dtype == np.dtype("|S16"):
            py_string = fort_string.decode().strip()
        elif fort_string.dtype == np.dtype("|S40"):
            py_string = fort_string.decode().strip()
        elif fort_string.dtype == np.dtype("|S80"):
            py_string = fort_string.decode().strip()
        else:
            raise TypeError(f"Unable to convert {fort_string} of type {fort_string.dtype} to string")

        return py_string

    # Utility functions
    def get_num_surfaces(self) -> int:
        """Get the number of surfaces in the geometry"""
        return self.get_avl_fort_arr("CASE_I", "NSURF")

    def get_num_sections(self, surf_name) -> int:
        """Get the number of surfaces in the geometry"""
        surf_names = self.get_surface_names()
        idx_surf = surf_names.index(surf_name)
        slice_idx_surf = (idx_surf,)
        return self.get_avl_fort_arr("SURF_GEOM_I", "NSEC", slicer=slice_idx_surf)

    def get_num_strips(self) -> int:
        """Get the number of strips in the mesh"""
        return self.get_avl_fort_arr("CASE_I", "NSTRIP")

    def get_num_control_surfs(self) -> int:
        return self.get_avl_fort_arr("CASE_I", "NCONTROL")

    def get_mesh_size(self) -> int:
        """Get the number of panels in the mesh"""
        return int(self.get_avl_fort_arr("CASE_I", "NVOR"))

    def _createFortranStringArray(self, strList, num_max_char):
        """Setting arrays of strings in Fortran can be kinda nasty. This
        takesa list of strings and returns the array"""

        arr = np.zeros((len(strList), num_max_char), dtype="str")
        arr[:] = " "
        for i, s in enumerate(strList):
            for j in range(len(s)):
                arr[i, j] = s[j]

        return arr

    def _convertFortranStringArrayToList(self, fortArray):
        """Undoes the _createFotranStringArray"""
        strList = []

        if fortArray.size == 1:
            # we must handle the 0-d array case sperately
            return self.__decodeFortranString(fortArray[()])

        for ii in range(fortArray.size):
            py_string = self.__decodeFortranString(fortArray[ii])
            if py_string != "":
                strList.append(py_string)

        return strList

    # Derivative routines
    def get_constraint_ad_seeds(self) -> Dict[str, float]:
        con_seeds = {}
        for con in self.con_var_to_fort_var:
            idx_con = self.conval_idx_dict[con]
            blk = "CASE_R" + self.ad_suffix
            var = "CONVAL" + self.ad_suffix
            slicer = (0, idx_con)

            fort_arr = self.get_avl_fort_arr(blk, var, slicer=slicer)
            con_seeds[con] = copy.deepcopy(fort_arr)

        return con_seeds

    def set_constraint_ad_seeds(self, con_seeds: Dict[str, Dict[str, float]], mode: str = "AD", scale=1.0) -> None:
        for con in con_seeds:
            # determine the proper index

            idx_con = self.conval_idx_dict[con]
            # determine the proper index

            con_seed_arr = con_seeds[con]

            if mode == "AD":
                # [0] because pyavl only supports 1 run case
                blk = "CASE_R" + self.ad_suffix
                var = "CONVAL" + self.ad_suffix
                val = con_seed_arr * scale
                slicer = (0, idx_con)

                self.set_avl_fort_arr(blk, var, val, slicer=slicer)

            elif mode == "FD":
                # reverse lookup in the con_var_to_fort_var dict

                val = self.get_case_constraint(con)

                val += con_seed_arr * scale

                # use the contraint API to adjust the value
                self.add_constraint(con, val)

    def get_geom_ad_seeds(self) -> Dict[str, Dict[str, float]]:
        geom_seeds = {}
        for surf_key in self.surf_geom_to_fort_var:
            geom_seeds[surf_key] = {}
            for geom_key in self.surf_geom_to_fort_var[surf_key]:
                blk, var, slicer = self.surf_geom_to_fort_var[surf_key][geom_key]

                blk += self.ad_suffix
                var += self.ad_suffix

                geom_seeds[surf_key][geom_key] = copy.deepcopy(self.get_avl_fort_arr(blk, var, slicer=slicer))

        return geom_seeds

    def set_geom_ad_seeds(self, geom_seeds: Dict[str, float], mode: str = "AD", scale=1.0) -> None:
        for surf_key in geom_seeds:
            for geom_key in geom_seeds[surf_key]:
                blk, var, slicer = self.surf_geom_to_fort_var[surf_key][geom_key]

                if mode == "AD":
                    blk += self.ad_suffix
                    var += self.ad_suffix
                    val = geom_seeds[surf_key][geom_key] * scale
                elif mode == "FD":
                    val = self.get_avl_fort_arr(blk, var, slicer=slicer)
                    val += geom_seeds[surf_key][geom_key] * scale
                # print(blk, var, val, slicer)
                self.set_avl_fort_arr(blk, var, val, slicer=slicer)

    def get_gamma_ad_seeds(self) -> np.ndarray:
        slicer = (slice(0, self.get_mesh_size()),)
        blk = "VRTX_R_DIFF"
        var = "GAM_DIFF"

        gamma_seeds = copy.deepcopy(self.get_avl_fort_arr(blk, var, slicer=slicer))
        return gamma_seeds

    def get_gamma_d_ad_seeds(self) -> np.ndarray:
        slicer = (slice(0, self.get_num_control_surfs()), slice(0, self.get_mesh_size()))
        blk = "VRTX_R_DIFF"
        var = "GAM_d_DIFF"

        gamma_d_seeds = copy.deepcopy(self.get_avl_fort_arr(blk, var, slicer=slicer))
        return gamma_d_seeds

    def set_gamma_ad_seeds(self, gamma_seeds: np.ndarray, mode: str = "AD", scale=1.0) -> None:
        slicer = (slice(0, self.get_mesh_size()),)
        if mode == "AD":
            blk = "VRTX_R_DIFF"
            var = "GAM_DIFF"
            val = gamma_seeds * scale
        elif mode == "FD":
            blk = "VRTX_R"
            var = "GAM"
            val = self.get_avl_fort_arr(blk, var, slicer=slicer)
            val += gamma_seeds * scale

        self.set_avl_fort_arr(blk, var, val, slicer=slicer)

    def set_gamma_d_ad_seeds(self, gamma_d_seeds: np.ndarray, mode: str = "AD", scale=1.0) -> None:
        slicer = (slice(0, self.get_num_control_surfs()), slice(0, self.get_mesh_size()))
        if mode == "AD":
            blk = "VRTX_R_DIFF"
            var = "GAM_D_DIFF"
            val = gamma_d_seeds * scale
        elif mode == "FD":
            blk = "VRTX_R"
            var = "GAM_D"
            val = self.get_avl_fort_arr(blk, var, slicer=slicer)
            val += gamma_d_seeds * scale

        self.set_avl_fort_arr(blk, var, val, slicer=slicer)

    def get_function_ad_seeds(self):
        func_seeds = {}
        for _var in self.case_var_to_fort_var:
            blk, var = self.case_var_to_fort_var[_var]
            blk += self.ad_suffix
            var += self.ad_suffix
            val = self.get_avl_fort_arr(blk, var)
            func_seeds[_var] = copy.deepcopy(val)

        return func_seeds

    def set_function_ad_seeds(self, func_seeds: Dict[str, float], scale=1.0):
        for _var in func_seeds:
            blk, var = self.case_var_to_fort_var[_var]
            blk += self.ad_suffix
            var += self.ad_suffix
            val = func_seeds[_var] * scale
            self.set_avl_fort_arr(blk, var, val)

    def get_consurf_derivs_ad_seeds(self):
        cs_deriv_seeds = {}
        consurf_names = self.get_control_names()
        num_consurf = len(consurf_names)

        for _var in self.case_derivs_to_fort_var:
            blk, var = self.case_derivs_to_fort_var[_var]
            slicer = (slice(0, num_consurf),)
            blk += self.ad_suffix
            var += self.ad_suffix
            val_arr = self.get_avl_fort_arr(blk, var, slicer=slicer)
            cs_deriv_seeds[_var] = {}

            for idx_control, val in enumerate(val_arr):
                control = consurf_names[idx_control]
                cs_deriv_seeds[_var][control] = val[()]

        return cs_deriv_seeds

    def set_consurf_derivs_ad_seeds(self, cs_deriv_seeds: Dict[str, float], scale=1.0):
        consurf_names = self.get_control_names()
        num_consurf = len(consurf_names)

        for _var in cs_deriv_seeds:
            val_arr = np.zeros((num_consurf))
            for cs_name in cs_deriv_seeds[_var]:
                idx_cs = consurf_names.index(cs_name)
                val_arr[idx_cs] = cs_deriv_seeds[_var][cs_name] * scale

            blk, var = self.case_derivs_to_fort_var[_var]
            slicer = (slice(0, num_consurf),)

            blk += self.ad_suffix
            var += self.ad_suffix

            self.set_avl_fort_arr(blk, var, val_arr, slicer=slicer)

    def get_residual_ad_seeds(self) -> np.ndarray:
        res_slice = (slice(0, self.get_mesh_size()),)
        res_seeds = copy.deepcopy(self.get_avl_fort_arr("VRTX_R_DIFF", "RES_DIFF", slicer=res_slice))
        return res_seeds

    def set_residual_ad_seeds(self, res_seeds: np.ndarray, scale=1.0) -> None:
        res_slice = (slice(0, self.get_mesh_size()),)
        self.set_avl_fort_arr("VRTX_R_DIFF", "RES_DIFF", res_seeds * scale, slicer=res_slice)
        return

    def get_residual_d_ad_seeds(self) -> np.ndarray:
        res_slice = (slice(0, self.get_num_control_surfs()), slice(0, self.get_mesh_size()))
        res_seeds = copy.deepcopy(self.get_avl_fort_arr("VRTX_R_DIFF", "RES_D_DIFF", slicer=res_slice))
        return res_seeds

    def set_residual_d_ad_seeds(self, res_d_seeds: np.ndarray, scale=1.0) -> None:
        res_slice = (slice(0, self.get_num_control_surfs()), slice(0, self.get_mesh_size()))

        self.set_avl_fort_arr("VRTX_R_DIFF", "RES_D_DIFF", res_d_seeds * scale, slicer=res_slice)
        return

    def clear_ad_seeds(self):
        for att in dir(self.avl):
            if att.endswith(self.ad_suffix):
                # loop over the attributes of the common block
                diff_blk = getattr(self.avl, att)
                for _var in dir(diff_blk):
                    if not (_var.startswith("__") and _var.endswith("__")):
                        val = getattr(diff_blk, _var)
                        setattr(diff_blk, _var, val * 0.0)

    def clear_ad_seeds_fast(self):
        # Only clear the seeds that are used in Make_tapenade file
        num_vor = self.get_mesh_size()
        num_vor_max = 6000  # HACK: hardcoded value from AVL.inc

        for att in dir(self.avl):
            if att.endswith(self.ad_suffix):
                # loop over the attributes of the common block
                diff_blk = getattr(self.avl, att)
                for _var in dir(diff_blk):
                    if not (_var.startswith("__") and _var.endswith("__")):

                        val = getattr(diff_blk, _var)

                        # trim sizes set to NVMAX to NVOR
                        shape = val.shape
                        slices = []
                        for idx_dim in range(len(shape)):
                            dim_size = shape[idx_dim]
                            if dim_size == num_vor_max:
                                dim_size = num_vor

                            slices.append(slice(0, dim_size))
                        slicer = tuple(slices)
                        val[slicer] = 0.0
                        # setattr(diff_blk, _var, val)
                        # print(diff_blk, _var, val.shape, slicer)
                        # mb_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
                        # print(f"    Memory usage: {mb_memory} MB")


    def print_ad_seeds(self, print_non_zero: bool = False):
        for att in dir(self.avl):
            if att.endswith(self.ad_suffix):
                # loop over the attributes of the common block
                diff_blk = getattr(self.avl, att)
                for _var in dir(diff_blk):
                    if not (_var.startswith("__") and _var.endswith("__")):
                        val = getattr(diff_blk, _var)
                        norm = np.linalg.norm(val)
                        if norm == 0.0 and print_non_zero:
                            continue

                        print(att, _var, norm)

    def execute_jac_vec_prod_fwd(
        self,
        con_seeds: Optional[Dict[str, float]] = None,
        geom_seeds: Optional[Dict[str, Dict[str, any]]] = None,
        gamma_seeds: Optional[np.ndarray] = None,
        gamma_d_seeds: Optional[np.ndarray] = None,
        mode="AD",
        step=1e-7,
    ) -> Tuple[Dict[str, float], np.ndarray]:
        # TODO: add better name for gamma_d, it is too confusing

        # TODO: add error if data is not initailzed properly
        #   The easiest fix is to run an analysis before hand

        mesh_size = self.get_mesh_size()
        num_control_surfs = self.get_num_control_surfs()

        if con_seeds is None:
            con_seeds = {}

        if geom_seeds is None:
            geom_seeds = {}

        if gamma_seeds is None:
            gamma_seeds = np.zeros(mesh_size)

        if gamma_d_seeds is None:
            gamma_d_seeds = np.zeros((num_control_surfs, mesh_size))

        res_slice = (slice(0, mesh_size),)
        res_d_slice = (slice(0, num_control_surfs), slice(0, mesh_size))

        if mode == "AD":
            # set derivative seeds
            # self.clear_ad_seeds()
            self.set_constraint_ad_seeds(con_seeds)
            self.set_geom_ad_seeds(geom_seeds)
            self.set_gamma_ad_seeds(gamma_seeds)
            self.set_gamma_d_ad_seeds(gamma_d_seeds)

            self.avl.update_surfaces_d()
            self.avl.get_res_d()
            self.avl.velsum_d()
            self.avl.aero_d()

            # extract derivatives seeds and set the output dict of functions
            func_seeds = self.get_function_ad_seeds()
            res_seeds = self.get_residual_ad_seeds()
            consurf_derivs_seeds = self.get_consurf_derivs_ad_seeds()
            res_d_seeds = self.get_residual_d_ad_seeds()

            self.set_constraint_ad_seeds(con_seeds, scale=0.0)
            self.set_geom_ad_seeds(geom_seeds, scale=0.0)
            self.set_gamma_ad_seeds(gamma_seeds, scale=0.0)

            self.set_avl_fort_arr("VRTX_R_DIFF", "GAM_DIFF", gamma_seeds * 0.0, slicer=res_slice)

        if mode == "FD":
            self.set_constraint_ad_seeds(con_seeds, mode="FD", scale=step)
            self.set_geom_ad_seeds(geom_seeds, mode="FD", scale=step)
            self.set_gamma_ad_seeds(gamma_seeds, mode="FD", scale=step)
            self.set_gamma_d_ad_seeds(gamma_d_seeds, mode="FD", scale=step)

            # propogate the seeds through without resolving
            self.avl.update_surfaces()
            self.avl.get_res()
            self.avl.velsum()
            self.avl.aero()

            coef_data_peturb = self.get_case_total_data()
            consurf_derivs_petrub = self.get_case_coef_derivs()

            res_peturbed = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES", slicer=res_slice))
            res_d_peturbed = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES_D", slicer=res_d_slice))
            self.set_constraint_ad_seeds(con_seeds, mode="FD", scale=-1 * step)
            self.set_geom_ad_seeds(geom_seeds, mode="FD", scale=-1 * step)
            self.set_gamma_ad_seeds(gamma_seeds, mode="FD", scale=-1 * step)
            self.set_gamma_d_ad_seeds(gamma_d_seeds, mode="FD", scale=-1 * step)

            self.avl.update_surfaces()
            self.avl.get_res()
            self.avl.velsum()
            self.avl.aero()

            coef_data = self.get_case_total_data()
            consurf_derivs = self.get_case_coef_derivs()
            res = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES", slicer=res_slice))
            res_d = copy.deepcopy(self.get_avl_fort_arr("VRTX_R", "RES_D", slicer=res_d_slice))

            func_seeds = {}
            for func_key in coef_data:
                func_seeds[func_key] = (coef_data_peturb[func_key] - coef_data[func_key]) / step

            consurf_derivs_seeds = {}
            for func_key in consurf_derivs:
                consurf_derivs_seeds[func_key] = {}
                for surf_key in consurf_derivs[func_key]:
                    consurf_derivs_seeds[func_key][surf_key] = (
                        consurf_derivs_petrub[func_key][surf_key] - consurf_derivs[func_key][surf_key]
                    ) / step

            res_seeds = (res_peturbed - res) / step
            res_d_seeds = (res_d_peturbed - res_d) / step

        return func_seeds, res_seeds, consurf_derivs_seeds, res_d_seeds

    def execute_jac_vec_prod_rev(
        self,
        func_seeds: Optional[Dict[str, float]] = None,
        res_seeds: Optional[np.ndarray] = None,
        consurf_derivs_seeds: Optional[Dict[str, Dict[str, float]]] = None,
        res_d_seeds: Optional[np.ndarray] = None,
        print_timings=False,
    ) -> Tuple[Dict[str, float], Dict[str, Dict[str, any]], np.ndarray]:
        # extract derivatives seeds and set the output dict of functions
        mesh_size = self.get_mesh_size()
        num_surf = self.get_num_control_surfs()

        if func_seeds is None:
            func_seeds = {}

        if res_seeds is None:
            res_seeds = np.zeros(mesh_size)

        if res_d_seeds is None:
            res_d_seeds = np.zeros((num_surf, mesh_size))

        if consurf_derivs_seeds is None:
            consurf_derivs_seeds = {}

        # set derivative seeds
        # self.clear_ad_seeds()
        time_last = time.time()
        self.set_function_ad_seeds(func_seeds)
        self.set_residual_ad_seeds(res_seeds)
        self.set_residual_d_ad_seeds(res_d_seeds)
        self.set_consurf_derivs_ad_seeds(consurf_derivs_seeds)
        if print_timings:
            print(f"    Time to set seeds: {time.time() - time_last}")
            time_last = time.time()

        # propogate the seeds through without resolveing
        self.avl.aero_b()
        self.avl.velsum_b()
        self.avl.get_res_b()
        self.avl.update_surfaces_b()
        if print_timings:
            print(f"    Time to propogate seeds: {time.time() - time_last}")
            time_last = time.time()

        # extract derivatives seeds and set the output dict of functions
        con_seeds = self.get_constraint_ad_seeds()
        geom_seeds = self.get_geom_ad_seeds()
        gamma_seeds = self.get_gamma_ad_seeds()
        gamma_d_seeds = self.get_gamma_d_ad_seeds()
        if print_timings:
            print(f"    Time to extract seeds: {time.time() - time_last}")
            time_last = time.time()

        self.set_function_ad_seeds(func_seeds, scale=0.0)
        self.set_residual_ad_seeds(res_seeds, scale=0.0)
        self.set_residual_d_ad_seeds(res_d_seeds, scale=0.0)
        self.set_consurf_derivs_ad_seeds(consurf_derivs_seeds, scale=0.0)
        if print_timings:
            print(f"    Time to clear seeds: {time.time() - time_last}")
            time_last = time.time()

        return con_seeds, geom_seeds, gamma_seeds, gamma_d_seeds

    def scipy_solve():
        """
        use scipy's linear solves to solve the system of equations.
        """
        pass

    def execute_adjoint_solve():
        pass

    def execute_direct_solve():
        pass

    def execute_run_sensitivies(self, funcs, consurf_derivs=None, print_timings=False):

        """
        only runs in adjoint mode
        funcs: list of functions that need derivatives
        """
        sens = {}
        
        if self.get_avl_fort_arr("CASE_L", "LTIMING"):
            print_timings = True
        
        # set up and solve the adjoint for each function
        for func in funcs:
            sens[func] = {}
            # get the RHS of the adjiont equation (pFpU)
            # TODO: remove seeds if it doesn't effect accuracy
            # self.clear_ad_seeds()
            time_last = time.time()
            _, _, pfpU, pf_pU_d = self.execute_jac_vec_prod_rev(func_seeds={func: 1.0})
            if print_timings:
                print(f"Time to get RHS: {time.time() - time_last}")
                time_last = time.time()

            # self.clear_ad_seeds()
            # u solver adjoint equation with RHS
            self.set_gamma_ad_seeds(-1 * pfpU)
            self.set_gamma_d_ad_seeds(-1 * pf_pU_d)
            self.avl.solve_adjoint()
            if print_timings:
                print(f"Time to solve adjoint: {time.time() - time_last}")
                time_last = time.time()
            # get the resulting adjoint vector (dfunc/dRes) from fortran
            dfdR = self.get_residual_ad_seeds()
            dfdR_d = self.get_residual_d_ad_seeds()
            # self.clear_ad_seeds()
            con_seeds, geom_seeds, _, _ = self.execute_jac_vec_prod_rev(
                func_seeds={func: 1.0}, res_seeds=dfdR, res_d_seeds=dfdR_d
            )
            if print_timings:
                print(f"Time to combine derivs: {time.time() - time_last}")
                time_last = time.time()

            sens[func].update(con_seeds)
            sens[func].update(geom_seeds)

        if consurf_derivs is not None:
            if print_timings:
                print("Running consurf derivs")
                time_last = time.time()

            cs_deriv_seeds = {}
            for func_key in consurf_derivs:
                cs_deriv_seeds[func_key] = {}
                if func_key not in sens:
                    sens[func_key] = {}
                for cs_key in consurf_derivs[func_key]:
                    cs_deriv_seeds[func_key][cs_key] = 1.0
                    sens[func_key][cs_key] = {}

                    # get the RHS of the adjiont equation (pFpU)
                    # TODO: remove seeds if it doesn't effect accuracy
                    _, _, pfpU, pf_pU_d = self.execute_jac_vec_prod_rev(consurf_derivs_seeds=cs_deriv_seeds)
                    if print_timings:
                        print(f"Time to get RHS: {time.time() - time_last}")
                        time_last = time.time()

                    # self.clear_ad_seeds()
                    # u solver adjoint equation with RHS
                    self.set_gamma_ad_seeds(-1 * pfpU)
                    self.set_gamma_d_ad_seeds(-1 * pf_pU_d)
                    self.avl.solve_adjoint()
                    if print_timings:
                        print(f"Time to solve adjoint: {time.time() - time_last}")
                        time_last = time.time()

                    # get the resulting adjoint vector (dfunc/dRes) from fortran
                    dfdR = self.get_residual_ad_seeds()
                    dfdR_d = self.get_residual_d_ad_seeds()
                    # self.clear_ad_seeds()
                    con_seeds, geom_seeds, _, _ = self.execute_jac_vec_prod_rev(
                        consurf_derivs_seeds=cs_deriv_seeds, res_seeds=dfdR, res_d_seeds=dfdR_d
                    )
                    if print_timings:
                        print(f"Time to combine : {time.time() - time_last}")
                        time_last = time.time()

                    sens[func_key][cs_key].update(con_seeds)
                    sens[func_key][cs_key].update(geom_seeds)
                    cs_deriv_seeds[func_key][cs_key] = 0.0

        return sens
