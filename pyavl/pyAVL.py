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
        "pb/2V": 2,
        "qc/2V": 3,
        "rb/2V": 4,
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

    # fmt: off
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

    ad_suffix = "_diff"

    def __init__(self, geo_file=None, mass_file=None, debug=False, timing=False):
        # TODO: fix MExt to work with pyavl.lib packages
        #   - copy dependecies to the temp directory
        # This only needs to happen once for each fake directory

        # This is important for creating multiple instances of the AVL solver that do not share memory
        # It is very gross, but I cannot figure out a better way.
        # increment this counter for the hours you wasted on trying to remove this monkey patch
        # 4 hours

        module_dir = os.path.dirname(os.path.realpath(__file__))
        module_name = os.path.basename(module_dir)
        avl_lib_so_file = glob.glob(os.path.join(module_dir, "libavl*.so"))[0]
        # # get just the file name
        avl_lib_so_file = os.path.basename(avl_lib_so_file)
        self.avl = MExt.MExt("libavl", module_name, lib_so_file=avl_lib_so_file, debug=debug)._module

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
                    f"Could not open the file {file} from python. This is usually an issue with the specified file path"
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
                "xyzles": ["SURF_GEOM_R", "xyzles", slice_surf_secs_all],
                "chords": ["SURF_GEOM_R", "chords", slice_surf_secs],
                "aincs": ["SURF_GEOM_R", "aincs", slice_surf_secs],
                "xasec": ["SURF_GEOM_R", "xasec", slice_surf_secs_nasec],
                "sasec": ["SURF_GEOM_R", "sasec", slice_surf_secs_nasec],
                "tasec": ["SURF_GEOM_R", "tasec", slice_surf_secs_nasec],
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
                f"specified variable `{var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.control_variables.items()]}."
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
            raise TypeError(f"contraint `val` must be a int or float. Got {type(val)}")

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
        # fmt: off
        # This dict has the following structure:
        # python key: [common block name, fortran varaiable name]
        var_to_fort_var = {
            # lift and drag from surface integration (wind frame)
            "CL": ["CASE_R", "CLTOT"],
            "CD": ["CASE_R", "CDTOT"],
            "CDv": ["CASE_R", "CDVTOT"], # viscous drag
            
            # lift and drag calculated from farfield integration
            "CLff": ["CASE_R", "CLFF"],
            "CYff": ["CASE_R", "CDFF"],
            "CDi": ["CASE_R", "CDFF"], # viscous drag
            
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
        # fmt: on

        total_data = {}

        for key, avl_key in var_to_fort_var.items():
            val = self.get_avl_fort_arr(*avl_key)
            # [()] because all the data is stored as a ndarray.
            # for scalars this results in a 0-d array.
            # It is easier to work with floats so we extract the value with [()]
            total_data[key] = val[()]

        return total_data

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
        common_block = getattr(self.avl, common_block.upper())

        # get the value of the variable from the common block
        if slicer is None:
            setattr(common_block, variable.upper(), val)
        else:
            # flip the order of the slicer to match the cordinates of the val
            new_slicer = slicer[::-1]

            original_val = getattr(common_block, variable.upper())
            original_val[new_slicer] = val
            setattr(common_block, variable.upper(), original_val)

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

    def set_case_parameter(self, param_key: str, param_val: float) -> None:
        """
        analogous to ruinont Modify parameters for the oper menu to view parameters.
        """
        # warn the user that alpha, beta,
        if param_key is ["alpha", "beta", "pb/2V", "qc/2V", "rb/2V", "CL"]:
            raise ValueError(
                "alpha, beta, pb/2V, qc/2V, rb/2V, and CL are not allowed to be set,\n\
                             they are calculated during each run based on the contraints. to specify\n\
                             one of these values use the set_contrain method."
            )

        parvals = self.get_avl_fort_arr("CASE_R", "PARVAL")
        # [0] because pyavl only supports 1 run case
        parvals[0][self.param_idx_dict[param_key]] = param_val

        self.set_avl_fort_arr("CASE_R", "PARVAL", parvals)

        # (1) here because we want to set the first runcase with fortran indexing (the only one)
        self.avl.set_params(1)

    # def get_condition_data(self) -> Dict[str, float]:
    #     """get the data defining the run condition from the fortran layer"""
    #     # fmt: off
    #     var_to_fort_var = {
    #         "alpha": ["CASE_R", "ALPHA"],
    #         "beta": ["CASE_R", "BETA"],
    #         "Vinf": ["CASE_R", "VINF"],
    #         "mach": ["CASE_R", "MACH"],
    #     }
    #     # fmt: on
    #     condition_data = {}

    #     for key, avl_key in var_to_fort_var.items():
    #         val = self.get_avl_fort_arr(*avl_key)
    #         # [()] because all the data is stored as a ndarray.
    #         # for scalars this results in a 0-d array.
    #         # It is easier to work with floats so we extract the value with [()]
    #         condition_data[key] = val[()]

    #     return condition_data

    def get_strip_data(self) -> Dict[str, Dict[str, np.ndarray]]:
        # fmt: off
        var_to_fort_var = {
            # geometric quantities
            "chord": ["STRP_R", "CHORD"],
            "width": ["STRP_R", "WSTRIP"],
            "XYZ LE": ["STRP_R", "RLE"], #control point leading edge coordinates
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

        return strip_data

    def get_surface_strip_indices(self, idx_surf):
        num_strips = np.trim_zeros(self.avl.surf_i.nj)
        idx_srp_beg = np.sum(num_strips[:idx_surf])
        idx_srp_end = np.sum(num_strips[: idx_surf + 1])

        return idx_srp_beg, idx_srp_end

    def executeRun(self):
        warnings.warn("executeRun is deprecated, use execute_run instead")
        self.execute_run()

    def execute_run(self):
        # run the analysis (equivalent to the avl command `x` in the oper menu)
        self.avl.oper()

        # needed for stability derivaties
        self.avl.calcst()

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
        # print(fort_var[0], fort_var[1], fort_var[2][idx_slice])
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
        # print(fort_var[0], fort_var[1], fort_var[2][idx_slice])
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
        self, include_geom: bool = True, include_panneling: bool = False, include_con_surf: bool = False
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

        return surf_data

    def set_surface_params(self, surf_data: Dict[str, Dict[str, any]]) -> None:
        """set the give surface data into the current avl object.
        ASSUMES THE CONTROL SURFACE DATA STAYS AT THE SAME LOCATION"""
        surf_names = self.get_surface_names()
        unique_surf_names = self.get_surface_names(remove_dublicated=True)

        for surf_name in unique_surf_names:
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

    # def write_geom_file(self, filename):
    #     """write the current avl geometry to a file"""
    #     surf_geom = self.get_surface_params()

    # def __write_surface(self, fid, surf_name, surf_data):
    #     """write a surface to a file"""
    #     fid.write("SURFACE

    # Utility functions
    def get_num_surfaces(self) -> int:
        """Get the number of surfaces in the geometry"""
        return self.get_avl_fort_arr("CASE_I", "NSURF")

    def get_num_strips(self) -> int:
        """Get the number of strips in the mesh"""
        return self.get_avl_fort_arr("CASE_I", "NSTRIP")

    def get_mesh_size(self) -> int:
        """Get the number of panels in the mesh"""
        return self.get_avl_fort_arr("CASE_I", "NVOR")

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
        for ii in range(len(fortArray)):
            if fortArray[ii].dtype == np.dtype("|S0"):
                # there is characters in the sting to add
                continue
            if fortArray[ii].dtype == np.dtype("|S1"):
                strList.append(b"".join(fortArray[ii]).decode().strip())
            elif fortArray[ii].dtype == np.dtype("<U1"):
                strList.append("".join(fortArray[ii]).strip())
            # TODO: need a more general solution for |S<variable> type
            elif fortArray[ii].dtype == np.dtype("|S16"):
                strList.append(fortArray[ii].decode().strip())
            elif fortArray[ii].dtype == np.dtype("|S40"):
                strList.append(fortArray[ii].decode().strip())
            else:
                raise TypeError(f"Unable to convert {fortArray[ii]} of type {fortArray[ii].dtype} to string")

        return strList

    # Derivative routines
    # def get_deriv(self):
    #     self.avl.case_r_d.alfad = 1.0

    #     self.avl.sfforc_d()

    #     print("cl_d", self.avl.case_r_d.cltotd)
    def set_constraint_ad_seeds(con_seeds: Dict[str, float]) -> None:
        for con in con_seeds:
            blk, var = self.con_var_to_fort_var[con]
            blk += self.ad_suffix
            var += self.ad_suffix
            self.set_avl_fort_arr(blk, var, con_seeds[con])

    def set_surface_ad_seeds(self, surf_name, param, val):
        if param in self.surf_geom_to_fort_var[surf_name].keys():
            fort_var = self.surf_geom_to_fort_var[surf_name][param]
        else:
            raise ValueError(
                f"param, {param}, not in found for {surf_name}, that has geom data {list(self.surf_geom_to_fort_var[surf_name].keys()) + list(self.surf_pannel_to_fort_var[surf_name].keys())}"
            )
        fort_var[0] += self.ad_suffix
        fort_var[1] += self.ad_suffix

        self.set_avl_fort_arr(fort_var[0], fort_var[1], val, slicer=fort_var[2])

    def get_function_ad_seeds(self):
        pass

    def clear_ad_seeds(self):
        for att in dir(self.avl):
            if att.endswith(self.ad_suffix):
                # loop over the attributes of the common block
                diff_blk = getattr(self.avl, att)
                for _var in dir(diff_blk):
                    val = getattr(diff_blk, _var)
                    setattr(diff_blk, _var, val * 0.0)

    # def set_geom_ad_seeds(surf_seeds: Dict[str, Dict[str, any]]) -> None:
    #     pass
    #     # set the seeds as you would the surface data

    def execute_jac_vec_prod_fwd(con_seeds: Dict[str, float], surf_seeds: Dict[str, Dict[str, any]]):
        # set derivative seeds
        self.clear_ad_seeds()
        self.set_constraint_ad_seeds(con_seeds)

        # propogate the seeds through without resolving
        self.avl_solver.avl.update_surfaces_d()
        self.avl_solver.avl.aero_d()

        # extract derivatives seeds and set the output dict of functions
        self.get_function_ad_seeds()

    def execute_jac_vec_prod_rev():
        pass

    def execute_adjoint_solve():
        pass

    def execute_direct_solve():
        pass

    def execute_run_sensitivies(funcs, con_var, geom_var):
        """
        funcs: dictionary of functions that need derivatives
        con_var: dictionary of contraints that need derivatives
        geom_var: dictionary of geometric variables the need derivatives

        """
        pass
