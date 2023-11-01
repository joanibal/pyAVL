# =============================================================================
# Extension modules
# =============================================================================
from pyavl import AVLSolver
import copy
import resource

# =============================================================================
# Standard Python Modules
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import unittest
import numpy as np


base_dir = os.path.dirname(os.path.abspath(__file__))  # Path to current folder
geom_file = os.path.join(base_dir, "aircraft.avl")
mass_file = os.path.join(base_dir, "aircraft.mass")
geom_mod_file = os.path.join(base_dir, "aircraft_mod.avl")


class TestFunctionPartials(unittest.TestCase):
    def setUp(self):
        # self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        # self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")
        self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.add_constraint("alpha", 25.0)
        self.avl_solver.add_constraint("beta", 5.0)
        self.avl_solver.execute_run()

    def test_fwd_aero_constraint(self):
        for con_key in self.avl_solver.con_var_to_fort_var:
            func_seeds, _, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, geom_seeds={})

            func_seeds_FD, _, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                con_seeds={con_key: 1.0}, geom_seeds={}, mode="FD", step=1e-7
            )

            for func_key in func_seeds:
                print(f"{func_key} wrt {con_key}", func_seeds[func_key], func_seeds_FD[func_key])
                tol = 1e-13
                if np.abs(func_seeds[func_key]) < tol or np.abs(func_seeds_FD[func_key]) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        func_seeds[func_key],
                        func_seeds_FD[func_key],
                        atol=1e-6,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        func_seeds[func_key],
                        func_seeds_FD[func_key],
                        rtol=1e-5,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )

    def test_rev_aero_constraint(self):
        for con_key in self.avl_solver.con_var_to_fort_var:
            self.avl_solver.clear_ad_seeds_fast()

            func_seeds_fwd, _, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, geom_seeds={})
            self.avl_solver.clear_ad_seeds_fast()

            for func_key in self.avl_solver.case_var_to_fort_var:
                con_seeds_rev, _, _, _ = self.avl_solver.execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})
                print(f"{func_key} wrt {con_key}", "fwd", func_seeds_fwd[func_key], "rev", con_seeds_rev[con_key])
                tol = 1e-14

                if np.abs(func_seeds_fwd[func_key]) < tol or np.abs(con_seeds_rev[con_key]) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        func_seeds_fwd[func_key],
                        con_seeds_rev[con_key],
                        atol=1e-14,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        func_seeds_fwd[func_key],
                        con_seeds_rev[con_key],
                        rtol=1e-12,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )

    def test_fwd_geom(self):
        np.random.seed(111)
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                func_seeds, _, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )

                func_seeds_FD, _, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-7
                )

                for func_key in func_seeds:
                    rel_error = np.linalg.norm(func_seeds[func_key] - func_seeds_FD[func_key]) / np.linalg.norm(
                        func_seeds_FD[func_key] + 1e-15
                    )

                    print(
                        f"{func_key:10} wrt {surf_key}:{geom_key} AD:{func_seeds[func_key]: .5e} FD:{func_seeds_FD[func_key]: .5e} rel_error:{rel_error: .3e}"
                    )

                    tol = 1e-13
                    if np.abs(func_seeds[func_key]) < tol or np.abs(func_seeds_FD[func_key]) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            func_seeds[func_key],
                            func_seeds_FD[func_key],
                            atol=1e-6,
                            err_msg=f"func_key {func_key} w.r.t. {geom_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            func_seeds[func_key],
                            func_seeds_FD[func_key],
                            rtol=5e-4,
                            err_msg=f"func_key {func_key} w.r.t. {geom_key}",
                        )

    def test_rev_geom(self):
        np.random.seed(111)

        sens_dict_rev = {}
        for func_key in self.avl_solver.case_var_to_fort_var:
            _, sens_dict_rev[func_key], _, _ = self.avl_solver.execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})
        self.avl_solver.clear_ad_seeds_fast()

        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                func_seeds_fwd, _, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )
                for func_key in func_seeds_fwd:
                    # use dot product test as design variables maybe arrays
                    rev_sum = np.sum(sens_dict_rev[func_key][surf_key][geom_key] * geom_seeds)
                    fwd_sum = np.sum(func_seeds_fwd[func_key])

                    # # print(geom_seeds_rev)
                    tol = 1e-13
                    print(f"{func_key} wrt {surf_key}:{geom_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(fwd_sum) < tol or np.abs(rev_sum) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            atol=1e-14,
                            err_msg=f"func_key {func_key} w.r.t. {surf_key}:{geom_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            rtol=1e-12,
                            err_msg=f"func_key {func_key} w.r.t. {surf_key}:{geom_key}",
                        )

    def test_fwd_gamma(self):
        num_gamma = self.avl_solver.get_mesh_size()
        gamma_seeds = np.random.rand(num_gamma)

        func_seeds, _, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds)
        func_seeds_FD, _, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds, mode="FD", step=1e-7)

        for func_key in func_seeds:
            rel_error = np.linalg.norm(func_seeds[func_key] - func_seeds_FD[func_key]) / np.linalg.norm(
                func_seeds_FD[func_key] + 1e-15
            )
            print(
                f"{func_key:10} AD:{func_seeds[func_key]: .5e} FD:{func_seeds_FD[func_key]: .5e} rel_error:{rel_error: .3e}"
            )

            tol = 1e-13
            if np.abs(func_seeds[func_key]) < tol or np.abs(func_seeds_FD[func_key]) < tol:
                # If either value is basically zero, use an absolute tolerance
                np.testing.assert_allclose(
                    func_seeds[func_key],
                    func_seeds_FD[func_key],
                    atol=1e-6,
                    err_msg=f"func_key {func_key} w.r.t. gamma",
                )
            else:
                np.testing.assert_allclose(
                    func_seeds[func_key],
                    func_seeds_FD[func_key],
                    rtol=5e-5,
                    err_msg=f"func_key {func_key} w.r.t. gamma",
                )

    def test_rev_gamma(self):
        num_gamma = self.avl_solver.get_mesh_size()
        gamma_seeds_fwd = np.random.rand(num_gamma)

        func_seeds_fwd, _, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds_fwd)
        self.avl_solver.clear_ad_seeds_fast()

        for func_key in self.avl_solver.case_var_to_fort_var:
            _, _, gamma_seeds_rev, _ = self.avl_solver.execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})

            rev_sum = np.sum(gamma_seeds_rev * gamma_seeds_fwd)
            fwd_sum = np.sum(func_seeds_fwd[func_key])

            np.testing.assert_allclose(
                fwd_sum,
                rev_sum,
                atol=1e-14,
                err_msg=f"func_key {func_key} w.r.t. gamma",
            )


class TestResidualPartials(unittest.TestCase):
    def setUp(self):
        # self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.add_constraint("alpha", 25.0)
        self.avl_solver.add_constraint("beta", 5.0)
        self.avl_solver.execute_run()

    def test_fwd_aero_constraint(self):
        for con_key in self.avl_solver.con_var_to_fort_var:
            _, res_seeds, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, geom_seeds={})

            _, res_seeds_FD, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                con_seeds={con_key: 1.0}, geom_seeds={}, mode="FD", step=1e-8
            )

            np.testing.assert_allclose(
                res_seeds,
                res_seeds_FD,
                rtol=1e-5,
            )

    def test_rev_aero_constraint(self):
        num_res = self.avl_solver.get_mesh_size()
        res_seeds_rev = np.random.rand(num_res)
        con_seeds_rev, _, _, _ = self.avl_solver.execute_jac_vec_prod_rev(res_seeds=res_seeds_rev)

        self.avl_solver.clear_ad_seeds_fast()

        for con_key in self.avl_solver.con_var_to_fort_var:
            _, res_seeds_fwd, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})

            # do dot product
            res_sum = np.sum(res_seeds_rev * res_seeds_fwd)
            con_sum = np.sum(con_seeds_rev[con_key])

            # print(f"res wrt {con_key}", "rev", con_sum, "fwd", res_sum)

    def test_fwd_geom(self):
        np.random.seed(111)
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                _, res_seeds, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )

                _, res_seeds_FD, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-8
                )
                abs_error = np.abs(res_seeds - res_seeds_FD)
                rel_error = (res_seeds - res_seeds_FD) / (res_seeds + 1e-15)
                idx_max_rel_error = np.argmax(np.abs(rel_error))
                idx_max_abs_error = np.argmax(np.abs(abs_error))

                # print(
                #     f"{surf_key:10} {geom_key:10} AD:{np.linalg.norm(res_seeds): .5e} FD:{np.linalg.norm(res_seeds_FD): .5e} max rel err:{(rel_error[idx_max_rel_error]): .3e} max abs err:{(np.max(abs_error)): .3e}"
                # )
                np.testing.assert_allclose(
                    res_seeds,
                    res_seeds_FD,
                    atol=1e-5,
                    err_msg=f"func_key res w.r.t. {surf_key}:{geom_key}",
                )

    def test_rev_geom(self):
        num_res = self.avl_solver.get_mesh_size()
        res_seeds_rev = np.random.seed(111)
        res_seeds_rev = np.random.rand(num_res)
        _, geom_seeds_rev, _, _ = self.avl_solver.execute_jac_vec_prod_rev(res_seeds=res_seeds_rev)

        self.avl_solver.clear_ad_seeds_fast()
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                _, res_seeds, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )

                # do dot product
                res_sum = np.sum(res_seeds_rev * res_seeds)
                geom_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds)

                # print(f"res wrt {surf_key}:{geom_key}", "rev", geom_sum, "fwd", res_sum)

                np.testing.assert_allclose(
                    res_sum,
                    geom_sum,
                    atol=1e-14,
                    err_msg=f"func_key res w.r.t. {surf_key}:{geom_key}",
                )

    def test_fwd_gamma(self):
        num_gamma = self.avl_solver.get_mesh_size()
        np.random.seed(111)
        gamma_seeds = np.random.rand(num_gamma)

        _, res_seeds, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds)
        _, res_seeds_FD, _, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds, mode="FD", step=1e-7)

        np.testing.assert_allclose(
            res_seeds,
            res_seeds_FD,
            rtol=1e-5,
        )


class TestResidualDPartials(unittest.TestCase):
    def setUp(self):
        # self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.add_constraint("alpha", 25.0)
        self.avl_solver.add_constraint("beta", 5.0)
        self.avl_solver.execute_run()

    # def tearDown(self) -> None:
    #     mb_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
    #     print("Memory usage: %s MB" % mb_memory)

    def test_fwd_aero_constraint(self):
        for con_key in self.avl_solver.con_var_to_fort_var:
            _, _, _, res_d_seeds = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})

            _, _, _, res_d_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(
                con_seeds={con_key: 1.0}, mode="FD", step=1e-5
            )

            np.testing.assert_allclose(
                res_d_seeds,
                res_d_seeds_FD,
                rtol=1e-5,
            )

    def test_rev_aero_constraint(self):
        # for con_key in self.avl_solver.con_var_to_fort_var:
        for con_key in ["beta", "beta"]:
            _, _, _, res_d_seeds = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})

            num_gamma = self.avl_solver.get_mesh_size()
            num_consurf = self.avl_solver.get_num_control_surfs()
            res_d_seeds_rev = np.random.rand(num_consurf, num_gamma)
            res_d_seeds_rev = np.ones_like(res_d_seeds_rev)

            self.avl_solver.clear_ad_seeds_fast()

            con_seeds, _, _, _ = self.avl_solver.execute_jac_vec_prod_rev(res_d_seeds=res_d_seeds_rev)
            self.avl_solver.clear_ad_seeds_fast()

            # do dot product
            res_sum = np.sum(res_d_seeds_rev * res_d_seeds)

            print(f"res wrt {con_key}", "fwd", res_sum, "rev", con_seeds[con_key])

            np.testing.assert_allclose(
                res_sum,
                con_seeds[con_key],
                atol=1e-14,
                err_msg=f"func_key res w.r.t. {con_key}",
            )

    def test_fwd_geom(self):
        np.random.seed(111)
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                _, _, _, res_d_seeds = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )

                _, _, _, res_d_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-8
                )
                abs_error = np.abs(res_d_seeds.flatten() - res_d_seeds_FD.flatten())
                rel_error = np.abs((res_d_seeds.flatten() - res_d_seeds_FD.flatten()) / (res_d_seeds.flatten() + 1e-15))

                idx_max_rel_error = np.argmax(rel_error)
                idx_max_abs_error = np.argmax(abs_error)
                print(
                    f"{surf_key:10} {geom_key:10} AD:{np.linalg.norm(res_d_seeds): .5e} FD:{np.linalg.norm(res_d_seeds_FD): .5e} max rel err:{(rel_error[idx_max_rel_error]): .3e} max abs err:{(np.max(abs_error)): .3e}"
                )
                np.testing.assert_allclose(
                    res_d_seeds,
                    res_d_seeds_FD,
                    atol=1e-4,
                    err_msg=f"func_key res w.r.t. {surf_key}:{geom_key}",
                )

    def test_rev_geom(self):
        np.random.seed(111)
        num_gamma = self.avl_solver.get_mesh_size()
        num_consurf = self.avl_solver.get_num_control_surfs()
        res_d_seeds_rev = np.random.rand(num_consurf, num_gamma)
        res_d_seeds_rev = np.ones_like(res_d_seeds_rev)

        self.avl_solver.clear_ad_seeds_fast()

        _, geom_seeds_rev, _, _ = self.avl_solver.execute_jac_vec_prod_rev(res_d_seeds=res_d_seeds_rev)
        self.avl_solver.clear_ad_seeds_fast()

        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                _, _, _, res_d_seeds_fwd = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )

                # do dot product
                res_sum = np.sum(res_d_seeds_rev * res_d_seeds_fwd)
                geom_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds)

                print(f"res wrt {surf_key}:{geom_key}", "rev", geom_sum, "fwd", res_sum)

                np.testing.assert_allclose(
                    res_sum,
                    geom_sum,
                    atol=1e-14,
                    err_msg=f"func_key res w.r.t. {surf_key}:{geom_key}",
                )
        self.avl_solver.clear_ad_seeds_fast()

    def test_fwd_gamma_d(self):
        num_gamma = self.avl_solver.get_mesh_size()
        num_consurf = self.avl_solver.get_num_control_surfs()
        gamma_d_seeds = np.random.rand(num_consurf, num_gamma)
        # gamma_d_seeds = np.array([[1],[0]])

        _, _, _, res_d_seeds = self.avl_solver.execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds)

        _, _, _, res_d_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(
            gamma_d_seeds=gamma_d_seeds, mode="FD", step=1e-0
        )
        # print("res_d_seeds", res_d_seeds)
        # print("res_d_seeds_FD", res_d_seeds_FD)

        np.testing.assert_allclose(
            res_d_seeds,
            res_d_seeds_FD,
            rtol=1e-5,
        )

    def test_rev_gamma_d(self):
        num_gamma = self.avl_solver.get_mesh_size()
        num_consurf = self.avl_solver.get_num_control_surfs()
        gamma_d_seeds_fwd = np.random.rand(num_consurf, num_gamma)

        res_d_seeds_rev = np.random.rand(num_consurf, num_gamma)

        _, _, _, res_d_seeds_fwd = self.avl_solver.execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds_fwd)
        self.avl_solver.clear_ad_seeds_fast()

        _, _, _, gamma_d_seeds_rev = self.avl_solver.execute_jac_vec_prod_rev(res_d_seeds=res_d_seeds_rev)

        gamma_sum = np.sum(gamma_d_seeds_rev * gamma_d_seeds_fwd)
        res_sum = np.sum(res_d_seeds_rev * res_d_seeds_fwd)

        print("fwd_sum", gamma_sum, "rev_sum", res_sum)
        np.testing.assert_allclose(
            gamma_sum,
            res_sum,
            atol=1e-14,
            err_msg=f"res w.r.t. gamma",
        )


class TestConSurfDerivsPartials(unittest.TestCase):
    def setUp(self):
        # self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.add_constraint("alpha", 45.0)
        self.avl_solver.add_constraint("beta", 45.0)
        self.avl_solver.execute_run()
        self.avl_solver.clear_ad_seeds_fast()

    def test_fwd_aero_constraint(self):
        for con_key in self.avl_solver.con_var_to_fort_var:
            _, _, cs_d, _ = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})

            _, _, cs_d_fd, _ = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, mode="FD", step=1e-8)

            for func_key in cs_d:
                for cs_key in cs_d[func_key]:
                    sens_label = f"d{func_key}/d{cs_key} wrt {con_key}"
                    np.testing.assert_allclose(
                        cs_d[func_key][cs_key],
                        cs_d_fd[func_key][cs_key],
                        rtol=1e-3,
                        err_msg=sens_label,
                    )

    def test_rev_aero_constraint(self):
        cs_names = self.avl_solver.get_control_names()

        cs_deriv_seeds = {}
        for func_key in self.avl_solver.case_derivs_to_fort_var:
            cs_deriv_seeds[func_key] = {}
            for cs_key in cs_names:
                cs_deriv_seeds[func_key][cs_key] = np.random.rand(1)[0]

        con_seeds_rev, _, _, _ = self.avl_solver.execute_jac_vec_prod_rev(consurf_derivs_seeds=cs_deriv_seeds)

        self.avl_solver.clear_ad_seeds_fast()

        for con_key in self.avl_solver.con_var_to_fort_var:
            _, _, cs_deriv_seeds_fwd, _ = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})

            cs_deriv_sum = 0.0
            for func_key in cs_deriv_seeds_fwd:
                for cs_key in cs_deriv_seeds_fwd[func_key]:
                    cs_deriv_sum += cs_deriv_seeds[func_key][cs_key] * cs_deriv_seeds_fwd[func_key][cs_key]

            # do dot product
            con_sum = np.sum(con_seeds_rev[con_key])

            print(f"cs_dervs wrt {con_key}", "rev", con_sum, "fwd", cs_deriv_sum)

            np.testing.assert_allclose(
                con_sum,
                cs_deriv_sum,
                atol=1e-14,
                err_msg=f"cs_dervs wrt {con_key}",
            )

    def test_fwd_geom(self):
        np.random.seed(111)
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                _, _, cs_d, _ = self.avl_solver.execute_jac_vec_prod_fwd(geom_seeds={surf_key: {geom_key: geom_seeds}})

                _, _, cs_d_fd, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-8
                )

                for func_key in cs_d:
                    for cs_key in cs_d[func_key]:
                        sens_label = f"d{func_key}/d{cs_key} wrt {surf_key}:{geom_key:5}"

                        # rel_err = np.linalg.norm(cs_d[func_key][cs_key] - cs_d_fd[func_key][cs_key]) / np.linalg.norm(cs_d[func_key][cs_key])

                        print(f"{sens_label} AD:{cs_d[func_key][cs_key]} FD:{cs_d_fd[func_key][cs_key]}")

                        #
                        # TODO: fix this: add absolute tolerance when true value is zero
                        # np.testing.assert_allclose(
                        #     cs_d[func_key][cs_key],
                        #     cs_d_fd[func_key][cs_key],
                        #     rtol=1e-3,
                        #     err_msg=sens_label,
                        # )

    def test_rev_geom(self):
        np.random.seed(111)
        cs_names = self.avl_solver.get_control_names()

        cs_d_rev = {}
        for func_key in self.avl_solver.case_derivs_to_fort_var:
            cs_d_rev[func_key] = {}
            for cs_key in cs_names:
                cs_d_rev[func_key][cs_key] = np.random.rand(1)[0]

        _, geom_seeds_rev, _, _ = self.avl_solver.execute_jac_vec_prod_rev(consurf_derivs_seeds=cs_d_rev)
        self.avl_solver.clear_ad_seeds_fast()

        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds_fwd = np.random.rand(*arr.shape)

                func_seeds_fwd, _, cs_d_fwd, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds_fwd}}
                )

                for func_key in func_seeds_fwd:
                    # use dot product test as design variables maybe arrays
                    rev_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds_fwd)

                    fwd_sum = 0.0
                    for func_key in cs_d_fwd:
                        for cs_key in cs_d_fwd[func_key]:
                            fwd_sum += cs_d_rev[func_key][cs_key] * cs_d_fwd[func_key][cs_key]

                    # fwd_sum = np.sum(func_seeds_fwd[func_key])

                    # # print(geom_seeds_rev)
                    tol = 1e-13
                    print(f"{func_key} wrt {surf_key}:{geom_key}", "fwd", fwd_sum, "rev", rev_sum)
                    if np.abs(fwd_sum) < tol or np.abs(rev_sum) < tol:
                        # If either value is basically zero, use an absolute tolerance
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            atol=1e-14,
                            err_msg=f"func_key {func_key} w.r.t. {surf_key}:{geom_key}",
                        )
                    else:
                        np.testing.assert_allclose(
                            fwd_sum,
                            rev_sum,
                            rtol=1e-12,
                            err_msg=f"func_key {func_key} w.r.t. {surf_key}:{geom_key}",
                        )

    def test_fwd_gamma_d(self):
        num_gamma = self.avl_solver.get_mesh_size()
        num_consurf = self.avl_solver.get_num_control_surfs()
        gamma_d_seeds = np.random.rand(num_consurf, num_gamma)

        _, _, cs_d, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds)
        _, _, cs_d_fd, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds, mode="FD", step=1e-7)

        for func_key in cs_d:
            for cs_key in cs_d[func_key]:
                sens_label = f"d{func_key}/d{cs_key} wrt gamma_d"
                print(sens_label, cs_d[func_key][cs_key], cs_d_fd[func_key][cs_key])
                np.testing.assert_allclose(
                    cs_d[func_key][cs_key],
                    cs_d_fd[func_key][cs_key],
                    rtol=1e-3,
                    err_msg=sens_label,
                )

    def test_rev_gamma_d(self):
        num_gamma = self.avl_solver.get_mesh_size()
        num_consurf = self.avl_solver.get_num_control_surfs()
        gamma_d_seeds_fwd = np.random.rand(num_consurf, num_gamma)

        _, _, cs_d_fwd, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_d_seeds=gamma_d_seeds_fwd)
        self.avl_solver.clear_ad_seeds_fast()

        for func_key in cs_d_fwd:
            for cs_key in cs_d_fwd[func_key]:
                cs_d_rev = {func_key: {cs_key: 1.0}}

                _, _, _, gamma_d_seeds_rev = self.avl_solver.execute_jac_vec_prod_rev(consurf_derivs_seeds=cs_d_rev)

                rev_sum = np.sum(gamma_d_seeds_rev * gamma_d_seeds_fwd)

                fwd_sum = np.sum(cs_d_fwd[func_key][cs_key])

                print("fwd_sum", fwd_sum, "rev_sum", rev_sum)
                np.testing.assert_allclose(
                    fwd_sum,
                    rev_sum,
                    atol=1e-14,
                    err_msg=f"func_key {func_key} w.r.t. gamma",
                )


if __name__ == "__main__":
    unittest.main()
