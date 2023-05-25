# =============================================================================
# Extension modules
# =============================================================================
from pyavl import AVLSolver
import copy

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


class TestNewSubroutines(unittest.TestCase):
    def setUp(self):
        self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.add_constraint("alpha", 25.0)

    def test_residual(self):
        self.avl_solver.avl.get_res()
        res = self.avl_solver.avl.VRTX_R.RES[:]
        rhs = self.avl_solver.avl.VRTX_R.RHS[:]
        np.testing.assert_allclose(
            res,
            -1 * rhs,
            atol=1e-15,
        )

        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.get_res()
        res = self.avl_solver.avl.VRTX_R.RES[:]

        np.testing.assert_allclose(
            res,
            np.zeros_like(res),
            atol=5e-14,
        )

    def test_new_solve(self):
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()
        gam_new = self.avl_solver.avl.VRTX_R.GAM[:]
        coef_data_new = self.avl_solver.get_case_total_data()

        self.avl_solver.execute_run()
        gam = self.avl_solver.avl.VRTX_R.GAM[:]
        coef_data = self.avl_solver.get_case_total_data()

        np.testing.assert_allclose(
            gam,
            gam_new,
            atol=1e-15,
        )
        for func_key in coef_data:
            np.testing.assert_allclose(
                coef_data[func_key],
                coef_data_new[func_key],
                atol=1e-14,
                err_msg=f"func_key {func_key}",
            )


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
            func_seeds, _ = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, geom_seeds={})

            func_seeds_FD, _ = self.avl_solver.execute_jac_vec_prod_fwd(
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
            self.avl_solver.clear_ad_seeds()

            func_seeds_fwd, _ = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, geom_seeds={})
            self.avl_solver.clear_ad_seeds()

            for func_key in self.avl_solver.case_var_to_fort_var:
                con_seeds_rev, _, _ = self.avl_solver.execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})
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

                func_seeds, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )

                func_seeds_FD, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-7
                )

                for func_key in ["CD"]:
                    # for func_key in func_seeds:
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
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                self.avl_solver.clear_ad_seeds()

                func_seeds_fwd, _ = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )
                self.avl_solver.clear_ad_seeds()
                
                #TODO: split the loops for speed. see total tests
                for func_key in func_seeds_fwd:
                    _, geom_seeds_rev, _ = self.avl_solver.execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})

                    # use dot product test as design variables maybe arrays
                    rev_sum = np.sum(geom_seeds_rev[surf_key][geom_key] * geom_seeds)
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

        func_seeds, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds)
        func_seeds_FD, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds, mode="FD", step=1e-7)

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

        func_seeds_fwd, _ = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds_fwd)
        self.avl_solver.clear_ad_seeds()

        for func_key in self.avl_solver.case_var_to_fort_var:
            _, _, gamma_seeds_rev = self.avl_solver.execute_jac_vec_prod_rev(func_seeds={func_key: 1.0})

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
            _, res_seeds = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0}, geom_seeds={})

            _, res_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(
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
        con_seeds_rev, _, _ = self.avl_solver.execute_jac_vec_prod_rev(res_seeds=res_seeds_rev)

        self.avl_solver.clear_ad_seeds()

        for con_key in self.avl_solver.con_var_to_fort_var:
            _, res_seeds_fwd = self.avl_solver.execute_jac_vec_prod_fwd(con_seeds={con_key: 1.0})

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

                _, res_seeds = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}
                )

                _, res_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(
                    con_seeds={}, geom_seeds={surf_key: {geom_key: geom_seeds}}, mode="FD", step=1e-10
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
        _, geom_seeds_rev, _ = self.avl_solver.execute_jac_vec_prod_rev(res_seeds=res_seeds_rev)

        self.avl_solver.clear_ad_seeds()
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            for geom_key in self.avl_solver.surf_geom_to_fort_var[surf_key]:
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                geom_seeds = np.random.rand(*arr.shape)

                _, res_seeds = self.avl_solver.execute_jac_vec_prod_fwd(
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

        _, res_seeds = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds)
        _, res_seeds_FD = self.avl_solver.execute_jac_vec_prod_fwd(gamma_seeds=gamma_seeds, mode="FD", step=1e-7)

        np.testing.assert_allclose(
            res_seeds,
            res_seeds_FD,
            rtol=1e-5,
        )


class TestTotals(unittest.TestCase):
    # TODO: beta derivatives likely wrong
    # TODO: CN SA and CR SA are not implemneted (calculated in calcST)

    def setUp(self):
        # self.avl_solver = AVLSolver(geo_file=geom_file, mass_file=mass_file)
        self.avl_solver = AVLSolver(geo_file="aircraft_L1.avl")
        # self.avl_solver = AVLSolver(geo_file="rect.avl")
        self.avl_solver.add_constraint("alpha", 25.0)
        self.avl_solver.add_constraint("beta", 5.0)
        self.avl_solver.execute_run()

    def finite_dif(self, con_list, geom_seeds, step=1e-7):
        con_seeds = {}

        for con in con_list:
            con_seeds[con] = 1.0

        self.avl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=step)
        self.avl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=step)

        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()
        # self.avl_solver.execute_run()
        coef_data_peturb = self.avl_solver.get_case_total_data()

        self.avl_solver.set_constraint_ad_seeds(con_seeds, mode="FD", scale=-step)
        self.avl_solver.set_geom_ad_seeds(geom_seeds, mode="FD", scale=-step)

        self.avl_solver.avl.update_surfaces()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.exec_rhs()
        self.avl_solver.avl.get_res()
        self.avl_solver.avl.velsum()
        self.avl_solver.avl.aero()
        # self.avl_solver.execute_run()

        coef_data = self.avl_solver.get_case_total_data()

        func_seeds = {}
        for func_key in coef_data:
            func_seeds[func_key] = (coef_data_peturb[func_key] - coef_data[func_key]) / step

        return func_seeds

    def test_aero_constraint(self):
        # compare the analytical gradients with finite difference for each constraint and function
        func_vars = self.avl_solver.case_var_to_fort_var
        sens = self.avl_solver.execute_run_sensitivies(func_vars)

        for con_key in self.avl_solver.con_var_to_fort_var:
            # for con_key in ['beta']:
            sens_FD = self.finite_dif([con_key], {}, step=1.0e-5)

            for func_key in func_vars:
                ad_dot = sens[func_key][con_key]
                fd_dot = sens_FD[func_key]

                # print(f"{func_key} wrt {con_key}", "AD", ad_dot, "FD", fd_dot)
                rel_err = np.abs((ad_dot - fd_dot) / fd_dot)

                print(f"{func_key:5} wrt {con_key:5} | AD:{ad_dot: 5e} FD:{fd_dot: 5e} rel err:{rel_err:.2e}")

                tol = 1e-13
                if np.abs(ad_dot) < tol or np.abs(fd_dot) < tol:
                    # If either value is basically zero, use an absolute tolerance
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        atol=1e-5,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )
                else:
                    np.testing.assert_allclose(
                        ad_dot,
                        fd_dot,
                        rtol=5e-5,
                        err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    )

    def test_geom(self):
        # compare the analytical gradients with finite difference for each
        # geometric variable and function

        surf_key = list(self.avl_solver.surf_geom_to_fort_var.keys())[0]
        geom_vars = self.avl_solver.surf_geom_to_fort_var[surf_key]
        func_vars = self.avl_solver.case_var_to_fort_var.keys()

        sens = self.avl_solver.execute_run_sensitivies(func_vars)

        # for con_key in self.avl_solver.con_var_to_fort_var:
        sens_FD = {}
        for surf_key in self.avl_solver.surf_geom_to_fort_var:
            sens_FD[surf_key] = {}
            for geom_key in geom_vars:
                # _, res_seeds = self.avl_solver.execute_jac_vec_prod_fwd(
                #     con_seeds={}, geom_seeds={surf_key: {geom_key: 1.0}}
                # )
                arr = self.avl_solver.get_surface_param(surf_key, geom_key)
                np.random.seed(arr.size)
                rand_arr = np.random.rand(*arr.shape)
                rand_arr /= np.linalg.norm(rand_arr)
                # geom_seeds[surf_key][geom_key] = rand_arr

                sens_FD[surf_key][geom_key] = {}

                func_seeds = self.finite_dif([], {surf_key: {geom_key: rand_arr}}, step=1.0e-3)

                for func_key in func_seeds:
                    sens_FD[surf_key][geom_key][func_key] = func_seeds[func_key]

                    geom_dot = np.sum(sens[func_key][surf_key][geom_key] * rand_arr)
                    func_dot = sens_FD[surf_key][geom_key][func_key]

                    rel_err = np.abs(geom_dot - func_dot) / np.abs(func_dot)

                    print(
                        f"{func_key:5} wrt {surf_key}:{geom_key:10} | AD:{geom_dot: 5e} FD:{func_dot: 5e} rel err:{rel_err:.2e}"
                    )

                    # tol = 1e-13
                    # if np.abs(sens[func_key][con_key]) < tol or np.abs(sens_FD[func_key]) < tol:
                    #     # If either value is basically zero, use an absolute tolerance
                    #     np.testing.assert_allclose(
                    #         sens[func_key][con_key],
                    #         sens_FD[func_key],
                    #         atol=1e-5,
                    #         err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    #     )
                    # else:
                    #     np.testing.assert_allclose(
                    #         sens[func_key][con_key],
                    #         sens_FD[func_key],
                    #         rtol=5e-5,
                    #         err_msg=f"func_key {func_key} w.r.t. {con_key}",
                    #     )


if __name__ == "__main__":
    unittest.main()
