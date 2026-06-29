

import numpy as np
from scipy.optimize import minimize

from .import_functions import data_import, PTES_geometry_import
from .H2O_prop import density_water, specific_heat_water

class STES:
    def __init__(self, n_layers, T_min=12, T_max=90, T_ref=10):
        self.n_layers = n_layers
        self.T_min = T_min
        self.T_max = T_max
        self.T_ref = T_ref

        # These must be defined by subclasses
        self.V_storage = None
        self.V_layer = None
        self.A_side = None
        self.A_top = None
        self.A_bottom = None
        self.A_side_layers = None

        self.Q_storage_max = None
        self.Q_storage_min = None

        self.T_curves = None
        self.T_curves_Q = None

        self.U_lid = None
        self.U_side = None
        self.U_bottom = None

        self.eta_self_discharge = None

    def compute_energy_bounds(self):
        """Common energy calculation"""

        self.Q_storage_max = (self.V_storage * density_water(self.T_max) * specific_heat_water(self.T_max) * (self.T_max - self.T_ref) / 3.6e9)
        self.Q_storage_min = (self.V_storage * density_water(self.T_min) * specific_heat_water(self.T_min) * (self.T_min - self.T_ref) / 3.6e9)
    
    def set_temperature_map(self, T_curves):
        """Set the temperature map for the storage, which can be used for efficiency calculations"""
        self.T_curves = T_curves
        # calculate volume per layer, density and specific heat of water for each temperature state in the storage
        rho_water = density_water(self.T_curves)
        cp_water = specific_heat_water(self.T_curves)

        # calculate Q_storage for each temperature state in the storage based on the synthetic temperature curves and the geometry of the storage
        self.T_curves_Q = self.V_layer @ (rho_water * cp_water * (self.T_curves - self.T_ref)) / 3.6e9 # MWh, energy content of the storage based on the synthetic temperature curves and the geometry of the storage
    
    def get_temperature_layers(self, Q_storage):

        if Q_storage >= self.Q_storage_min and Q_storage <= self.Q_storage_max:    
            
            Q_1_idx = np.where(self.T_curves_Q <= Q_storage)[0][np.argmax(self.T_curves_Q[self.T_curves_Q <= Q_storage])]
            Q_2_idx = np.where(self.T_curves_Q >= Q_storage)[0][np.argmin(self.T_curves_Q[self.T_curves_Q >= Q_storage])]

            # interpolate between the two closest points to get the temperature profile corresponding to the given Q_storage
            if self.T_curves_Q[Q_2_idx] != self.T_curves_Q[Q_1_idx]:
                T_layer = self.T_curves[:,Q_1_idx] + (self.T_curves[:,Q_2_idx] - self.T_curves[:,Q_1_idx]) * (Q_storage - self.T_curves_Q[Q_1_idx]) / (self.T_curves_Q[Q_2_idx] - self.T_curves_Q[Q_1_idx])
            else:
                T_layer = self.T_curves[:,Q_1_idx]
        
        else:
            print("WARNING: Q_storage is not between Q_storage_min and Q_storage_max")
            if Q_storage <= self.Q_storage_min:
                T_layer = np.full(
                    self.T_curves.shape[0],
                    Q_storage * 3.6e9 / self.V_storage / density_water(self.T_min) / specific_heat_water(self.T_min) + self.T_ref
                )
            else:
                T_layer = np.full(
                    self.T_curves.shape[0],
                    Q_storage * 3.6e9 / self.V_storage / density_water(self.T_max) / specific_heat_water(self.T_max) + self.T_ref
                )
        
        return T_layer.T
    
    def set_U_values(self, U_lid, U_side, U_bottom):
        self.U_lid = U_lid
        self.U_side = U_side
        self.U_bottom = U_bottom


# =====================================================================
# STANDALONE SIMPLE STORAGE MODEL
# =====================================================================
#
# Extracted from PTES.simulate_PTES so it can be reused on its own,
# without running the full (slow) detailed simulation first. Behavior
# is unchanged from the original nested function.


def simulate_storage_simple(eta, Q_charge, Q_discharge, Q_storage_start):
    """
    Simple self-discharge storage model.

    Q_simple[t] = (1 - eta) * Q_simple[t-1] + Q_charge[t] - Q_discharge[t]

    Parameters
    ----------
    eta : float
        Self-discharge rate per timestep.
    Q_charge : array-like
        Charging energy per timestep.
    Q_discharge : array-like
        Discharging energy per timestep.
    Q_storage_start : float
        Storage energy content used to seed Q_simple[-1] (i.e. the
        value the recursion wraps around from at t=0), matching the
        convention used in simulate_PTES.

    Returns
    -------
    Q_simple : np.ndarray
        Simulated storage energy content, same length as Q_charge.
    """

    Q_charge = np.asarray(Q_charge)
    Q_discharge = np.asarray(Q_discharge)

    Q_simple = np.zeros_like(Q_charge, dtype=float)
    Q_simple[-1] = Q_storage_start

    for t in range(len(Q_charge)):
        Q_simple[t] = (1 - eta) * Q_simple[t - 1] + Q_charge[t] - Q_discharge[t]

    return Q_simple


class PTES(STES):
    def __init__(self, h, a, b, c, d, n_layers, T_min=12, T_max=90, T_ref=10):
        super().__init__(n_layers, T_min, T_max, T_ref)

        self.h, self.a, self.b, self.c, self.d = h, a, b, c, d

        self.V_storage = self.volume_truncated_pyramid(h, a, b, c, d)
        self.V_layer = np.flip(self.volume_per_layer_truncated_pyramid(h, a, b, c, d, n_layers))

        self.A_side, self.A_bottom, self.A_top = self.surface_area_truncated_pyramid(h, a, b, c, d)

        self.A_side_layers = np.flip(self.surface_area_per_layer_truncated_pyramid(h, a, b, c, d, n_layers))

        self.compute_energy_bounds()

    def volume_truncated_pyramid(self, h, a, b, c, d):
        '''(https://doi.org/10.1016/j.energy.2018.03.152)'''
        return (h/6)*((2*a+c)*b + (2*c+a)*d)

    def volume_per_layer_truncated_pyramid(self, h, a, b, c, d, n):
        volume_per_layer = []

        diff_a_c = (a - c) / n
        diff_b_d = (b - d) / n

        h_layer = h / n
        cx, dx = c, d

        for _ in range(n):
            ax = cx + diff_a_c
            bx = dx + diff_b_d
            V = self.volume_truncated_pyramid(h_layer, ax, bx, cx, dx)
            volume_per_layer.append(float(V))
            cx, dx = ax, bx

        return volume_per_layer

    def surface_area_truncated_pyramid(self, h, a, b, c, d):
        '''(https://doi.org/10.1016/j.energy.2018.03.152)'''
        A_side = (a+c)*((h**2)+((a-c)/2)**2)**0.5 + (b+d)*((h**2)+((b-d)/2)**2)**0.5

        return A_side, c*d, a*b

    def surface_area_per_layer_truncated_pyramid(self, h, a, b, c, d, n):
        surface_area_per_layer = []

        diff_a_c = (a - c) / n
        diff_b_d = (b - d) / n

        h_layer = h / n
        cx, dx = c, d

        for _ in range(n):
            ax = cx + diff_a_c
            bx = dx + diff_b_d
            A = self.surface_area_truncated_pyramid(h_layer, ax, bx, cx, dx)[0]
            surface_area_per_layer.append(float(A))
            cx, dx = ax, bx

        return surface_area_per_layer
    
    def simulate_PTES(self, file_path, Q_storage_start, sim_start=None, sim_end=None):
        # Extracting data from the DataFrame
        data = data_import(file_path)

        time = data.index
        time = time[sim_start:sim_end]
        Q_charge = data['Q_charge'].to_numpy()
        Q_charge = Q_charge[sim_start:sim_end]
        Q_discharge = data['Q_discharge'].to_numpy()
        Q_discharge = Q_discharge[sim_start:sim_end]
        T_air = data['T_air'].to_numpy()
        T_air = T_air[sim_start:sim_end]
        T_soil = data['T_soil'].to_numpy()
        T_soil = T_soil[sim_start:sim_end]

        hours = time.view('int64') / 1e6 / 3600
        time_step = (hours[-1] - hours[0]) / time.shape[0] # hours, time step of the data

        Q_loss_sim = np.zeros_like(Q_charge)
        Q_storage_sim = np.zeros_like(Q_charge)
        Q_storage_sim[-1] = Q_storage_start

        for i in range(len(Q_storage_sim)):
            
            T_layers = self.get_temperature_layers(Q_storage_sim[i-1])

            T_side_average = np.nansum(T_layers * self.A_side_layers / self.A_side)
            T_bottom = T_layers[-1] # Assuming the first column corresponds to the bottom layer
            T_top = T_layers[0] # Assuming the last column corresponds to the top layer

            Q_loss_sim[i] = self.U_lid * self.A_top * (T_top - T_air[i]) * time_step / 1e6 + self.U_side * self.A_side * (T_side_average - T_soil[i]) * time_step / 1e6 + self.U_bottom * self.A_bottom * (T_bottom - T_soil[i]) * time_step / 1e6
            Q_storage_sim[i] = Q_storage_sim[i-1] + Q_charge[i] - Q_discharge[i] - Q_loss_sim[i]

        def loss_function(eta):
            eta = float(eta[0])
            Q_simple = simulate_storage_simple(eta, Q_charge, Q_discharge, Q_storage_start)

            return np.nanmean((Q_simple - Q_storage_sim)**2)

        result = minimize(loss_function, x0=[0.01], bounds=[(0, 1)])
        eta_self_discharge = result.x[0]
        Q_storage_simple = simulate_storage_simple(eta_self_discharge, Q_charge, Q_discharge, Q_storage_start)

        self.eta_self_discharge = eta_self_discharge
        
        return time, Q_storage_sim, Q_loss_sim, Q_storage_simple

    def calculate_U_values_PTES(self, file_path, Q_storage_start, Q_storage_end, share_loss_lid=0.56, share_loss_side=0.41, share_loss_bottom=0.03, start_idx=None, end_idx=None):

        """
        Iterative determination of U-values based on:
        - initial storage energy
        - final storage energy
        - charge/discharge profiles
        - ambient temperatures

        The U-values are optimized such that the simulated final
        storage energy matches Q_storage_end.
        """

        # ---------------------------------------------------------
        # Normalize loss shares
        # ---------------------------------------------------------

        total_share = share_loss_lid + share_loss_side + share_loss_bottom

        share_loss_lid /= total_share
        share_loss_side /= total_share
        share_loss_bottom /= total_share

        # ---------------------------------------------------------
        # Import data
        # ---------------------------------------------------------

        data = data_import(file_path)

        time = data.index[start_idx:end_idx]

        T_air = data['T_air'].to_numpy()[start_idx:end_idx]
        T_soil = data['T_soil'].to_numpy()[start_idx:end_idx]

        Q_charge = data['Q_charge'].to_numpy()[start_idx:end_idx]
        Q_discharge = data['Q_discharge'].to_numpy()[start_idx:end_idx]

        # ---------------------------------------------------------
        # Time step
        # ---------------------------------------------------------

        hours = time.view('int64') / 1e6 / 3600
        time_step = (hours[-1] - hours[0]) / time.shape[0]

        n = len(time)

        # ---------------------------------------------------------
        # Initial guesses for U-values
        # ---------------------------------------------------------

        # Typical ranges for large PTES systems:
        #
        # lid:    0.2 - 0.6 W/m²K
        # side:   0.1 - 0.4 W/m²K
        # bottom: 0.1 - 0.3 W/m²K
        #
        # These are realistic starting points.

        k0 = np.array([
            0.35,   # U_lid
            0.20,   # U_side
            0.15    # U_bottom
        ])

        # ---------------------------------------------------------
        # Bounds to keep optimization physical
        # ---------------------------------------------------------

        bounds = [
            (0.01, 2.0),   # lid
            (0.01, 2.0),   # side
            (0.01, 2.0)    # bottom
        ]

        # ---------------------------------------------------------
        # Simulation function
        # ---------------------------------------------------------

        def simulate_storage(U_values):

            U_lid, U_side, U_bottom = U_values

            Q_storage_sim = np.zeros(n)
            Q_storage_sim[0] = Q_storage_start

            Q_loss = np.zeros(n)
            Q_loss_lid_array = np.zeros(n)
            Q_loss_side_array = np.zeros(n)
            Q_loss_bottom_array = np.zeros(n)

            for i in range(1, n):

                # Temperatures from current storage state
                T_layers = self.get_temperature_layers(
                    Q_storage_sim[i - 1]
                )

                T_side_average = np.nansum(
                    T_layers * self.A_side_layers / self.A_side
                )

                T_bottom = T_layers[-1]
                T_top = T_layers[0]

                # Heat losses
                Q_loss_lid = U_lid * self.A_top * (T_top - T_air[i]) * time_step / 1e6
                Q_loss_side = U_side * self.A_side * (T_side_average - T_soil[i]) * time_step / 1e6
                Q_loss_bottom = U_bottom * self.A_bottom * (T_bottom - T_soil[i]) * time_step / 1e6

                Q_loss_lid_array[i] = Q_loss_lid
                Q_loss_side_array[i] = Q_loss_side
                Q_loss_bottom_array[i] = Q_loss_bottom

                Q_loss[i] = (Q_loss_lid + Q_loss_side + Q_loss_bottom)

                # Storage balance
                Q_storage_sim[i] = Q_storage_sim[i - 1] + Q_charge[i] - Q_discharge[i] - Q_loss[i]

            return (Q_storage_sim, Q_loss, Q_loss_lid_array, Q_loss_side_array, Q_loss_bottom_array)

        # ---------------------------------------------------------
        # Objective function
        # ---------------------------------------------------------

        def objective(U_values):

            (Q_storage_sim, Q_loss, Q_loss_lid, Q_loss_side, Q_loss_bottom) = simulate_storage(U_values)

            # ---------------------------------------------
            # Final energy matching
            # ---------------------------------------------

            final_error = Q_storage_sim[-1] - Q_storage_end

            # ---------------------------------------------
            # Physical heat loss shares
            # ---------------------------------------------

            total_loss = (np.sum(Q_loss_lid) + np.sum(Q_loss_side) + np.sum(Q_loss_bottom))

            simulated_shares = np.array([np.sum(Q_loss_lid) / total_loss, np.sum(Q_loss_side) / total_loss, np.sum(Q_loss_bottom) / total_loss])

            target_shares = np.array([share_loss_lid, share_loss_side, share_loss_bottom])

            share_error = np.sum((simulated_shares - target_shares)**2)

            # ---------------------------------------------
            # Combined objective
            # ---------------------------------------------

            return (final_error**2 + 1e6 * share_error)

        # ---------------------------------------------------------
        # Optimization
        # ---------------------------------------------------------

        result = minimize(objective, x0=k0, bounds=bounds, method='L-BFGS-B')

        # ---------------------------------------------------------
        # Store results
        # ---------------------------------------------------------

        self.U_lid, self.U_side, self.U_bottom = result.x

        Q_storage_sim, Q_loss, _, _, _ = simulate_storage(result.x)

        return {
            "U_lid": self.U_lid,
            "U_side": self.U_side,
            "U_bottom": self.U_bottom,
            "Q_storage_sim": Q_storage_sim,
            "Q_loss": Q_loss,
            "optimization_success": result.success,
            "optimization_message": result.message
        }

    def calibrate_self_discharge_yearly(self, file_path, Q_storage_start_by_year, Q_storage_end_by_year):
        """
        Calibrate a self-discharge rate (eta) separately for each
        year, using only the simple storage model (no detailed
        simulation), so this runs fast even for many plants/years.

        For each year, the simple model is run starting from
        Q_storage_start_by_year[year] and eta is optimized so that
        the simple model's value at the last timestep of that year
        matches Q_storage_end_by_year[year] as closely as possible
        (single end-of-year point, per your earlier yearly-error
        approach).

        Parameters
        ----------
        file_path : str
            Path to the data file (same format as used by
            simulate_PTES / calculate_U_values_PTES).
        Q_storage_start_by_year : dict
            {year: Q_storage_start} - storage energy content at the
            start of each year to calibrate.
        Q_storage_end_by_year : dict
            {year: Q_storage_end} - measured/target storage energy
            content at the end of each year to calibrate against.
            Must have the same keys as Q_storage_start_by_year.

        Returns
        -------
        results : dict
            {
                year: {
                    "eta_self_discharge": float,
                    "time": <time index for that year>,
                    "Q_storage_simple": np.ndarray,
                    "optimization_success": bool,
                    "optimization_message": str,
                }
            }
        """

        data = data_import(file_path)

        years = np.array(data.index.year)

        results = {}

        for year, Q_storage_start in Q_storage_start_by_year.items():

            if year not in Q_storage_end_by_year:
                continue

            Q_storage_end = Q_storage_end_by_year[year]

            # ---------------------------------------------------
            # Slice data to this year only
            # ---------------------------------------------------

            idx = np.where(years == year)[0]

            if len(idx) == 0:
                continue

            time_year = data.index[idx]
            Q_charge_year = data['Q_charge'].to_numpy()[idx]
            Q_discharge_year = data['Q_discharge'].to_numpy()[idx]

            # ---------------------------------------------------
            # Objective: end-of-year simple-model value vs target
            # ---------------------------------------------------

            def loss_function(eta, Q_charge_year=Q_charge_year, Q_discharge_year=Q_discharge_year,
                               Q_storage_start=Q_storage_start, Q_storage_end=Q_storage_end):

                eta = float(eta[0])

                Q_simple = simulate_storage_simple(eta, Q_charge_year, Q_discharge_year, Q_storage_start)

                rel_error = (Q_simple[-1] - Q_storage_end) / Q_storage_end

                return rel_error**2

            result = minimize(loss_function, x0=[0.01], bounds=[(0, 1)])

            eta_opt = float(result.x[0])

            Q_storage_simple_year = simulate_storage_simple(eta_opt, Q_charge_year, Q_discharge_year, Q_storage_start)

            results[year] = {
                "eta_self_discharge": eta_opt,
                "time": time_year,
                "Q_storage_simple": Q_storage_simple_year,
                "optimization_success": result.success,
                "optimization_message": result.message,
            }

        return results


class TTES(STES):
    def __init__(self, h, r, n_layers, T_min=12, T_max=90, T_ref=10):
        super().__init__(n_layers, T_min, T_max, T_ref)

        self.h = h
        self.r = r

        self.V_storage = self.volume_cylinder(h, r)
        self.V_layer = self.volume_per_layer(h, r, n_layers)

        self.A_side, self.A_top, self.A_bottom = self.surface_area_cylinder(h, r)

        self.A_side_layers = self.surface_area_per_layer(h, r, n_layers)

        self.compute_energy_bounds()

    def volume_cylinder(self, h, r):
        return np.pi * r**2 * h

    def volume_per_layer(self, h, r, n):
        h_layer = h / n
        V_layer = np.pi * r**2 * h_layer
        return [V_layer] * n

    def surface_area_cylinder(self, h, r):
        A_side = 2 * np.pi * r * h
        A_top = np.pi * r**2
        A_bottom = np.pi * r**2
        return A_side, A_top, A_bottom

    def surface_area_per_layer(self, h, r, n):
        h_layer = h / n
        A_layer = 2 * np.pi * r * h_layer
        return [A_layer] * n