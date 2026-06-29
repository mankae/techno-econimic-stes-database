from dataclasses import dataclass
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pathlib import Path

# define path for data
BASE_DIR = Path(__file__).resolve().parent
OPEX_DATA_PATH = BASE_DIR / "data" / "STES_OPEX_data.csv"
CAPEX_DATA_PATH = BASE_DIR / "data" / "STES_CAPEX_data.csv"

from .H2O_prop import density_water, specific_heat_water


# ------------------------------------------------------------------
# Dataclass configuration
# ------------------------------------------------------------------

@dataclass(frozen=True)
class TechnologyConfig:
    name: str
    energy_density: float | None = None


TECHNOLOGIES = {
    "PTES": TechnologyConfig(
        name="PTES",
        energy_density=0.070,
    ),
    "TTES": TechnologyConfig(
        name="TTES",
        energy_density=0.070,
    ),
    "BTES": TechnologyConfig(
        name="BTES",
        energy_density=0.023,
    ),
    "ATES": TechnologyConfig(
        name="ATES",
        energy_density=0.035
    ),
}


# ------------------------------------------------------------------
# OPEX
# ------------------------------------------------------------------

def OPEX_STES(technology):

    if technology not in TECHNOLOGIES:
        raise ValueError(
            "Invalid technology. Please choose from "
            "'PTES', 'TTES', 'BTES', or 'ATES'."
        )

    config = TECHNOLOGIES[technology]

    OPEX_DATA = pd.read_csv(OPEX_DATA_PATH)

    return OPEX_DATA[f"OPEX_{config.name}"][0]


# ------------------------------------------------------------------
# CAPEX
# ------------------------------------------------------------------

def CAPEX_STES(
    technology="PTES",
    unit="per_volume",
    capacity=10000,
    T_min=None,
    T_max=None
):

    '''
    This function return the CAPEX for a given type of STES (PTES, TTES, BTES, ATES)
    based on the capacity and the unit of the capacity.
    
    The units for the capacity can be:
    - "per_volume": CAPEX is calculated per unit of volume [CHF/m³] (PTES, TTES, BTES, ATES)
    - "per_energy": CAPEX is calculated per unit of energy capacity [CHF/MWh] (PTES, TTES, BTES, ATES)
    - "per_borehole_length": CAPEX is calculated per total length of the boreholes [CHF/m] (BTES, ATES)
    
    heat density if no temperature range is given (https://www.euroheat.org/dhc/knowledge-hub/large-storage-systems-for-dhc-networks):
    - PTES: 0.070 MWh/m³
    - TTES: 0.070 MWh/m³
    - BTES: 0.023 MWh/m³
    - ATES: 0.035 MWh/m³

    Inputs:
    technology: PTES, TTES, BTES, ATES
    unit: "per_volume", "per_energy", "per_borehole_length"
    capacity: the capacity for which to calculate the CAPEX (in MWh for "per_energy", in m³ for "per_volume", in m for "per_borehole_length")
              if capacity =0, the points for the piecewise linear CAPEX curve will be returned
    T_min: minimum temperature of the storage (in °C, only meaningful for "per_energy" unit)
    T_max: maximum temperature of the storage (in °C, only meaningful for "per_energy" unit)

    Outputs:
    CAPEX: the calculated CAPEX for the given capacity and unit (in CHF/m³ for "per_volume", in CHF/MWh for "per_energy", in CHF/m for "per_borehole_length")
    capacity_lin: the capacity values for the piecewise linear CAPEX curve (in m³ for "per_volume", in MWh for "per_energy", in m for "per_borehole_length")
    CAPEX_lin: the CAPEX values for the piecewise linear CAPEX curve (in CHF/m³ for "per_volume", in CHF/MWh for "per_energy", in CHF/m for "per_borehole_length")
    '''

    # Return error if technology is invalid
    if technology not in TECHNOLOGIES:
        raise ValueError(
            "Invalid technology. Please choose from "
            "'PTES', 'TTES', 'BTES', or 'ATES'."
        )

    # choose configuration based on technology
    config = TECHNOLOGIES[technology]

    # Read CAPEX data from CSV file
    CAPEX_DATA = pd.read_csv(CAPEX_DATA_PATH)

    # initialize variables
    CAPEX_lin = None
    capacity_lin = None
    CAPEX = None

    # Check if the unit is supported by the technology
    if f'CAPEX_{config.name}_{unit}' not in CAPEX_DATA.columns:
        raise ValueError(
            f"{technology} does not support unit '{unit}'."
        )
    
    # if yes, read the linearized CAPEX data for the given technology and unit
    else:

        # Read linearized CAPEX data for the given technology and unit
        CAPEX_lin = CAPEX_DATA[f"CAPEX_{config.name}_{unit}"].to_numpy()
        
        # Read capacity data for the given technology and unit
        if unit == "per_energy" or unit == "per_volume":
            capacity_lin = CAPEX_DATA[f"volume_{config.name}"].to_numpy()
        
        # if unit is per energy, convert capacity from energy to volume using the energy density
        if unit == "per_energy":
            # if no temperature range is given, use the default energy density for the technology, otherwise calculate the energy density based on the given temperature range
            if T_min is None or T_max is None:
                capacity_lin *= config.energy_density
                CAPEX_lin /= config.energy_density
            else:
                energy_density = ((T_max - T_min) * specific_heat_water((T_max + T_min) / 2) * density_water((T_max + T_min) / 2) / 3.6e9)
                capacity_lin *= energy_density
                CAPEX_lin /= energy_density

        # if unit is per borehole length, read the capacity data for the borehole length
        if unit == "per_borehole_length":
            capacity_lin = CAPEX_DATA[f"{config.name}_borehole_length"].to_numpy()

        # Check if the given capacity is within the range of the linearized data, if not, print a warning and use extrapolation
        if capacity < min(capacity_lin) and capacity != 0:
            print(f"Warning: capacity {capacity} is below the minimum capacity {min(capacity_lin)} for {technology} with unit '{unit}'. Extrapolation will be used.")
        
        if capacity > max(capacity_lin):
            print(f"Warning: capacity {capacity} is above the maximum capacity {max(capacity_lin)} for {technology} with unit '{unit}'. Extrapolation will be used.")

        # Interpolate CAPEX for the given capacity using the linearized data
        CAPEX = np.interp(capacity, capacity_lin, CAPEX_lin)

        if capacity != 0:
            return CAPEX
        else:
            return capacity_lin, CAPEX_lin