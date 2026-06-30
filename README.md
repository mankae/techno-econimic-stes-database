# techno-economic-stes-database

This repository was created during a semester project "Seasonal Thermal Energy Storage: Recent Techno-Economic Developments and Modeling Approaches" at Eidgenössische Technische Hochschule Zürich (ETHZ). The motivation of the project was the modeling of Seasonal Thermal Energy Storage (STES) in optimization-based energy system models (ESMs)

The repository contains CAPEX and OPEX data for PTES, TTES, BTES, and ATES. There are functions (OPEX_STES() and CAPEX_STES()) to get the CAPEX and OPEX value of these technologies depending on the type of technology and its capacity.

There is also a simulation framework for simple heat-loss modeling of PTES and TTES. The goal of the simulation is, to describe the heat loss with a parameter called self-discharge rate ($\eta_{\mathrm{self}}$). This parameter can be used in the simple storage model (SSM) which is often used in ESM:

$$Q_{\mathrm{sto,t+1}}=\eta_{\mathrm{self}} \cdot Q_{\mathrm{sto,t}}+\left(\eta_{\mathrm{ch}} \cdot \dot{Q}_{\mathrm{ch}} -\frac{\dot{Q}_{\mathrm{disch}}}{\eta_{\mathrm{disch}}}\right) \cdot \Delta t$$

The python package stes-tools contains the cost data function and also the simulation framework to simulate heat loss of PTES and TTES. The following functions are contained within the package:
- density_water(T): Density (rho) of water in kg/m^3 based on fluid temperature (T) nearest the flow meter in degrees Celsius
- specific_heat_water(T): Specific heat (cp) of water in J/(kg K) based on mean fluid temperature (T) in degrees Celsius
- CAPEX_STES(technology, unit, capacity, T_min, T_max): CAPEX for a given type of STES (PTES, TTES, BTES, ATES) based on the capacity and the unit of the capacity.
- OPEX_STES(technology): OPEX for a given type of STES (PTES, TTES, BTES, ATES)
- data_import(file_path): imports data for heat loss simulation
- PTES_geometry_import(file_path): imports the geometry data for heat loss simulation (specifically for PTES: truncated pyramid)
- STES: Is the class which is used to define the simulated storages. It is divided in the following two subclasses:
  - PTES(h, a, b, c, d, n_layers, T_min, T_max, T_ref): Define a PTES storage by defining the geometry (h, a, b, c, d), the number of layers, the temperature range (T_min to T_max) and the reference temperature
  - TTES(h, r, n_layers, T_min, T_max, T_ref): Define a TTES storage by defining the geometry (h, r), the number of layers, the temperature range (T_min to T_max) and the reference temperature
- temperature_map

## Examples
- [Cost function data treatment and example](notebooks/CAPEX_OPEX_database_STES.ipynb)
- [Heat loss simulation of PTES](notebooks/Heat_Loss_Simulation_of_PTES.ipynb)

```python
# density function
>>> rho = density_water(25)
>>> rho
997.0680068359376

# heat capacity function
>>> c_p = specific_heat_water(25)
>>> c_p
4181.562794921874
```
```python
# get the CAPEX and OPEX value of a PTES plant
>>> CAPEX = st.CAPEX_STES('PTES', 'per_volume', 70000, T_min=45, T_max=85) * 70000
>>> OPEX = CAPEX * st.OPEX_STES('PTES')
>>> print("CAPEX of a PTES with a volume of 70000 m^3 and temperature range from 45°C to 85°C:", round(CAPEX), "CHF, OPEX of the same PTES:", round(OPEX), "CHF/a")
CAPEX of a PTES with a volume of 70000 m^3 and temperature range from 45°C to 85°C: 4522103 CHF, OPEX of the same PTES: 44560 CHF/a
```

