# techno-economic-stes-database

This repository was created during a semester project "Seasonal Thermal Energy Storage: Recent Techno-Economic Developments and Modeling Approaches" at Eidgenössische Technische Hochschule Zürich (ETHZ). The motivation of the project was the modeling of Seasonal Thermal Energy Storage (STES) in optimization-based energy system models (ESMs)

The repository contains CAPEX and OPEX data for PTES, TTES, BTES, and ATES. There are functions to get the CAPEX and OPEX value of these technologies depending on the type of technology and its capacity. 

There is also a simulation framework for simple heat-loss modeling of PTES and TTES. The goal of the simulation is, to describe the heat loss with a parameter called self-discharge rate ($\eta_{\mathrm{self}}$). This parameter can be used in the simple storage model (SSM):

$$Q_{\mathrm{sto,t+1}}=\eta_{\mathrm{self}} \cdot Q_{\mathrm{sto,t}}+\left(\eta_{\mathrm{ch}} \cdot \dot{Q}_{\mathrm{ch}} -\frac{\dot{Q}_{\mathrm{disch}}}{\eta_{\mathrm{disch}}}\right) \cdot \Delta t$$

## Examples
- [Cost function derivation](notebooks/CAPEX_OPEX_database_STES.ipynb)
- [Heat loss simulations of PTES](notebooks/Heat_Loss_Simulation_of_PTES.ipynb)
