# Cabin-Heatup-3D-Simulations
Contains the Java macro file that was used to run cabin heatup simulation of volvo xc40. With that it has the flow chart of the whole simulation how each step is being done. 

Initialization– The procedure start with initialization of initial and boundary condition such as initial temperatures of
domain, vehicle velocity, mass flow rate and temperature of air at cabin inlets and final physical time for the cabin
heatup.

Running Aero– Simulation is run on the outside of vehicle only for flow with no temperature data. It uses the script
AeroFlow. The simulation is run for 500 steps (1st order) then for 1000 steps (2nd order).

Running Airstream– Then for 4000 steps, simulation is being run with temperature solver added. At the end,
temperature data is being mapped from solids/airstream. CHT interfaces are setup between solids in contact with external
air (front/rear windshield, sunroof glass). At the end of this step, HTC data computed on external boundaries is being
mapped as convection boundary condition for solids (from solids/airstream (front/rear windshield, sunroof glass) to
solids/cabinair.

Run Defroster Flow- Simulation is being run in defroster ducts for 500 steps and at the end, velocity data is mapped at
the inlet of defroster.

Run cabin steady flow- Simulation is being run in cabin for 5500 steps. Only flow solver in the cabin to give us a rough
estimate of flow within a cabin. This will serve as a starting point for unsteady flow in cabin. The internal volume of air
in the cabin is modelled as a gas with constant density during domain initialization. During heating, internal volume is
modelled as an ideal gas. In order to save computational time, the solver is set to steady state first and then switched to
unsteady. Done only for cabin air without solids

Gravity- its added to the solver before cabin unsteady flow to account for the effects of buoyancy during heating while
the effects of solar radiation is neglected because of their low impact in cold weather.

Run Cabin Unsteady flow- The solver is switched to unsteady flow for both cabin air and solids, with a timestep of 0.02
secs. Once condition of steady state is being reached at monitor volumes for velocity a loop is being made.

First condition– is a time based stopping criteria to terminate the simulation if physical time is being reached (2400
secs)

Second condition- is to recompute the external flow and HTC on external interfaces every five minutes. In the
meantime we freeze the flow for 9 seconds with a timestep of 0.25sec and then unfreeze flow for a 1 seconds with a
timestep of 0.02secs. After 5 mins has reached we update the value of HTC on external interfaces if there is any change
and the loop continues until we achieve a final physical time of 2400secs. 
