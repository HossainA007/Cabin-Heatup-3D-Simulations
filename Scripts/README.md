# XC40 Cabin Heat-Up Simulation Optimization

This repository contains files, scripts, and documentation for optimizing cabin heat-up simulations for the Volvo XC40 in low-temperature environments. The project focused on reducing simulation time while maintaining accuracy, leveraging advanced techniques and solver optimizations.

---

## **Overview**

The project builds upon a prior study, focusing on reducing computational time from 115 hours to 23 hours while ensuring temperature deviations remain under ±1°C. Key improvements include:
1. Optimized mesh settings and solver parameters.
2. Advanced post-processing techniques for data analysis.
3. Trunk and vent optimization studies to minimize energy consumption.

---

## **Key Improvements**

1. **Larger Monitor Volumes**:
   - Replaced 10x10 cm volumes with 40x40 cm volumes for accurate temperature readings.
2. **Mean of Last N-Samples**:
   - Reduced noise in plots by using mean temperature over the last 2500 timesteps.
3. **Freezing and Unfreezing Flow**:
   - Reduced simulation time by freezing flow for 9 seconds and using dynamic timesteps.
4. **Improved Mesh**:
   - Refined prism layers and reduced total cells, ensuring better boundary layer resolution.

---

## **Workflow**

### **Simulation Framework**
1. **Initialization**:
   - Set boundary conditions: initial temperatures, mass flow rates, and vehicle velocity.
2. **Aero Flow**:
   - External airflow simulation for 500 steps.
3. **Air Stream**:
   - Simulated air interactions for 4000 steps.
4. **Defroster Flow**:
   - Ran defroster simulations for 500 steps.
5. **Cabin Steady Flow**:
   - Simulated cabin airflow for 5500 steps.
6. **Cabin Unsteady Flow**:
   - Dynamic unsteady simulation for 2400 seconds using timestep control.

### **Studies Conducted**
1. **Timestep Study**:
   - Adjusted unfreeze timestep to reduce simulation time by 30 hours with minimal deviation.
2. **Iteration Study**:
   - Reduced sub-iterations per timestep to save 31 hours.
3. **Mesh Study**:
   - Optimized surface sizes and reduced total simulation time to 23 hours.
4. **Prism Layer Study**:
   - Analyzed boundary resolution but discarded due to inaccuracies.
5. **Climate Refinement**:
   - Optimized vent jet profiles using cone and box refinement.
6. **Trunk Study**:
   - Evaluated energy consumption impacts of trunk size and luggage covers.

---

## **Results**
1. **Reduced Simulation Time**:
   - From 115 hours to 23 hours using 920 cores.
   - Further reduction possible to 13 hours with 1920 cores.
2. **Temperature Accuracy**:
   - Achieved deviations of less than 1°C across key regions.
3. **Energy Efficiency**:
   - Showed minimal heat loss through the trunk with a luggage cover.

---

## **Files and Scripts**

1. **Simulation Scripts**:
   - **`corrMesh_case12`**: Used for mesh optimization.
   - **`cabinHeatup_case12`**: Original case from Anandh's study.
   - **`Trunk1`**, **`Trunk2`**: Scripts for trunk analysis.

2. **Excel Data**:
   - Includes temperature deviations, solver timings, and mesh comparisons.

---

## **How to Use**

1. Clone the repository:
   ```bash
   git clone https://github.com/username/XC40-CabinHeatUp.git

Load simulation scripts in Simcenter STAR-CCM+.

Configure boundary conditions and solver settings as per your study.

Execute simulations and analyze results using provided scripts and datasets.

Acknowledgments
Special thanks to Randi Franzke and Emil Wileson for their guidance during this project. This study utilized computational resources at Volvo Cars for achieving optimal results.