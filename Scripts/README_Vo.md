RunSimulation.java: STAR-CCM+ Simulation Macro
This repository contains the RunSimulation.java script for automating cabin airflow and heat-up simulations in Simcenter STAR-CCM+. The macro is designed to streamline the simulation process, from initialization to post-processing, ensuring efficiency and accuracy in computational fluid dynamics (CFD) analysis.

Overview
The script is tailored for simulations involving:

Cabin Environment Studies:
Analyze thermal and airflow performance in vehicle cabins.
Solver Optimization:
Dynamic time-stepping for transient and steady-state phases.
Post-Processing Automation:
Extract meaningful results and generate reports.
Features
1. Mesh Quality Validation
Checks for minimum cell volume, quality, and face area connectivity.
Ensures mesh readiness for accurate simulation results.

2. Dynamic Solver Control
Automates adjustments to solver settings for transient and steady-state phases.
Improves simulation efficiency while maintaining accuracy.

3. Parameterized Simulation
Customizable inputs for:
Boundary conditions (e.g., mass flow rates, temperatures).
Simulation steps for different phases.

4. Post-Processing and Results Export
Automates the extraction of key metrics (e.g., temperature, velocity).
Saves processed data for external analysis.
Workflow

1. Package and Imports
Package:
package cabinStudy;
Groups the macro under a structured package for better organization.

Imports:
import star.common.*;
import java.io.*;
Includes STAR-CCM+ modules and Java utilities for simulation control and file handling.

2. Initialization
Fetches the active simulation instance:
simulation_0 = getActiveSimulation();
Sets up helper classes for modular operations:

Function function = new Function();
Collect collect = new Collect();

3. Simulation Configuration
Solver Settings:
double defrosterSteps = 500;
double airStreamSteps = 4000;
double cabSteadySteps = 5000;
Defines solver steps for transient and steady-state phases.

Initial and Boundary Conditions:
double U_ref = 50 / 3.6; // Reference velocity in m/s
double cabInitialTemp = 268.55; // Initial cabin temperature in Kelvin
Configures starting conditions for cabin and airflow properties.

Mass Flow Rates:
double defrosterMassFlow = 0.011416; // Defroster inlet flow rate in kg/s
Specifies HVAC component settings.

4. Mesh Validation
Defines criteria for mesh quality:
String invalidCriteria = " {MinCellVolume < 1e-14}";
Ensures the mesh meets quality standards for numerical accuracy.

5. Solver Execution
Implements dynamic time-stepping:
double freezeTimeStep = 0.25;
double unfreezeTimeStep = 0.04;
Adjusts solver settings for different simulation phases to balance accuracy and computational cost.

6. Post-Processing
Extracts key metrics:
Cabin temperature and velocity fields.
Pressure drops across HVAC components.
Generates visualizations:
Plots and 3D visualizations of airflow and temperature distributions.

7. Results Export
Saves processed data: 
simulation_0.exportData(filePath, format);
Allows seamless integration with external tools for further analysis.


How to Use
1. Import the Macro
Open Simcenter STAR-CCM+.
Navigate to the macros section and load RunSimulation.java.

2. Configure Parameters
Update boundary and initial conditions in the script to match your test case.
Modify solver settings for specific simulation requirements.

3. Execute the Macro
Run the macro within STAR-CCM+.
Monitor the simulation phases (e.g., defroster setup, steady-state).

4. Analyze Results
Review generated plots and reports.
Export data for external visualization or documentation.
Applications
Cabin Heat-Up Studies:
Analyze thermal comfort and energy efficiency in vehicle cabins.
Dynamic Solver Optimization:
Optimize transient and steady-state simulations.
Post-Processing Automation:
Automate data extraction and result generation.
Key Parameters
Parameter	Description
defrosterSteps	Solver steps for defroster setup phase.
airStreamSteps	Solver steps for airflow simulation.
cabSteadySteps	Solver steps for steady-state phase.
U_ref	Reference velocity in m/s.
cabInitialTemp	Initial cabin temperature in Kelvin.
defrosterMassFlow	Mass flow rate for defrosters.


## **Workflow**

### **1. Initialization**
- Retrieves the active STAR-CCM+ simulation instance.
- Loads custom helper classes (`Function`, `Collect`) for reusable operations.

### **2. Mesh Validation**
- Validates mesh quality using predefined criteria:
  - **Minimum Cell Volume**: Ensures adequate mesh resolution.
  - **Cell Quality**: Checks numerical quality for accurate results.
  - **Face Area Connectivity**: Validates surface continuity.

### **3. Simulation Configuration**
- **Solver Settings**:
  - Configures transient and steady-state solver steps:
    ```java
    double defrosterSteps = 500;
    double airStreamSteps = 4000;
    double cabSteadySteps = 5000;
    ```
- **Initial and Boundary Conditions**:
  - Sets initial temperatures, velocities, and mass flow rates:
    ```java
    double U_ref = 50 / 3.6; // Reference velocity (m/s)
    double defrosterMassFlow = 0.011416; // Defroster inlet mass flow (kg/s)
    ```

### **4. File Management**
- Reads external configuration files for flexible initialization:
  ```java
  File initialConditions = new File(resolvePath("/path/to/InitialConditions.csv"));
5. Execution
Automates simulation execution across multiple phases:
Defroster Setup
Air Stream Analysis
Cabin Steady-State Simulation
6. Post-Processing
Extracts and visualizes results:
Generates temperature, velocity, and pressure distributions.
Produces plots and reports for analysis.
7. Exporting Results
Saves processed data and visualizations for external use:
java
 
simulation_0.exportData(filePath, format);
Key Parameters
Parameter	Description
defrosterSteps	Solver steps for defroster setup phase.
airStreamSteps	Solver steps for airflow simulation.
cabSteadySteps	Solver steps for steady-state phase.
U_ref	Reference velocity in m/s.
cabInitialTemp	Initial cabin temperature in Kelvin.
defrosterMassFlow	Mass flow rate for defrosters.
How to Use
Import the Macro:

Load the RunSimulation.java file in Simcenter STAR-CCM+ under the macros section.
Configure Parameters:

Update boundary and initial conditions in the script to match your scenario.
Run the Script:

Execute the macro in STAR-CCM+ to start the automated workflow.
Analyze Results:

Review generated plots, reports, and exported data.
Applications
Cabin Heat-Up Analysis:
Study airflow and thermal performance in vehicle cabins.
Dynamic Solver Optimization:
Improve simulation efficiency with time-stepping control.
Post-Processing Automation:
Streamline result generation and analysis.
Acknowledgments
This macro leverages the Simcenter STAR-CCM+ API for CFD automation. It is tailored for cabin airflow studies but can be adapted to other simulation scenarios.