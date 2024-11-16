// Simcenter STAR-CCM+ macro: RunSimulation.java
// Written by Simcenter STAR-CCM+ 16.04.007
package cabinStudy;

import classes.Function;
import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.metrics.*;
import java.io.*;
import star.mapping.*;
import star.segregatedflow.*;
import star.keturb.*;
import star.flow.*;
import classes.*;
import star.base.report.*;
import star.segregatedenergy.*;
import star.energy.*;
import star.material.*;
import star.post.*;
import star.vis.*;
import star.base.query.*;
import star.meshing.*;
import star.turbulence.*;

public class RunSimulation extends StarMacro {
Simulation simulation_0;
Function function = new Function();
Collect collect = new Collect();
String invalidCriteria = " \'functionEnabled\': true,"
            + " \'minimumDiscontiguousCellsEnabled\': true, \'minimumDiscontiguousCells\': 10,"
            + " \'minimumCellVolumeEnabled\': true, \'minimumCellVolume\': 0.0,"
            + " \'minimumVolumeChangeEnabled\': true, \'minimumVolumeChange\': 1.0E-4,"
            + " \'minimumCellQualityEnabled\': true, \'minimumCellQuality\': 1.0E-5,"
            + " \'minimumContiguousFaceAreaEnabled\': true,  \'minimumContiguousFaceArea\': 5.0E-10,"
            + " \'minimumFaceValidityEnabled\': true, \'minimumFaceValidity\': 0.9}";

String invalidCriteriaCab = "{ \'functionEnabled\': false,"
            + " \'minimumDiscontiguousCellsEnabled\': true, \'minimumDiscontiguousCells\': 1000,"
            + " \'minimumCellVolumeEnabled\': true, \'minimumCellVolume\': 0.0,"
            + " \'minimumVolumeChangeEnabled\': true, \'minimumVolumeChange\': 1.0E-3,"
            + " \'minimumCellQualityEnabled\': true, \'minimumCellQuality\': 1.0E-4,"
            + " \'minimumContiguousFaceAreaEnabled\': true,  \'minimumContiguousFaceArea\': 5.0E-10,"
            + " \'minimumFaceValidityEnabled\': true, \'minimumFaceValidity\': 0.98}";

String invalidCriteriaSolids = "{ \'functionEnabled\': false,"
            + " \'minimumDiscontiguousCellsEnabled\': true, \'minimumDiscontiguousCells\': 10,"
            + " \'minimumCellVolumeEnabled\': true, \'minimumCellVolume\': 0.0,"
            + " \'minimumVolumeChangeEnabled\': true, \'minimumVolumeChange\': 1.0E-3,"
            + " \'minimumCellQualityEnabled\': true, \'minimumCellQuality\': 1.0E-4,"
            + " \'minimumContiguousFaceAreaEnabled\': true,  \'minimumContiguousFaceArea\': 5.0E-10,"
            + " \'minimumFaceValidityEnabled\': true, \'minimumFaceValidity\': 0.95}";

@Override
@SuppressWarnings("empty-statement")
  public void execute() {
    simulation_0 = getActiveSimulation();

	// Simulation settings:
    double defrosterSteps = 500;
    double airStreamSteps = 4000;
    double cabSteadySteps = 5000;
    double totalTime = 2400.0;
    double tsUnsteady = 10.0;
    double freezeTime = 9.0;
    double unfreezeTime = 1.0;
    double freezeTimeStep = 0.25;
    double unfreezeTimeStep = 0.04;
    String solverType = "Second";
    String saveFile = "java"; //"sim" or "java"
    boolean gravity = true;
    boolean simHistory = true;
    boolean srh =true;
        // Boundary & Initial conditions:
    double U_ref = 50/3.6; //kph
    double T_ref = 266.92; //K
    double cabInitialTemp = 268.55; //K
    double solidsInitialTemp = 267.15; //K
    double tempInlet_HVAC = 303.15; //K
    
        // Mass flow rates
    double defrosterMassFlow = 0.011416; //kg/s
    double frontVentInlets = 0.009472; 
    double frontFloorInlets = 0.015138; 
    double rearFloorInlets = 0.009361;
    double rearVentInlets = 0.001361;
    boolean tableInputs=true;
    boolean insulation=true;
    String outletCondition = "pressure"; // or "outlet"
    
    
    String table = "";
    if (tableInputs) {
        table = "BaselineTest_2";
        File initialConditions = new File(resolvePath("/vcc/cae/nobackup/cfdcar/thermo/methods/Cabin-heat-up/Anandh/cabinInitial/measurements/InitialConditions.csv"));
        try{            
            BufferedReader in = new BufferedReader(new FileReader(initialConditions));
            simulation_0.println("==Reading Initial conditions file==");
            in.readLine();
            String row;
            int linecount = 0;
            while ((row = in.readLine())!=null){
                linecount++;
                String[] line = row.split(",");
                if(line[0].equals(table)){
                    simulation_0.println("==Setting initial conditions==");
                    cabInitialTemp = Double.parseDouble(line[1]);
                    solidsInitialTemp = Double.parseDouble(line[2]);
                    T_ref = Double.parseDouble(line[3]);
                }
            }
        }catch (IOException | NumberFormatException ex)
        {
             simulation_0.println(ex.getMessage());
        }
    }
    if(insulation){
        String insulationType = "BaselineTest_2";
        double R_th = 0.82;
        double[] insulations={0.0,0.0,0.0,0.0,0.0,0.0}; // {dashboard_tunnel, doorPanels, roof, seats, walls, uninsulated}
        switch (insulationType) {
            case "FullInsulation":
                insulations[0]=R_th;
                insulations[1]=R_th;
                insulations[2]=R_th;
                insulations[3]=R_th;
                insulations[4]=R_th;
                break;
            case "RoofInsulation":
                insulations[2] = R_th;
                break;
            case "DoorInsulation":
                insulations[1] = R_th;
                break;
            case "RearSeatInsulation": 
                insulations[3] = R_th;
                break;
            case "IPInsulation":
                insulations[0] = R_th;
                break;
            case "WallInsulation":
                insulations[4] = R_th;
                break;
            default:
                break;
        }
            
        insulateCabin(insulation,insulations);
    }
    
    //     Solver settings
    double VelURF = 0.7;
    double PreURF = 0.3;
    double TurbURF = 0.8;
    double VelTol = 0.1;
    double PreTol = 0.1;
    double PreSweeps = 2;
    String pAcc = "BiCGStab";
    double SolidURF = 0.99;
    double FluidURF = 0.9;
    
            
    // Running Script
    
    RemoveInvalidCells();
    SetInputs(U_ref, T_ref, defrosterMassFlow, frontFloorInlets, frontVentInlets, rearFloorInlets, rearVentInlets, tempInlet_HVAC, cabInitialTemp, solidsInitialTemp, tableInputs, table,outletCondition);
    UpdateInterfaces();  
    
    //External Flow
    RunAero();  //Flow Only
    RunAirStream(airStreamSteps); //Temperature solver added
    setSolverSettings(VelURF,PreURF,TurbURF,VelTol,PreTol,PreSweeps,pAcc,SolidURF,FluidURF);
    RunDefroster(defrosterSteps);   
    RunSteadyCabin(cabSteadySteps);
	
	/*SaveSim1(saveFile);*/
    
	SetGravity(gravity,outletCondition);
    simHistorySetup(simHistory,saveFile);
	
	/*SaveSim2(saveFile);*/
    
	RunUnsteady(tsUnsteady,unfreezeTimeStep,solverType,srh);
    RunAirStream(1000.0);
    RunFreezingFlow(totalTime,freezeTime,freezeTimeStep,unfreezeTime,unfreezeTimeStep,solverType);
   
	SaveSim2(saveFile);
    
   
  }
/////////////// //////////////////////////////////////////////// Running Script////////////////////////////////////////////////////////////////////
	  
//     1. Remove Invalid cells
//     2. Set the inputs
//     3. Update all interfaces
//     4. Run Aero
//     5. Run Airstream 
//     6. Set solver settings
//     7. Run defroster steps
//     8. Run Steady cabin
//    8a. Save simulation 
//     9. Set gravity
//    10. Create Simulation history
//   10a. Save Simulation
//    11. Run Unsteady
//    12. Run Air stream
//    13. Run Freezing Flow
//    14. Save simulation
	  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	  
//////////////////////////////////////////////////////////// STEP 1 - REMOVE INVALID CELLS /////////////////////////////////////////////////////////
    

    private void RemoveInvalidCells() {
       Collection<Region> regions = collect.collectRegions(simulation_0);
       Collection<Region> regionCabin = collect.collectRegions(regions, "cabinAir", true);
       Collection<Region> regionDoor = collect.collectRegions(regions, "doorAir", true);
       Collection<Region> regionSolids = collect.collectRegions(regions, "solids", true);
       simulation_0.println("============");
       simulation_0.println("Removing invalid cells");
       simulation_0.println("============");
       
       for (int i=1;i<=3;i++){
           simulation_0.getMeshManager().removeInvalidCells(new NeoObjectVector(regionCabin.toArray()), NeoProperty.fromString(invalidCriteriaCab));
           simulation_0.getMeshManager().removeInvalidCells(new NeoObjectVector(regionDoor.toArray()), NeoProperty.fromString(invalidCriteriaCab));
           simulation_0.getMeshManager().removeInvalidCells(new NeoObjectVector(regionSolids.toArray()), NeoProperty.fromString(invalidCriteriaSolids));
       }
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	  
/////////////////////////////////////////////////////////////////// STEP 2 - SET INPUTS /////////////////////////////////////////////////////////////////////////
	

    private void SetInputs(double U_ref, double T_ref, double defrosterMassFlow, double frontFloorInlets, double frontVentInlets, double rearFloorInlets, double rearVentInlets, double tempInlet_HVAC, double cabInitialTemp,double solidsInitialTemp, boolean tableInputs, String table,String outletCondition) {
        simulation_0.println("============");
        simulation_0.println("Setting inputs");
        simulation_0.println("============");
        
        Units units_0 = ((Units) simulation_0.getUnitsManager().getObject(""));
        Units units_1 = ((Units) simulation_0.getUnitsManager().getObject("K"));
        Units units_2 = ((Units) simulation_0.getUnitsManager().getObject("kg/s"));
        
        ScalarGlobalParameter scalarGlobalParameter_0 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("U_ref"));
        scalarGlobalParameter_0.getQuantity().setValue(U_ref);
        scalarGlobalParameter_0.getQuantity().setUnits(units_0);
        
        ScalarGlobalParameter scalarGlobalParameter_1 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("T_ref"));
        scalarGlobalParameter_1.getQuantity().setValue(T_ref);
        scalarGlobalParameter_1.getQuantity().setUnits(units_0);
        
        ScalarGlobalParameter scalarGlobalParameter_2 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("defrosterMassFlow"));
        scalarGlobalParameter_2.getQuantity().setValue(defrosterMassFlow);
        scalarGlobalParameter_2.getQuantity().setUnits(units_2);
        
        ScalarGlobalParameter scalarGlobalParameter_3 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("frontFloorInlets"));
        scalarGlobalParameter_3.getQuantity().setValue(frontFloorInlets);
        scalarGlobalParameter_3.getQuantity().setUnits(units_2);
        
        ScalarGlobalParameter scalarGlobalParameter_4 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("frontVentInlets"));
        scalarGlobalParameter_4.getQuantity().setValue(frontVentInlets);
        scalarGlobalParameter_4.getQuantity().setUnits(units_2);
        
        ScalarGlobalParameter scalarGlobalParameter_5 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("rearVentInlets"));
        scalarGlobalParameter_5.getQuantity().setValue(rearVentInlets);
        scalarGlobalParameter_5.getQuantity().setUnits(units_2);
        
        ScalarGlobalParameter scalarGlobalParameter_6 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("rearFloorInlets"));
        scalarGlobalParameter_6.getQuantity().setValue(rearFloorInlets);
        scalarGlobalParameter_6.getQuantity().setUnits(units_2);
        
        ScalarGlobalParameter scalarGlobalParameter_7 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("tempInlet_HVAC"));
        scalarGlobalParameter_7.getQuantity().setValue(tempInlet_HVAC);
        scalarGlobalParameter_7.getQuantity().setUnits(units_1);
        
        simulation_0.get(GlobalParameterManager.class).createGlobalParameter(ScalarGlobalParameter.class, "cabInitialTemp");
        ScalarGlobalParameter scalarGlobalParameter_9 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("cabInitialTemp"));
        scalarGlobalParameter_9.getQuantity().setValue(cabInitialTemp);
        scalarGlobalParameter_9.getQuantity().setUnits(units_0);
        
        simulation_0.get(GlobalParameterManager.class).createGlobalParameter(ScalarGlobalParameter.class, "solidsInitialTemp");
        ScalarGlobalParameter scalarGlobalParameter_10 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("solidsInitialTemp"));
        scalarGlobalParameter_10.getQuantity().setValue(solidsInitialTemp);
        scalarGlobalParameter_10.getQuantity().setUnits(units_0);
        
        PhysicsContinuum physicsContinuum_0 = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("cabinAir"));
        StaticTemperatureProfile staticTemperatureProfile_0 = physicsContinuum_0.getInitialConditions().get(StaticTemperatureProfile.class);
        staticTemperatureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("$cabInitialTemp");
        staticTemperatureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_0);
        
        PhysicsContinuum physicsContinuum_1 = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Solids"));
        StaticTemperatureProfile staticTemperatureProfile_1 = physicsContinuum_1.getInitialConditions().get(StaticTemperatureProfile.class);
        staticTemperatureProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("$solidsInitialTemp");
        staticTemperatureProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_0);
        
        
        if (outletCondition=="outlet"){
            PhysicsContinuum cabinAir = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("cabinAir")); 
            Region CabinRegion = ((Region) simulation_0.getRegionManager().getRegion("cabinAir"));
            Boundary Outlet = CabinRegion.getBoundaryManager().getBoundary("outlets");
            Outlet.setBoundaryType(OutletBoundary.class);
            
        }
        
        
        if(tableInputs){
            FileTable fileTable_0 = (FileTable) simulation_0.getTableManager().createFromFile(resolvePath("/vcc/cae/nobackup/cfdcar/thermo/methods/Cabin-heat-up/Anandh/cabinInitial/measurements/"+table+".csv"));
            simulation_0.println("Loading Table data ===> "+table);
            Region cabinAir = simulation_0.getRegionManager().getRegion("cabinAir");

            //Defroster Inlet:
            Boundary defrosterInlet = cabinAir.getBoundaryManager().getBoundary("DefrosterInlet");
            //MassFlowBoundary massFlowBoundary_0 = ((MassFlowBoundary) simulation_0.get(ConditionTypeManager.class).get(MassFlowBoundary.class));
            //defrosterInlet.setBoundaryType(massFlowBoundary_0);
            
            //Setting up subgroups
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("front")).setName("frontLeft");
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("frontLeft")).setQuery(null);
            MeshOperationPart meshOperationPart_0 = ((MeshOperationPart) simulation_0.get(SimulationPartManager.class).getPart("cabinAir"));
            PartSurface partSurface_0 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.ClosingSurfaces.defroster-Frontleft"));
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("frontLeft")).setObjects(partSurface_0);
            ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 1")).setName("frontRight");
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("frontRight")).setQuery(null);
            PartSurface partSurface_1 =  ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.ClosingSurfaces.defroster-Frontright"));
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("frontRight")).setObjects(partSurface_1);
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("side")).setName("sideLeft");
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("sideLeft")).setQuery(null);
            PartSurface partSurface_2 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.ClosingSurfaces.defrosterSideLeft-inlet"));
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("sideLeft")).setObjects(partSurface_2);
            ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 1")).setName("sideRight");
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("sideRight")).setQuery(null);
            PartSurface partSurface_3 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.ClosingSurfaces.defroster-SideRight-inlet"));
            ((NamedPartGroup) ((PartGrouping) defrosterInlet.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("sideRight")).setObjects(partSurface_3);

            
            //Setup Mass Flows:
            /*MassFlowRateProfile massFlowRateProfile_0 = defrosterInlet.getValues().get(MassFlowRateProfile.class);

            ProxyProfile proxyProfile_0 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("frontLeft"));
            proxyProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.36*${defrosterMassFlow}");
            proxyProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_1 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("frontRight"));
            proxyProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.36*$defrosterMassFlow");
            proxyProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_2 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("sideLeft"));
            proxyProfile_2.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.05*$defrosterMassFlow");
            proxyProfile_2.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_3 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("sideRight"));
            proxyProfile_3.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.05*$defrosterMassFlow");
            proxyProfile_3.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ScalarGlobalParameter scalarGlobalParameter_8 = ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("windshieldFlow"));
            scalarGlobalParameter_8.getQuantity().setDefinition("0.18*${defrosterMassFlow}");
            scalarGlobalParameter_8.getQuantity().setUnits(units_0);

            ProxyProfile proxyProfile_4 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_1"));
            proxyProfile_4.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.15*${windshieldFlow}");
            proxyProfile_4.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_5 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_2"));
            proxyProfile_5.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.16*${windshieldFlow}");
            proxyProfile_5.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_6 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_3"));
            proxyProfile_6.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.10*${windshieldFlow}");
            proxyProfile_6.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_7 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_4"));
            proxyProfile_7.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.11*${windshieldFlow}");
            proxyProfile_7.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_8 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_5"));
            proxyProfile_8.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.10*${windshieldFlow}");
            proxyProfile_8.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_9 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_6"));
            proxyProfile_9.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.09*${windshieldFlow}");
            proxyProfile_9.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_10 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_7"));
            proxyProfile_10.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.08*${windshieldFlow}");
            proxyProfile_10.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_11 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_8"));
            proxyProfile_11.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.07*${windshieldFlow}");
            proxyProfile_11.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_12 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_9"));
            proxyProfile_12.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.07*${windshieldFlow}");
            proxyProfile_12.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_13 = ((ProxyProfile) massFlowRateProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_10"));
            proxyProfile_13.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.07*${windshieldFlow}");
            proxyProfile_13.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            //Setup Components
            UserFieldFunction userFieldFunction_1 = simulation_0.getFieldFunctionManager().createFieldFunction();
            userFieldFunction_1.getTypeOption().setSelected(FieldFunctionTypeOption.Type.VECTOR);
            userFieldFunction_1.setPresentationName("defrosterFlowComponents");
            userFieldFunction_1.setDefinition("$${MappedVertexVelocity}/mag($${MappedVertexVelocity})");
            userFieldFunction_1.setFunctionName("defrosterFlowComponents");
            
            defrosterInlet.getConditions().get(FlowDirectionOption.class).setSelected(FlowDirectionOption.Type.COMPONENTS);
            FlowDirectionProfile flowDirectionProfile_0 = defrosterInlet.getValues().get(FlowDirectionProfile.class);
            flowDirectionProfile_0.setMethod(FunctionVectorProfileMethod.class);
            //PrimitiveFieldFunction primitiveFieldFunction_0 = ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("MappedVertexVelocity"));
            flowDirectionProfile_0.getMethod(FunctionVectorProfileMethod.class).setFieldFunction(userFieldFunction_1);*/

            //Setup Temperature
            StaticTemperatureProfile totalTemperatureProfile_0 = defrosterInlet.getValues().get(StaticTemperatureProfile.class);
            totalTemperatureProfile_0.setMethod(ByPartProfileMethod.class);
            
            ProxyProfile proxyProfile_14 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("frontLeft"));
            proxyProfile_14.setMethod(TimeTabularScalarProfileMethod.class);
            proxyProfile_14.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_14.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterOutL1");

            ProxyProfile proxyProfile_15 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("frontRight"));
            proxyProfile_15.setMethod(TimeTabularScalarProfileMethod.class);
            proxyProfile_15.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_15.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterOutR1");

            ProxyProfile proxyProfile_16 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("sideLeft"));
            proxyProfile_16.setMethod(TimeTabularScalarProfileMethod.class);
            proxyProfile_16.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_16.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");

            ProxyProfile proxyProfile_17 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("sideRight"));
            proxyProfile_17.setMethod(TimeTabularScalarProfileMethod.class);
            proxyProfile_17.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_17.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutR");

            ProxyProfile proxyProfile_18 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_1"));
            proxyProfile_18.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_19 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_2"));
            proxyProfile_19.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_20 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_3"));
            proxyProfile_20.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_21 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_4"));
            proxyProfile_21.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_22 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_5"));
            proxyProfile_22.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_23 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_6"));
            proxyProfile_23.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_24 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_7"));
            proxyProfile_24.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_25 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_8"));
            proxyProfile_25.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_26 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_9"));
            proxyProfile_26.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_27 = ((ProxyProfile) totalTemperatureProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Ws_10"));
            proxyProfile_27.setMethod(TimeTabularScalarProfileMethod.class);

            proxyProfile_18.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_19.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_20.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_21.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_22.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_23.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_24.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_25.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_26.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_27.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);

            proxyProfile_18.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            proxyProfile_19.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            proxyProfile_20.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            proxyProfile_21.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            proxyProfile_22.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            proxyProfile_23.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            proxyProfile_24.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            proxyProfile_25.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            proxyProfile_26.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            proxyProfile_27.getMethod(TimeTabularScalarProfileMethod.class).setData("DefrosterSideOutL");
            
            //Front Floor inlets
            Boundary boundary_1 = cabinAir.getBoundaryManager().getBoundary("frontFloorInlets");
            boundary_1.setAllowPerPartValues(true);

            ((PartGrouping) boundary_1.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();
            ((PartGrouping) boundary_1.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();
            
            ((NamedPartGroup) ((PartGrouping) boundary_1.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setQuery(null);            
            PartSurface partSurface_4 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.Closing_surfaces.Driver-foot-inlet"));
            ((NamedPartGroup) ((PartGrouping) boundary_1.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setObjects(partSurface_4);

            ((NamedPartGroup) ((PartGrouping) boundary_1.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 3")).setQuery(null);
            PartSurface partSurface_5 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.Closing_surfaces.Passenger-foot-inlet"));
            ((NamedPartGroup) ((PartGrouping) boundary_1.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 3")).setObjects(partSurface_5);

            MassFlowRateProfile massFlowRateProfile_1 = boundary_1.getValues().get(MassFlowRateProfile.class);
            massFlowRateProfile_1.setMethod(ByPartProfileMethod.class);

            ProxyProfile proxyProfile_28 = ((ProxyProfile) massFlowRateProfile_1.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 2"));
            proxyProfile_28.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.4972*$frontFloorInlets");
            proxyProfile_28.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_29 = ((ProxyProfile) massFlowRateProfile_1.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 3"));
            proxyProfile_29.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.5028*$frontFloorInlets");
            proxyProfile_29.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            TotalTemperatureProfile totalTemperatureProfile_1 = boundary_1.getValues().get(TotalTemperatureProfile.class);
            totalTemperatureProfile_1.setMethod(ByPartProfileMethod.class);

            ProxyProfile proxyProfile_30 = ((ProxyProfile) totalTemperatureProfile_1.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 2"));
            proxyProfile_30.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_31 = ((ProxyProfile) totalTemperatureProfile_1.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 3"));
            proxyProfile_31.setMethod(TimeTabularScalarProfileMethod.class);

            proxyProfile_30.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_31.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);

            proxyProfile_30.getMethod(TimeTabularScalarProfileMethod.class).setData("FloorOut1L1");
            proxyProfile_31.getMethod(TimeTabularScalarProfileMethod.class).setData("FloorOut1R1");

            //Front vent inlets
            Boundary boundary_2 = cabinAir.getBoundaryManager().getBoundary("frontVentInlet");
            boundary_2.setAllowPerPartValues(true);

            ((PartGrouping) boundary_2.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();
            ((PartGrouping) boundary_2.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();
            ((PartGrouping) boundary_2.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();

            ((NamedPartGroup) ((PartGrouping) boundary_2.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setQuery(null);
            PartSurface partSurface_6 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.Closing_surfaces.Front-left-inlet"));
            ((NamedPartGroup) ((PartGrouping) boundary_2.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setObjects(partSurface_6);

            ((NamedPartGroup) ((PartGrouping) boundary_2.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 3")).setQuery(null);
            PartSurface partSurface_7 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.Closing_surfaces.Front-right-inlet"));
            ((NamedPartGroup) ((PartGrouping) boundary_2.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 3")).setObjects(partSurface_7);

            ((NamedPartGroup) ((PartGrouping) boundary_2.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 4")).setQuery(null);
            PartSurface partSurface_8 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.Closing_surfaces.Front-center-inlet"));
            ((NamedPartGroup) ((PartGrouping) boundary_2.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 4")).setObjects(partSurface_8);

            MassFlowRateProfile massFlowRateProfile_2 = boundary_2.getValues().get(MassFlowRateProfile.class);
            massFlowRateProfile_2.setMethod(ByPartProfileMethod.class);

            ProxyProfile proxyProfile_32 = ((ProxyProfile) massFlowRateProfile_2.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 2"));
            proxyProfile_32.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.4105*$frontVentInlets");
            proxyProfile_32.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_33 = ((ProxyProfile) massFlowRateProfile_2.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 3"));
            proxyProfile_33.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.4223*$frontVentInlets");
            proxyProfile_33.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_34 = ((ProxyProfile) massFlowRateProfile_2.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 4"));
            proxyProfile_34.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.1672*$frontVentInlets");
            proxyProfile_34.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            TotalTemperatureProfile totalTemperatureProfile_2 = boundary_2.getValues().get(TotalTemperatureProfile.class);
            totalTemperatureProfile_2.setMethod(ByPartProfileMethod.class);

            ProxyProfile proxyProfile_35 = ((ProxyProfile) totalTemperatureProfile_2.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 2"));
            proxyProfile_35.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_36 = ((ProxyProfile) totalTemperatureProfile_2.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 3"));
            proxyProfile_36.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_37 = ((ProxyProfile) totalTemperatureProfile_2.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 4"));
            proxyProfile_37.setMethod(TimeTabularScalarProfileMethod.class);

            proxyProfile_35.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_36.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_37.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);

            proxyProfile_35.getMethod(TimeTabularScalarProfileMethod.class).setData("VentOut1LO");
            proxyProfile_36.getMethod(TimeTabularScalarProfileMethod.class).setData("VentOut1RO");
            proxyProfile_37.getMethod(TimeTabularScalarProfileMethod.class).setData("VentOut1LC");

            //Rear Floor Inlets
            Boundary boundary_3 = cabinAir.getBoundaryManager().getBoundary("rearFootInlets");
            boundary_3.setAllowPerPartValues(true);

            ((PartGrouping) boundary_3.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();
            ((PartGrouping) boundary_3.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();

            ((NamedPartGroup) ((PartGrouping) boundary_3.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setQuery(null);
            PartSurface partSurface_9 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.Closing_surfaces.Rear-left-inlet"));
            ((NamedPartGroup) ((PartGrouping) boundary_3.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setObjects(partSurface_9);

            ((NamedPartGroup) ((PartGrouping) boundary_3.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 3")).setQuery(null);
            PartSurface partSurface_10 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Subtracts.Closing_surfaces.Rear-right-inlet"));
            ((NamedPartGroup) ((PartGrouping) boundary_3.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 3")).setObjects(partSurface_10);

            MassFlowRateProfile massFlowRateProfile_3 = boundary_3.getValues().get(MassFlowRateProfile.class);
            massFlowRateProfile_3.setMethod(ByPartProfileMethod.class);

            ProxyProfile proxyProfile_38 = ((ProxyProfile) massFlowRateProfile_3.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 2"));
            proxyProfile_38.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.4713*$rearFloorInlets");//0.4713
            proxyProfile_38.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            ProxyProfile proxyProfile_39 = ((ProxyProfile) massFlowRateProfile_3.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 3"));
            proxyProfile_39.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("0.5287*$rearFloorInlets");//0.5287
            proxyProfile_39.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_2);

            TotalTemperatureProfile totalTemperatureProfile_3 = boundary_3.getValues().get(TotalTemperatureProfile.class);
            totalTemperatureProfile_3.setMethod(ByPartProfileMethod.class);

            ProxyProfile proxyProfile_40 = ((ProxyProfile) totalTemperatureProfile_3.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 2"));
            proxyProfile_40.setMethod(TimeTabularScalarProfileMethod.class);

            ProxyProfile proxyProfile_41 = ((ProxyProfile) totalTemperatureProfile_3.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 3"));
            proxyProfile_41.setMethod(TimeTabularScalarProfileMethod.class);

            proxyProfile_40.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            proxyProfile_41.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);

            proxyProfile_40.getMethod(TimeTabularScalarProfileMethod.class).setData("FloorOut2L1");
            proxyProfile_41.getMethod(TimeTabularScalarProfileMethod.class).setData("FloorOut2R1");

            //Rear Vent Inlets - Single boundary
            Boundary boundary_4 = cabinAir.getBoundaryManager().getBoundary("rearVentInlets");
            TotalTemperatureProfile totalTemperatureProfile_4 = boundary_4.getValues().get(TotalTemperatureProfile.class);
            totalTemperatureProfile_4.setMethod(TimeTabularScalarProfileMethod.class);
            totalTemperatureProfile_4.getMethod(TimeTabularScalarProfileMethod.class).setTable(fileTable_0);
            totalTemperatureProfile_4.getMethod(TimeTabularScalarProfileMethod.class).setData("VentOut2LC");
        }
	}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////// STEP 3 - UPDATE INTERFACES ///////////////////////////////////////////////////////////////////////////////////
		
	
        private void UpdateInterfaces(){
        simulation_0.println("============");
        simulation_0.println("Updating Interface tolerance");
        simulation_0.println("============");
        BoundaryInterface car_AirStream = ((BoundaryInterface) simulation_0.getInterfaceManager().getInterface("solids/AirStream"));
        MappedInterfaceToleranceCondition mappedInterfaceToleranceCondition_0 = car_AirStream.getValues().get(MappedInterfaceToleranceCondition.class);
        MappedInterfaceToleranceConditionLeaf mappedInterfaceToleranceConditionLeaf_0 = mappedInterfaceToleranceCondition_0.getModelPartValue();
        mappedInterfaceToleranceConditionLeaf_0.setProximityTolerance(0.05);
    }
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
///////////////////////////////////////////////////////// STEP 4 - RUN AERO FLOW ////////////////////////////////////////////////////////////////////////////////////////////////
    
	
	private void RunAero() {
        File runAeroMethod = new File("/vcc/cae/nobackup/cfdcar/thermo/methods/Cabin-heat-up/Anandh/cabinInitial/STARCFD/src/cabinStudy/RunAero.java");
        simulation_0.println("============");
        simulation_0.println("Running Aero method");
        simulation_0.println("============");
        //Load physics
        PhysicsContinuum cabinAir = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("cabinAir"));
        PhysicsContinuum Solids = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Solids"));
        PhysicsContinuum defrosterFlow = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("defrosterFlow"));
        PhysicsContinuum AirStream = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));
        
        
        simulation_0.println("Running Flow only");
        SegregatedFluidTemperatureModel segregatedFluidTemperatureModel_0 =  AirStream.getModelManager().getModel(SegregatedFluidTemperatureModel.class);
        AirStream.disableModel(segregatedFluidTemperatureModel_0);
        
        //Load stopping criteria 
        StepStoppingCriterion steadySteps = ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));
        MonitorIterationStoppingCriterion tempCriterion = ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridTemperatureSum Monitor Criterion"));
        MonitorIterationStoppingCriterion velCriterion =  ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridVelocitySum Monitor Criterion"));
        
        //Set physics
        cabinAir.setIsActive(false);
        Solids.setIsActive(false);
        defrosterFlow.setIsActive(false);
        AirStream.setIsActive(true);
        
        // Aero method only steady
        steadySteps.setIsUsed(true);
        tempCriterion.setIsUsed(false);
        velCriterion.setIsUsed(false);
        
        //Run method
        try{
        new StarScript(simulation_0, runAeroMethod).play();
        } catch (NullPointerException e) {
            simulation_0.println("Can't Run Aero method.. *Sad smiley");
        }
        //Turn Off Physics 1
        AirStream.setIsActive(false);
        steadySteps.setIsUsed(false);
    
    
    }

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	///////////////////////////////////////////////////////////////////// STEP 5 - RUN AIRSTREAM ///////////////////////////////////////////////////////////////////////////////////////
		
                                                                        //Temperature solver added
	
	    private void RunAirStream(double airStreamSteps) {
        simulation_0.println("============");
        simulation_0.println("Running AirStream Steady for "+airStreamSteps+" iterations");
        simulation_0.println("============");

        //Load physics
        PhysicsContinuum cabinAir = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("cabinAir"));
        PhysicsContinuum Solids = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Solids"));
        PhysicsContinuum defrosterFlow = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("defrosterFlow"));
        PhysicsContinuum AirStream = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));
        
        simulation_0.println("Running Second order temperature model");
        AirStream.enable(SegregatedFluidTemperatureModel.class);
        SegregatedFluidTemperatureModel segregatedFluidTemperatureModel_0 =  AirStream.getModelManager().getModel(SegregatedFluidTemperatureModel.class);
        segregatedFluidTemperatureModel_0.getConvectiveFluxOption().setSelected(ConvectiveFluxOption.Type.SECOND_ORDER);
        
        //MonitorIterationStoppingCriterion monitorUMax = (MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("UMax Criterion");
        Collection<Region> regions = collect.collectRegions(simulation_0);
        Collection<Region> regionNoCore = collect.collectRegions(regions, "^Core.*", false);
        
        //Load stopping criteria 
        StepStoppingCriterion steadySteps = ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));
        MonitorIterationStoppingCriterion tempCriterion = ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridTemperatureSum Monitor Criterion"));
        MonitorIterationStoppingCriterion velCriterion =  ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridVelocitySum Monitor Criterion"));
        MonitorIterationStoppingCriterion monitorUMax = (MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("UMax Criterion");
        MonitorIterationStoppingCriterion monitorTVRMax = (MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("TVRMax Criterion");
        
        //Set physics
        cabinAir.setIsActive(false);
        Solids.setIsActive(false);
        defrosterFlow.setIsActive(false);
        AirStream.setIsActive(true);
        
        // Aero method only steady
        steadySteps.setIsUsed(true);
        tempCriterion.setIsUsed(false);
        velCriterion.setIsUsed(false);
        
        //Enable UMax and TVRMax monitors
        monitorUMax.setIsUsed(true);
        monitorTVRMax.setIsUsed(true);
        
        steadySteps.setMaximumNumberSteps((int) (simulation_0.getSimulationIterator().getCurrentIteration() + airStreamSteps));
        while(simulation_0.getSimulationIterator().getCurrentIteration() < steadySteps.getMaximumNumberSteps()){
            if(monitorUMax.getIsSatisfied()) {
                for (int i = 1; i < 5; i++) {
                    simulation_0.getMeshManager().removeInvalidCells(new NeoObjectVector(regionNoCore.toArray()), NeoProperty.fromString("{\'function\': \'VelocityMagnitude\', \'functionValue\': " + ((MonitorIterationStoppingCriterionMaxLimitType) monitorUMax.getCriterionType()).getLimit().evaluate()+ ", \'functionOperator\': 1," + invalidCriteria));
                }
                }
            if(monitorTVRMax.getIsSatisfied()) {
                for (int i = 1; i < 5; i++) {
                    simulation_0.getMeshManager().removeInvalidCells(new NeoObjectVector(regionNoCore.toArray()), NeoProperty.fromString("{\'function\': \'TurbulentViscosityRatio\', \'functionValue\': " + ((MonitorIterationStoppingCriterionMaxLimitType) monitorTVRMax.getCriterionType()).getLimit().evaluate()+ ", \'functionOperator\': 1," + invalidCriteria));
                }
                }
			
			
			simulation_0.getSimulationIterator().run();
            
            
        } 
        
        //Turn Off Physics 1
        AirStream.setIsActive(false);
        steadySteps.setIsUsed(false);
        monitorUMax.setIsUsed(false);
        monitorTVRMax.setIsUsed(false);
    }
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	//////////////////////////////////////////////////////////////////////////////// STEP 5 - SET SOLVER SETTINGS ///////////////////////////////////////////////////////////////////////////////////
    
    
	 private void setSolverSettings(double VelURF, double PreURF, double TurbURF, double VelTol, double PreTol, double PreSweeps, String pAcc, double SolidURF, double FluidURF) {
		 
        SegregatedFlowSolver segregatedFlowSolver_0 = ((SegregatedFlowSolver) simulation_0.getSolverManager().getSolver(SegregatedFlowSolver.class));
        VelocitySolver velocitySolver_0 = segregatedFlowSolver_0.getVelocitySolver();
        AMGLinearSolver aMGLinearSolver_1 = velocitySolver_0.getAMGLinearSolver();
        PressureSolver pressureSolver_0 = segregatedFlowSolver_0.getPressureSolver();
        AMGLinearSolver aMGLinearSolver_2 = pressureSolver_0.getAMGLinearSolver();
        SegregatedEnergySolver segregatedEnergySolver_0 = ((SegregatedEnergySolver) simulation_0.getSolverManager().getSolver(SegregatedEnergySolver.class));
        KeTurbSolver keTurbSolver_0 = ((KeTurbSolver) simulation_0.getSolverManager().getSolver(KeTurbSolver.class));
        Units units_2 = ((Units) simulation_0.getUnitsManager().getObject(""));
        
        velocitySolver_0.getUrfQuantity().setValue(VelURF);
        velocitySolver_0.getUrfQuantity().setUnits(units_2);

        aMGLinearSolver_1.setConvergeTol(VelTol);

        pressureSolver_0.getUrfQuantity().setValue(PreURF);
        pressureSolver_0.getUrfQuantity().setUnits(units_2);
 
        aMGLinearSolver_2.getAccelerationOption().setSelected(AMGAccelerationOption.Type.valueOf(pAcc));
        aMGLinearSolver_2.setConvergeTol(PreTol);
        
        ((AMGVCycle) aMGLinearSolver_2.getCycleType()).setPreSweeps((int) PreSweeps);
        ((AMGVCycle) aMGLinearSolver_2.getCycleType()).setPostSweeps((int) PreSweeps);

        segregatedEnergySolver_0.getFluidUrfQuantity().setValue(FluidURF);
        segregatedEnergySolver_0.getFluidUrfQuantity().setUnits(units_2);

        segregatedEnergySolver_0.getSolidUrfQuantity().setValue(SolidURF);
        segregatedEnergySolver_0.getSolidUrfQuantity().setUnits(units_2);

        keTurbSolver_0.getUrfQuantity().setValue(TurbURF);
        keTurbSolver_0.getUrfQuantity().setUnits(units_2);
    }
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
		///////////////////////////////////////////////////////////////////////////// STEP 6 - RUN DEFROSTER ///////////////////////////////////////////////////////////////////////////////////////
		  
    
	private void RunDefroster(double defrosterSteps) {
        simulation_0.println("============");
        simulation_0.println("Running defroster");
        simulation_0.println("============");
        
        //Load physics
        PhysicsContinuum cabinAir = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("cabinAir"));
        PhysicsContinuum Solids = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Solids"));
        PhysicsContinuum defrosterFlow = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("defrosterFlow"));
        PhysicsContinuum AirStream = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));
        
        //Load stopping criteria 
        StepStoppingCriterion steadySteps = ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));
        MonitorIterationStoppingCriterion tempCriterion = ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridTemperatureSum Monitor Criterion"));
        MonitorIterationStoppingCriterion velCriterion =  ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridVelocitySum Monitor Criterion"));
        
        //Set physics
        cabinAir.setIsActive(false);
        Solids.setIsActive(false);
        defrosterFlow.setIsActive(true);
        AirStream.setIsActive(false);
        
        //Only Steady
        steadySteps.setIsUsed(true);
        tempCriterion.setIsUsed(false);
        velCriterion.setIsUsed(false);
        
        steadySteps.setMaximumNumberSteps((int) (simulation_0.getSimulationIterator().getCurrentIteration() + defrosterSteps));
        simulation_0.getSimulationIterator().run();
        
        //Load surface mapper
        SurfaceDataMapper surfaceDataMapper_0 = ((SurfaceDataMapper) simulation_0.get(DataMapperManager.class).getObject("Surface Data Mapper 1"));
        surfaceDataMapper_0.mapData();
        
        //Turn Off defroster
        defrosterFlow.setIsActive(false);
        steadySteps.setIsUsed(false);      
        

	}

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
		/////////////////////////////////////////////////////////////// STEP 7 - RUN STEADY CABIN //////////////////////////////////////////////////////////////////////////////////////////////////
		
    
		 private void RunSteadyCabin(double cabSteadySteps) {
        simulation_0.println("============");
        simulation_0.println("Running Cab Steady: Flow Only");
        simulation_0.println("============");
        //Load physics
        PhysicsContinuum cabinAir = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("cabinAir"));
        PhysicsContinuum Solids = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Solids"));
        PhysicsContinuum defrosterFlow = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("defrosterFlow"));
        PhysicsContinuum AirStream = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));
        
        
        simulation_0.println("Running Flow only");
        //Disable gravity and temperature
        SegregatedFluidTemperatureModel segregatedFluidTemperatureModel_0 =  cabinAir.getModelManager().getModel(SegregatedFluidTemperatureModel.class);
        cabinAir.disableModel(segregatedFluidTemperatureModel_0);
        
        IdealGasModel idealGasModel_0 = cabinAir.getModelManager().getModel(IdealGasModel.class);
        cabinAir.disableModel(idealGasModel_0);
        cabinAir.enable(ConstantDensityModel.class);
        ImplicitUnsteadyModel implicitUnsteadyModel_0 = cabinAir.getModelManager().getModel(ImplicitUnsteadyModel.class);
        cabinAir.disableModel(implicitUnsteadyModel_0);
        cabinAir.enable(SteadyModel.class);
        
        Units units_0 = simulation_0.getUnitsManager().getInternalUnits(new IntVector(new int[] {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0}));
        SingleComponentGasModel singleComponentGasModel_0 = cabinAir.getModelManager().getModel(SingleComponentGasModel.class);

        Gas gas_0 = ((Gas) singleComponentGasModel_0.getMaterial());
        ConstantMaterialPropertyMethod constantMaterialPropertyMethod_0 = ((ConstantMaterialPropertyMethod) gas_0.getMaterialProperties().getMaterialProperty(ConstantDensityProperty.class).getMethod());

        constantMaterialPropertyMethod_0.getQuantity().setDefinition("101325/287.058/$cabInitialTemp");
        constantMaterialPropertyMethod_0.getQuantity().setUnits(units_0);
        
        //Load stopping criteria 
        StepStoppingCriterion steadySteps = ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));
        MonitorIterationStoppingCriterion tempCriterion = ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridTemperatureSum Monitor Criterion"));
        MonitorIterationStoppingCriterion velCriterion =  ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridVelocitySum Monitor Criterion"));
        
        //Set physics
        cabinAir.setIsActive(true);
        Solids.setIsActive(false);
        defrosterFlow.setIsActive(false);
        AirStream.setIsActive(false);
        
        
        steadySteps.setIsUsed(true);
        tempCriterion.setIsUsed(false);
        velCriterion.setIsUsed(false);
        //Run steady state method
        
        MonitorIterationStoppingCriterion monitorUMax = (MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("UMax Criterion");
        MonitorIterationStoppingCriterion monitorTVRMax = (MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("TVRMax Criterion");
        Collection<Region> regions=collect.collectRegions(simulation_0);
        Collection<Region> regionNoCore = collect.collectRegions(regions, "cabinAir", true);
        monitorUMax.setIsUsed(true);
        steadySteps.setMaximumNumberSteps((int) (simulation_0.getSimulationIterator().getCurrentIteration() + cabSteadySteps));
        while(simulation_0.getSimulationIterator().getCurrentIteration() < steadySteps.getMaximumNumberSteps()){
            if(monitorUMax.getIsSatisfied()) {
                for (int i = 1; i < 5; i++) {
                    simulation_0.getMeshManager().removeInvalidCells(new NeoObjectVector(regionNoCore.toArray()), NeoProperty.fromString("{\'function\': \'VelocityMagnitude\', \'functionValue\': " + ((MonitorIterationStoppingCriterionMaxLimitType) monitorUMax.getCriterionType()).getLimit().evaluate()+ ", \'functionOperator\': 1," + invalidCriteria));
                }
                }
        
        simulation_0.getSimulationIterator().run();
        }
        monitorUMax.setIsUsed(false);
        //Add temperature solvers: Disable contant density model and add ideal gas seg. temp model
        ConstantDensityModel constantDensityModel_0 = cabinAir.getModelManager().getModel(ConstantDensityModel.class);
        cabinAir.disableModel(constantDensityModel_0);
        cabinAir.enable(IdealGasModel.class);
        cabinAir.enable(SegregatedFluidTemperatureModel.class);
        SteadyModel steadyModel_0 = cabinAir.getModelManager().getModel(SteadyModel.class);
        cabinAir.disableModel(steadyModel_0);
        cabinAir.enable(ImplicitUnsteadyModel.class);

        //Set standard state temperature to T_ref
        Units units_1 = simulation_0.getUnitsManager().getInternalUnits(new IntVector(new int[] {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
        ConstantMaterialPropertyMethod constantMaterialPropertyMethod_1 = ((ConstantMaterialPropertyMethod) gas_0.getMaterialProperties().getMaterialProperty(StandardStateTemperatureProperty.class).getMethod());
        constantMaterialPropertyMethod_1.getQuantity().setDefinition("$T_ref");
        constantMaterialPropertyMethod_1.getQuantity().setUnits(units_1);
        
    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
////////////////////////////////////////////////////////////////////////////////////// STEP 7a - SAVE FILE ///////////////////////////////////////////////////////////////////////////////////////
			
     private void SaveSim1(String saveFile) {
        String dir = simulation_0.getSessionDir(); //get the name of the directory
        String sep = System.getProperty("file.separator"); //get the right separator for your operative system
        String filename = simulation_0.getPresentationName(); //get the name of the current sim file
        switch(saveFile){
            case "sim":
                simulation_0.saveState(dir + sep + filename + "_beforeGravity.sim"); //save the current sim file as| name_of_the_file@mod.sim
                break;
            case "java":
                simulation_0.saveState(resolvePath(filename + "_beforeGravity.sim")); //save the current sim file as| name_of_the_file@mod.sim
                break;
        }
        
    }
			
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	
/////////////////////////////////////////////////////////////////////////////////// STEP 8 - SETTING GRAVITY //////////////////////////////////////////////////////////////////////////////////////
		
    
	 private void SetGravity(boolean gravity, String outletCondition) {
		 
        simulation_0.println("============");
        simulation_0.println("Setting gravity");
        simulation_0.println("============");
        
        
        PhysicsContinuum cabinAir = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("cabinAir")); 
        Region CabinRegion = ((Region) simulation_0.getRegionManager().getRegion("cabinAir"));
        Boundary Outlet = CabinRegion.getBoundaryManager().getBoundary("outlets");
        
        if(gravity){
        
        simulation_0.println("Setting Reference Altitude and Density");
        
        cabinAir.enable(GravityModel.class);
        cabinAir.getReferenceValues().get(ReferenceAltitude.class).setDefinition("[0, 0, ${ground}]");
        Units units_0 = ((Units) simulation_0.getUnitsManager().getObject("m"));
        cabinAir.getReferenceValues().get(ReferenceAltitude.class).setUnits(units_0);
        
        cabinAir.getReferenceValues().get(ReferenceDensity.class).setDefinition("101325/(287.058*$T_ref)");
        
        function.createScalarFieldFunction(simulation_0, "hydrostaticPressure", "((-101325/(287.058*$T_ref)*9.81*($${Centroid}[2]-${ground})))");
        function.createScalarFieldFunction(simulation_0,"outletPressure","(101325*(exp(-9.81*($${Centroid}[2]-${ground})/287.058/$T_ref)-1)+9.81*101325/(287.058*$T_ref)*($${Centroid}[2]-${ground}))");
        
        UserFieldFunction hydroStaticPressure = ((UserFieldFunction) simulation_0.getFieldFunctionManager().getFunction("hydrostaticPressure"));
        UserFieldFunction outletPressure = ((UserFieldFunction) simulation_0.getFieldFunctionManager().getFunction("outletPressure")); 
        
        
        simulation_0.println("Setting hydrostatic pressure as initial condition");
        
        InitialPressureProfile initialPressure = cabinAir.getInitialConditions().get(InitialPressureProfile.class);
        initialPressure.setMethod(FunctionScalarProfileMethod.class);
        initialPressure.getMethod(FunctionScalarProfileMethod.class).setFieldFunction(outletPressure);
        
        if ("pressure".equals(outletCondition)){
        simulation_0.println("Setting temperature and altitude adjusted pressure outlet");
        simulation_0.println("============");
        StaticPressureProfile staticPressureOutlet = Outlet.getValues().get(StaticPressureProfile.class);
        staticPressureOutlet.setMethod(FunctionScalarProfileMethod.class);
        staticPressureOutlet.getMethod(FunctionScalarProfileMethod.class).setFieldFunction(outletPressure);
        }
        
        
        } else{
        simulation_0.println("===Gravity Disabled===");
        //GravityModel gravityModel_0 = cabinAir.getModelManager().getModel(GravityModel.class);
        //cabinAir.disableModel(gravityModel_0);
        }
    }
    

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
////////////////////////////////////////////////////////////////////////////////////// STEP 9 - SIMULATION HISTORY ////////////////////////////////////////////////////////////////////////////////
	
   
	private void simHistorySetup(boolean simHistory, String saveFile) {
		
        if(simHistory){
            simulation_0.println("============");
            simulation_0.println("Setting up Solution History");
            simulation_0.println("============");
            
            SolutionHistory solutionHistory_0 = null;
            switch(saveFile){
                case "sim":
                    String dir = simulation_0.getSessionDir(); //get the name of the directory
                    String sep = System.getProperty("file.separator"); //get the right separator for your operative system
                    solutionHistory_0 = simulation_0.get(SolutionHistoryManager.class).createForFile(dir+sep+"simulationHistory.simh", false,true);
                    break;
                case "java":
                    solutionHistory_0 = simulation_0.get(SolutionHistoryManager.class).createForFile(resolvePath("simulationHistory.simh"), false,true);
                    break;
            }
            

            solutionHistory_0.setAddReportFieldFunctions(true);
            
            solutionHistory_0.getInputs().setQuery(new Query(new CompoundPredicate(CompoundOperator.And, Arrays.<QueryPredicate>asList(new TypePredicate(TypeOperator.Is, Part.class), new CompoundPredicate(CompoundOperator.Or, Arrays.<QueryPredicate>asList(new NamePredicate(NameOperator.Contains, "X"), new NamePredicate(NameOperator.Contains, "Y"), new NamePredicate(NameOperator.Contains, "Z"))))), Query.STANDARD_MODIFIERS));

            HtcUserYPlusFunction htcUserYPlusFunction_0 = ((HtcUserYPlusFunction) simulation_0.getFieldFunctionManager().getFunction("HeatTransferCoefficientUserYPlus"));
            PrimitiveFieldFunction primitiveFieldFunction_0 = ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("Temperature"));
            PrimitiveFieldFunction primitiveFieldFunction_2 = ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("KeTurbBuoyancyC3"));
            PrimitiveFieldFunction primitiveFieldFunction_1 = ((PrimitiveFieldFunction) simulation_0.getFieldFunctionManager().getFunction("KeTurbBuoyancySource"));
            
            UserFieldFunction userFieldFunction_0 = ((UserFieldFunction) simulation_0.getFieldFunctionManager().getFunction("UserFieldFunction_1"));
            VectorMagnitudeFieldFunction vectorMagnitudeFieldFunction_0 = ((VectorMagnitudeFieldFunction) userFieldFunction_0.getMagnitudeFunction());
            VectorComponentFieldFunction vectorComponentFieldFunction_0 = ((VectorComponentFieldFunction) userFieldFunction_0.getComponentFunction(0));
            VectorComponentFieldFunction vectorComponentFieldFunction_1 = ((VectorComponentFieldFunction) userFieldFunction_0.getComponentFunction(1));
            VectorComponentFieldFunction vectorComponentFieldFunction_2 = ((VectorComponentFieldFunction) userFieldFunction_0.getComponentFunction(2));

            

            solutionHistory_0.setFunctions(new NeoObjectVector(new Object[] {primitiveFieldFunction_0,primitiveFieldFunction_1,primitiveFieldFunction_2}));

            StarUpdate starUpdate_6 = solutionHistory_0.getUpdate();
            starUpdate_6.getUpdateModeOption().setSelected(StarUpdateModeOption.Type.DELTATIME);

            DeltaTimeUpdateFrequency deltaTimeUpdateFrequency_0 = starUpdate_6.getDeltaTimeUpdateFrequency();
            Units units_1 = ((Units) simulation_0.getUnitsManager().getObject("s"));
            deltaTimeUpdateFrequency_0.setDeltaTime("30.0", units_1);

            solutionHistory_0.setExportAtAllStencils(true);
            solutionHistory_0.setAutoRescan(false);
        
        } else {
            simulation_0.println("============");
            simulation_0.println("Solution History Disabled");
            simulation_0.println("============");
        }
    }
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	////////////////////////////////////////////////////////////////////////////////////// STEP 10a - SAVE FILE ///////////////////////////////////////////////////////////////////////////////////////
			
     private void SaveSim2(String saveFile) {
        String dir = simulation_0.getSessionDir(); //get the name of the directory
        String sep = System.getProperty("file.separator"); //get the right separator for your operative system
        String filename = simulation_0.getPresentationName(); //get the name of the current sim file
        switch(saveFile){
            case "sim":
                simulation_0.saveState(dir + sep + filename + "_cabinsteady.sim"); //save the current sim file as| name_of_the_file@mod.sim
                break;
            case "java":
                simulation_0.saveState(resolvePath(filename + "_cabinsteady.sim")); //save the current sim file as| name_of_the_file@mod.sim
                break;
        }
        
    }
			
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
/////////////////////////////////////////////////////////////////////////////// STEP 10 - RUN UNSTEADY CABIN ///////////////////////////////////////////////////////////////////////////////////////
    
   
		 private void RunUnsteady(double tsUnsteady, double unfreezeTimeStep, String solverType, boolean srh) {
			 
        simulation_0.println("============");
        simulation_0.println("Running unsteady cabin");
        simulation_0.println("============");

        //Load physics
        PhysicsContinuum cabinAir = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("cabinAir")); 
        PhysicsContinuum Solids = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Solids"));
        ImplicitUnsteadySolver tS = ((ImplicitUnsteadySolver) simulation_0.getSolverManager().getSolver(ImplicitUnsteadySolver.class));
        PhysicalTimeStoppingCriterion stopTime = ((PhysicalTimeStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Physical Time"));
        Units units_0 = ((Units) simulation_0.getUnitsManager().getObject("s"));
        MonitorIterationStoppingCriterion monitorUMax = (MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("UMax Criterion");
        MonitorIterationStoppingCriterion monitorTVRMax = (MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("TVRMax Criterion");
        Collection<Region> regions = collect.collectRegions(simulation_0);
        Collection<Region> regionNoCore = collect.collectRegions(regions, "^Core.*", false);
        
        monitorUMax.setIsUsed(false);
        monitorTVRMax.setIsUsed(false);
        
        //Set physics
        cabinAir.setIsActive(true);
        Solids.setIsActive(true);
        
        MonitorIterationStoppingCriterion tempCriterion = ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridTemperatureSum Monitor Criterion"));
        MonitorIterationStoppingCriterion velCriterion =  ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridVelocitySum Monitor Criterion"));
        tempCriterion.setIsUsed(true);
        velCriterion.setIsUsed(true);
        
        ReportMonitor iterationsPerStep = ((ReportMonitor) simulation_0.getMonitorManager().getMonitor("IterationsPerTimeStep Monitor"));
        ReportMonitor iterationsAtStep = ((ReportMonitor) simulation_0.getMonitorManager().getMonitor("IterationAtTimeStep Monitor"));
        
        StarUpdate starUpdate_0 = iterationsPerStep.getStarUpdate();
        StarUpdate starUpdate_1 = iterationsAtStep.getStarUpdate();

        starUpdate_0.getUpdateModeOption().setSelected(StarUpdateModeOption.Type.TIMESTEP);
        starUpdate_1.getUpdateModeOption().setSelected(StarUpdateModeOption.Type.TIMESTEP);
        
        switch (solverType) {
          case "First":
              tS.getTimeDiscretizationOption().setSelected(TimeDiscretizationOption.Type.FIRST_ORDER);
              simulation_0.println("Running First Order..");
              break;
          case "Second":
              tS.getTimeDiscretizationOption().setSelected(TimeDiscretizationOption.Type.SECOND_ORDER);
              simulation_0.println("Running Second Order..");
              break;
          default:
            simulation_0.println("Running First Order..");
            break;
      }
        if (srh){
            cabinAir.enable(SrhModel.class);
        }
        
      stopTime.getMaximumTime().setValue(0.5);
      stopTime.getMaximumTime().setUnits(units_0);
      tS.getTimeStep().setValue(unfreezeTimeStep);
      tS.getTimeStep().setUnits(units_0);
	
      
	simulation_0.getSimulationIterator().run();
      
      //monitorUMax.setIsUsed(true);
      //monitorUMax.setInnerIterationCriterion(false);
      //monitorTVRMax.setIsUsed(true);
      
      stopTime.getMaximumTime().setValue(tsUnsteady);
      stopTime.getMaximumTime().setUnits(units_0);

             
      simulation_0.getSimulationIterator().run();
        
    }
	


 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  
	//////////////////////////////////////////////////////////////////////////////////////////////// STEP 12 -FREEZING FLOW ////////////////////////////////////////////////////////////////////////////////////////////
    
     
	  private void RunFreezingFlow(double totalTime, double freezeTime, double freezeTimeStep, double unfreezeTime, double unfreezeTimeStep, String solverType) {
      
		simulation_0.println("============");
        simulation_0.println("Running cabin unsteady with freezing flow");
        simulation_0.println("============");
        
        PhysicsContinuum cabinAir = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("cabinAir")); 
        PhysicsContinuum Solids = ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Solids"));
        ImplicitUnsteadySolver tS = ((ImplicitUnsteadySolver) simulation_0.getSolverManager().getSolver(ImplicitUnsteadySolver.class));
        PhysicalTimeStoppingCriterion stopTime = ((PhysicalTimeStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Physical Time"));
        
        SegregatedFlowSolver segregatedFlowSolver_0 = ((SegregatedFlowSolver) simulation_0.getSolverManager().getSolver(SegregatedFlowSolver.class));
        KeTurbSolver keTurbSolver_0 = ((KeTurbSolver) simulation_0.getSolverManager().getSolver(KeTurbSolver.class));
        KeTurbViscositySolver keTurbViscositySolver_0 = ((KeTurbViscositySolver) simulation_0.getSolverManager().getSolver(KeTurbViscositySolver.class));
        Units units_0 = ((Units) simulation_0.getUnitsManager().getObject("s"));
        
        MonitorIterationStoppingCriterion monitorUMax = (MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("UMax Criterion");
        Collection<Region> regions = collect.collectRegions(simulation_0);
        Collection<Region> regionNoCore = collect.collectRegions(regions, "^Core.*", false);
        
        //Set physics
        cabinAir.setIsActive(true);
        Solids.setIsActive(true);
        
        MonitorIterationStoppingCriterion tempCriterion = ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridTemperatureSum Monitor Criterion"));
        MonitorIterationStoppingCriterion velCriterion =  ((MonitorIterationStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("PresentationGridVelocitySum Monitor Criterion"));
        tempCriterion.setIsUsed(true);
        velCriterion.setIsUsed(true);
        
        switch (solverType) {
          case "First":
              tS.getTimeDiscretizationOption().setSelected(TimeDiscretizationOption.Type.FIRST_ORDER);
              simulation_0.println("Running First Order..");
              break;
          case "Second":
              tS.getTimeDiscretizationOption().setSelected(TimeDiscretizationOption.Type.SECOND_ORDER);
              simulation_0.println("Running Second Order..");
              break;
          default:
            simulation_0.println("Running First Order..");
            break;
        }
        
        FixedPhysicalTimeStoppingCriterion fixedPhysicalTimeStoppingCriterion_0;
        

        try {
        fixedPhysicalTimeStoppingCriterion_0 = ((FixedPhysicalTimeStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Fixed Physical Time"));    
        } 
        catch (Exception e) {
        fixedPhysicalTimeStoppingCriterion_0 = simulation_0.getSolverStoppingCriterionManager().createSolverStoppingCriterion(FixedPhysicalTimeStoppingCriterion.class);
        }
        
        stopTime.getMaximumTime().setValue(totalTime);

        double startTime;
        while(simulation_0.getSolution().getPhysicalTime()<totalTime) {

            //Run Frozen solver for ...s @ ... dT
            startTime = simulation_0.getSolution().getPhysicalTime();
            fixedPhysicalTimeStoppingCriterion_0.getFixedPhysicalTime().setValue(freezeTime);
            fixedPhysicalTimeStoppingCriterion_0.getFixedPhysicalTime().setUnits(units_0);
            tS.getTimeStep().setValue(freezeTimeStep);
            tS.getTimeStep().setUnits(units_0);
            segregatedFlowSolver_0.setFreezeFlow(true);    
            keTurbSolver_0.setFrozen(true);
            keTurbViscositySolver_0.setFrozen(true);
			
            
			simulation_0.getSimulationIterator().run();
                
            
            //Unfreeze flow and run for ...s @ 0.02s dT
            startTime = simulation_0.getSolution().getPhysicalTime();
            tS.getTimeStep().setValue(unfreezeTimeStep);
            tS.getTimeStep().setUnits(units_0);
            segregatedFlowSolver_0.setFreezeFlow(false);
            keTurbSolver_0.setFrozen(false);
            keTurbViscositySolver_0.setFrozen(false);
            fixedPhysicalTimeStoppingCriterion_0.getFixedPhysicalTime().setValue(unfreezeTime);
            fixedPhysicalTimeStoppingCriterion_0.getFixedPhysicalTime().setUnits(units_0);
            //while(simulation_0.getSolution().getPhysicalTime()<startTime+unfreezeTime){               
            simulation_0.getSimulationIterator().run();                
            //}   
            
            if (totalTime<simulation_0.getSolution().getPhysicalTime()) {
                break;
            }
            
            if ((simulation_0.getSolution().getPhysicalTime()%300.0==0.0) && (simulation_0.getSolution().getPhysicalTime()<1200.0)) {
                fixedPhysicalTimeStoppingCriterion_0.setIsUsed(false);
                RunAirStream(500.0);
                cabinAir.setIsActive(true);
                Solids.setIsActive(true);
                tempCriterion.setIsUsed(true);
                velCriterion.setIsUsed(true);
                fixedPhysicalTimeStoppingCriterion_0.setIsUsed(true);
                
                
            }
            
    }
    }
	  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
////////////////////////////////////////////////////////////////////////////////////// STEP 13 - SAVE FILE ///////////////////////////////////////////////////////////////////////////////////////
			

     private void SaveSim3(String saveFile) {
        String dir = simulation_0.getSessionDir(); //get the name of the directory
        String sep = System.getProperty("file.separator"); //get the right separator for your operative system
        String filename = simulation_0.getPresentationName(); //get the name of the current sim file
        switch(saveFile){
            case "sim":
                simulation_0.saveState(dir + sep + filename + "_Run.sim"); //save the current sim file as| name_of_the_file@mod.sim
                break;
            case "java":
                simulation_0.saveState(resolvePath(filename + "_Run.sim")); //save the current sim file as| name_of_the_file@mod.sim
                break;
        }
        
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
        
    private void insulateCabin(boolean insulation, double[] insulations){
        
        if(insulation){
            simulation_0.println("======================");
            simulation_0.println("Cabin Mode: Insulated");
            simulation_0.println("======================");
            
            BoundaryInterface boundaryInterface_0 = ((BoundaryInterface) simulation_0.getInterfaceManager().getInterface("solids/cabinAir"));

            boundaryInterface_0.setAllowPerPartValues(true);

            ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setQuery(null);

            CompositePart compositePart_1 = ((CompositePart) simulation_0.get(SimulationPartManager.class).getPart("Solids"));

            CompositePart compositePart_0 = ((CompositePart) compositePart_1.getChildParts().getPart("Door"));

            MeshPart meshPart_0 = ((MeshPart) compositePart_0.getChildParts().getPart("Armrest-inner"));

            MeshOperationPart meshOperationPart_0 = ((MeshOperationPart) simulation_0.get(SimulationPartManager.class).getPart("cabinAir"));

            PartContact partContact_0 = simulation_0.get(PartContactManager.class).getPartContact(meshPart_0,meshOperationPart_0);

            PartSurface partSurface_0 = ((PartSurface) meshPart_0.getPartSurfaceManager().getPartSurface("front-left-armrest"));

            PartSurface partSurface_1 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Armrest-inner.front-left-armrest"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_0 = ((InPlacePartSurfaceContact) partContact_0.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_0,partSurface_1));

            PartSurface partSurface_2 = ((PartSurface) meshPart_0.getPartSurfaceManager().getPartSurface("front-right-armrest"));

            PartSurface partSurface_3 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Armrest-inner.front-right-armrest"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_1 = ((InPlacePartSurfaceContact) partContact_0.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_2,partSurface_3));

            PartSurface partSurface_4 = ((PartSurface) meshPart_0.getPartSurfaceManager().getPartSurface("rear-left-armrest"));

            PartSurface partSurface_5 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Armrest-inner.rear-left-armrest"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_2 = ((InPlacePartSurfaceContact) partContact_0.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_4,partSurface_5));

            PartSurface partSurface_6 = ((PartSurface) meshPart_0.getPartSurfaceManager().getPartSurface("rear-right-armrest"));

            PartSurface partSurface_7 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Armrest-inner.rear-right-armrest"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_3 = ((InPlacePartSurfaceContact) partContact_0.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_6,partSurface_7));

            MeshPart meshPart_1 = ((MeshPart) compositePart_0.getChildParts().getPart("Door-inner-front-left"));

            PartContact partContact_1 = simulation_0.get(PartContactManager.class).getPartContact(meshPart_1,meshOperationPart_0);

            PartSurface partSurface_8 = ((PartSurface) meshPart_1.getPartSurfaceManager().getPartSurface("front-left-door-inner-1"));

            PartSurface partSurface_9 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Door-inner-front-left.front-left-door-inner-1"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_4 = ((InPlacePartSurfaceContact) partContact_1.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_8,partSurface_9));

            PartSurface partSurface_10 = ((PartSurface) meshPart_1.getPartSurfaceManager().getPartSurface("left-front-door-inner-2"));

            PartSurface partSurface_11 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Door-inner-front-left.left-front-door-inner-2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_5 = ((InPlacePartSurfaceContact) partContact_1.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_10,partSurface_11));

            MeshPart meshPart_2 = ((MeshPart) compositePart_0.getChildParts().getPart("Door-inner-front-right"));

            PartContact partContact_2 = simulation_0.get(PartContactManager.class).getPartContact(meshPart_2,meshOperationPart_0);

            PartSurface partSurface_12 = ((PartSurface) meshPart_2.getPartSurfaceManager().getPartSurface("front-right-door-inner-1"));

            PartSurface partSurface_13 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Door-inner-front-right.front-right-door-inner-1"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_6 = ((InPlacePartSurfaceContact) partContact_2.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_12,partSurface_13));

            PartSurface partSurface_14 = ((PartSurface) meshPart_2.getPartSurfaceManager().getPartSurface("right-front-door-inner-2"));

            PartSurface partSurface_15 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Door-inner-front-right.right-front-door-inner-2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_7 = ((InPlacePartSurfaceContact) partContact_2.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_14,partSurface_15));

            MeshPart meshPart_3 = ((MeshPart) compositePart_0.getChildParts().getPart("Door-inner-rear-left"));

            PartContact partContact_3 = simulation_0.get(PartContactManager.class).getPartContact(meshPart_3,meshOperationPart_0);

            PartSurface partSurface_16 = ((PartSurface) meshPart_3.getPartSurfaceManager().getPartSurface("left-rear-door-inner-2"));

            PartSurface partSurface_17 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Door-inner-rear-left.left-rear-door-inner-2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_8 = ((InPlacePartSurfaceContact) partContact_3.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_16,partSurface_17));

            PartSurface partSurface_18 = ((PartSurface) meshPart_3.getPartSurfaceManager().getPartSurface("rear-left-door-inner-1"));

            PartSurface partSurface_19 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Door-inner-rear-left.rear-left-door-inner-1"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_9 = ((InPlacePartSurfaceContact) partContact_3.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_18,partSurface_19));

            MeshPart meshPart_4 = ((MeshPart) compositePart_0.getChildParts().getPart("Door-inner-rear-right"));

            PartContact partContact_4 = simulation_0.get(PartContactManager.class).getPartContact(meshPart_4,meshOperationPart_0);

            PartSurface partSurface_20 = ((PartSurface) meshPart_4.getPartSurfaceManager().getPartSurface("rear-right-door-inner-1"));

            PartSurface partSurface_21 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Door-inner-rear-right.rear-right-door-inner-1"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_10 = ((InPlacePartSurfaceContact) partContact_4.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_20,partSurface_21));

            PartSurface partSurface_22 = ((PartSurface) meshPart_4.getPartSurfaceManager().getPartSurface("right-rear-door-inner-2"));

            PartSurface partSurface_23 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.Door-inner-rear-right.right-rear-door-inner-2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_11 = ((InPlacePartSurfaceContact) partContact_4.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_22,partSurface_23));

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setObjects(inPlacePartSurfaceContact_0, inPlacePartSurfaceContact_1, inPlacePartSurfaceContact_2, inPlacePartSurfaceContact_3, inPlacePartSurfaceContact_4, inPlacePartSurfaceContact_5, inPlacePartSurfaceContact_6, inPlacePartSurfaceContact_7, inPlacePartSurfaceContact_8, inPlacePartSurfaceContact_9, inPlacePartSurfaceContact_10, inPlacePartSurfaceContact_11);

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setName("DoorPanels");

            ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setQuery(null);

            CompositePart compositePart_3 = ((CompositePart) compositePart_1.getChildParts().getPart("cabinInterior"));

            CompositePart compositePart_2 = ((CompositePart) compositePart_3.getChildParts().getPart("InteriorWalls"));

            MeshPart meshPart_5 = ((MeshPart) compositePart_2.getChildParts().getPart("A-pillar-floor"));

            PartContact partContact_5 = simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_5);

            PartSurface partSurface_24 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.A-pillar-floor.A-pillar-floor-left 2"));

            PartSurface partSurface_25 = ((PartSurface) meshPart_5.getPartSurfaceManager().getPartSurface("A-pillar-floor-left 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_12 = ((InPlacePartSurfaceContact) partContact_5.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_24,partSurface_25));

            PartSurface partSurface_26 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.A-pillar-floor.A-pillar-floor-right 2"));

            PartSurface partSurface_27 = ((PartSurface) meshPart_5.getPartSurfaceManager().getPartSurface("A-pillar-floor-right 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_13 = ((InPlacePartSurfaceContact) partContact_5.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_26,partSurface_27));

            MeshPart meshPart_6 = ((MeshPart) compositePart_2.getChildParts().getPart("A-pillar-top"));

            PartContact partContact_6 = simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_6);

            PartSurface partSurface_28 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.A-pillar-top.A-pillar-left 2"));

            PartSurface partSurface_29 = ((PartSurface) meshPart_6.getPartSurfaceManager().getPartSurface("A-pillar-left 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_14 = ((InPlacePartSurfaceContact) partContact_6.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_28,partSurface_29));

            PartSurface partSurface_30 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.A-pillar-top.A-pillar-right 2"));

            PartSurface partSurface_31 = ((PartSurface) meshPart_6.getPartSurfaceManager().getPartSurface("A-pillar-right 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_15 = ((InPlacePartSurfaceContact) partContact_6.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_30,partSurface_31));

            MeshPart meshPart_7 = ((MeshPart) compositePart_2.getChildParts().getPart("B-pillar"));

            PartContact partContact_7 = simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_7);

            PartSurface partSurface_32 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.B-pillar.B-pillar-left 2"));

            PartSurface partSurface_33 = ((PartSurface) meshPart_7.getPartSurfaceManager().getPartSurface("B-pillar-left 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_16 = ((InPlacePartSurfaceContact) partContact_7.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_32,partSurface_33));

            PartSurface partSurface_34 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.B-pillar.B-pillar-right 2"));

            PartSurface partSurface_35 = ((PartSurface) meshPart_7.getPartSurfaceManager().getPartSurface("B-pillar-right 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_17 = 
              ((InPlacePartSurfaceContact) partContact_7.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_34,partSurface_35));

            MeshPart meshPart_8 = ((MeshPart) compositePart_2.getChildParts().getPart("Boot-sidewall"));

            PartContact partContact_8 = simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_8);

            PartSurface partSurface_36 = ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.Boot-sidewall.boot-sidewall-left"));

            PartSurface partSurface_37 = 
              ((PartSurface) meshPart_8.getPartSurfaceManager().getPartSurface("boot-sidewall-left"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_18 = 
              ((InPlacePartSurfaceContact) partContact_8.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_36,partSurface_37));

            PartSurface partSurface_38 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.Boot-sidewall.boot-sidewall-right"));

            PartSurface partSurface_39 = 
              ((PartSurface) meshPart_8.getPartSurfaceManager().getPartSurface("boot-sidewall-right"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_19 = 
              ((InPlacePartSurfaceContact) partContact_8.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_38,partSurface_39));

            MeshPart meshPart_9 = 
              ((MeshPart) compositePart_2.getChildParts().getPart("C-pillar"));

            PartContact partContact_9 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_9);

            PartSurface partSurface_40 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.C-pillar.C-pillar-left 2"));

            PartSurface partSurface_41 = 
              ((PartSurface) meshPart_9.getPartSurfaceManager().getPartSurface("C-pillar-left 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_20 = 
              ((InPlacePartSurfaceContact) partContact_9.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_40,partSurface_41));

            PartSurface partSurface_42 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.C-pillar.C-pillar-right 2"));

            PartSurface partSurface_43 = 
              ((PartSurface) meshPart_9.getPartSurfaceManager().getPartSurface("C-pillar-right 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_21 = 
              ((InPlacePartSurfaceContact) partContact_9.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_42,partSurface_43));

            MeshPart meshPart_10 = 
              ((MeshPart) compositePart_2.getChildParts().getPart("C-pillar-floor"));

            PartContact partContact_10 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_10);

            PartSurface partSurface_44 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.C-pillar-floor.C-pillar-floor-left 2"));

            PartSurface partSurface_45 = 
              ((PartSurface) meshPart_10.getPartSurfaceManager().getPartSurface("C-pillar-floor-left 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_22 = 
              ((InPlacePartSurfaceContact) partContact_10.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_44,partSurface_45));

            PartSurface partSurface_46 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.C-pillar-floor.C-pillar-floor-right 2"));

            PartSurface partSurface_47 = 
              ((PartSurface) meshPart_10.getPartSurfaceManager().getPartSurface("C-pillar-floor-right 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_23 = 
              ((InPlacePartSurfaceContact) partContact_10.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_46,partSurface_47));

            MeshPart meshPart_11 = 
              ((MeshPart) compositePart_0.getChildParts().getPart("windowframe-inner"));

            PartContact partContact_11 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshPart_11,meshOperationPart_0);

            PartSurface partSurface_48 = 
              ((PartSurface) meshPart_11.getPartSurfaceManager().getPartSurface("front-left-windowframe-inner"));

            PartSurface partSurface_49 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.windowframe-inner.front-left-windowframe-inner"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_24 = 
              ((InPlacePartSurfaceContact) partContact_11.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_48,partSurface_49));

            PartSurface partSurface_50 = 
              ((PartSurface) meshPart_11.getPartSurfaceManager().getPartSurface("front-right-windowframe-inner"));

            PartSurface partSurface_51 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.windowframe-inner.front-right-windowframe-inner"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_25 = 
              ((InPlacePartSurfaceContact) partContact_11.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_50,partSurface_51));

            PartSurface partSurface_52 = 
              ((PartSurface) meshPart_11.getPartSurfaceManager().getPartSurface("rear-left-windowframe-inner"));

            PartSurface partSurface_53 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.windowframe-inner.rear-left-windowframe-inner"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_26 = 
              ((InPlacePartSurfaceContact) partContact_11.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_52,partSurface_53));

            PartSurface partSurface_54 = 
              ((PartSurface) meshPart_11.getPartSurfaceManager().getPartSurface("rear-right-windowframe-inner"));

            PartSurface partSurface_55 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Door.windowframe-inner.rear-right-windowframe-inner"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_27 = 
              ((InPlacePartSurfaceContact) partContact_11.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_54,partSurface_55));

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setObjects(inPlacePartSurfaceContact_12, inPlacePartSurfaceContact_13, inPlacePartSurfaceContact_14, inPlacePartSurfaceContact_15, inPlacePartSurfaceContact_16, inPlacePartSurfaceContact_17, inPlacePartSurfaceContact_18, inPlacePartSurfaceContact_19, inPlacePartSurfaceContact_20, inPlacePartSurfaceContact_21, inPlacePartSurfaceContact_22, inPlacePartSurfaceContact_23, inPlacePartSurfaceContact_24, inPlacePartSurfaceContact_25, inPlacePartSurfaceContact_26, inPlacePartSurfaceContact_27);

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setName("Walls");

            ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setQuery(null);

            CompositePart compositePart_5 = 
              ((CompositePart) compositePart_1.getChildParts().getPart("Dashboard-CentreConsole-Steering"));

            CompositePart compositePart_4 = 
              ((CompositePart) compositePart_5.getChildParts().getPart("Dashboard"));

            MeshPart meshPart_12 = 
              ((MeshPart) compositePart_4.getChildParts().getPart("Dashboard-bottom"));

            PartContact partContact_12 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_12);

            PartSurface partSurface_56 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Dashboard-CentreConsole-Steering.Dashboard.Dashboard-bottom.Dashboard-bottom"));

            PartSurface partSurface_57 = 
              ((PartSurface) meshPart_12.getPartSurfaceManager().getPartSurface("Dashboard-bottom"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_28 = 
              ((InPlacePartSurfaceContact) partContact_12.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_56,partSurface_57));

            MeshPart meshPart_13 = 
              ((MeshPart) compositePart_4.getChildParts().getPart("Dashboard-leftcover"));

            PartContact partContact_13 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_13);

            PartSurface partSurface_58 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Dashboard-CentreConsole-Steering.Dashboard.Dashboard-leftcover.Dashboard-leftcover"));

            PartSurface partSurface_59 = 
              ((PartSurface) meshPart_13.getPartSurfaceManager().getPartSurface("Dashboard-leftcover"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_29 = 
              ((InPlacePartSurfaceContact) partContact_13.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_58,partSurface_59));

            MeshPart meshPart_14 = 
              ((MeshPart) compositePart_4.getChildParts().getPart("Dashboard-rightcover"));

            PartContact partContact_14 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_14);

            PartSurface partSurface_60 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Dashboard-CentreConsole-Steering.Dashboard.Dashboard-rightcover.Dashboard-rightcover"));

            PartSurface partSurface_61 = 
              ((PartSurface) meshPart_14.getPartSurfaceManager().getPartSurface("Dashboard-rightcover"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_30 = 
              ((InPlacePartSurfaceContact) partContact_14.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_60,partSurface_61));

            CompositePart compositePart_6 = 
              ((CompositePart) compositePart_4.getChildParts().getPart("Dashboard-top"));

            MeshPart meshPart_15 = 
              ((MeshPart) compositePart_6.getChildParts().getPart("DashboardTop"));

            PartContact partContact_15 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_15);

            PartSurface partSurface_62 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Dashboard-CentreConsole-Steering.Dashboard.Dashboard-top.DashboardTop.Dashboard-top"));

            PartSurface partSurface_63 = 
              ((PartSurface) meshPart_15.getPartSurfaceManager().getPartSurface("Dashboard-top"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_31 = 
              ((InPlacePartSurfaceContact) partContact_15.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_62,partSurface_63));

            PartSurface partSurface_64 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Dashboard-CentreConsole-Steering.Dashboard.Dashboard-top.DashboardTop.DriverInfo"));

            PartSurface partSurface_65 = 
              ((PartSurface) meshPart_15.getPartSurfaceManager().getPartSurface("DriverInfo"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_32 = 
              ((InPlacePartSurfaceContact) partContact_15.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_64,partSurface_65));

            MeshPart meshPart_16 = 
              ((MeshPart) compositePart_5.getChildParts().getPart("Centre-console"));

            PartContact partContact_16 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshPart_16,meshOperationPart_0);

            PartSurface partSurface_66 = 
              ((PartSurface) meshPart_16.getPartSurfaceManager().getPartSurface("Ccon-left-ext"));

            PartSurface partSurface_67 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Dashboard-CentreConsole-Steering.Centre-console.Ccon-left-ext"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_33 = 
              ((InPlacePartSurfaceContact) partContact_16.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_66,partSurface_67));

            PartSurface partSurface_68 = 
              ((PartSurface) meshPart_16.getPartSurfaceManager().getPartSurface("Ccon-middle-ext"));

            PartSurface partSurface_69 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Dashboard-CentreConsole-Steering.Centre-console.Ccon-middle-ext"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_34 = 
              ((InPlacePartSurfaceContact) partContact_16.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_68,partSurface_69));

            PartSurface partSurface_70 = 
              ((PartSurface) meshPart_16.getPartSurfaceManager().getPartSurface("Ccon-right-ext"));

            PartSurface partSurface_71 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Dashboard-CentreConsole-Steering.Centre-console.Ccon-right-ext"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_35 = 
              ((InPlacePartSurfaceContact) partContact_16.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_70,partSurface_71));

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setObjects(inPlacePartSurfaceContact_28, inPlacePartSurfaceContact_29, inPlacePartSurfaceContact_30, inPlacePartSurfaceContact_31, inPlacePartSurfaceContact_32, inPlacePartSurfaceContact_33, inPlacePartSurfaceContact_34, inPlacePartSurfaceContact_35);

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setName("DashBoard_Tunnel");

            ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setQuery(null);

            MeshPart meshPart_17 = 
              ((MeshPart) compositePart_2.getChildParts().getPart("CabinRoof"));

            PartContact partContact_17 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshOperationPart_0,meshPart_17);

            PartSurface partSurface_72 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.InteriorWalls.CabinRoof.roof-inner 2"));

            PartSurface partSurface_73 = 
              ((PartSurface) meshPart_17.getPartSurfaceManager().getPartSurface("roof-inner 2"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_36 = 
              ((InPlacePartSurfaceContact) partContact_17.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_72,partSurface_73));

            MeshPart meshPart_18 = 
              ((MeshPart) compositePart_3.getChildParts().getPart("Sunroof-glassPlanels"));

            PartContact partContact_18 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshPart_18,meshOperationPart_0);

            PartSurface partSurface_74 = 
              ((PartSurface) meshPart_18.getPartSurfaceManager().getPartSurface("GlassFront"));

            PartSurface partSurface_75 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.Sunroof-glassPlanels.GlassFront"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_37 = 
              ((InPlacePartSurfaceContact) partContact_18.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_74,partSurface_75));

            PartSurface partSurface_76 = 
              ((PartSurface) meshPart_18.getPartSurfaceManager().getPartSurface("GlassRear"));

            PartSurface partSurface_77 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.cabinInterior.Sunroof-glassPlanels.GlassRear"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_38 = 
              ((InPlacePartSurfaceContact) partContact_18.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_76,partSurface_77));

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setObjects(inPlacePartSurfaceContact_36, inPlacePartSurfaceContact_37, inPlacePartSurfaceContact_38);

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setName("Roof");

            ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).createNewGroup();

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setQuery(null);

            CompositePart compositePart_7 = 
              ((CompositePart) compositePart_1.getChildParts().getPart("Seats"));

            MeshPart meshPart_19 = 
              ((MeshPart) compositePart_7.getChildParts().getPart("DriverSeat"));

            PartContact partContact_19 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshPart_19,meshOperationPart_0);

            PartSurface partSurface_78 = 
              ((PartSurface) meshPart_19.getPartSurfaceManager().getPartSurface("Driver-seats"));

            PartSurface partSurface_79 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Seats.DriverSeat.Driver-seats"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_39 = 
              ((InPlacePartSurfaceContact) partContact_19.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_78,partSurface_79));

            MeshPart meshPart_20 = 
              ((MeshPart) compositePart_7.getChildParts().getPart("PassengerSeat"));

            PartContact partContact_20 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshPart_20,meshOperationPart_0);

            PartSurface partSurface_80 = 
              ((PartSurface) meshPart_20.getPartSurfaceManager().getPartSurface("Passenger-seat"));

            PartSurface partSurface_81 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Seats.PassengerSeat.Passenger-seat"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_40 = 
              ((InPlacePartSurfaceContact) partContact_20.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_80,partSurface_81));

            MeshPart meshPart_21 = 
              ((MeshPart) compositePart_7.getChildParts().getPart("RearSeats"));

            PartContact partContact_21 = 
              simulation_0.get(PartContactManager.class).getPartContact(meshPart_21,meshOperationPart_0);

            PartSurface partSurface_82 = 
              ((PartSurface) meshPart_21.getPartSurfaceManager().getPartSurface("Rear-seats"));

            PartSurface partSurface_83 = 
              ((PartSurface) meshOperationPart_0.getPartSurfaceManager().getPartSurface("Solids.Seats.RearSeats.Rear-seats"));

            InPlacePartSurfaceContact inPlacePartSurfaceContact_41 = 
              ((InPlacePartSurfaceContact) partContact_21.get(PartSurfaceContactManager.class).getPartSurfaceContact(partSurface_82,partSurface_83));

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setObjects(inPlacePartSurfaceContact_39, inPlacePartSurfaceContact_40, inPlacePartSurfaceContact_41);

            ((NamedPartGroup) ((PartGrouping) boundaryInterface_0.get(PartGroupingManager.class).getObject("Subgrouping 1")).getObject("Subgroup 2")).setName("Seats");

            ThermalContactResistanceProfile thermalContactResistanceProfile_0 = 
              boundaryInterface_0.getValues().get(ThermalContactResistanceProfile.class);

            thermalContactResistanceProfile_0.setMethod(ByPartProfileMethod.class);

            
            ProxyProfile proxyProfile_0 = 
              ((ProxyProfile) thermalContactResistanceProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("DashBoard_Tunnel"));

            proxyProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(insulations[0]);

            Units units_0 = 
              ((Units) simulation_0.getUnitsManager().getObject("m^2-K/W"));

            proxyProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_0);

            ProxyProfile proxyProfile_1 = 
              ((ProxyProfile) thermalContactResistanceProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("DoorPanels"));

            proxyProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(insulations[1]);

            proxyProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_0);

            ProxyProfile proxyProfile_2 = 
              ((ProxyProfile) thermalContactResistanceProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Roof"));

            proxyProfile_2.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(insulations[2]);

            proxyProfile_2.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_0);

            ProxyProfile proxyProfile_3 = 
              ((ProxyProfile) thermalContactResistanceProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Seats"));

            proxyProfile_3.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(insulations[3]);

            proxyProfile_3.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_0);

            ProxyProfile proxyProfile_4 = 
              ((ProxyProfile) thermalContactResistanceProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Subgroup 1"));

            proxyProfile_4.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(insulations[5]);

            proxyProfile_4.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_0);

            ProxyProfile proxyProfile_5 = 
              ((ProxyProfile) thermalContactResistanceProfile_0.getMethod(ByPartProfileMethod.class).getProfileManager().getObject("Walls"));

            proxyProfile_5.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(insulations[4]);

            proxyProfile_5.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_0);
    }
        else{
            simulation_0.println("======================");
            simulation_0.println("Cabin Baseline-Mode");
            simulation_0.println("======================");
        }
    }
}

   

   



   
   

    

    
