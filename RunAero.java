/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cabinStudy;

/**
 *
 * @author thobeika
 */
import java.util.*;
import star.base.neo.*;
import star.base.report.*;
import star.common.*;
import star.turbulence.*;
import star.kwturb.*;
import star.segregatedflow.*;
import star.metrics.*;
import star.motion.*;
import java.io.*;
import star.flow.*;
import star.meshing.LatestMeshProxyRepresentation;
import star.segregatedenergy.*;
import classes.*;
import java.io.BufferedReader;
import java.text.DecimalFormat;
import star.coupledflow.AutomaticCourantNumberControl;
import star.coupledflow.CoupledEnergyModel;
import star.coupledflow.CoupledFlowModel;
import star.coupledflow.CoupledImplicitSolver;
import star.energy.MaximumAllowableTemperature;
import star.energy.MinimumAllowableTemperature;
import star.energy.StaticTemperatureProfile;
import star.energy.WallThermalOption;
import star.meshing.MeshOperationPart;
import star.meshing.SimpleBlockPart;
import star.post.SolutionHistory;
import star.post.SolutionHistoryManager;
import star.walldistance.WallDistanceSolver;
import star.vis.*;

// The main class
public class RunAero extends StarMacro {

    // global variables and instances
    Simulation sim = null;
    Collection<Region> regions = null;
    Collection<Region> regionsNoCores = null;

    DecimalFormat dfTime = new DecimalFormat("#####0.0");

    Collect collect = new Collect();
    Analyze analyze = new Analyze();
    Solve solve = new Solve();
    Function function = new Function();
    General general = new General();
    PVT pvt = new PVT();

    String invalidCriteria = " \'functionEnabled\': true,"
            + " \'minimumDiscontiguousCellsEnabled\': true, \'minimumDiscontiguousCells\': 10,"
            + " \'minimumCellVolumeEnabled\': true, \'minimumCellVolume\': 0.0,"
            + " \'minimumVolumeChangeEnabled\': true, \'minimumVolumeChange\': 1.0E-4,"
            + " \'minimumCellQualityEnabled\': true, \'minimumCellQuality\': 1.0E-5,"
            + " \'minimumContiguousFaceAreaEnabled\': true,  \'minimumContiguousFaceArea\': 1.0E-9,"
            + " \'minimumFaceValidityEnabled\': true, \'minimumFaceValidity\': 0.9}";

    String invalidCriteriaAir = " \'functionEnabled\': true,"
            + " \'minimumDiscontiguousCellsEnabled\': true, \'minimumDiscontiguousCells\': 1000000,"
            + " \'minimumCellVolumeEnabled\': true, \'minimumCellVolume\': 0.0,"
            + " \'minimumVolumeChangeEnabled\': true, \'minimumVolumeChange\': 1.0E-3,"
            + " \'minimumCellQualityEnabled\': true, \'minimumCellQuality\': 1.0E-5,"
            + " \'minimumContiguousFaceAreaEnabled\': true,  \'minimumContiguousFaceArea\': 1.0E-8,"
            + " \'minimumFaceValidityEnabled\': true, \'minimumFaceValidity\': 0.98}";

    boolean PVT;
    Boolean LagEB = false;
    boolean detailedMethodEvaluation;
    final Collection<String> wheels = Arrays.asList("front-left,front-right,rear-left,rear-right".split(","));

    @Override
    public void execute() {
        long startTime = System.nanoTime();
        sim = getActiveSimulation();
        boolean backup = ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("backup")).getQuantity().getDefinition().equals("1.0");
        regions = collect.collectRegions(sim);
        regionsNoCores = collect.collectRegions(regions, "^Core_.*", false);
        PVT = ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("Domain")).getQuantity().getDefinition().equals("1.0");
        detailedMethodEvaluation = ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("DetailedMethodEvaluation")).getQuantity().getDefinition().equals("1.0");   // Saves data for evaluating properties if CFD setup (resolution etc.)
        String volname = sim.getPresentationName();
        // Parameters not typically changed
        boolean localDynURF = false; // local dynamic under relaxation for velocity // ref is false // default is false        
        boolean gravity = false;
        boolean bcdAll = false;
        int compIter = 00;
        String formulation = "IDDES"; // ref is IDDES
        String pAccelerator = "BiCGStab"; // ref is BiCGStab // default is CG 
        String vAccelerator = "NONE"; // ref is BiCGStab // default is NONE 
        String kwAccelerator = "NONE"; // ref is BiCGStab // default is NONE

        double unsteadyPressureAMGtol = 0.05;  // ref is 0.05 // default is 0.1        
        double unsteadyVelAMGtol = 0.05;  // ref is 0.05 // default is 0.1

        double unsteadyTurbAMGtol = 0.1; // ref is 0.05 // default is 0.1
        double magicValue = 1.0E-20; // ref is 1.0E-20 // default is 1.0
        double U_ref = 140*5/18;
        double T_ref = 273.15;
        double T_wall = 278.15;
        
        //Default settings and parameters
        // Parameters
        double timeStep = 0.02;   // ref is 2.5e-4
        int maxfirst = 500;
        int maxSteadyIterations = 1000;
        int autoSaveInterval = 2000000000;
        double averagingStartTime = 60.0; // ref is 2.7
        double maxPhysicalTime = 180.0;  // new ref is 3.7 
        int innerIterations = 6; // ref is 6
        double vURF = 0.8;   // ref is 0.9
        double pURF = 0.4;   // ref is 0.4
        String HOTDO = "NONE"; // ref is NONE // default is NONE BDF2_5 highest order
        boolean iterate = true; // Select to try to run iteration by iteration with control over velocity

        // Define the timestep ramp
        List<RampInstance> timeStepRamp = new ArrayList<>();
        timeStepRamp.add(new RampInstance(1e-1, 5e-3, 20, 0.7, 0.25, 0.1, 0.1)); //  400 iterations
        timeStepRamp.add(new RampInstance(2.0, 5e-3, 10, 0.7, 0.3, 0.1, 0.1)); //  3800 iterations
        timeStepRamp.add(new RampInstance(2.4, 1e-3, 8, 0.8, 0.4, 0.1, 0.1)); //  3200 iterations
        //timeStepRamp.add(new RampInstance(2.45, timeStep, 6, 0.9, 0.4, 0.1, 0.1)); //  1200 iterations
        timeStepRamp.add(new RampInstance(maxPhysicalTime, timeStep, innerIterations, vURF, pURF, unsteadyVelAMGtol, unsteadyPressureAMGtol)); // 33600 iterations

        if (PVT) { // modify the solver parameters to work with the PVT Domain setup
            // Solver parameters
            maxfirst = 1000;
            maxSteadyIterations = 3000;
            autoSaveInterval = 2000000000;
            averagingStartTime = 2.7; // ref is 2.7
            maxPhysicalTime = 3.7;  // new ref is 3.7 
            innerIterations = 10; // ref is 6
            vURF = 0.8;   // ref is 0.9
            pURF = 0.4;   // ref is 0.4

            timeStepRamp.clear();
            timeStepRamp.add(new RampInstance(1e-1, 5e-3, 20, 0.7, 0.25, 0.1, 0.1)); //  400 iterations
            timeStepRamp.add(new RampInstance(2.0, 5e-3, 12, 0.7, 0.3, 0.1, 0.1)); //  3800 iterations
            timeStepRamp.add(new RampInstance(2.4, 1e-3, 10, 0.8, 0.4, 0.1, 0.1)); //  3200 iterations
            timeStepRamp.add(new RampInstance(maxPhysicalTime, timeStep, innerIterations, vURF, pURF, unsteadyVelAMGtol, unsteadyPressureAMGtol)); // 33600 iterations
        }

        sim.getInterfaceManager().getDirectIntersectorOption().setSelected(DirectIntersectorOption.Type.GEOMETRY_BASED);

//======================================================================
        //
        // DO NOT TOUCH FROM HERE! // Rage Against the Machine- Killing in the Name
        //
        //======================================================================
        sim.getUnitsManager().getObject("m").setPreferred(true);
        // Turn on auto save
        AutoSave autoSave = sim.getSimulationIterator().getAutoSave();
        StarUpdate starUpdate = autoSave.getStarUpdate();
        starUpdate.setEnabled(true);
        starUpdate.getIterationUpdateFrequency().setIterations(autoSaveInterval);

        solve.defineRemFF(sim);
        removeBadAndInvalidCells();
        // enable bad cell detection model
        solve.enableBCD(sim, bcdAll);
        //solve.setVelocityChangeMonitoring(sim);

        //==========================
        // Setting up reports, field function, etc... 
        //==========================      
        if (gravity) {
            solve.setGravity(sim);
        }
        setGeneralSolverSettings(); // to remove done in Setup
        solve.turnOffNormalization(sim); // to remove done in Setup
        setFieldFunctionRef(sim);

        if (detailedMethodEvaluation) {
            setDetailedME();
        }
        solve.createWaypointEvent(sim, "AveragingStart", "Physical Time", averagingStartTime);
        solve.createWaypointEvent(sim, "First500Iterations", "Iteration", 500);
        createWDUReports(PVT);
        createReportsMonitors();
        solve.setMonitorUpdate(sim, averagingStartTime, timeStep, false);

        //==========================
        // Running Steady 
        //==========================
        set_steady_state_settings(magicValue, pAccelerator, vAccelerator, kwAccelerator); // to remove done in Setup
        //set_gradient_limits();
        //setCabinConditions(U_ref,T_ref,T_wall);
        run_1st_order(maxfirst);
        
        

        run_2nd_order(maxSteadyIterations, compIter, gravity);
        //export_plots_monitors("Steady", PVT, startTime);

        // in the future write a boolean to save the steady state based on boolean in settings.toml
        //general.saveState(sim, volname.replace("_vol", "_steady") + ".sim", backup);
        //==========================
        // Running Unsteady Part
        //==========================
        //solve.setMonitorUpdate(sim, averagingStartTime, timeStep, false);
        //set_unsteady_solver_settings(formulation, pAccelerator, vAccelerator, kwAccelerator, unsteadyVelAMGtol, unsteadyPressureAMGtol, unsteadyTurbAMGtol, localDynURF, HOTDO);
        //solve.setMonitorUpdate(sim, averagingStartTime, timeStep, true); // Ran another time just to update the STET. should later be replaced with a neater solution to update in a method locally
        //run_unsteady(timeStepRamp, timeStep, iterate, HOTDO, averagingStartTime);
        //if (sim.getSolution().getPhysicalTime() < maxPhysicalTime) {
        //    general.saveState(sim, volname.replace("_vol", "_@it-" + sim.getSimulationIterator().getCurrentIteration()) + ".sim", backup);
        //} else {
        //general.saveState(sim, volname.replace("_vol", "_Finished") + ".sim", backup);
        /*    //==========================
            // Posting
            //==========================
            export_plots_monitors(".", PVT, startTime);
            export_XYZ_csv();
            exportWduImages(PVT);
            python();
            general.bashRun(sim, "grep -n Removed log/run.log");
            general.bashRun(sim, "grep -n BCD log/run.log");
            sim.println("=====> Total time to run = " + dfTime.format(general.deltaTime(startTime) / 3600) + " hours <=====");
        }*/
    }

    private void setTemperature() {

        PhysicsContinuum continuum = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));

        //if (((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("BrakediscSetup")).getQuantity().getDefinition().equals("1.0")) {
        continuum.disableModel(continuum.getModelManager().getModel(ConstantDensityModel.class));

        //solve.setGravity(sim);
        //continuum.getReferenceValues().get(ReferenceDensity.class).setDefinition("${AirDensity}");
        continuum.enable(IdealGasModel.class);
        continuum.getModelManager().getModel(IdealGasModel.class).setIncompressible(true);
        continuum.enable(SegregatedFluidTemperatureModel.class);

        ((SegregatedEnergySolver) sim.getSolverManager().getSolver(SegregatedEnergySolver.class)).getFluidUrfQuantity().setValue(0.7);

        continuum.getReferenceValues().get(MinimumAllowableTemperature.class).setDefinition("${T_ref}-3");
        continuum.getReferenceValues().get(MaximumAllowableTemperature.class).setDefinition("${BrakeTemp}+3");
        continuum.getInitialConditions().get(StaticTemperatureProfile.class).getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("${T_ref}");
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.inlet").getValues().get(StaticTemperatureProfile.class).getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("${T_ref}");
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.outlet").getValues().get(StaticTemperatureProfile.class).getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("${T_ref}");

        for (Boundary boundary : collect.collectBoundaries(sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundaries(), ".*Brakedisc-.*")) {
            sim.println("setting temperature for " + boundary.getPresentationName() + " in region AirStream ");
            boundary.getConditions().get(WallThermalOption.class).setSelected(WallThermalOption.Type.TEMPERATURE);
            boundary.getValues().get(StaticTemperatureProfile.class).getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("${BrakeTemp}");
        }

        for (Region region : collect.collectRegions(collect.collectRegions(sim), "^Brakedisc_.*", true)) {
            for (Boundary boundary : collect.collectBoundaries(region.getBoundaryManager().getBoundaries(), ".*Brakedisc-.*")) {
                sim.println("setting temperature for " + boundary.getPresentationName() + " in region " + region.getPresentationName());
                boundary.getConditions().get(WallThermalOption.class).setSelected(WallThermalOption.Type.TEMPERATURE);
                boundary.getValues().get(StaticTemperatureProfile.class).getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("${BrakeTemp}");
            }

        }

        Collection<Boundary> boundariesCarAndDomain = collect.collectCarBoundariesAndDomain(sim);

        function.createMeanMonitor(sim, "HeatTransferCoefficientUserYPlus", "Mean Specified Y+ Heat Transfer Coefficient", "Scalar", boundariesCarAndDomain);
        function.createMeanMonitor(sim, "HeatTransferCoefficient", "Mean Heat Transfer Coefficient", "Scalar", boundariesCarAndDomain);

        function.createSurfAvgReport(sim, "YPlusAvgHTC-front-left", collect.collectBoundaries(boundariesCarAndDomain, ".*Brakedisc-front-left$"), "HeatTransferCoefficientUserYPlus", true);
        function.createSurfAvgReport(sim, "YPlusAvgHTC-front-right", collect.collectBoundaries(boundariesCarAndDomain, ".*Brakedisc-front-right$"), "HeatTransferCoefficientUserYPlus", true);
        function.createSurfAvgReport(sim, "YPlusAvgHTC-rear-left", collect.collectBoundaries(boundariesCarAndDomain, ".*Brakedisc-rear-left$"), "HeatTransferCoefficientUserYPlus", true);
        function.createSurfAvgReport(sim, "YPlusAvgHTC-rear-right", collect.collectBoundaries(boundariesCarAndDomain, ".*Brakedisc-rear-right$"), "HeatTransferCoefficientUserYPlus", true);

        function.createSurfAvgReport(sim, "AvgHTC-front-left", collect.collectBoundaries(boundariesCarAndDomain, ".*Brakedisc-front-left$"), "HeatTransferCoefficient", true);
        function.createSurfAvgReport(sim, "AvgHTC-front-right", collect.collectBoundaries(boundariesCarAndDomain, ".*Brakedisc-front-right$"), "HeatTransferCoefficient", true);
        function.createSurfAvgReport(sim, "AvgHTC-rear-left", collect.collectBoundaries(boundariesCarAndDomain, ".*Brakedisc-rear-left$"), "HeatTransferCoefficient", true);
        function.createSurfAvgReport(sim, "AvgHTC-rear-right", collect.collectBoundaries(boundariesCarAndDomain, ".*Brakedisc-rear-right$"), "HeatTransferCoefficient", true);

        function.createStatisticsReport(sim, "MeanYPlusAvgHTC-front-left", "YPlusAvgHTC-front-left Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanYPlusAvgHTC-front-right", "YPlusAvgHTC-front-right Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanYPlusAvgHTC-rear-left", "YPlusAvgHTC-rear-left Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanYPlusAvgHTC-rear-right", "YPlusAvgHTC-rear-right Monitor", "Mean", "AveragingStart");

        function.createStatisticsReport(sim, "MeanAvgHTC-front-left", "AvgHTC-front-left Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanAvgHTC-front-right", "AvgHTC-front-right Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanAvgHTC-rear-left", "AvgHTC-rear-left Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanAvgHTC-rear-right", "AvgHTC-rear-right Monitor", "Mean", "AveragingStart");

        // Create plots
        function.createPlot(sim, "YPlusHTC Plot", "YPlusAvgHTC-front-left Monitor,YPlusAvgHTC-front-right Monitor,YPlusAvgHTC-rear-left Monitor,YPlusAvgHTC-rear-right Monitor,"
                + "MeanYPlusAvgHTC-front-left Monitor,MeanYPlusAvgHTC-front-right Monitor,MeanYPlusAvgHTC-rear-left Monitor,MeanYPlusAvgHTC-rear-right Monitor");
        function.createPlot(sim, "_YPlusHTC Plot", "YPlusAvgHTC-front-left Monitor,YPlusAvgHTC-front-right Monitor,YPlusAvgHTC-rear-left Monitor,YPlusAvgHTC-rear-right Monitor,"
                + "MeanYPlusAvgHTC-front-left Monitor,MeanYPlusAvgHTC-front-right Monitor,MeanYPlusAvgHTC-rear-left Monitor,MeanYPlusAvgHTC-rear-right Monitor", true);

        function.createPlot(sim, "HTC Plot", "AvgHTC-front-left Monitor,AvgHTC-front-right Monitor,AvgHTC-rear-left Monitor,AvgHTC-rear-right Monitor,"
                + "MeanAvgHTC-front-left Monitor,MeanAvgHTC-front-right Monitor,MeanAvgHTC-rear-left Monitor,MeanAvgHTC-rear-right Monitor");
        function.createPlot(sim, "_HTC Plot", "AvgHTC-front-left Monitor,AvgHTC-front-right Monitor,AvgHTC-rear-left Monitor,AvgHTC-rear-right Monitor,"
                + "MeanAvgHTC-front-left Monitor,MeanAvgHTC-front-right Monitor,MeanAvgHTC-rear-left Monitor,MeanAvgHTC-rear-right Monitor", true);

        function.createMaxReport(sim, "rhoMax", regions, "Density", true);
        function.createMinReport(sim, "rhoMin", regions, "Density", true);
        function.createMaxReport(sim, "TMax", regions, "Temperature", true);
        function.createMinReport(sim, "TMin", regions, "Temperature", true);

        sim.println("left");
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.left").setBoundaryType((WallBoundary) sim.get(ConditionTypeManager.class).get(WallBoundary.class));
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.left").getConditions().get(WallShearStressOption.class).setSelected(WallShearStressOption.Type.SLIP);
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.left").getConditions().get(WallThermalOption.class).setSelected(WallThermalOption.Type.ADIABATIC);
        sim.println("right");
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.right").setBoundaryType((WallBoundary) sim.get(ConditionTypeManager.class).get(WallBoundary.class));
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.right").getConditions().get(WallShearStressOption.class).setSelected(WallShearStressOption.Type.SLIP);
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.right").getConditions().get(WallThermalOption.class).setSelected(WallThermalOption.Type.ADIABATIC);
        sim.println("top");
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.top").setBoundaryType((WallBoundary) sim.get(ConditionTypeManager.class).get(WallBoundary.class));
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.top").getConditions().get(WallShearStressOption.class).setSelected(WallShearStressOption.Type.SLIP);
        sim.getRegionManager().getRegion("AirStream").getBoundaryManager().getBoundary("DOMAIN_OR.top").getConditions().get(WallThermalOption.class).setSelected(WallThermalOption.Type.ADIABATIC);
        //}
    }

    private void set_steady_state_settings(double magicValue, String pAccelerator, String vAccelerator, String kwAccelerator) {

        PhysicsContinuum continuum = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));
        SegregatedFlowSolver segregatedFlowSolver = ((SegregatedFlowSolver) sim.getSolverManager().getSolver(SegregatedFlowSolver.class));
        VelocitySolver velocitySolver = segregatedFlowSolver.getVelocitySolver();
        PressureSolver pressureSolver = segregatedFlowSolver.getPressureSolver();
        velocitySolver.getUrfQuantity().setValue(0.7);
        pressureSolver.getUrfQuantity().setValue(0.3);
        segregatedFlowSolver.setBadCellRobustnessTreatment(false);
        //segregatedFlowSolver.setContinuityInitialization(true);
        //segregatedFlowSolver_0.getMaximumUnlimitedVelocity().setValue(100.0);
        //segregatedFlowSolver_0.setMinPressureCorrectionScale(0.8);
        //segregatedFlowSolver_0.setAcceptableVolumeChange(1.0E-3);
        //segregatedFlowSolver_0.setAcceptableVelocityChangeFactor(0.15);

        //segregatedFlowSolver.setContinuityInitialization(true);
        if (!LagEB) {
            KwTurbSolver kwTurbSolver = ((KwTurbSolver) sim.getSolverManager().getSolver(KwTurbSolver.class));
            KwTurbViscositySolver kwTurbViscositySolver = ((KwTurbViscositySolver) sim.getSolverManager().getSolver(KwTurbViscositySolver.class));

            kwTurbSolver.getUrfQuantity().setValue(0.8);
            kwTurbViscositySolver.setViscosityUrf(1.0);

            //kwTurbSolver.setTurbulenceInitialization(true);
        }
        //pressureSolver.getAMGLinearSolver().getAccelerationOption().setSelected(AMGAccelerationOption.Type.BiCGStab);
        //((AMGVCycle) pressureSolver.getAMGLinearSolver().getCycleType()).setPreSweeps(2);
        //((AMGVCycle) pressureSolver.getAMGLinearSolver().getCycleType()).setPostSweeps(3);
        // Stability improving settings
        continuum.getModelManager().getModel(GradientsModel.class).setNormalizedCurvatureFactor(magicValue);
        continuum.getModelManager().getModel(GradientsModel.class).setMaxSafeSkewAngle(75.0);
        continuum.getModelManager().getModel(GradientsModel.class).setMinUnsafeSkewAngle(88.0);

        // test initializing the domain with high visvocity ratio
        continuum.getInitialConditions().get(TurbulentViscosityRatioProfile.class).getMethod(ConstantScalarProfileMethod.class).getQuantity().setValue(200.0);

        //segregatedFlowSolver.getVelocitySolver().getAMGLinearSolver().getAccelerationOption().setSelected(AMGAccelerationOption.Type.valueOf(pAccelerator));
        //segregatedFlowSolver.getVelocitySolver().getAMGLinearSolver().getAccelerationOption().setSelected(AMGAccelerationOption.Type.valueOf(vAccelerator));
        //kwTurbSolver.getAMGLinearSolver().getAccelerationOption().setSelected(AMGAccelerationOption.Type.valueOf(kwAccelerator));
    }

    private void set_gradient_limits() {

        sim.println("=======================================");
        sim.println("====>    Setting Gradients");
        sim.println("=======================================");

        PhysicsContinuum continuum = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));

        GradientsModel gradients = continuum.getModelManager().getModel(GradientsModel.class);
        //gradients.getLimiterMethodOption().setSelected(LimiterMethodOption.Type.MINMOD);
        gradients.setUseTVBGradientLimiting(true);
        gradients.setAcceptableFieldVariationFactor(0.1);
    }

    private void run_2nd_order(int maxSteadyIterations, int compIter, boolean gravity) {
        sim.println("=======================================");
        sim.println("====>    Running 2nd order for " + maxSteadyIterations);
        sim.println("=======================================");

        // create stopping criteria for iteration monitor
        MeshManager meshManager = sim.getMeshManager();

        StepStoppingCriterion stepStoppingCriterion = ((StepStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));

        IterationMonitor iterationMonitor = ((IterationMonitor) sim.getMonitorManager().getMonitor("Iteration"));
        MonitorIterationStoppingCriterion monitorIterationStoppingCriterionMin = sim.getSolverStoppingCriterionManager().createIterationStoppingCriterion(iterationMonitor);
        ((MonitorIterationStoppingCriterionOption) monitorIterationStoppingCriterionMin.getCriterionOption()).setSelected(MonitorIterationStoppingCriterionOption.Type.MAXIMUM);
        MonitorIterationStoppingCriterionMaxLimitType monitorIterationStoppingCriterionMaxLimitType_0 = ((MonitorIterationStoppingCriterionMaxLimitType) monitorIterationStoppingCriterionMin.getCriterionType());
        monitorIterationStoppingCriterionMin.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.Type.AND);
        monitorIterationStoppingCriterionMaxLimitType_0.getLimit().setValue(500);
        monitorIterationStoppingCriterionMin.setPresentationName("Min Criterion");
        monitorIterationStoppingCriterionMin.setIsUsed(true);
        monitorIterationStoppingCriterionMin.setInnerIterationCriterion(false);

        // create stopping criteria for turbulent viscosity  
        ReportMonitor reportMonitorTVR = ((ReportMonitor) sim.getMonitorManager().getMonitor("TVRMax Monitor"));
        MonitorIterationStoppingCriterion monitorIterationStoppingCriterionTVR = sim.getSolverStoppingCriterionManager().createIterationStoppingCriterion(reportMonitorTVR);
        monitorIterationStoppingCriterionTVR.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.Type.OR);
        ((MonitorIterationStoppingCriterionOption) monitorIterationStoppingCriterionTVR.getCriterionOption()).setSelected(MonitorIterationStoppingCriterionOption.Type.MAXIMUM);
        MonitorIterationStoppingCriterionMaxLimitType monitorIterationStoppingCriterionMaxLimitType_1 = ((MonitorIterationStoppingCriterionMaxLimitType) monitorIterationStoppingCriterionTVR.getCriterionType());
        monitorIterationStoppingCriterionMaxLimitType_1.getLimit().setValue(99999.0);
        monitorIterationStoppingCriterionTVR.setPresentationName("TVRMax Criterion");
        monitorIterationStoppingCriterionTVR.setIsUsed(true);
        monitorIterationStoppingCriterionTVR.setInnerIterationCriterion(true);

        MonitorIterationStoppingCriterion monitorIterationStoppingCriterionU = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("UMax Criterion"));

        // run the simulation with a smart loop to control deleting high turb viscosity cells
        stepStoppingCriterion.setMaximumNumberSteps(maxSteadyIterations);
        while (sim.getSimulationIterator().getCurrentIteration() < stepStoppingCriterion.getMaximumNumberSteps()) {

            if (monitorIterationStoppingCriterionTVR.getIsSatisfied()) {
                for (int i = 1; i < 5; i++) {
                    meshManager.removeInvalidCells(new NeoObjectVector(regionsNoCores.toArray()), NeoProperty.fromString("{\'function\': \'TurbulentViscosityRatio\', \'functionValue\': 99000.0, \'functionOperator\': 1," + invalidCriteria));
                }
            }
            if (monitorIterationStoppingCriterionU.getIsSatisfied()) {
                for (int i = 1; i < 5; i++) {
                    meshManager.removeInvalidCells(new NeoObjectVector(regionsNoCores.toArray()), NeoProperty.fromString("{\'function\': \'VelocityMagnitude\', \'functionValue\': " + ((MonitorIterationStoppingCriterionMaxLimitType) monitorIterationStoppingCriterionU.getCriterionType()).getLimit().evaluate() + ", \'functionOperator\': 1," + invalidCriteria));
                }
            }
            for (int i = 1; i < 5; i++) {
                //meshManager.removeInvalidCells(new NeoObjectVector(regionsNoCores.toArray()), NeoProperty.fromString("{\'function\': \'Cells2Remove\', \'functionValue\': 0.999999, \'functionOperator\': 1," + invalidCriteria));
            }

            sim.getSimulationIterator().run();
            monitorIterationStoppingCriterionMaxLimitType_0.getLimit().setValue(sim.getSimulationIterator().getCurrentIteration() + 20);

            if (compIter > 0 && sim.getSimulationIterator().getCurrentIteration() == maxSteadyIterations) {
                set_compressible();
                stepStoppingCriterion.setMaximumNumberSteps(maxSteadyIterations + compIter);
            }

            if (((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("BrakeTemp")).getQuantity().evaluate() > 0 && sim.getSimulationIterator().getCurrentIteration() == maxSteadyIterations) {
                setTemperature();
                stepStoppingCriterion.setMaximumNumberSteps(maxSteadyIterations + 1000);
            }
        }
        monitorIterationStoppingCriterionMin.setIsUsed(false);
    }

    private void removeBadAndInvalidCells() {

        Collection<Region> regionAir = collect.collectRegions(regions, "^AirStream.*", true);
        Collection<Region> regionWheels = collect.collectRegions(regions, "^Wheel.*", true);
        MeshManager meshManager = sim.getMeshManager();
        sim.getInterfaceManager().initialize();
        for (int i = 1; i < 8; i++) {
            meshManager.removeInvalidCells(new NeoObjectVector(regionAir.toArray()), NeoProperty.fromString("{\'function\': \'Cells2Remove\', \'functionValue\': 0.999999, \'functionOperator\': 1," + invalidCriteriaAir));
        }
        for (int i = 1; i < 8; i++) {
            meshManager.removeInvalidCells(new NeoObjectVector(regionWheels.toArray()), NeoProperty.fromString("{\'function\': \'Cells2Remove\', \'functionValue\': 0.999999, \'functionOperator\': 1," + invalidCriteria));
        }

        //solve.internal2wall(sim);
    }

    private void run_1st_order(int maxfirst) {

        sim.println("=======================================");
        sim.println("====>    Running 1st order for " + maxfirst + " iterations");
        sim.println("=======================================");

        MeshManager meshManager = sim.getMeshManager();
        PhysicsContinuum continuum = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));
        StepStoppingCriterion stepStoppingCriterion = ((StepStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));
        SegregatedFlowModel segregatedFlowModel = continuum.getModelManager().getModel(SegregatedFlowModel.class);

        if (!LagEB) {
            segregatedFlowModel.getUpwindOption().setSelected(FlowUpwindOption.Type.FIRST_ORDER);
            SstKwTurbModel sstKwTurbModel = continuum.getModelManager().getModel(SstKwTurbModel.class);
            sstKwTurbModel.getUpwindOption().setSelected(UpwindOption.Type.FIRST_ORDER);
        }
        stepStoppingCriterion.setMaximumNumberSteps(50);
        sim.getSimulationIterator().run();

        stepStoppingCriterion.setMaximumNumberSteps(maxfirst);

        // create stopping criteria for velocity  
        ReportMonitor reportMonitorU = ((ReportMonitor) sim.getMonitorManager().getMonitor("UMax Monitor"));
        MonitorIterationStoppingCriterion monitorIterationStoppingCriterionU = sim.getSolverStoppingCriterionManager().createIterationStoppingCriterion(reportMonitorU);
        monitorIterationStoppingCriterionU.getLogicalOption().setSelected(SolverStoppingCriterionLogicalOption.Type.OR);
        ((MonitorIterationStoppingCriterionOption) monitorIterationStoppingCriterionU.getCriterionOption()).setSelected(MonitorIterationStoppingCriterionOption.Type.MAXIMUM);
        MonitorIterationStoppingCriterionMaxLimitType monitorIterationStoppingCriterionMaxLimitTypeU = ((MonitorIterationStoppingCriterionMaxLimitType) monitorIterationStoppingCriterionU.getCriterionType());
        monitorIterationStoppingCriterionMaxLimitTypeU.getLimit().setDefinition("${U_ref}*4");
        monitorIterationStoppingCriterionU.setPresentationName("UMax Criterion");
        monitorIterationStoppingCriterionU.setIsUsed(true);
        monitorIterationStoppingCriterionU.setInnerIterationCriterion(true);

        //run 1st order
        while (sim.getSimulationIterator().getCurrentIteration() < stepStoppingCriterion.getMaximumNumberSteps()) {

            if (monitorIterationStoppingCriterionU.getIsSatisfied()) {
                for (int i = 1; i < 5; i++) {
                    meshManager.removeInvalidCells(new NeoObjectVector(regions.toArray()), NeoProperty.fromString("{\'function\': \'VelocityMagnitude\', \'functionValue\': " + ((MonitorIterationStoppingCriterionMaxLimitType) monitorIterationStoppingCriterionU.getCriterionType()).getLimit().evaluate() + ", \'functionOperator\': 1," + invalidCriteria));
                }
            }

            sim.getSimulationIterator().run();

        }

        if (!LagEB) {
            sim.println("Setting second order");
            segregatedFlowModel.getUpwindOption().setSelected(FlowUpwindOption.Type.SECOND_ORDER);
            SstKwTurbModel sstKwTurbModel = continuum.getModelManager().getModel(SstKwTurbModel.class);
            sstKwTurbModel.getUpwindOption().setSelected(UpwindOption.Type.SECOND_ORDER);
        }
    }

    private void export_plots_monitors(String dir, boolean PVT, long Time) {
        long startTime = System.nanoTime();
        sim.println("=======================================");
        sim.println("====>    Exporting plots and monitors");
        sim.println("=======================================");

        // exporting all plots
        Collection<StarPlot> Plots = sim.getPlotManager().getPlots();
        for (StarPlot plot : Plots) {
            String name = plot.getPresentationName().replaceAll(" ", "_");
            sim.println("====> Exporting plot " + plot.getPresentationName() + " as " + name);
            if (plot.getPresentationName().startsWith("_")) {
                function.autoSetMinMaxYaxis(sim, plot, "AveragingStart");
            } else if (!plot.getPresentationName().matches("Residuals")) {
                function.autoSetMinMaxYaxis(sim, plot, "First500Iterations");
            }
            plot.encode(dir + "/Plots/" + name + ".png", "png", 1600, 1200);
        }

        Collection<Monitor> monitors = sim.getMonitorManager().getMonitors();
        Collection<ReportMonitor> reportMonitors = new ArrayList();
        for (Monitor monitor : monitors) {
            if (monitor instanceof ReportMonitor) {
                sim.println("Adding monitor " + monitor.getPresentationName() + " to Monitors.csv");
                reportMonitors.add((ReportMonitor) monitor);
            }
        }
        sim.println("====> Exporting Monitors.csv");
        sim.getMonitorManager().export(dir + "/Monitors.csv", ",", new NeoObjectVector(reportMonitors.toArray()));

        DecimalFormat df = new DecimalFormat("###0.0000");
        Collection<Report> reports = sim.getReportManager().getObjects();
        Map<String, String> map = new HashMap();
        for (Report report : reports) {
            if (report instanceof StatisticsReport) {
                map.put(report.getPresentationName(), df.format(((StatisticsReport) report).getValue()));
                sim.println("Adding report " + report.getPresentationName() + " to the Means.csv output");
            }
        }
        map.put("tag", sim.get(CommentManager.class).getCommentFor((ClientServerObject) ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("tag"))).toString());
        map.put("repo", sim.get(CommentManager.class).getCommentFor((ClientServerObject) ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("repo"))).toString());
        map.put("Runtime (h)", dfTime.format(general.deltaTime(Time) / 3600));
        map.put("Frontal area", df.format(sim.getReportManager().getObject("VehicleFrontalArea").getReportMonitorValue()));

        sim.println("====> Exporting Means.csv");
        general.createFile(sim, dir + "/Means.csv", map.keySet().toString().replaceAll("^.|.$", "").replaceAll(" ", "") + "\n");
        general.appendFile(sim, dir + "/Means.csv", map.values().toString().replaceAll("^.|.$", "").replaceAll(" ", "") + "\n");

        if (PVT) {

            pvt.exportPVTcsv(sim, dir);

        }
        sim.println("=====> Total time to export plots and monitors = " + dfTime.format(general.deltaTime(startTime) / 60) + " minutes");
    }

    private void set_unsteady_solver_settings(String formulation, String pAccelerator, String vAccelerator, String kwAccelerator, double unsteadyVelAMGtol, double unsteadyPressureAMGtol, double unsteadyTurbAMGtol, boolean localDynURF, String HOTDO) {

        sim.println("=======================================");
        sim.println("====>    Setting up unsteady solver settings");
        sim.println("=======================================");

        PhysicsContinuum continuum = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));
        if (!LagEB) {
            // Change to unsteady settings
            // Disable steady models
            continuum.getModelManager().getModel(KwAllYplusWallTreatment.class).disable();
            continuum.getModelManager().getModel(SstKwTurbModel.class).disable();
            continuum.getModelManager().getModel(KOmegaTurbulence.class).disable();
            continuum.getModelManager().getModel(SteadyModel.class).disable();
            continuum.getModelManager().getModel(RansTurbulenceModel.class).disable();

            // Enable unsteady models
            continuum.enable(ImplicitUnsteadyModel.class);
            continuum.enable(DesTurbulenceModel.class);
            continuum.enable(SstKwTurbDesModel.class);
            continuum.enable(KwAllYplusWallTreatment.class);

            //continuum.getModelManager().getModel(SstKwTurbDesModel.class).getRealizableTimeParameter().setRealizableTimeCoefficient(1.0);

            // Solver settings    
            sim.getSolverManager().getSolver(ImplicitUnsteadySolver.class).getTimeDiscretizationOption().setSelected(TimeDiscretizationOption.Type.SECOND_ORDER);
            sim.getSolverManager().getSolver(SegregatedFlowSolver.class).getVelocitySolver().setUseDynamicLocalUrf(localDynURF);

            sim.getSolverManager().getSolver(SegregatedFlowSolver.class).getHighOrderTemporalDiscretizationOption().setSelected(HighOrderTemporalDiscretizationOption.Type.valueOf(HOTDO));

            SegregatedFlowSolver segregatedFlowSolver = ((SegregatedFlowSolver) sim.getSolverManager().getSolver(SegregatedFlowSolver.class));
            KwTurbSolver kwTurbSolver = ((KwTurbSolver) sim.getSolverManager().getSolver(KwTurbSolver.class));

            continuum.getModelManager().getModel(SstKwTurbDesModel.class).getDesFormulationOption().setSelected(DesFormulationOption.Type.valueOf(formulation));
            //continuum.getModelManager().getModel(SstKwTurbDesModel.class).getRealizableTimeParameter().setRealizableTimeCoefficient(1.0); // ref ans def is 0.6
            //aMGLinearSolver_0.getAccelerationOption().setSelected(AMGAccelerationOption.Type.NONE);
            segregatedFlowSolver.getPressureSolver().getAMGLinearSolver().getAccelerationOption().setSelected(AMGAccelerationOption.Type.valueOf(pAccelerator));
            segregatedFlowSolver.getVelocitySolver().getAMGLinearSolver().getAccelerationOption().setSelected(AMGAccelerationOption.Type.valueOf(vAccelerator));
            kwTurbSolver.getAMGLinearSolver().getAccelerationOption().setSelected(AMGAccelerationOption.Type.valueOf(kwAccelerator));

            segregatedFlowSolver.getVelocitySolver().getAMGLinearSolver().setConvergeTol(unsteadyVelAMGtol);
            segregatedFlowSolver.getPressureSolver().getAMGLinearSolver().setConvergeTol(unsteadyPressureAMGtol);
            kwTurbSolver.getAMGLinearSolver().setConvergeTol(unsteadyTurbAMGtol);
        } else {

            continuum.getModelManager().getModel(SteadyModel.class).disable();
            continuum.enable(ImplicitUnsteadyModel.class);
            continuum.enable(SrhModel.class);

            sim.getSolverManager().getSolver(ImplicitUnsteadySolver.class).getTimeDiscretizationOption().setSelected(TimeDiscretizationOption.Type.SECOND_ORDER);

            sim.getSolverManager().getSolver(SegregatedFlowSolver.class).getHighOrderTemporalDiscretizationOption().setSelected(HighOrderTemporalDiscretizationOption.Type.valueOf(HOTDO));
        }

        function.createSTET(sim);
        function.createPlot(sim, "STET Plot", "STET Monitor");
        function.createStatisticsReport(sim, "MeanSTET", "STET Monitor", "Mean", "AveragingStart");

    }

    private void setCoreToFirstOrder() {
        sim.println("=======================================");
        sim.println("====>    Setting Core regions to first order");
        sim.println("=======================================");

        PhysicsContinuum continuum = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));
        // split the Cores into a first order continuum for enhanced stability
        PhysicsContinuum continuum2 = sim.getContinuumManager().createContinuum(PhysicsContinuum.class);
        continuum2.setPresentationName("Physics 2");
        continuum2.copyProperties(continuum);

        ((SstKwTurbModel) continuum2.getModelManager().getModel(SstKwTurbModel.class)).getUpwindOption().setSelected(UpwindOption.Type.FIRST_ORDER);
        ((GradientsModel) continuum2.getModelManager().getModel(GradientsModel.class)).setCustomAccuracyLevelSelector(1.0);
        for (Region region : collect.collectRegions(regions, "^Core_.*", true)) {
            continuum2.add(region);
        }
    }

    private void run_unsteady(List<RampInstance> timeStepRamp, double timeStep, boolean iterate, String HOTDO, double averagingStartTime) {

        ImplicitUnsteadySolver solver = ((ImplicitUnsteadySolver) sim.getSolverManager().getSolver(ImplicitUnsteadySolver.class));
        MeshManager meshManager = sim.getMeshManager();

        // Stopping criteria
        ((StepStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps")).setIsUsed(false);
        InnerIterationStoppingCriterion innerStopping = ((InnerIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Inner Iterations"));
        PhysicalTimeStoppingCriterion physicalTimeStoppingCriterion = ((PhysicalTimeStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Physical Time"));

        MonitorIterationStoppingCriterion minCriterion = (MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Min Criterion");
        minCriterion.setIsUsed(true);
        minCriterion.setInnerIterationCriterion(false);

        // create dummy criterion to prevent mincriterion from stopping solution every timestep
        MonitorIterationStoppingCriterion minCriterionDummy = sim.getSolverStoppingCriterionManager().createSolverStoppingCriterion(MonitorIterationStoppingCriterion.class, "Min Criterion");
        minCriterionDummy.copyProperties(minCriterion);
        ((MonitorIterationStoppingCriterionMaxLimitType) minCriterionDummy.getCriterionType()).getLimit().setDefinition("1.0E9");

        //disable inner iteration stopping criterion for TVRMax
        MonitorIterationStoppingCriterion monitorIterationStoppingCriterionTVR = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("TVRMax Criterion"));
        monitorIterationStoppingCriterionTVR.setInnerIterationCriterion(false);

        //disable inner iteration stopping criterion for UMax
        MonitorIterationStoppingCriterion monitorIterationStoppingCriterionU = ((MonitorIterationStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("UMax Criterion"));
        ((MonitorIterationStoppingCriterionMaxLimitType) monitorIterationStoppingCriterionU.getCriterionType()).getLimit().setDefinition("${U_ref}*4");
        monitorIterationStoppingCriterionU.setInnerIterationCriterion(false);

        // Enable mean IDDES blendning function
        if (detailedMethodEvaluation) {
            sim.getSolverManager().getSolver(KwTurbSolver.class).setLeaveTemporaryStorage(true);

            Collection<Boundary> boundariesCarAndDomain = collect.collectCarBoundariesAndDomain(sim);
            Collection<NamedObject> regions_boundaries = new ArrayList();
            regions_boundaries.addAll(regions);
            regions_boundaries.addAll(boundariesCarAndDomain);
            function.createMeanMonitor(sim, "SaTurbDesFd", "Mean IDDES Blendning Function", "Scalar", regions_boundaries);

            function.createMeanMonitor(sim, "CourantNumber", "Mean Convective Courant Number", "Scalar", regions_boundaries);

            solve.setMonitorUpdate(sim, averagingStartTime, timeStep, false);
        }

        // Run through the defined timestep ramp        
        for (RampInstance rampInstance : timeStepRamp) {
            // Retrieve values from the ramp instance and run until the next instance
            innerStopping.setMaximumNumberInnerIterations(rampInstance.getInnerIterations());
            physicalTimeStoppingCriterion.getMaximumTime().setValue(rampInstance.getMaxTime());
            sim.getSolverManager().getSolver(SegregatedFlowSolver.class).getVelocitySolver().getUrfQuantity().setValue(rampInstance.getV_URF());
            sim.getSolverManager().getSolver(SegregatedFlowSolver.class).getPressureSolver().getUrfQuantity().setValue(rampInstance.getP_URF());
            solver.getTimeStep().setValue(rampInstance.getTimeStep());
            sim.getSolverManager().getSolver(SegregatedFlowSolver.class).getVelocitySolver().getAMGLinearSolver().setConvergeTol(rampInstance.getVAMGtol());
            sim.getSolverManager().getSolver(SegregatedFlowSolver.class).getPressureSolver().getAMGLinearSolver().setConvergeTol(rampInstance.getPAMGtol());

            if (rampInstance.getTimeStep() <= timeStep && ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("WheelSetup")).getQuantity().getDefinition().equals("2.0")) {

                Collection<Boundary> SM_boundaries = new ArrayList();
                for (Region region : regionsNoCores) {
                    if (region.getPresentationName().matches("^Wheel_.*")) {
                        region.getValues().get(MotionSpecification.class).setMotion(sim.get(MotionManager.class).getObject(region.getPresentationName()));
                        SM_boundaries.addAll(region.getBoundaryManager().getBoundaries());
                    }
                }
                for (Boundary boundary : SM_boundaries) {
                    if (boundary.getBoundaryType().toString().equals("Wall")) {
                        boundary.getConditions().get(WallSlidingOption.class).setSelected(WallSlidingOption.Type.FIXED);
                    }
                }
                sim.getSolverManager().getSolver(SegregatedFlowSolver.class).getHighOrderTemporalDiscretizationOption().setSelected(HighOrderTemporalDiscretizationOption.Type.valueOf(HOTDO));
            }

            while (sim.getSolution().getPhysicalTime() < rampInstance.getMaxTime()) {

                if (monitorIterationStoppingCriterionTVR.getIsSatisfied() && minCriterion.getIsSatisfied()) {
                    sim.println("=====>TVR satisfied");                    
                    
                    if (sim.getSolution().getPhysicalTime() > averagingStartTime) {
                        //calling BCD
                        sim.println("=====>Applying new function BCD");
                        //long startTime = System.nanoTime();
                        solve.applyBCD(sim, (((MonitorIterationStoppingCriterionMaxLimitType) monitorIterationStoppingCriterionU.getCriterionType()).getLimit().evaluate()));
                        //sim.println("=====> Total time to BCD = " + dfTime.format(general.deltaTime(startTime) / 60) + " minutes");
                    } else {
                        sim.println("=====>Removing cells");
                        for (int i = 1; i < 5; i++) {
                            meshManager.removeInvalidCells(new NeoObjectVector(regionsNoCores.toArray()), NeoProperty.fromString("{\'function\': \'TurbulentViscosityRatio\', \'functionValue\': 99000.0, \'functionOperator\': 1," + invalidCriteria));
                        }
                    }
                    ((MonitorIterationStoppingCriterionMaxLimitType) minCriterion.getCriterionType()).getLimit().setValue(sim.getSimulationIterator().getCurrentIteration() + 2);
                }

                if (monitorIterationStoppingCriterionU.getIsSatisfied() && minCriterion.getIsSatisfied()) {
                    sim.println("=====>Umax satisfied");

                    if (sim.getSolution().getPhysicalTime() > averagingStartTime) {
                        //calling BCD
                        sim.println("=====>Applying new function BCD");
                        //long startTime = System.nanoTime();
                        solve.applyBCD(sim, (((MonitorIterationStoppingCriterionMaxLimitType) monitorIterationStoppingCriterionU.getCriterionType()).getLimit().evaluate()));
                        //sim.println("=====> Total time to BCD = " + dfTime.format(general.deltaTime(startTime) / 60) + " minutes");
                    } else {
                        sim.println("=====>Removing cells");
                        for (int i = 1; i < 5; i++) {
                            meshManager.removeInvalidCells(new NeoObjectVector(regionsNoCores.toArray()), NeoProperty.fromString("{\'function\': \'VelocityMagnitude\', \'functionValue\': " + ((MonitorIterationStoppingCriterionMaxLimitType) monitorIterationStoppingCriterionU.getCriterionType()).getLimit().evaluate() + ", \'functionOperator\': 1," + invalidCriteria));
                        }
                    }
                    ((MonitorIterationStoppingCriterionMaxLimitType) minCriterion.getCriterionType()).getLimit().setValue(sim.getSimulationIterator().getCurrentIteration() + 2);
                }

                sim.getSimulationIterator().run();
                
                if (((AbortFileStoppingCriterion) sim.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Stop File")).getIsSatisfied()) {
                    sim.println("=====> Abort File detected, exiting run_unsteady method");
                    //find command to stop the simualtion here skip all other post processing
                    return;

                }
            }
            sim.println("=====> Exiting RampInstances For loop");
        }

    }

    private void python() {
        long startTime = System.nanoTime();
        sim.println("=======================================");
        sim.println("====>    Running python PostRun.py");
        sim.println("=======================================");

        try {
            String[] command = {"/bin/bash", "-c", resolvePath("../../python/aero/PostRun.py")};
            Process proc = new ProcessBuilder(command).start();

            InputStream is = proc.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line;

            while ((line = br.readLine()) != null) {
                sim.println(line);
            }
            if (line == null) {
                sim.println("The output of the command is null");
            }
            sim.println("Python script finished");
        } catch (IOException ex) {
            sim.println("Error running python script");
            sim.println("Error Message:" + ex);
        }
        sim.println("=====> Total time to run python = " + dfTime.format(general.deltaTime(startTime) / 60) + " minutes");
    }

    private void export_XYZ_csv() {
        long startTime = System.nanoTime();
        sim.println("=======================================");
        sim.println("====>    Export Cell data for Mean CdA and ClA");
        sim.println("=======================================");

        // Attempt to find and delte if previous table existed
        try {
            XyzInternalTable xyzInternalTable_0 = ((XyzInternalTable) sim.getTableManager().getTable("XYZ Cell Data"));
            sim.getTableManager().remove(xyzInternalTable_0);
        } catch (Exception e) {
        }

        // Create XYZ table for surface parts
        XyzInternalTable xyzInternalTable = sim.getTableManager().createTable(XyzInternalTable.class);
        xyzInternalTable.setPresentationName("XYZ Cell Data");

        Collection<Boundary> car_boundaries = collect.collectCarBoundaries(sim);

        FieldFunction Cd = sim.getFieldFunctionManager().getFunction("MeanCdASurfaceMonitor");
        FieldFunction Cl = sim.getFieldFunctionManager().getFunction("MeanClASurfaceMonitor");

        xyzInternalTable.setFieldFunctions(new NeoObjectVector(new Object[]{Cd, Cl}));
        xyzInternalTable.getParts().setObjects(car_boundaries);

        xyzInternalTable.extract();
        xyzInternalTable.export("Cell_Data.csv", ",");
        sim.println("=====> Total time to export XYZ csv = " + dfTime.format(general.deltaTime(startTime) / 60) + " minutes");
    }

    private void createReportsMonitors() {

        sim.println("=======================================");
        sim.println("====>    Create Mean Monitors, Reports, and Plots");
        sim.println("=======================================");

        // Update Paramters from reports
        Collection<Boundary> boundariesCar = collect.collectCarBoundaries(sim);
        Collection<Boundary> boundariesCarAndDomain = collect.collectCarBoundariesAndDomain(sim);
        //Collection<PartSurface> partSurfacesCarAndDomain =collect.collectPartSurfacesCarAndDomain(sim);
        //Collection<PartSurface> partSurfacesCar = collect.collectPartSurfacesFromBoundaries(sim,boundariesCar);

        Collection<NamedObject> regions_boundaries = new ArrayList();
        regions_boundaries.addAll(regions);
        regions_boundaries.addAll(boundariesCarAndDomain);

        Collection<NamedObject> boundariesCar_wdu = new ArrayList();
        boundariesCar_wdu.addAll(boundariesCar);
        boundariesCar_wdu.addAll(collect.collectParts(sim, "^wdu_.*"));
/*
        // Create mean monitors
        function.createScalarFieldFunction(sim, "CdA Surface", "((${Pressure}-${P_ref})*$${Area}[0] +alternateValue($${WallShearStress}[0], 0.0)*mag($${Area}))/(0.5*${Rho_ref}*${U_ref}*${U_ref})");
        function.createScalarFieldFunction(sim, "ClA Surface", "((${Pressure}-${P_ref})*$${Area}[2] +alternateValue($${WallShearStress}[2], 0.0)*mag($${Area}))/(0.5*${Rho_ref}*${U_ref}*${U_ref})");
        function.createScalarFieldFunction(sim, "CyA Surface", "((${Pressure}-${P_ref})*$${Area}[1] +alternateValue($${WallShearStress}[1], 0.0)*mag($${Area}))/(0.5*${Rho_ref}*${U_ref}*${U_ref})");

        function.createMeanMonitor(sim, "CdA Surface", "Mean CdA Surface", "Scalar", boundariesCarAndDomain);
        function.createMeanMonitor(sim, "ClA Surface", "Mean ClA Surface", "Scalar", boundariesCarAndDomain);
        function.createMeanMonitor(sim, "CyA Surface", "Mean CyA Surface", "Scalar", boundariesCarAndDomain);
        
        function.createMeanMonitor(sim, "Pressure", "Mean Pressure", "Scalar", regions_boundaries);
        function.createMeanMonitor(sim, "Total Pressure", "Mean Total Pressure", "Scalar", regions_boundaries);
        function.createMeanMonitor(sim, "Skin Friction Coefficient", "Mean Skin Friction", "Scalar", boundariesCarAndDomain);*/

        function.createMeanMonitor(sim, "Velocity", "Mean Velocity X", "VectorX", regions_boundaries);
        function.createMeanMonitor(sim, "Velocity", "Mean Velocity Y", "VectorY", regions_boundaries);
        function.createMeanMonitor(sim, "Velocity", "Mean Velocity Z", "VectorZ", regions_boundaries);
/*
        function.createScalarFieldFunction(sim, "Pressure Squared", "${Pressure}*${Pressure}");
        function.createMeanMonitor(sim, "Pressure Squared", "Mean Pressure Squared", "Scalar", regions_boundaries);
        function.createScalarFieldFunction(sim, "Pressure RMS", "sqrt(${MeanPressureSquaredMonitor}-${MeanPressureMonitor}*${MeanPressureMonitor})");*/
        
        // TKE
        function.createScalarFieldFunction(sim, "Velocity X Squared", "$${Velocity}[0]*$${Velocity}[0]");
        function.createScalarFieldFunction(sim, "Velocity Y Squared", "$${Velocity}[1]*$${Velocity}[1]");
        function.createScalarFieldFunction(sim, "Velocity Z Squared", "$${Velocity}[2]*$${Velocity}[2]");
        function.createMeanMonitor(sim, "VelocityXSquared", "Mean Velocity X Squared", "Scalar", regions_boundaries);
        function.createMeanMonitor(sim, "VelocityYSquared", "Mean Velocity Y Squared", "Scalar", regions_boundaries);
        function.createMeanMonitor(sim, "VelocityZSquared", "Mean Velocity Z Squared", "Scalar", regions_boundaries);
        function.createMeanMonitor(sim, "Turbulent Kinetic Energy", "Mean TKE", "Scalar", regions_boundaries);
        function.createScalarFieldFunction(sim, "Mean TKE", "${MeanTKEMonitor} + 0.5*((${MeanVelocityXSquaredMonitor}-$${MeanVelocity}[0]*$${MeanVelocity}[0]) + (${MeanVelocityYSquaredMonitor}-$${MeanVelocity}[1]*$${MeanVelocity}[1]) + (${MeanVelocityZSquaredMonitor}-$${MeanVelocity}[2]*$${MeanVelocity}[2]))");
        
        // Q-criterion (built from averaged fields)
        function.createVectorFieldFunction(sim, "gradVi", "grad($${MeanVelocity}[0])");
        function.createVectorFieldFunction(sim, "gradVj", "grad($${MeanVelocity}[1])");
        function.createVectorFieldFunction(sim, "gradVk", "grad($${MeanVelocity}[2])");
        function.createScalarFieldFunction(sim, "Mean Q-criterion", "-0.5 * ($${gradVi}[0]*$${gradVi}[0] + $${gradVj}[1]*$${gradVj}[1] + $${gradVk}[2]*$${gradVk}[2]) -$${gradVi}[1]*$${gradVj}[0] - $${gradVi}[2]*$${gradVk}[0] -$${gradVj}[2]*$${gradVk}[1]");
        function.createScalarFieldFunction(sim, "Mean Q-criterion Iso", "$WallDistance > 0.005 ? ${MeanQ-criterion} : -9999.9");
        /*
        // Surface force plot values
        function.createScalarFieldFunction(sim, "Mean ClrA Surface", "(${MeanClASurfaceMonitor}*($${Position}-[$${Wheel_front-left_origin}[0], 0 , ${ground}])[0] - ${MeanCdASurfaceMonitor}*($${Position}-[$${Wheel_front-left_origin}[0], 0 , ${ground}])[2]) / ${WheelBase}");
        function.createScalarFieldFunction(sim, "Mean ClfA Surface", "(-${MeanClASurfaceMonitor}*($${Position}-[$${Wheel_rear-left_origin}[0], 0 , ${ground}])[0] + ${MeanCdASurfaceMonitor}*($${Position}-[$${Wheel_rear-left_origin}[0], 0 , ${ground}])[2]) / ${WheelBase}");
        function.createScalarFieldFunction(sim, "Mean CdA Plot", "${MeanCdASurfaceMonitor}/mag($${Area})");
        function.createScalarFieldFunction(sim, "Mean ClA Plot", "${MeanClASurfaceMonitor}/mag($${Area})");
        function.createScalarFieldFunction(sim, "Mean CyA Plot", "${MeanCyASurfaceMonitor}/mag($${Area})");
        function.createScalarFieldFunction(sim, "Mean ClrA Plot", "${MeanClrASurface}/mag($${Area})");
        function.createScalarFieldFunction(sim, "Mean ClfA Plot", "${MeanClfASurface}/mag($${Area})");
        */
               
        //function.createMeanMonitors(sim, "Wall Shear Stress", "Mean Wall Shear X", "VectorX", false);
        //function.createMeanMonitors(sim, "Wall Shear Stress", "Mean Wall Shear Y", "VectorY", false);
        //function.createMeanMonitors(sim, "Wall Shear Stress", "Mean Wall Shear Z", "VectorZ", false);

        //function.createMeanMonitors(sim, "Lambda2", "Mean Lambda2", "Scalar", true);
        //function.createMeanMonitors(sim, "Qcriterion", "Mean Qcriterion", "Scalar", true);
        // create additional field functions
        function.createVectorFieldFunction(sim, "Mean Velocity", "[${MeanVelocityXMonitor},${MeanVelocityYMonitor},${MeanVelocityZMonitor}]");
        //function.createVectorFieldFunction(sim, "Mean Wall Shear Stress", "[${MeanWallShearXMonitor},${MeanWallShearYMonitor},${MeanWallShearZMonitor}]");
/*
        function.createScalarFieldFunction(sim, "Mean Cp", "${MeanPressureMonitor}/(0.5*${Rho_ref}*${U_ref}*${U_ref})");
        function.createScalarFieldFunction(sim, "Mean Cp Tot", "${MeanTotalPressureMonitor}/(0.5*${Rho_ref}*${U_ref}*${U_ref})");
        function.createScalarFieldFunction(sim, "Mean Cp Tot Iso", "$WallDistance > 0.001 ? $MeanCpTot : 0.2");
        
        // Crossflow
        function.createVectorFieldFunction(sim, "Mean Crossflow", "[0,${MeanVelocityYMonitor},${MeanVelocityZMonitor}]");
        function.createVectorFieldFunction(sim, "Mean Crossflow Iso", "$MeanCpTot < .9 ? $${MeanCrossflow} : [0.0, 0.0, 0.0]");
        
        // CdA
        function.createScalarFieldFunction(sim, "Mean CdA", "-${MeanCp}-2*((${MeanVelocityXSquaredMonitor}/(${U_ref}*${U_ref}))-(${MeanVelocityXMonitor}/(${U_ref})))");
*/
        function.createScalarFieldFunction(sim, "Mean Vorticity X", "alternateValue(grad(${MeanVelocityZMonitor})[1] - grad(${MeanVelocityYMonitor})[2] , 0.0)");
        function.createScalarFieldFunction(sim, "Mean Vorticity Y", "alternateValue(grad(${MeanVelocityXMonitor})[2] - grad(${MeanVelocityZMonitor})[0] , 0.0)");
        function.createScalarFieldFunction(sim, "Mean Vorticity Z", "alternateValue(grad(${MeanVelocityYMonitor})[0] - grad(${MeanVelocityXMonitor})[1] , 0.0)");
        function.createVectorFieldFunction(sim, "Mean Vorticity", "[${MeanVorticityX},${MeanVorticityY},${MeanVorticityZ}]");

        function.createScalarFieldFunction(sim, "Face Angle X", "alternateValue(asin(-$${Normal}[0])*180/3.141592653589794 , 0.0)");
        function.createScalarFieldFunction(sim, "Face Angle Y", "alternateValue(asin(-$${Normal}[1])*180/3.141592653589794 , 0.0)");
        function.createScalarFieldFunction(sim, "Face Angle Z", "alternateValue(asin(-$${Normal}[2])*180/3.141592653589794 , 0.0)");
        function.createVectorFieldFunction(sim, "Face Angle", "[${FaceAngleX},${FaceAngleY},${FaceAngleZ}]");
/*
        function.createScalarFieldFunction(sim, "Force X", "alternateValue(${MeanCdASurfaceMonitor}*(0.5*${Rho_ref}*${U_ref}*${U_ref}) , 0.0)");
        function.createScalarFieldFunction(sim, "Force Y", "alternateValue(${MeanCyASurfaceMonitor}*(0.5*${Rho_ref}*${U_ref}*${U_ref}) , 0.0)");
        function.createScalarFieldFunction(sim, "Force Z", "alternateValue(${MeanClASurfaceMonitor}*(0.5*${Rho_ref}*${U_ref}*${U_ref}) , 0.0)");
        function.createVectorFieldFunction(sim, "Force", "[${ForceX},${ForceY},${ForceZ}]");

        // create reports             
        function.createForceCoefficientReport(sim, "CdA", boundariesCar, "[1,0,0]", "${A_ref}", true);
        function.createForceCoefficientReport(sim, "Cd", boundariesCar, "[1,0,0]", "${FrontalArea}", true);
        function.createForceCoefficientReport(sim, "ClA", boundariesCar, "[0,0,1]", "${A_ref}", true);
        function.createForceCoefficientReport(sim, "Cl", boundariesCar, "[0,0,1]", "${FrontalArea}", true);
        function.createMomentCoefficientReport(sim, "Clf", boundariesCar, "[0,1,0]", "[$${Wheel_rear-left_origin}[0], 0 , ${ground}]", "${FrontalArea}", true);
        function.createMomentCoefficientReport(sim, "Clr", boundariesCar, "[0,-1,0]", "[$${Wheel_front-left_origin}[0], 0, ${ground}]", "${FrontalArea}", true);
        function.createMomentCoefficientReport(sim, "Clf_wdu", boundariesCar_wdu, "[0,1,0]", "[$${Wheel_rear-left_origin}[0], 0 , ${ground}]", "${FrontalArea}", true);
        function.createMomentCoefficientReport(sim, "Clr_wdu", boundariesCar_wdu, "[0,-1,0]", "[$${Wheel_front-left_origin}[0], 0, ${ground}]", "${FrontalArea}", true);*/

        function.createMaxReport(sim, "TVRMax", regions, "TurbulentViscosityRatio", true);
        function.createMaxReport(sim, "UMax", regions, "Velocity", true);
        function.createMinReport(sim, "VolMin", regions, "Volume", false);

        //function.createSIET(sim);

        //function.createFrontalAreaReport(sim, "VehicleFrontalArea", boundariesCar);

        //function.createSumReport(sim, "MeanCdASurfaceSum", boundariesCar, "MeanCdASurfaceMonitor");

        // Monitors and reports for ventilation moment are slowing down the simualtion, they will be replaced by a report from mean surface pressures instead
        /*
        String forPlot="";
        for (String wheel : wheels) {
            Collection<Boundary> boundariesWheel = collect.collectBoundaries(boundariesCar, ".*" + wheel + ".*");
            
            LabCoordinateSystem labCSys = sim.getCoordinateSystemManager().getLabCoordinateSystem();
            CoordinateSystem cSys = ((CoordinateSystem) labCSys.getLocalCoordinateSystemManager().getObject("Wheel_" + wheel));
            function.createMomentReport(sim, "Ventilation Moment " + wheel, boundariesWheel, "[1,0,0]", "[0,0,0]", cSys, false);
            
            double ground = ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("ground")).getQuantity().getInternalValue();
            double distance = cSys.getOriginVector().getComponent(2) - ground;
            function.createMomentToForceReport(sim, "Ventilation Force " + wheel, "Ventilation Moment " + wheel, String.valueOf(distance), false);
            function.createForceToCoefficientReport(sim, "Ventilation CdA " + wheel, "Ventilation Force " + wheel, "1", true);
        
            function.createStatisticsReport(sim, "MeanVentilation CdA " + wheel, "Ventilation CdA " + wheel + " Monitor", "Mean", "AveragingStart");
            forPlot = forPlot + "Ventilation CdA " + wheel + " Monitor,MeanVentilation CdA " + wheel + " Monitor,";  
        }
        function.createPlot(sim, "Ventilation Force CdAs Plot", forPlot);
         
        // final reports for ventilation moments obtaines from surface averaged pressure
        for (String wheel : wheels) {
            Collection<Boundary> boundariesWheel = collect.collectBoundaries(boundariesCar, ".*" + wheel + ".*");
            CoordinateSystem cSys = (((CoordinateSystem) sim.getCoordinateSystemManager().getLabCoordinateSystem()).getLocalCoordinateSystemManager().getObject("Wheel_" + wheel));
            // Create expression reports here to calcualte ventilation moment from final mean Force vector
        }
*/
        if (detailedMethodEvaluation) {
            function.createMeanMonitor(sim, "Turbulent Viscosity Ratio", "Mean Turbulent Viscosity Ratio", "Scalar", regions_boundaries);
            function.createMeanMonitor(sim, "Turbulent Kinetic Energy", "Mean Modelled Turbulent Kinetic Energy", "Scalar", regions_boundaries);
            function.createMeanMonitor(sim, "WallYplus", "Mean Wall Y+", "Scalar", boundariesCarAndDomain);

            //function.createForceCoefficientReport(sim, "Cd Iteration", boundariesCar, "[1,0,0]", "${FrontalArea}", true);
            //function.createMomentCoefficientReport(sim, "Clf Iteration", boundariesCar, "[0,1,0]", "[$${Wheel_rear-left_origin}[0], 0 , ${ground}]", "${FrontalArea}", true);
            //function.createMomentCoefficientReport(sim, "Clr Iteration", boundariesCar, "[0,-1,0]", "[$${Wheel_front-left_origin}[0], 0, ${ground}]", "${FrontalArea}", true);
        }

        // Set latest representation for reports
        LatestMeshProxyRepresentation latestMeshProxyRepresentation = ((LatestMeshProxyRepresentation) sim.getRepresentationManager().getObject("Latest Surface/Volume"));

        Collection<Report> reportCollection = sim.getReportManager().getObjects();
        for (Report thisReport : reportCollection) {
            thisReport.getReportManager().applyRepresentation(latestMeshProxyRepresentation);
        }

        //update the Frontal Area
        //((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("FrontalArea")).getQuantity().setValue(sim.getReportManager().getObject("VehicleFrontalArea").getReportMonitorValue());
/*
        // create statistics report
        function.createStatisticsReport(sim, "MeanCdA", "CdA Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanCd", "Cd Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanClA", "ClA Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanCl", "Cl Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanClf", "Clf Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanClr", "Clr Monitor", "Mean", "AveragingStart");

        function.createStatisticsReport(sim, "MeanClf_wdu", "Clf_wdu Monitor", "Mean", "AveragingStart");
        function.createStatisticsReport(sim, "MeanClr_wdu", "Clr_wdu Monitor", "Mean", "AveragingStart");

        function.createStatisticsReport(sim, "MeanSIET", "SIET Monitor", "Mean", "AveragingStart");

        // Create plots
        function.createPlot(sim, "CdA Plot", "CdA Monitor,MeanCdA Monitor");
        function.createPlot(sim, "Cd Plot", "Cd Monitor,MeanCd Monitor");
        function.createPlot(sim, "Cl Plot", "Clr Monitor,MeanClr Monitor,Clf Monitor,MeanClf Monitor,Cl Monitor,MeanCl Monitor");
        function.createPlot(sim, "Cl_wdu Plot", "Clr_wdu Monitor,MeanClr_wdu Monitor,Clf_wdu Monitor,MeanClf_wdu Monitor");*/
        function.createPlot(sim, "TVRMax Plot", "TVRMax Monitor");
        function.createPlot(sim, "UMax Plot", "UMax Monitor");
        //function.createPlot(sim, "SIET Plot", "SIET Monitor");

        //function.createPlot(sim, "_CdA Plot", "CdA Monitor,MeanCdA Monitor", true);
        //function.createPlot(sim, "_Cd Plot", "Cd Monitor,MeanCd Monitor", true);
        //function.createPlot(sim, "_Cl Plot", "Clr Monitor,MeanClr Monitor,Clf Monitor,MeanClf Monitor", true);
        //function.createPlot(sim, "_Cl_wdu Plot", "Clr_wdu Monitor,MeanClr_wdu Monitor,Clf_wdu Monitor,MeanClf_wdu Monitor", true);

        if (PVT) {
            //pvt.createPVTFieldFunctions(sim); The pvt field functions are created in prep (before pvt boundaries) in order to set the boundary conditions for suction etc...
            pvt.createPVTReports(sim);
            pvt.createPVTPlots(sim);
        }
    }

    private void setGeneralSolverSettings() {
        ((WallDistanceSolver) sim.getSolverManager().getSolver(WallDistanceSolver.class)).setFrozen(true);
        ((PartitioningSolver) sim.getSolverManager().getSolver(PartitioningSolver.class)).getPartitioningOption().setSelected(PartitioningOption.Type.PER_CONTINUUM);
    }

    private void setFieldFunctionRef(Simulation sim) {

        ((SkinFrictionCoefficientFunction) sim.getFieldFunctionManager().getFunction("SkinFrictionCoefficient")).getReferenceDensity().setDefinition("${Rho_ref}");
        ((SkinFrictionCoefficientFunction) sim.getFieldFunctionManager().getFunction("SkinFrictionCoefficient")).getReferenceVelocity().setDefinition("${U_ref}");

        ((TotalPressureCoefficientFunction) sim.getFieldFunctionManager().getFunction("TotalPressureCoefficient")).getReferenceDensity().setDefinition("${Rho_ref}");
        ((TotalPressureCoefficientFunction) sim.getFieldFunctionManager().getFunction("TotalPressureCoefficient")).getReferenceVelocity().setDefinition("${U_ref}");
        ((TotalPressureCoefficientFunction) sim.getFieldFunctionManager().getFunction("TotalPressureCoefficient")).getReferencePressure().setDefinition("${P_ref}");

        ((PressureCoefficientFunction) sim.getFieldFunctionManager().getFunction("PressureCoefficient")).getReferenceDensity().setDefinition("${Rho_ref}");
        ((PressureCoefficientFunction) sim.getFieldFunctionManager().getFunction("PressureCoefficient")).getReferenceVelocity().setDefinition("${U_ref}");
        ((PressureCoefficientFunction) sim.getFieldFunctionManager().getFunction("PressureCoefficient")).getReferencePressure().setDefinition("${P_ref}");

    }

    private void set_compressible() {

        sim.println("=======================================");
        sim.println("====>    Setting Compressible Flow");
        sim.println("=======================================");

        PhysicsContinuum continuum = ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));
        continuum.getModelManager().getModel(ConstantDensityModel.class).disable();
        continuum.enable(IdealGasModel.class);
        continuum.enable(SegregatedFluidIsothermalModel.class);
        SegregatedFluidIsothermalModel segregatedFluidIsothermalModel = continuum.getModelManager().getModel(SegregatedFluidIsothermalModel.class);
        segregatedFluidIsothermalModel.getContinuumTemperature().setValue(293.0);
        function.createMaxReport(sim, "rhoMax", regions, "Density");
        function.createMinReport(sim, "rhoMin", regions, "Density");
    }

    private void createSolutionHistory(double[] historyPlanesX, double[] historyPlanesY, double[] historyPlanesZ, double[] planeLimits, int solutionHistoryFrequency, int solutionHistoryStart) {
        sim.println("=======================================");
        sim.println("====>    Create Solution History");
        sim.println("=======================================");

        //SolutionHistory solutionHistory = sim.get(SolutionHistoryManager.class).createForFile(resolvePath("solution_history_planes.simh"), false, false);
        SolutionHistory solutionHistory = sim.get(SolutionHistoryManager.class).createForFile("solution_history_planes.simh", false, false);

        PrimitiveFieldFunction velocityFF = ((PrimitiveFieldFunction) sim.getFieldFunctionManager().getFunction("Velocity"));
        VectorComponentFieldFunction velocityXFF = ((VectorComponentFieldFunction) velocityFF.getComponentFunction(0));
        VectorComponentFieldFunction velocityYFF = ((VectorComponentFieldFunction) velocityFF.getComponentFunction(1));
        VectorComponentFieldFunction velocityZFF = ((VectorComponentFieldFunction) velocityFF.getComponentFunction(2));
        solutionHistory.setFunctions(new NeoObjectVector(new Object[]{velocityXFF, velocityYFF, velocityZFF}));

        StarUpdate starUpdate = solutionHistory.getUpdate();

        solutionHistory.getUpdate().getUpdateModeOption().setSelected(StarUpdateModeOption.Type.TIMESTEP);
        solutionHistory.getUpdate().getTimeStepUpdateFrequency().setStart(solutionHistoryStart);
        solutionHistory.getUpdate().getIterationUpdateFrequency().setIterations(solutionHistoryFrequency);

        for (double plane : historyPlanesX) {
            DoubleVector origin = new DoubleVector(new double[]{plane, 0.0, 0.0});
            DoubleVector xaxis = new DoubleVector(new double[]{1.0, 0.0, 0.0});
            double[] loop = {plane, planeLimits[3], planeLimits[4],
                plane, planeLimits[3], planeLimits[5],
                plane, planeLimits[2], planeLimits[5],
                plane, planeLimits[2], planeLimits[4]};

            ConstrainedPlaneSection constrainedPlaneSection = analyze.createConstrainedPlane(sim, "CPS_x=" + plane + "m", new NeoObjectVector(regions.toArray()), origin, xaxis, new DoubleVector(loop));
            solutionHistory.getInputs().addObjects(constrainedPlaneSection);
        }

        for (double plane : historyPlanesY) {
            DoubleVector origin = new DoubleVector(new double[]{0.0, plane, 0.0});
            DoubleVector yaxis = new DoubleVector(new double[]{0.0, 1.0, 0.0});
            double[] loop = {planeLimits[0], plane, planeLimits[4],
                planeLimits[0], plane, planeLimits[5],
                planeLimits[1], plane, planeLimits[5],
                planeLimits[1], plane, planeLimits[4]};

            ConstrainedPlaneSection constrainedPlaneSection = analyze.createConstrainedPlane(sim, "CPS_y=" + plane + "m", new NeoObjectVector(regions.toArray()), origin, yaxis, new DoubleVector(loop));
            solutionHistory.getInputs().addObjects(constrainedPlaneSection);
            solutionHistory.getRegions().setObjects();  // Do not include any regions
        }

        for (double plane : historyPlanesZ) {
            DoubleVector origin = new DoubleVector(new double[]{0.0, 0.0, plane});
            DoubleVector zaxis = new DoubleVector(new double[]{0.0, 0.0, 1.0});
            double[] loop = {planeLimits[0], planeLimits[2], plane,
                planeLimits[0], planeLimits[3], plane,
                planeLimits[1], planeLimits[3], plane,
                planeLimits[1], planeLimits[2], plane};

            ConstrainedPlaneSection constrainedPlaneSection = analyze.createConstrainedPlane(sim, "CPS_z=" + plane + "m", new NeoObjectVector(regions.toArray()), origin, zaxis, new DoubleVector(loop));
            solutionHistory.getInputs().addObjects(constrainedPlaneSection);
        }
    }

    private void createPointMonitors(double[] pointMonitorsXCoords, double[] pointMonitorsYCoords, double[] pointMonitorsZCoords) {
        sim.println("=======================================");
        sim.println("====>    Create Point Monitors");
        sim.println("=======================================");

        if (pointMonitorsXCoords.length == pointMonitorsYCoords.length && pointMonitorsXCoords.length == pointMonitorsZCoords.length) {
            for (int i = 0; i < pointMonitorsXCoords.length; i++) {
                DoubleVector coordinate = new DoubleVector(new double[]{pointMonitorsXCoords[i], pointMonitorsYCoords[i], pointMonitorsZCoords[i]});
                PointPart point = analyze.createPoint(sim, "Point", new NeoObjectVector(regions.toArray()), coordinate);

                Collection<NamedObject> pointCollection = new ArrayList();
                pointCollection.add(point);
                String ff = "Velocity";
                String name = ff + " " + point.getPresentationName() + " Iteration";    // Add iteration to trigger monitor on iteration instead of time step
                function.createMaxReport(sim, name, pointCollection, ff);
            }
        } else {
            sim.println("The lengths of the coordinate vectors are not equal");
        }
    }

    private void createWDUReports(boolean PVT) {

        String forPlot = "";

        if (PVT) {
            for (String wheel : wheels) {
                String match = ".*wdu_" + wheel + ".*belt.*";
                Collection<PartSurface> WDUPart = collect.collectPartSurfacesFromParts(sim, ".*Subtract.*", match);
                function.createForceCoefficientReport(sim, "Cl WDU " + wheel, WDUPart, "[0,0,1]", "${FrontalArea}", true);
                function.createStatisticsReport(sim, "MeanCl WDU " + wheel, "Cl WDU " + wheel + " Monitor", "Mean", "AveragingStart");
                function.createPlot(sim, "Cl WDU Plot " + wheel, "Cl WDU "+wheel+" Monitor,MeanCl WDU "+wheel+" Monitor");
                function.createPlot(sim, "_Cl WDU Plot " + wheel, "Cl WDU "+wheel+" Monitor,MeanCl WDU "+wheel+" Monitor", true);

                /*    // the update of mesh representation is done later after all other reports are created hence it is not needed here
                LatestMeshProxyRepresentation latestMeshProxyRepresentation = ((LatestMeshProxyRepresentation) sim.getRepresentationManager().getObject("Latest Surface/Volume"));
                report.getReportManager().applyRepresentation(latestMeshProxyRepresentation);
                meanReport.getReportManager().applyRepresentation(latestMeshProxyRepresentation);
                 */
                
            }
        } else {

            String fieldFunctionName = "ZThresholdWheel";
            function.createScalarFieldFunction(sim, fieldFunctionName, "${ground}<$${Position}[2] && (${ground}+0.1)>$${Position}[2]");

            for (String wheel : wheels) {

                //double[] centerCoordinate =((VectorGlobalParameter) sim.get(GlobalParameterManager.class).getObject("ContactPatchCenter_"+wheel)).getQuantity().evaluate().toDoubleArray(); 
                double[] corner1 = ((SimpleBlockPart) sim.get(SimulationPartManager.class).getPart("ContactPatchRefinement_" + wheel)).getCorner1().evaluate().toDoubleArray();
                double[] corner2 = ((SimpleBlockPart) sim.get(SimulationPartManager.class).getPart("ContactPatchRefinement_" + wheel)).getCorner2().evaluate().toDoubleArray();
                function.createScalarFieldFunction(sim, "BoxFunction_WDU_" + wheel, corner1[0] + "<$${Position}[0] && $${Position}[0]<" + corner2[0] + " && " + corner1[1] + "<$${Position}[1] && $${Position}[1]<" + corner2[1]);
                Collection<Boundary> groundBoundary = collect.collectBoundaries(collect.collectCarBoundariesAndDomain(sim), "(?i).*ground.*");  // (?i) for case insensativity
                ThresholdPart groundThreshold = analyze.createThresholdAbove(sim, "wdu_" + wheel, new NeoObjectVector(groundBoundary.toArray()), "BoxFunction_WDU_" + wheel, 0.0, 0.9);
                /*
                function.createForceCoefficientReport(sim, "Cl WDU " + wheel, Arrays.asList(groundThreshold), "[0,0,1]", "${FrontalArea}",true);
                Report report = sim.getReportManager().getObject("Cl WDU " + wheel);
                function.createStatisticsReport(sim, "Mean Cl WDU " + wheel, "Cl WDU " + wheel + " Monitor", "Mean", "AveragingStart");
                Report meanReport = sim.getReportManager().getObject("Mean Cl WDU " + wheel);

                LatestMeshProxyRepresentation latestMeshProxyRepresentation = ((LatestMeshProxyRepresentation) sim.getRepresentationManager().getObject("Latest Surface/Volume"));
                report.getReportManager().applyRepresentation(latestMeshProxyRepresentation);
                meanReport.getReportManager().applyRepresentation(latestMeshProxyRepresentation);

                forPlot = "Cl WDU " + wheel + " Monitor,Mean Cl WDU " + wheel + " Monitor";
                function.createPlot(sim, "Cl WDU Plot " + wheel, forPlot);*/
            }

            sim.getFieldFunctionManager().remove(sim.getFieldFunctionManager().getFunction(fieldFunctionName));
        }
    }

    private void exportWduImages(boolean PVT) {
        sim.println("Exporting WDU pics");
        Scene scene = analyze.createScene(sim, "WDU Surface");

        ScalarDisplayer scalarDisplayer = scene.getDisplayerManager().createScalarDisplayer("Scalar");
        scalarDisplayer.setPresentationName("Scalar");
        scalarDisplayer.setRepresentation(sim.getRepresentationManager().getObject("Latest Surface/Volume"));
        Legend legendScalar = scalarDisplayer.getLegend();
        legendScalar.setLevels(11);
        legendScalar.setLabelFormat("%-2.2g");
        legendScalar.setLookupTable(sim.get(LookupTableManager.class).getObject("cool-warm"));
        legendScalar.updateLayout(new DoubleVector(new double[]{0.9, 0.1}), 0.02, 0.6, 1);

        FieldFunction ff = sim.getFieldFunctionManager().getFunction("MeanCp");
        scalarDisplayer.getScalarDisplayQuantity().setFieldFunction(ff);

        scalarDisplayer.getScalarDisplayQuantity().setClip(ClipMode.NONE);
        scalarDisplayer.getScalarDisplayQuantity().setAutoRange(AutoRangeMode.NONE);
        scalarDisplayer.getScalarDisplayQuantity().setRangeMin(-1.0);
        scalarDisplayer.getScalarDisplayQuantity().setRangeMax(1.0);
        
        scalarDisplayer.setFillMode(ScalarFillMode.NODE_FILLED);

        DoubleVector viewFocal = new DoubleVector(new double[]{0, 0, 0});
        DoubleVector viewPosition = new DoubleVector(new double[]{-1, 0, 0});
        DoubleVector viewZ = new DoubleVector(new double[]{0, 0, 1});
        double zoom = 1.0;
        double extrazoom = 0.0;

        CurrentView currentView = scene.getCurrentView();
        currentView.setInput(viewFocal, viewPosition, viewZ, zoom, 1, 30.0);

        for (String wheel : wheels) {
            if (PVT) {
                String match = ".*wdu_" + wheel + ".*belt.*";
                Collection<PartSurface> WDUPart = collect.collectPartSurfacesFromParts(sim, ".*Subtract.*", match);
                scalarDisplayer.getInputParts().setObjects(WDUPart);
            } else {
                ThresholdPart WDUThreshold = (ThresholdPart) sim.getPartManager().getObject("wdu_" + wheel);
                scalarDisplayer.getInputParts().setObjects(WDUThreshold);
            }

            viewFocal = new DoubleVector(new double[]{4, 0, 0});
            viewPosition = new DoubleVector(new double[]{4, 0, 1});
            viewZ = new DoubleVector(new double[]{0, 1, 0});
            zoom = 2.0;

            currentView.setInput(viewFocal, viewPosition, viewZ, zoom, 1, 30.0);
            scene.resetCamera();

            viewFocal = currentView.getFocalPoint();
            viewPosition = currentView.getPosition();
            ParallelScale ps = currentView.getParallelScale();
            double zooom = ps.getValue();
            zoom = zooom + extrazoom;
            currentView.setInput(viewFocal, viewPosition, viewZ, zoom, 1, 30.0);

            scene.printAndWait(("Check/wdu_" + wheel + ".png"), 1, 3840, 2160, true, false);
        }

        scene.close(true);
        sim.getSceneManager().deleteScene(scene);
        sim.println("WDU pics done");
    }

    private void setDetailedME() {
        double[] historyPlanesX = {};
        double[] historyPlanesY = {0.0};
        double[] historyPlanesZ = {1.0};
        double[] planeLimits = {0.0, 8.0, -1.5, 1.5, 0.0, 2.0}; // Limits for constrained planes minX, maxX, minY, maxY, minZ, maxZ
        int solutionHistoryFrequency = 1;
        int solutionHistoryStart = 1994;   // Time step to start collection solution history
        double[] pointMonitorsXCoords = {6.0, 6.0, 6.0, 6.0, 6.0, 6.0};
        double[] pointMonitorsYCoords = {0.0, 0.0, 0.0, 0.0, -0.5, 0.5};
        double[] pointMonitorsZCoords = {0.5, 0.75, 1.0, 1.25, 1.0, 1.0};
        createSolutionHistory(historyPlanesX, historyPlanesY, historyPlanesZ, planeLimits, solutionHistoryFrequency, solutionHistoryStart);
        createPointMonitors(pointMonitorsXCoords, pointMonitorsYCoords, pointMonitorsZCoords);
    }

    private void setCabinConditions(double U_ref, double T_ref, double T_wall) {
        PhysicsContinuum physicsContinuum_0 = 
      ((PhysicsContinuum) sim.getContinuumManager().getContinuum("Physics 1"));

    physicsContinuum_0.enable(SegregatedFluidTemperatureModel.class);

    physicsContinuum_0.enable(BoussinesqModel.class);

    physicsContinuum_0.enable(GravityModel.class);

    ScalarGlobalParameter scalarGlobalParameter_0 = 
      ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("U_ref"));

    scalarGlobalParameter_0.getQuantity().setValue(U_ref);

    Units units_0 = 
      ((Units) sim.getUnitsManager().getObject(""));

    scalarGlobalParameter_0.getQuantity().setUnits(units_0);

    ScalarGlobalParameter scalarGlobalParameter_1 = 
      ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("T_ref"));

    scalarGlobalParameter_1.getQuantity().setValue(T_ref);

    scalarGlobalParameter_1.getQuantity().setUnits(units_0);

    ScalarGlobalParameter scalarGlobalParameter_2 = 
      ((ScalarGlobalParameter) sim.get(GlobalParameterManager.class).getObject("T_wall"));

    scalarGlobalParameter_2.getQuantity().setValue(T_wall);

    Units units_1 = 
      ((Units) sim.getUnitsManager().getObject("K"));

    scalarGlobalParameter_2.getQuantity().setUnits(units_1);

    Region region_0 = 
      sim.getRegionManager().getRegion("AirStream");

    Boundary boundary_0 = 
      region_0.getBoundaryManager().getBoundary("DOMAIN_OR.outlet");

    StaticTemperatureProfile staticTemperatureProfile_0 = 
      boundary_0.getValues().get(StaticTemperatureProfile.class);

    staticTemperatureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("$T_ref");

    staticTemperatureProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_1);

    Boundary boundary_1 = 
      region_0.getBoundaryManager().getBoundary("DOMAIN_OR.inlet");

    StaticTemperatureProfile staticTemperatureProfile_1 = 
      boundary_1.getValues().get(StaticTemperatureProfile.class);

    staticTemperatureProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition("$T_ref");

    staticTemperatureProfile_1.getMethod(ConstantScalarProfileMethod.class).getQuantity().setUnits(units_1);
    }
}
