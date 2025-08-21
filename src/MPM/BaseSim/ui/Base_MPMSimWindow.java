package MPM.BaseSim.ui;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import MPM.BaseSim.sim.Base_MPMSim;
import MPM.BaseSim.sim.SimResetProcess;
import MPM.BaseSim.utils.MPM_SimUpdateFromUIData;
import base_Math_Objects.vectorObjs.doubles.myPoint;
import base_Render_Interface.IGraphicsAppInterface;
import base_UI_Objects.GUI_AppManager;
import base_UI_Objects.baseApp.GUI_AppUIFlags;
import base_UI_Objects.windowUI.base.Base_DispWindow;
import base_UI_Objects.windowUI.base.GUI_AppWinVals;
import base_UI_Objects.windowUI.drawnTrajectories.DrawnSimpleTraj;
import base_UI_Objects.windowUI.uiData.UIDataUpdater;
import base_UI_Objects.windowUI.uiObjs.base.GUIObj_Params;
import base_Utils_Objects.io.messaging.MsgCodes;

public abstract class Base_MPMSimWindow extends Base_DispWindow {

    //simulator world within which simulation executes
    //TODO support multiple sim worlds possibly within the same window, with different configurations
    protected Base_MPMSim currSim;
        
    //////////////////////////////////
    //initial values of simulation variables
    //ints
    protected final int initNumGridCellsPerDim = 100;
    protected final int initSimStepsPerFrame = 10;
    protected final int initNumBalls = 2;
    protected final int initNumParts = 100000;
    protected final int initDrawPtIncr = 1;
    //floats
    protected final float initDeltaT = 4e-4f;
    protected final float initCellSize = .10f;
    protected final float initVel = 30.0f;
    
    //snow density varies from 50 to ~800
    //powdery snow is about 100
    //wet firm compacted snow is about 600
    protected final float initParticleDensity = 100.0f;
    protected final float initParticleMass = initCellSize *initCellSize *initCellSize * initParticleDensity; 
    protected final float initWallFric = 1.0f;
    protected final float initColFric = 1.0f;
    protected final float initDrawnVecScale = .01f;
    
    //initial values of material quantities
    protected final float
        //Assuming the snow is isotropic
        //Units are kg / (m * s^2) 
        init_initYoungMod              = 4.8e4f,//1.4e5f,
        init_poissonRatio              = 0.2f,
        //Governs strain hardening, where the material gets stronger past a certain amount of stress in plastic deformation
        init_hardeningCoeff          = 15.0f, 
        init_criticalCompression      = 0.040f, 
        init_criticalStretch          = 0.0075f, 
        init_alphaPicFlip              = 0.95f;
    
    
    //////////////////////////////////////
    //ui vals
    public final static int
        gIDX_TimeStep                 = 0,        
        gIDX_SimStepsPerFrame        = 1,
        gIDX_NumParticles            = 2,
        gIDX_NumSnowballs             = 3,
        gIDX_InitVel                = 4,
        gIDX_PartMass                = 5,
        gIDX_GridCellSize            = 6,
        gIDX_GridCount                = 7,
        gIDX_InitYoungMod             = 8,
        gIDX_PoissonRatio             = 9,
        gIDX_HardeningCoeff         = 10,
        gIDX_CriticalCompression     = 11,
        gIDX_CriticalStretch         = 12,
        gIDX_AlphaPicFlip             = 13,
        gIDX_wallFricCoeff             = 14,
        gIDX_CollFricCoeff            = 15,
        gIDX_DrawnValScale            = 16,
        gIDX_DrawPointIncr            = 17;
    
    public final static int numBaseMPMWinUIGuiObjs = 18;
    
    //////////////////////////////////////
    //custom debug/function ui button names -empty will do nothing

    //////////////////////////////////////
    //private child-class flags - window specific
    public static final int 
        //idx 0 is debug in flags structure
        resetSimIDX                = 1,            //whether or not to reset sim    
        rebuildSimIDX            = 2,            //reset of sim is required for best results
        showLocColors            = 3,            //display particles by color of their initial location
        showCollider            = 4,            //show collider cylinder
        showParticles            = 5,
        showParticleVelArrows     = 6,            //plot velocity arrows for each particle    
        showGrid                = 7,             //plot the computational grid
        showGridVelArrows         = 8,            //plot velocity arrows for each gridNode
        showGridAccelArrows     = 9,            //plot acceleration arrows for each gridNode
        showGridMass              = 10,            //plot variable sized spheres proportional to gridnode mass
        showActiveNodes          = 11,                //show the grid nodes influenced by each particle
        showExecTime            = 12;            //show time of execution for each step in simulation

    public static final int numBaseMPMWinUIPrivFlags = 12;
        
    public Base_MPMSimWindow(IGraphicsAppInterface _p, GUI_AppManager _AppMgr, int _winIdx) {
        super(_p, _AppMgr, _winIdx);
    }//Base_MPMSimWindow

    /**
     * Initialize any UI control flags appropriate for window application
     * @param appUIFlags Snapshot of the initial flags structure for the application. 
     * Will not reflect future changes, so should not be retained
     */
    @Override
    protected final void initDispFlags(GUI_AppUIFlags appUIFlags) {
        //this window is runnable
        dispFlags.setIsRunnable(true);
        //this window uses a customizable camera
        dispFlags.setUseCustCam(true);
        // capable of using right side menu
        dispFlags.setHasRtSideInfoDisp(true);
        //any app-specific disp flags to set
        initDispFlags_Indiv(appUIFlags);
    }
    /**
     * Initialize any UI control flags appropriate for specific instanced MPM window
     * @param appUIFlags Snapshot of the initial flags structure for the application. 
     * Will not reflect future changes, so should not be retained
     */
    protected abstract void initDispFlags_Indiv(GUI_AppUIFlags appUIFlags);
    
    @Override
    protected void initMe() {//all ui objects set by here
        //init simulation construct here
        msgObj.dispInfoMessage(className,"initMe","Start building simulation now.");
        //build sim(s) here
        currSim = buildSim();        
        //instance-class-specific init
        initMe_Indiv();
    }//initMe    
    
    /**
     * Build sim used by instancing window/application
     * @return
     */
    protected abstract Base_MPMSim buildSim();
    
    /**
     * instancing class-specific functionality
     */
    protected abstract void initMe_Indiv();
    
    /**
     * This function implements the instantiation of a child window owned by this window, if such exists.
     * The implementation should be similar to how the main windows are implemented in GUI_AppManager::initAllDispWindows.
     * If no child window exists, this implementation of this function can be empty
     * If a child window is instantiated, it MUST have its init called (childWin.initThisWin(false))
     * 
     * @param GUI_AppWinVals the window control values for the child window.
     */
    @Override
    protected final void buildAndSetChildWindow_Indiv(GUI_AppWinVals _appVals) {} 
    
    @Override
    protected int[] getFlagIDXsToInitToTrue() {
        return new int[] {showLocColors, showParticles};
    }

    //call this to initialize or reinitialize simulation (on reset)
    protected void reinitSim(SimResetProcess process) {    
        currSim.resetSim(process);
    }

    /**
     * UI code-level Debug mode functionality. Called only from flags structure
     * @param val
     */
    @Override
    protected final void handleDispFlagsDebugMode_Indiv(boolean val) {}
    
    /**
     * Application-specific Debug mode functionality (application-specific). Called only from privflags structure
     * @param val
     */
    @Override
    protected final void handlePrivFlagsDebugMode_Indiv(boolean val) {    }
    
    /**
     * Handle application-specific flag setting
     */
    @Override
    public void handlePrivFlags_Indiv(int idx, boolean val, boolean oldVal){
        switch(idx){
            case resetSimIDX            : {
                if(val) {
                    reinitSim(SimResetProcess.RemakeKernel);
                    addPrivSwitchToClear(resetSimIDX);
                }
                break;}                        
            case rebuildSimIDX            : {//Rebuild sim environment
                if(val) {
                    reinitSim(SimResetProcess.RebuildSim);
                    addPrivSwitchToClear(rebuildSimIDX);
                }
                break;}    
            case showLocColors            : {
                currSim.simFlags.setShowLocColors(val);
                break;}
            case showCollider            : {//show collider
                currSim.simFlags.setShowCollider(val);
                break;}
            case showParticles            : {
                currSim.simFlags.setShowParticles(val);                
                break;}
            case showParticleVelArrows    : {
                currSim.simFlags.setShowPartVels(val);
                break;} 
            case showGrid                : {
                currSim.simFlags.setShowGrid(val);
                break;}                 
            case showGridVelArrows         : {
                currSim.simFlags.setShowGridVel(val);
                break;}         
            case showGridAccelArrows    : {
                currSim.simFlags.setShowGridAccel(val);
                break;}     
            case showGridMass              : {
                currSim.simFlags.setShowGridMass(val);
                break;}         
            case showActiveNodes     : {
                currSim.simFlags.setShowActiveNodes(val);
                break;}     
            case showExecTime            : { 
                break;}
            default: {            setPrivFlags_Indiv(idx, val);}
        }        
    }//setPrivFlags    
    
    /**
     * set values for instancing class-specific boolean flags
     * 
     * @param idx
     * @param val
     */
    protected abstract void setPrivFlags_Indiv(int idx, boolean val);
        
    /**
     * Build all UI objects to be shown in left side bar menu for this window. This is the first child class function called by initThisWin
     * @param tmpUIObjMap : map of GUIObj_Params, keyed by unique string, with values describing the UI object
     *             - The object IDX                   
     *          - A double array of min/max/mod values                                                   
     *          - The starting value                                                                      
     *          - The label for object                                                                       
     *          - The object type (GUIObj_Type enum)
     *          - A boolean array of behavior configuration values : (unspecified values default to false)
     *               idx 0: value is sent to owning window,  
     *               idx 1: value is sent on any modifications (while being modified, not just on release), 
     *               idx 2: changes to value must be explicitly sent to consumer (are not automatically sent),
     *          - A boolean array of renderer format values :(unspecified values default to false) - Behavior Boolean array must also be provided!
     *                 - Should be multiline
     *                 - One object per row in UI space (i.e. default for multi-line and btn objects is false, single line non-buttons is true)
     *                 - Force this object to be on a new row/line (For side-by-side layouts)
     *                 - Text should be centered (default is false)
     *                 - Object should be rendered with outline (default for btns is true, for non-buttons is false)
     *                 - Should have ornament
     *                 - Ornament color should match label color  
     */
    @Override
    protected final void setupGUIObjsAras(LinkedHashMap<String, GUIObj_Params> tmpUIObjMap){        //keyed by object idx (uiXXXIDX), entries are lists of values to use for list select ui objects            
//        //keyed by object idx (uiXXXIDX), entries are lists of values to use for list select ui objects
        tmpUIObjMap.put("gIDX_TimeStep ", uiMgr.uiObjInitAra_Float(gIDX_TimeStep, new double[]{.00005f, .0010f, .00005f}, 1.0*initDeltaT, "Sim Time Step"));//delta T for simulation init  MPM_ABS_Sim.deltaT = 1e-3f;
        tmpUIObjMap.put("gIDX_SimStepsPerFrame", uiMgr.uiObjInitAra_Int(gIDX_SimStepsPerFrame, new double[]{1, 20, 1}, 1.0*initSimStepsPerFrame, "Sim Steps per Drawn Frame"));//gIDX_simStepsPerFrame  init 5
        tmpUIObjMap.put("gIDX_NumParticles", uiMgr.uiObjInitAra_Int(gIDX_NumParticles, getMinMaxModParts(), getInitNumParts(), "# of Particles"));//number of particles
        tmpUIObjMap.put("gIDX_NumSnowballs", uiMgr.uiObjInitAra_Int(gIDX_NumSnowballs, new double[]{2, 20, 1}, 1.0*initNumBalls, "# of Snowballs"));//number of snowballs
        tmpUIObjMap.put("gIDX_InitVel", uiMgr.uiObjInitAra_Float(gIDX_InitVel, new double[]{1, 40, .1}, 1.0*initVel, "Initial Speed"));//initial speed of collisions
        tmpUIObjMap.put("gIDX_PartMass", uiMgr.uiObjInitAra_Float(gIDX_PartMass, new double[]{.0005, 5.00, .0005}, 1.0*initParticleMass, "Particle Mass"));//particle mass
        tmpUIObjMap.put("gIDX_GridCellSize", uiMgr.uiObjInitAra_Float(gIDX_GridCellSize, new double[]{.001, .5, .001}, 1.0*initCellSize, "Grid Cell Size"));//grid cell size
        tmpUIObjMap.put("gIDX_GridCount", uiMgr.uiObjInitAra_Int(gIDX_GridCount, new double[]{50.0f, 300.0f, 10.0f}, 1.0*initNumGridCellsPerDim, "Grid Cell Count Per Side")); //# of grid cells per side
        tmpUIObjMap.put("gIDX_InitYoungMod", uiMgr.uiObjInitAra_Float(gIDX_InitYoungMod, new double[]{1000.0f, 200000.0f, 100.0f}, 1.0*init_initYoungMod, "Initial Young's Modulus"));//gIDX_InitYoungMod init 4.8e4f, 
        tmpUIObjMap.put("gIDX_PoissonRatio", uiMgr.uiObjInitAra_Float(gIDX_PoissonRatio, new double[]{.01f, 0.6f, .01f}, 1.0*init_poissonRatio, "Poisson Ratio"));//gIDX_PoissonRatio init 0.2f, 
        tmpUIObjMap.put("gIDX_HardeningCoeff ", uiMgr.uiObjInitAra_Float(gIDX_HardeningCoeff, new double[]{1.0f, 20.0f, 1.0f}, 1.0*init_hardeningCoeff, "Hardening Coefficient"));//gIDX_HardeningCoeff init 15.0f, 
        tmpUIObjMap.put("gIDX_CriticalCompression", uiMgr.uiObjInitAra_Float(gIDX_CriticalCompression, new double[]{0.001f, 0.1f, 0.001f}, 1.0*init_criticalCompression, "Critical Compression"));//gIDX_CriticalCompression  init .019f, 
        tmpUIObjMap.put("gIDX_CriticalStretch ", uiMgr.uiObjInitAra_Float(gIDX_CriticalStretch, new double[]{0.0005f, 0.01f, 0.0005f}, 1.0*init_criticalStretch, "Critical Stretch"));//gIDX_CriticalStretch init .0075f, 
        tmpUIObjMap.put("gIDX_AlphaPicFlip", uiMgr.uiObjInitAra_Float(gIDX_AlphaPicFlip, new double[]{0.0f, 1.0f, 0.01f}, 1.0*init_alphaPicFlip, "PIC/FLIP Mix Slider (0==PIC, 1== FLIP)"));//gIDX_AlphaPicFlip init 0.95f, 
        tmpUIObjMap.put("gIDX_wallFricCoeff", uiMgr.uiObjInitAra_Float(gIDX_wallFricCoeff, new double[]{0.01f, 1.0f, 0.01f}, 1.0*initWallFric, "Wall Friction Coefficient"));//gIDX_wallfricCoeffinit 1.0f  
        tmpUIObjMap.put("gIDX_CollFricCoeff", uiMgr.uiObjInitAra_Float(gIDX_CollFricCoeff, new double[]{0.01f, 1.0f, 0.01f}, 1.0*initColFric, "Collider Friction Coefficient"));//gIDX_CollfricCoeffinit 1.0f  
        tmpUIObjMap.put("gIDX_DrawnValScale", uiMgr.uiObjInitAra_Float(gIDX_DrawnValScale, new double[]{0.01f, 1.0f, 0.01f}, 1.0*initDrawnVecScale, "Scale Drawn Vectors"));//gIDX_CollfricCoeffinit 1.0f  
        tmpUIObjMap.put("gIDX_DrawPointIncr", uiMgr.uiObjInitAra_Int(gIDX_DrawPointIncr, new double[]{1, 50, 1}, 1.0*initDrawPtIncr, "Draw Every x'th Point"));//every x'th point to draw
        setupGUIObjsAras_Indiv(tmpUIObjMap);
    }//setupGUIObjsAras
    
    /**
     * Build UI button objects to be shown in left side bar menu for this window.  This is the first child class function called by initThisWin
     * @param firstIdx : the first index to use in the map/as the objIdx
     * @param tmpUIBoolSwitchObjMap : map of GUIObj_Params to be built containing all flag-backed boolean switch definitions, keyed by sequential value == objId
     *                 the first element is the object index
     *                 the second element is true label
     *                 the third element is false label
     *                 the final element is integer flag idx 
     */
    @Override
    protected final void setupGUIBoolSwitchAras(int firstIdx, LinkedHashMap<String, GUIObj_Params> tmpUIBoolSwitchObjMap) {        
        //add an entry for each button, in the order they are wished to be displayed
        int idx=firstIdx;
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.buildDebugButton(idx++,"Debugging", "Enable Debug"));        
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Resetting Sim Env", "Reset Sim Environment", resetSimIDX));      
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Remaking Simulation", "Remake Simulation", rebuildSimIDX));
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Init Loc Clr", "Show Init Loc Clr", showLocColors));          
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Collider", "Show Collider", showCollider));          
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Particles", "Show Particles", showParticles));  
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Particle Vel", "Show Particle Vel", showParticleVelArrows));  
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Grid", "Show Grid", showGrid));           
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Grid Vel", "Show Grid Vel", showGridVelArrows));     
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Grid Accel", "Show Grid Accel", showGridAccelArrows));    
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Grid Mass", "Show Grid Mass", showGridMass));         
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Active Nodes", "Show Active Nodes", showActiveNodes));     
        tmpUIBoolSwitchObjMap.put("Button_"+idx, uiMgr.uiObjInitAra_Switch(idx, "button_"+idx++, "Showing Execution Time", "Show Execution Time",showExecTime));
        
        // populate instancing application objects
        setupGUIBoolSwitchAras_Indiv(idx, tmpUIBoolSwitchObjMap);
    }//setupGUIBoolSwitchAras
    
    /**
     * Build all UI objects to be shown in left side bar menu for this window. This is the first child class function called by initThisWin
     * @param tmpUIObjMap : map of GUIObj_Params, keyed by unique string, with values describing the UI object
     *             - The object IDX                   
     *          - A double array of min/max/mod values                                                   
     *          - The starting value                                                                      
     *          - The label for object                                                                       
     *          - The object type (GUIObj_Type enum)
     *          - A boolean array of behavior configuration values : (unspecified values default to false)
     *               idx 0: value is sent to owning window,  
     *               idx 1: value is sent on any modifications (while being modified, not just on release), 
     *               idx 2: changes to value must be explicitly sent to consumer (are not automatically sent),
     *          - A boolean array of renderer format values :(unspecified values default to false) - Behavior Boolean array must also be provided!
     *                 - Should be multiline
     *                 - One object per row in UI space (i.e. default for multi-line and btn objects is false, single line non-buttons is true)
     *                 - Force this object to be on a new row/line (For side-by-side layouts)
     *                 - Text should be centered (default is false)
     *                 - Object should be rendered with outline (default for btns is true, for non-buttons is false)
     *                 - Should have ornament
     *                 - Ornament color should match label color r 
     */
    protected abstract void setupGUIObjsAras_Indiv(LinkedHashMap<String, GUIObj_Params> tmpUIObjMap);

    /**
     * Build all UI buttons to be shown in left side bar menu for this window. This is for instancing windows to add to button region
     * @param firstIdx : the first index to use in the map/as the objIdx
     * @param tmpUIBoolSwitchObjMap : map of GUIObj_Params to be built containing all flag-backed boolean switch definitions, keyed by sequential value == objId
     *                 the first element is the object index
     *                 the second element is true label
     *                 the third element is false label
     *                 the final element is integer flag idx 
     */
    protected abstract void setupGUIBoolSwitchAras_Indiv(int firstIdx, LinkedHashMap<String, GUIObj_Params> tmpUIBoolSwitchObjMap);
    
    /**
     * Return an array holding [min, max, mod] for particle count. This is simulation dependent due
     * to the vast disparity in capabilities between the cpu and cuda implementations.
     * @return
     */
    protected abstract double[] getMinMaxModParts();  
    
    /**
     * Return the initial number of particles to build. This is simulation dependent due 
     * to the vast disparity in capabilities between the cpu and cuda implementations.
     * @return
     */
    protected abstract double getInitNumParts();
    
    /**
     * This function provides an instance of the override class for base_UpdateFromUIData, which would
     * be used to communicate changes in UI settings directly to the value consumers.
     */
    @Override
    protected UIDataUpdater buildUIDataUpdateObject() {
        return new MPM_SimUpdateFromUIData(this);
    }
    /**
     * This function is called on ui value update, to pass new ui values on to window-owned consumers
     */
    @Override
    protected final void updateCalcObjUIVals() {
        var uiUpdateData = ((MPM_SimUpdateFromUIData) getUIDataUpdater());
        currSim.updateSimVals_FromUI(uiUpdateData);
    }//updateCalcObjUIVals

    /**
     * Called if int-handling guiObjs_Numeric[UIidx] (int or list) has new data which updated UI adapter. 
     * Intended to support custom per-object handling by owning window.
     * Only called if data changed!
     * @param UIidx Index of gui obj with new data
     * @param ival integer value of new data
     * @param oldVal integer value of old data in UIUpdater
     */
    @Override
    protected final void setUI_IntValsCustom(int UIidx, int ival, int oldVal) {
        switch(UIidx){        
            case gIDX_NumParticles            :{break;}
            case gIDX_NumSnowballs            :{break;}
            case gIDX_GridCount                :{break;}
            case gIDX_SimStepsPerFrame         :{break;}
            case gIDX_DrawPointIncr            :{break;}
            default : {
                boolean found = setUI_IntValsCustom_Indiv( UIidx, ival, oldVal);
                if (!found) {
                    msgObj.dispWarningMessage(className, "setUI_IntValsCustom", "No int-defined gui object mapped to idx :"+UIidx);
                }
                break;}
        }                
    }//setUI_IntValsCustom
    
    protected abstract boolean setUI_IntValsCustom_Indiv(int UIidx, int ival, int oldVal);

    /**
     * Called if float-handling guiObjs_Numeric[UIidx] has new data which updated UI adapter.  
     * Intended to support custom per-object handling by owning window.
     * Only called if data changed!
     * @param UIidx Index of gui obj with new data
     * @param val float value of new data
     * @param oldVal float value of old data in UIUpdater
     */
    @Override
    protected final void setUI_FloatValsCustom(int UIidx, float val, float oldVal) {
        switch(UIidx){        
            case gIDX_TimeStep                 :{break;}
            case gIDX_InitVel                :{break;}
            case gIDX_PartMass                 :{break;}
            case gIDX_GridCellSize            :{break;}
            case gIDX_InitYoungMod             :{break;}
            case gIDX_PoissonRatio             :{break;}
            case gIDX_HardeningCoeff         :{break;}
            case gIDX_CriticalCompression     :{break;}
            case gIDX_CriticalStretch         :{break;}
            case gIDX_AlphaPicFlip             :{break;}
            case gIDX_wallFricCoeff         :{break;}
            case gIDX_CollFricCoeff         :{break;}
            case gIDX_DrawnValScale            :{break;}
            default : {
                boolean found = setUI_FloatValsCustom_Indiv( UIidx, val, oldVal);
                if (!found) {
                    msgObj.dispWarningMessage(className, "setUI_FloatValsCustom", "No float-defined gui object mapped to idx :"+UIidx);
                }
                break;}
        }    
    }//setUI_FloatValsCustom
    
    protected abstract boolean setUI_FloatValsCustom_Indiv(int UIidx, float ival, float oldVal);
        
    @Override
    public void initDrwnTraj_Indiv(){}

    @Override
    protected void setCamera_Indiv(float[] camVals) {
        // No custom camera handling
        setCameraBase(camVals);
    }//setCameraIndiv
    
    @Override
    //modAmtMillis is time passed per frame in milliseconds
    protected final boolean simMe(float modAmtMillis) {//run simulation
        boolean done = false;
        if (uiMgr.getPrivFlag(showExecTime)){
            done = currSim.simMeDebug(modAmtMillis);
        } else {
            done = currSim.simMe(modAmtMillis);    
        }
        return done;    
    }//simMe        
    
    @Override
    //animTimeMod is in seconds.
    protected void drawMe(float animTimeMod, boolean isGlblAppDebug) {
        currSim.drawMe(animTimeMod, isGlblAppDebug);
    }//drawMe    
        
    //draw custom 2d constructs below interactive component of menu
    @Override
    public void drawCustMenuObjs(float animTimeMod, boolean isGlblAppDebug){    }//drawCustMenuObjs

    //things to do when swapping this window out for another window - release objects that take up a lot of memory, for example.
    @Override
    protected void closeMe() {}    
    @Override
    //stopping simulation
    protected void stopMe() {    msgObj.dispInfoMessage(className,"stopMe","Simulation Finished");    }    

    @Override
    protected boolean hndlMouseMove_Indiv(int mouseX, int mouseY){                    return false;    }
    @Override
    protected boolean hndlMouseClick_Indiv(int mouseX, int mouseY, int mseBtn) {    return false;}//hndlMouseClickIndiv
    @Override
    protected boolean hndlMouseDrag_Indiv(int mouseX, int mouseY, int pmouseX, int pmouseY, int mseBtn) {return false;}    
    @Override
    protected boolean handleMouseWheel_Indiv(int ticks, float mult) {        return false;    }
    
    @Override
    protected void hndlMouseRel_Indiv() {    }
    @Override
    protected void endShiftKey_Indiv() {}
    @Override
    protected void endAltKey_Indiv() {}
    @Override
    protected void endCntlKey_Indiv() {}
    @Override
    protected void snapMouseLocs(int oldMouseX, int oldMouseY, int[] newMouseLoc) {}    
    @Override
    protected void addSScrToWin_Indiv(int newWinKey){}
    @Override
    protected void addTrajToScr_Indiv(int subScrKey, String newTrajKey){}
    @Override
    protected void delSScrToWin_Indiv(int idx) {}    
    @Override
    protected void delTrajToScr_Indiv(int subScrKey, String newTrajKey) {}
    //resize drawn all trajectories
    @Override
    protected void resizeMe(float scale) { }
    @Override
    protected void showMe() {}
    /**
     * type is row of buttons (1st idx in curCustBtn array) 2nd idx is btn
     * @param funcRow idx for button row
     * @param btn idx for button within row (column)
     * @param label label for this button (for display purposes)
     */
    @Override
    protected final void launchMenuBtnHndlr(int funcRow, int btn, String label){
        switch (funcRow) {
            case 0: {// row 1 of menu side bar buttons
                switch (btn) {
                    case 0: {
                        //Reset all UI vals to be initial values
                        uiMgr.resetUIVals(true);
                        resetButtonState();
                        break;
                    }
                    case 1: {
                        int[] winDims = AppMgr.getWindowLoc();
                        msgObj.dispInfoMessage(className, "launchMenuBtnHndlr", "Window dims : X:"+winDims[0]+" | Y:"+winDims[1]+" | width:"+winDims[2]+" | height :"+winDims[3]);
                        resetButtonState();
                        break;
                    }
                    case 2: {
                        resetButtonState();
                        break;
                    }
                    default: {
                        msgObj.dispMessage(className,"launchMenuBtnHndlr", "Unknown Functions 1 btn : " + btn, MsgCodes.warning2);
                        break;
                    }
                }
                break;
            } // row 1 of menu side bar buttons
    
            case 1: {// row 2 of menu side bar buttons
                switch (btn) {
                    case 0: {    
                        resetButtonState();
                        break;
                    }
                    case 1: {    
                        resetButtonState();
                        break;
                    }
                    case 2: {
                        resetButtonState();
                        break;
                    }
                    case 3: {
                        resetButtonState();
                        break;
                    }
                    default: {
                        msgObj.dispMessage(className,"launchMenuBtnHndlr", "Unknown Functions 2 btn : " + btn, MsgCodes.warning2);
                        resetButtonState();
                        break;
                    }
                }
                break;
            } // row 2 of menu side bar buttons
            case 2: {// row 3 of menu side bar buttons
                switch (btn) {
                    case 0: {
                        resetButtonState();
                        break;
                    }
                    case 1: {
                        resetButtonState();
                        break;
                    }
                    case 2: {
                        resetButtonState();
                        break;
                    }
                    case 3: {
                        resetButtonState();
                        break;
                    }
                    default: {
                        msgObj.dispMessage(className,"launchMenuBtnHndlr", "Unknown Functions 3 btn : " + btn,
                                MsgCodes.warning2);
                        resetButtonState();
                        break;
                    }
                }
                break;
            } // row 3 of menu side bar buttons
            case 3: {// row 3 of menu side bar buttons
                switch (btn) {
                    case 0: {
                        resetButtonState();
                        break;
                    }
                    case 1: {
                        resetButtonState();
                        break;
                    }
                    case 2: {
                        resetButtonState();
                        break;
                    }
                    case 3: {
                        resetButtonState();
                        break;
                    }
                    default: {
                        msgObj.dispMessage(className, "launchMenuBtnHndlr", "Unknown Functions 4 btn : " + btn, MsgCodes.warning2);
                        resetButtonState();
                        break;
                    }
                }
                break;
            } // row 3 of menu side bar buttons
            default : {
                msgObj.dispWarningMessage(className,"launchMenuBtnHndlr","Clicked Unknown Btn row : " + funcRow +" | Btn : " + btn);
                break;
            }            
        }
    }
    @Override
    public void handleSideMenuMseOvrDispSel(int btn, boolean val) {}
    @Override
    protected final void handleSideMenuDebugSelEnable(int btn) {
        switch (btn) {
            case 0: {                break;            }
            case 1: {                break;            }
            case 2: {                break;            }
            case 3: {                break;            }
            case 4: {                break;            }
            case 5: {                break;            }
            default: {
                msgObj.dispMessage(className, "handleSideMenuDebugSelEnable", "Unknown Debug btn : " + btn,MsgCodes.warning2);
                break;
            }
        }
    }
    
    @Override
    protected final void handleSideMenuDebugSelDisable(int btn) {
        switch (btn) {
            case 0: {                break;            }
            case 1: {                break;            }
            case 2: {                break;            }
            case 3: {                break;            }
            case 4: {                break;            }
            case 5: {                break;            }
        default: {
            msgObj.dispMessage(className, "handleSideMenuDebugSelDisable", "Unknown Debug btn : " + btn,MsgCodes.warning2);
            break;
            }
        }
    }
    
    @Override
    protected String[] getSaveFileDirNamesPriv() {        return null;    }
    @Override
    protected myPoint getMsePtAs3DPt(myPoint mseLoc) {        return new myPoint(mseLoc.x,mseLoc.y,0);    }
    @Override
    protected void setVisScreenDimsPriv() {    }
    @Override
    protected void setCustMenuBtnLabels() {                }
    @Override
    public void hndlFileLoad(File file, String[] vals, int[] stIdx) {    }
    @Override
    public ArrayList<String> hndlFileSave(File file) {        return null;    }
    @Override
    protected void drawOnScreenStuffPriv(float modAmtMillis, boolean isGlblAppDebug) {    }
    @Override
    protected void drawRightSideInfoBarPriv(float modAmtMillis, boolean isGlblAppDebug) {    }
    @Override
    public void processTraj_Indiv(DrawnSimpleTraj drawnTraj) {    }

}//DESSimWindow

