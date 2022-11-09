package MPM_CudaSim.ui;

import java.io.File;
import java.util.ArrayList;
import java.util.TreeMap;

import MPM_CudaSim.sim.MPM_CudaBalls;
import MPM_CudaSim.sim.SimResetProcess;
import MPM_CudaSim.sim.base.base_MPMCudaSim;
import MPM_CudaSim.utils.MPM_SimUpdateFromUIData;
import base_JavaProjTools_IRender.base_Render_Interface.IRenderInterface;
import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.doubles.myPoint;
import base_Math_Objects.vectorObjs.doubles.myVector;
import base_UI_Objects.GUI_AppManager;
import base_UI_Objects.windowUI.base.Base_DispWindow;
import base_UI_Objects.windowUI.drawnObjs.DrawnSimpleTraj;
import base_UI_Objects.windowUI.uiData.UIDataUpdater;
import base_UI_Objects.windowUI.uiObjs.base.GUIObj_Type;
import base_Utils_Objects.io.messaging.MsgCodes;
import base_Utils_Objects.tools.myTools;

public class MPM_SimWindow extends Base_DispWindow {

	//simulator world within which simulation executes
	//TODO support multiple sim worlds possibly, with different configurations
	public base_MPMCudaSim currSim;
		
	//////////////////////////////////
	//initial values of simulation variables
	//ints
	protected final int initNumGridCellsPerDim = 200;
	protected final int initSimStepsPerFrame = 5;
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
	protected final float initParticleDensity = 50.0f;
	
	protected final float initParticleMass = initCellSize *initCellSize *initCellSize * initParticleDensity; 
	protected final float initWallFric = 1.0f;
	protected final float initColFric = 1.0f;
	protected final float initDrawnVecScale = .01f;
	
	//initial values of material quantities
	protected final float	
		init_initYoungMod 			 = 1.4e5f,
		init_poissonRatio 			 = 0.2f,
		init_hardeningCoeff 		 = 15.0f, 
		init_criticalCompression 	 = 0.010f, 
		init_criticalStretch 		 = 0.0025f, 
		init_alphaPicFlip 			 = 0.95f;
	
	
	//////////////////////////////////////
	//ui vals
	public final static int
		gIDX_TimeStep		 		= 0,		
		gIDX_SimStepsPerFrame		= 1,
		gIDX_NumParticles			= 2,
		gIDX_NumSnowballs 			= 3,
		gIDX_InitVel				= 4,
		gIDX_PartMass				= 5,
		gIDX_GridCellSize			= 6,
		gIDX_GridCount				= 7,
		gIDX_InitYoungMod 			= 8,
		gIDX_PoissonRatio 			= 9,
		gIDX_HardeningCoeff 		= 10,
		gIDX_CriticalCompression 	= 11,
		gIDX_CriticalStretch 		= 12,
		gIDX_AlphaPicFlip 			= 13,
		gIDX_wallFricCoeff 			= 14,
		gIDX_CollFricCoeff			= 15,
		gIDX_DrawnValScale			= 16,
		gIDX_DrawPointIncr			= 17;
	
	public final static int numGuiObjs = 18;
	
	//////////////////////////////////////
	//custom debug/function ui button names -empty will do nothing

	//////////////////////////////////////
	//private child-class flags - window specific
	public static final int 
		debugAnimIDX 			= 0,			//debug
		resetSimIDX				= 1,			//whether or not to reset sim	
		rebuildSimIDX			= 2,			//reset of sim is required for best results
		showLocColors			= 3,			//display particles by color of their initial location
		showCollider			= 4,			//show collider cylinder
		showParticles			= 5,
		showParticleVelArrows 	= 6,			//plot velocity arrows for each particle	
		showGrid				= 7,  			//plot the computational grid
		showGridVelArrows 		= 8,			//plot velocity arrows for each gridNode
		showGridAccelArrows 	= 9,			//plot acceleration arrows for each gridNode
		showGridMass  			= 10,			//plot variable sized spheres proportional to gridnode mass
		showActiveNodes  		= 11,				//show the grid nodes influenced by each particle
		showExecTime			= 12;			//show time of execution for each step in simulation

	public static final int numPrivFlags = 12;
		
	public MPM_SimWindow(IRenderInterface _p, GUI_AppManager _AppMgr, int _winIdx, int _flagIdx) {
		super(_p, _AppMgr, _winIdx, _flagIdx);
		super.initThisWin(false);
	}//DancingBallWin
	
	@Override
	//initialize all private-flag based UI buttons here - called by base class
	public int initAllPrivBtns(ArrayList<Object[]> tmpBtnNamesArray){	
		tmpBtnNamesArray.add(new Object[]{"Visualization Debug",    "Enable Debug",       debugAnimIDX});          
		tmpBtnNamesArray.add(new Object[]{"Resetting Sim Env",   	"Reset Sim Environment",  resetSimIDX});      
		tmpBtnNamesArray.add(new Object[]{"Remaking Simulation",    "Remake Simulation",   rebuildSimIDX});
		tmpBtnNamesArray.add(new Object[]{"Showing Init Loc Clr",   "Show Init Loc Clr",  showLocColors});          
		tmpBtnNamesArray.add(new Object[]{"Showing Collider",       "Show Collider",      showCollider});          
		tmpBtnNamesArray.add(new Object[]{"Showing Particles",      "Show Particles",     showParticles});  
		tmpBtnNamesArray.add(new Object[]{"Showing Particle Vel",   "Show Particle Vel",  showParticleVelArrows});  
		tmpBtnNamesArray.add(new Object[]{"Showing Grid",           "Show Grid",          showGrid});           
		tmpBtnNamesArray.add(new Object[]{"Showing Grid Vel",       "Show Grid Vel",      showGridVelArrows});     
		tmpBtnNamesArray.add(new Object[]{"Showing Grid Accel",     "Show Grid Accel",    showGridAccelArrows});    
		tmpBtnNamesArray.add(new Object[]{"Showing Grid Mass",      "Show Grid Mass",     showGridMass});         
		tmpBtnNamesArray.add(new Object[]{"Showing Active Nodes",   "Show Active Nodes",  showActiveNodes});     
		tmpBtnNamesArray.add(new Object[]{"Showing Execution Time", "Show Execution Time",showExecTime});        
		return numPrivFlags;
	}//initAllPrivBtns
	//set labels of boolean buttons 

	
	@Override
	protected void initMe() {//all ui objects set by here
		//this window is runnable
		setFlags(isRunnable, true);
		//this window uses a customizable camera
		setFlags(useCustCam, true);
		// capable of using right side menu
		setFlags(drawRightSideMenu, true);
		
		//init simulation construct here
		msgObj.dispInfoMessage(className+"(MPM_SimWindow)","initMe","Start building simulation now.");
		//build sim(s) here
		currSim = new MPM_CudaBalls(pa, this, (MPM_SimUpdateFromUIData) uiUpdateData);		
		//initialize simulation here to simple world sim
		custMenuOffset = uiClkCoords[3];	//495	
		setPrivFlags(showParticles, true);
		setPrivFlags(showLocColors, true);

	}//initMe	

	@Override
	protected int[] getFlagIDXsToInitToTrue() {
		//Does not pass value to button handler, just sets value in flag array to true
		return null;// new int[] {showLocColors, showParticles};
	}

	//call this to initialize or reinitialize simulation (on reset)
	protected void reinitSim(SimResetProcess process) {	
		currSim.resetSim(process);
	}
		
	@Override
	//set flag values and execute special functionality for this sequencer
	//skipKnown will allow settings to be reset if passed redundantly
	public void setPrivFlags(int idx, boolean val){	
		boolean curVal = getPrivFlags(idx);
		if(val == curVal){return;}
		int flIDX = idx/32, mask = 1<<(idx%32);
		privFlags[flIDX] = (val ?  privFlags[flIDX] | mask : privFlags[flIDX] & ~mask);
		switch(idx){
			case debugAnimIDX 			: {
				break;}
			case resetSimIDX			: {
				if(val) {
					reinitSim(SimResetProcess.RemakeKernel);
					addPrivBtnToClear(resetSimIDX);
				}
				break;}						
			case rebuildSimIDX			: {//Rebuild sim environment
				if(val) {
					reinitSim(SimResetProcess.RebuildSim);
					addPrivBtnToClear(rebuildSimIDX);
				}
				break;}	
			case showLocColors			: {
				currSim.setSimFlags(base_MPMCudaSim.showLocColors, val);
				break;}
			case showCollider			: {//show collider
				currSim.setSimFlags(base_MPMCudaSim.showCollider, val);
				break;}
			case showParticles			: {
				currSim.setSimFlags(base_MPMCudaSim.showParticles, val);				
				break;}
			case showParticleVelArrows	: {
				currSim.setSimFlags(base_MPMCudaSim.showParticleVelArrows, val);
				break;} 
			case showGrid				: {
				currSim.setSimFlags(base_MPMCudaSim.showGrid, val);
				break;} 				
			case showGridVelArrows 		: {
				currSim.setSimFlags(base_MPMCudaSim.showGridVelArrows, val);
				break;} 		
			case showGridAccelArrows	: {
				currSim.setSimFlags(base_MPMCudaSim.showGridAccelArrows, val);
				break;} 	
			case showGridMass  			: {
				currSim.setSimFlags(base_MPMCudaSim.showGridMass, val);
				break;} 		
			case showActiveNodes 	: {
				currSim.setSimFlags(base_MPMCudaSim.showActiveNodes, val);
				break;} 	
			case showExecTime			: { 
				break;}
			default:					{}
		}		
	}//setPrivFlags	
		
	/**
	 * Build all UI objects to be shown in left side bar menu for this window.  This is the first child class function called by initThisWin
	 * @param tmpUIObjArray : map of object data, keyed by UI object idx, with array values being :                    
	 *           the first element double array of min/max/mod values                                                   
	 *           the 2nd element is starting value                                                                      
	 *           the 3rd elem is label for object                                                                       
	 *           the 4th element is object type (GUIObj_Type enum)
	 *           the 5th element is boolean array of : (unspecified values default to false)
	 *           	{value is sent to owning window, 
	 *           	value is sent on any modifications (while being modified, not just on release), 
	 *           	changes to value must be explicitly sent to consumer (are not automatically sent)}    
	 * @param tmpListObjVals : map of list object possible selection values
	 */
	@Override
	protected final void setupGUIObjsAras(TreeMap<Integer, Object[]> tmpUIObjArray, TreeMap<Integer, String[]> tmpListObjVals){		
		tmpUIObjArray.put(gIDX_TimeStep , new Object[] {new double[]{.00005f, .0010f, .00005f}, 1.0*initDeltaT,  "Sim Time Step", GUIObj_Type.FloatVal, new boolean[]{true}});//delta T for simulation init  MPM_ABS_Sim.deltaT = 1e-3f;
		tmpUIObjArray.put(gIDX_SimStepsPerFrame, new Object[] {new double[]{1, 20, 1}, 1.0*initSimStepsPerFrame, "Sim Steps per Drawn Frame", GUIObj_Type.IntVal, new boolean[]{true}});//gIDX_simStepsPerFrame  init 5
		tmpUIObjArray.put(gIDX_NumParticles, new Object[] {new double[]{1000, 1000000, 1000}, 1.0*initNumParts, "# of Particles", GUIObj_Type.IntVal, new boolean[]{true}});//number of particles
		tmpUIObjArray.put(gIDX_NumSnowballs, new Object[] {new double[]{2, 10, 1}, 1.0*initNumBalls, "# of Snowballs", GUIObj_Type.IntVal, new boolean[]{true}});//number of snowballs
		tmpUIObjArray.put(gIDX_InitVel, new Object[] {new double[]{1, 40, .1}, 1.0*initVel, "Initial Speed", GUIObj_Type.FloatVal, new boolean[]{true}});//initial speed of collisions
		tmpUIObjArray.put(gIDX_PartMass, new Object[] {new double[]{.0005, 5.00, .0005}, 1.0*initParticleMass, "Particle Mass", GUIObj_Type.FloatVal, new boolean[]{true}});//particle mass
		tmpUIObjArray.put(gIDX_GridCellSize, new Object[] {new double[]{.001, .5, .001}, 1.0*initCellSize, "Grid Cell Size", GUIObj_Type.FloatVal, new boolean[]{true}});//grid cell size
		tmpUIObjArray.put(gIDX_GridCount, new Object[] {new double[]{50, 300, 1}, 1.0*initNumGridCellsPerDim,  "Grid Cell Count Per Side", GUIObj_Type.IntVal, new boolean[]{true}}); //# of grid cells per side
		tmpUIObjArray.put(gIDX_InitYoungMod, new Object[] {new double[]{1000.0f, 200000.0f, 100.0f}, 1.0*init_initYoungMod, "Young's Modulus", GUIObj_Type.FloatVal, new boolean[]{true}});//gIDX_InitYoungMod init 4.8e4f, 
		tmpUIObjArray.put(gIDX_PoissonRatio, new Object[] {new double[]{.01f, 0.6f, .01f}, 1.0*init_poissonRatio, "Poisson Ratio", GUIObj_Type.FloatVal, new boolean[]{true}});//gIDX_PoissonRatio init 0.2f,  
		tmpUIObjArray.put(gIDX_HardeningCoeff , new Object[] {new double[]{1.0f, 20.0f, 1.0f}, 1.0*init_hardeningCoeff, "Hardening Coefficient", GUIObj_Type.FloatVal, new boolean[]{true}});//gIDX_HardeningCoeff init 15.0f, 
		tmpUIObjArray.put(gIDX_CriticalCompression, new Object[] {new double[]{0.001f, 0.1f, 0.001f}, 1.0*init_criticalCompression, "Critical Compression", GUIObj_Type.FloatVal, new boolean[]{true}});//gIDX_CriticalCompression  init .019f, 
		tmpUIObjArray.put(gIDX_CriticalStretch , new Object[] {new double[]{0.0005f, 0.01f, 0.0005f}, 1.0*init_criticalStretch, "Critical Stretch", GUIObj_Type.FloatVal, new boolean[]{true}});//gIDX_CriticalStretch init .0075f, 
		tmpUIObjArray.put(gIDX_AlphaPicFlip, new Object[] {new double[]{0.0f, 1.0f, 0.01f}, 1.0*init_alphaPicFlip, "Particle PIC/FLIP Vel Ratio", GUIObj_Type.FloatVal, new boolean[]{true}});//gIDX_AlphaPicFlip init 0.95f, 
		tmpUIObjArray.put(gIDX_wallFricCoeff, new Object[] {new double[]{0.01f, 1.0f, 0.01f}, 1.0*initWallFric, "Wall Friction Coefficient", GUIObj_Type.FloatVal, new boolean[]{true}});//gIDX_wallfricCoeffinit 1.0f  
		tmpUIObjArray.put(gIDX_CollFricCoeff, new Object[] {new double[]{0.01f, 1.0f, 0.01f}, 1.0*initColFric, "Collider Friction Coefficient", GUIObj_Type.FloatVal, new boolean[]{true}});//gIDX_CollfricCoeffinit 1.0f  
		tmpUIObjArray.put(gIDX_DrawnValScale, new Object[] {new double[]{0.01f, 1.0f, 0.01f}, 1.0*initDrawnVecScale, "Scale Drawn Vectors", GUIObj_Type.FloatVal, new boolean[]{true}});//gIDX_CollfricCoeffinit 1.0f  
		tmpUIObjArray.put(gIDX_DrawPointIncr, new Object[] {new double[]{1, 50, 1}, 1.0*initDrawPtIncr, "Draw Every x'th Point", GUIObj_Type.IntVal, new boolean[]{true}});//every x'th point to draw
	}//setupGUIObjsAras
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
		currSim.updateSimVals_FromUI((MPM_SimUpdateFromUIData) uiUpdateData);
	}//updateCalcObjUIVals

	/**
	 * Called if int-handling guiObjs[UIidx] (int or list) has new data which updated UI adapter. 
	 * Intended to support custom per-object handling by owning window.
	 * Only called if data changed!
	 * @param UIidx Index of gui obj with new data
	 * @param ival integer value of new data
	 * @param oldVal integer value of old data in UIUpdater
	 */
	@Override
	protected final void setUI_IntValsCustom(int UIidx, int ival, int oldVal) {
		switch(UIidx){		
			case gIDX_NumParticles			:{break;}
			case gIDX_NumSnowballs			:{break;}
			case gIDX_GridCount				:{break;}
			case gIDX_SimStepsPerFrame 		:{break;}
			case gIDX_DrawPointIncr			:{break;}
			default : {
				msgObj.dispWarningMessage(className+"(MPM_SimWindow)", "setUI_IntValsCustom", "No int-defined gui object mapped to idx :"+UIidx);
				break;}
		}				
	}//setUI_IntValsCustom

	/**
	 * Called if float-handling guiObjs[UIidx] has new data which updated UI adapter.  
	 * Intended to support custom per-object handling by owning window.
	 * Only called if data changed!
	 * @param UIidx Index of gui obj with new data
	 * @param val float value of new data
	 * @param oldVal float value of old data in UIUpdater
	 */
	@Override
	protected final void setUI_FloatValsCustom(int UIidx, float val, float oldVal) {
		switch(UIidx){		
			case gIDX_TimeStep 				:{break;}
			case gIDX_InitVel				:{break;}
			case gIDX_PartMass 				:{break;}
			case gIDX_GridCellSize			:{break;}
			case gIDX_InitYoungMod 			:{break;}
			case gIDX_PoissonRatio 			:{break;}
			case gIDX_HardeningCoeff 		:{break;}
			case gIDX_CriticalCompression 	:{break;}
			case gIDX_CriticalStretch 		:{break;}
			case gIDX_AlphaPicFlip 			:{break;}
			case gIDX_wallFricCoeff 		:{break;}
			case gIDX_CollFricCoeff 		:{break;}
			case gIDX_DrawnValScale			:{break;}
			default : {
				msgObj.dispWarningMessage(className+"(MPM_SimWindow)", "setUI_FloatValsCustom", "No float-defined gui object mapped to idx :"+UIidx);
				break;}
		}				
	}//setUI_FloatValsCustom	
	
	
	@Override
	public void initDrwnTrajIndiv(){}
	
	//overrides function in base class mseClkDisp
	@Override
	public void drawTraj3D(float animTimeMod,myPoint trans){}//drawTraj3D	
	//set camera to either be global or from pov of one of the boids
	@Override
	protected void setCameraIndiv(float[] camVals){		
		//, float rx, float ry, float dz are now member variables of every window
		pa.setCameraWinVals(camVals);//(camVals[0],camVals[1],camVals[2],camVals[3],camVals[4],camVals[5],camVals[6],camVals[7],camVals[8]);      
		// puts origin of all drawn objects at screen center and moves forward/away by dz
		pa.translate(camVals[0],camVals[1],(float)dz); 
	    setCamOrient();	
	}//setCameraIndiv
	
	@Override
	//modAmtMillis is time passed per frame in milliseconds
	protected boolean simMe(float modAmtMillis) {//run simulation
		boolean done = false;
		if (getPrivFlags(showExecTime)){
			done = currSim.simMeDebug(modAmtMillis);
		} else {
			done = currSim.simMe(modAmtMillis);	
		}
		
		return done;	
	}//simMe		
	
	@Override
	//animTimeMod is in seconds.
	protected void drawMe(float animTimeMod) {
		currSim.drawMe(animTimeMod);
	}//drawMe	
		
	//draw custom 2d constructs below interactive component of menu
	@Override
	public void drawCustMenuObjs(){
		pa.pushMatState();	
		//all sub menu drawing within push mat call
		pa.translate(0,custMenuOffset+yOff);		
		//draw any custom menu stuff here
		pa.popMatState();
	}//drawCustMenuObjs

	//things to do when swapping this window out for another window - release objects that take up a lot of memory, for example.
	@Override
	protected void closeMe() {}	
	@Override
	//stopping simulation
	protected void stopMe() {	msgObj.dispInfoMessage(className+"(MPM_SimWindow)","stopMe","Simulation Finished");	}	

	@Override
	protected boolean hndlMouseMoveIndiv(int mouseX, int mouseY, myPoint mseClckInWorld){
		return false;
	}
	
	//alt key pressed handles trajectory
	
	//cntl key pressed handles unfocus of spherey
	@Override
	protected boolean hndlMouseClickIndiv(int mouseX, int mouseY, myPoint mseClckInWorld, int mseBtn) {	
		return false;}//hndlMouseClickIndiv
	@Override
	protected boolean hndlMouseDragIndiv(int mouseX, int mouseY, int pmouseX, int pmouseY, myPoint mouseClickIn3D, myVector mseDragInWorld, int mseBtn) {
		boolean res = false;
		return res;}	
	@Override
	protected void hndlMouseRelIndiv() {	}
	@Override
	protected void endShiftKeyI() {}
	@Override
	protected void endAltKeyI() {}
	@Override
	protected void endCntlKeyI() {}
	@Override
	protected void snapMouseLocs(int oldMouseX, int oldMouseY, int[] newMouseLoc) {}	
	@Override
	protected void addSScrToWinIndiv(int newWinKey){}
	@Override
	protected void addTrajToScrIndiv(int subScrKey, String newTrajKey){}
	@Override
	protected void delSScrToWinIndiv(int idx) {}	
	@Override
	protected void delTrajToScrIndiv(int subScrKey, String newTrajKey) {}
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
		msgObj.dispMessage(className+"(MPM_SimWindow)", "launchMenuBtnHndlr", "Begin requested action : Click '" + label +"' (Row:"+(funcRow+1)+"|Col:"+btn+") in " + name, MsgCodes.info4);
		switch (funcRow) {
			case 0: {// row 1 of menu side bar buttons
				switch (btn) {
					case 0: {
						//Reset all UI vals to be initial values
						resetUIVals(false);
						resetButtonState();
						break;
					}
					case 1: {
						int[] winDims = AppMgr.getWindowLoc();
						msgObj.dispInfoMessage(className+"(MPM_SimWindow)", "launchMenuBtnHndlr", "Window dims : X:"+winDims[0]+" | Y:"+winDims[1]+" | width:"+winDims[2]+" | height :"+winDims[3]);
						resetButtonState();
						break;
					}
					case 2: {
						resetButtonState();
						break;
					}
					default: {
						msgObj.dispMessage(className+"(MPM_SimWindow)","launchMenuBtnHndlr", "Unknown Functions 1 btn : " + btn, MsgCodes.warning2);
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
						msgObj.dispMessage(className+"(MPM_SimWindow)","launchMenuBtnHndlr", "Unknown Functions 2 btn : " + btn, MsgCodes.warning2);
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
						msgObj.dispMessage(className+"(MPM_SimWindow)","launchMenuBtnHndlr", "Unknown Functions 3 btn : " + btn,
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
						msgObj.dispMessage(className+"(MPM_SimWindow)", "launchMenuBtnHndlr", "Unknown Functions 4 btn : " + btn, MsgCodes.warning2);
						resetButtonState();
						break;
					}
				}
				break;
			} // row 3 of menu side bar buttons
		}
		msgObj.dispMessage(className+"(MPM_SimWindow)","launchMenuBtnHndlr", "End requested action (multithreaded actions may still be working) : Click Functions "+(funcRow+1)+" in " + name + " : btn : " + btn, MsgCodes.info4);
	}
	@Override
	public void handleSideMenuMseOvrDispSel(int btn, boolean val) {}
	@Override
	public final void handleSideMenuDebugSelEnable(int btn) {
		msgObj.dispMessage(className+"(MPM_SimWindow)", "handleSideMenuDebugSelEnable","Click Debug functionality on in " + name + " : btn : " + btn, MsgCodes.info4);
		switch (btn) {
			case 0: {				
				//TEMP Factorial test
				msgObj.dispInfoMessage(className+"(MPM_SimWindow)", "handleSideMenuDebugSelEnable","Factorials 1->1000");
				for(int i=1;i<1000;++i) {
					byte[] ans = MyMathUtils.bigFact(i);
					String str = myTools.digitAraToString(ans);
					msgObj.dispInfoMessage(className+"(MPM_SimWindow)", "handleSideMenuDebugSelEnable","\t"+i+" ("+ans.length +") : "+str);
					
				}
				break;			}
			case 1: {				break;			}
			case 2: {				break;			}
			case 3: {				break;			}
			case 4: {				break;			}
			case 5: {				break;			}
			default: {
				msgObj.dispMessage(className+"(MPM_SimWindow)", "handleSideMenuDebugSelEnable", "Unknown Debug btn : " + btn,MsgCodes.warning2);
				break;
			}
		}
		msgObj.dispMessage(className+"(MPM_SimWindow)", "handleSideMenuDebugSelEnable", "End Debug functionality on selection.",MsgCodes.info4);
	}
	
	@Override
	public final void handleSideMenuDebugSelDisable(int btn) {
		msgObj.dispMessage(className+"(MPM_SimWindow)", "handleSideMenuDebugSelDisable","Click Debug functionality off in " + name + " : btn : " + btn, MsgCodes.info4);
		switch (btn) {
			case 0: {				break;			}
			case 1: {				break;			}
			case 2: {				break;			}
			case 3: {				break;			}
			case 4: {				break;			}
			case 5: {				break;			}
		default: {
			msgObj.dispMessage(className+"(MPM_SimWindow)", "handleSideMenuDebugSelDisable", "Unknown Debug btn : " + btn,MsgCodes.warning2);
			break;
			}
		}
		msgObj.dispMessage(className+"(MPM_SimWindow)", "handleSideMenuDebugSelDisable", "End Debug functionality off selection.",MsgCodes.info4);
	}
	
	@Override
	protected String[] getSaveFileDirNamesPriv() {		return null;	}
	@Override
	protected myPoint getMsePtAs3DPt(myPoint mseLoc) {		return new myPoint(mseLoc.x,mseLoc.y,0);	}
	@Override
	protected void setVisScreenDimsPriv() {	}
	@Override
	protected void setCustMenuBtnLabels() {				}
	@Override
	public void hndlFileLoad(File file, String[] vals, int[] stIdx) {	}
	@Override
	public ArrayList<String> hndlFileSave(File file) {		return null;	}
	@Override
	protected void drawOnScreenStuffPriv(float modAmtMillis) {	}
	@Override
	protected void drawRightSideInfoBarPriv(float modAmtMillis) {	}
	@Override
	public void processTrajIndiv(DrawnSimpleTraj drawnTraj) {	}

}//DESSimWindow

