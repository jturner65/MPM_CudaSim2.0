package MPM_CudaSim;

import java.io.File;
import java.util.*;

import base_JavaProjTools_IRender.base_Render_Interface.IRenderInterface;
import base_UI_Objects.GUI_AppManager;
import base_UI_Objects.windowUI.base.base_UpdateFromUIData;
import base_UI_Objects.windowUI.base.myDispWindow;
import base_UI_Objects.windowUI.drawnObjs.myDrawnSmplTraj;
import base_Math_Objects.vectorObjs.doubles.myPoint;
import base_Math_Objects.vectorObjs.doubles.myVector;

public class MPM_SimWindow extends myDispWindow {

	//simulator world within which simulation executes
	//TODO support multiple sim worlds possibly, with different configurations
	public MPM_Abs_CUDASim currSim;
	//grid variables
	private int numGridCells = MPM_Abs_CUDASim.numGridCells;
	private float cellSize = MPM_Abs_CUDASim.cellSize;
	private int numParts = MPM_Abs_CUDASim.numPartsUI_Init;
	private float csCube = cellSize*cellSize*cellSize;
	
	//motion
	///////////
	//ui vals
	public final static int
		gIDX_timeStep		 		= 0,
		gIDX_simStepsPerFrame		= 1,
		gIDX_numParticles			= 2,
		gIDX_partMass				= 3,
		gIDX_gridCellSize			= 4,
		gIDX_gridCount				= 5,
		gIDX_initYoungMod 			= 6,
		gIDX_poissonRatio 			= 7,
		gIDX_hardeningCoeff 		= 8,
		gIDX_criticalCompression 	= 9,
		gIDX_criticalStretch 		= 10,
		gIDX_alphaPicFlip 			= 11,
		gIDX_wallFricCoeff 			= 12,
		gIDX_collFricCoeff			= 13;
	
	//dims x density
	public float partMass = csCube*50.0f;

	//display variables
	//private float[] UIrectBox;	//box holding x,y,w,h values of black rectangle to hold UI sim display values
	
	/////////
	//custom debug/function ui button names -empty will do nothing
	public String[] menuDbgBtnNames = new String[] {};//must have literals for every button or this is ignored
	public String[] menuFuncBtnNames = new String[] {"Func 1", "Func 2", "Func 3", "Func 4", "Func 5"};//must have literals for every button or ignored
	public String[][] menuBtnNames = new String[][] {	//each must have literals for every button defined in side bar menu, or ignored
		{"Func 00", "Func 01", "Func 02"},				//row 1
		{"Func 10", "Func 11", "Func 12", "Func 13"},	//row 2
		{"Func 10", "Func 11", "Func 12", "Func 13"},	//row 3
		{"Func 20", "Func 21", "Func 22", "Func 23"},	//row 4
		{"Func 30", "Func 31", "Func 32", "Func 33","Func 34"}	//dbg
	};
	
	//private child-class flags - window specific
	public static final int 
		debugAnimIDX 			= 0,			//debug
		resetSimIDX				= 1,			//whether or not to reset sim	
		resetReqIDX				= 2,			//reset of sim is required for best results
		showCollider			= 3,			//show collider cylinder
		showParticleVelArrows 	= 4,			//plot velocity arrows for each particle	
		showGrid				= 5,  			//plot the computational grid
		showGridVelArrows 		= 6,			//plot velocity arrows for each gridNode
		showGridAccelArrows 	= 7,			//plot acceleration arrows for each gridNode
		showGridMass  			= 8,			//plot variable sized spheres proportional to gridnode mass
		showActiveNodes  		= 9,				//show the grid nodes influenced by each particle
		showExecTime			= 10,			//show time of execution for each step in simulation
		runMultiThdSim			= 11;			//run the multi-threaded version of the simulation

	public static final int numPrivFlags = 12;
		
	public MPM_SimWindow(IRenderInterface _p, GUI_AppManager _AppMgr, String _n, int _flagIdx, int[] fc, int[] sc, float[] rd, float[] rdClosed, String _winTxt) {
		super(_p, _AppMgr, _n, _flagIdx, fc, sc, rd, rdClosed, _winTxt);
		super.initThisWin(false);
	}//DancingBallWin
	
	@Override
	//initialize all private-flag based UI buttons here - called by base class
	public int initAllPrivBtns(ArrayList<Object[]> tmpBtnNamesArray){	
		tmpBtnNamesArray.add(new Object[]{"Visualization Debug",    "Enable Debug",       debugAnimIDX});          
		tmpBtnNamesArray.add(new Object[]{"Resetting Simulation",   "Reset Simulation",   resetSimIDX});           
		tmpBtnNamesArray.add(new Object[]{"Showing Collider",       "Show Collider",      showCollider});          
		tmpBtnNamesArray.add(new Object[]{"Showing Particle Vel",   "Show Particle Vel",  showParticleVelArrows});  
		tmpBtnNamesArray.add(new Object[]{"Showing Grid",           "Show Grid",          showGrid});           
		tmpBtnNamesArray.add(new Object[]{"Showing Grid Vel",       "Show Grid Vel",      showGridVelArrows});     
		tmpBtnNamesArray.add(new Object[]{"Showing Grid Accel",     "Show Grid Accel",    showGridAccelArrows});    
		tmpBtnNamesArray.add(new Object[]{"Showing Grid Mass",      "Show Grid Mass",     showGridMass});         
		tmpBtnNamesArray.add(new Object[]{"Showing Active Nodes",   "Show Active Nodes",  showActiveNodes});     
		tmpBtnNamesArray.add(new Object[]{"Showing Execution Time", "Show Execution Time",showExecTime});         
		tmpBtnNamesArray.add(new Object[]{"Running Multi-Thd Sim",   "Run Multi-Thd Sim",   runMultiThdSim });         
		return numPrivFlags;
	}//initAllPrivBtns
	//set labels of boolean buttons 

	
	@Override
	protected void initMe() {//all ui objects set by here
		//this window is runnable
		setFlags(isRunnable, true);
		//this window uses a customizable camera
		setFlags(useCustCam, true);
		//this window uses right side info window
		//setFlags(drawRightSideMenu, true);
		//called once
		//initPrivFlags(numPrivFlags);
		
//		//init simulation construct here
		//init simulation construct here
		msgObj.dispInfoMessage("MPM_SimWindow","initMe","Start building simulation now.");
		currSim = new MPM_Cuda2Balls(pa,numGridCells, cellSize,numParts);		
		//initialize simulation here to simple world sim
		custMenuOffset = uiClkCoords[3];	//495	
		
		AppMgr.setAllMenuBtnNames(menuBtnNames);	

	}//initMe	

	@Override
	protected int[] getFlagIDXsToInitToTrue() {
		// TODO Auto-generated method stub
		return new int[] {showCollider};
	}

	//call this to initialize or reinitialize simulation (on reset)
	protected void reinitSim() {
		AppMgr.setSimIsRunning( false);		
		currSim.resetSim(true);
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
				//TODO set debug state for simulation here
				break;}
			case resetSimIDX			: {
				if(val) {
					reinitSim();
					addPrivBtnToClear(resetSimIDX);
				}
				break;}						
			case resetReqIDX			: {//simulation reset not forced, but required for good results
				
				break;}	
			case showCollider			: {//show collider
				currSim.setSimFlags(MPM_Abs_CUDASim.showCollider, val);
				break;}
			case showParticleVelArrows	: {
				currSim.setSimFlags(MPM_Abs_CUDASim.showParticleVelArrows, val);
				break;} 
			case showGrid				: {
				currSim.setSimFlags(MPM_Abs_CUDASim.showGrid, val);
				break;} 				
			case showGridVelArrows 		: {
				currSim.setSimFlags(MPM_Abs_CUDASim.showGridVelArrows, val);
				break;} 		
			case showGridAccelArrows	: {
				currSim.setSimFlags(MPM_Abs_CUDASim.showGridAccelArrows, val);
				break;} 	
			case showGridMass  			: {
				currSim.setSimFlags(MPM_Abs_CUDASim.showGridMass, val);
				break;} 		
			case showActiveNodes 	: {
				currSim.setSimFlags(MPM_Abs_CUDASim.showActiveNodes, val);
				break;} 	
			case showExecTime			: { 
				break;}
			default:					{}
		}		
	}//setPrivFlags	
		
	//initialize structure to hold modifiable menu regions
	@Override
	protected void setupGUIObjsAras(TreeMap<Integer, Object[]> tmpUIObjArray, TreeMap<Integer, String[]> tmpListObjVals){	
		
		tmpUIObjArray.put(gIDX_timeStep , new Object[] {new double[]{.00005f, .0010f, .00005f}, 1.0*MPM_Abs_CUDASim.getDeltaT(),  "Sim Time Step", new boolean[]{false, false, true}});//delta T for simulation init  MPM_ABS_Sim.deltaT = 1e-3f;
		tmpUIObjArray.put(gIDX_simStepsPerFrame, new Object[] {new double[]{1, 20, 1}, 1.0*MPM_Abs_CUDASim.simStepsPerFrame, "Sim Steps per Drawn Frame",  new boolean[]{true, false, true}});//gIDX_simStepsPerFrame  init 5
		tmpUIObjArray.put(gIDX_numParticles, new Object[] {new double[]{100, 1000000, 100}, 1.0*numParts, "# of Particles", new boolean[]{true, false, true}});//number of particles
		tmpUIObjArray.put(gIDX_partMass, new Object[] {new double[]{.00005, 5.00, .00005}, 1.0*partMass, "Particle Mass", new boolean[]{false, false, true}});//particle mass
		tmpUIObjArray.put(gIDX_gridCellSize, new Object[] {new double[]{.001, .1, .001}, 1.0*cellSize, "Grid Cell Size", new boolean[]{false, false, true}});//grid cell size
		tmpUIObjArray.put(gIDX_gridCount, new Object[] {new double[]{50, 300, 1}, 1.0*numGridCells,  "Grid Cell Count Per Side", new boolean[]{true, false, true}}); //# of grid cells per side
		tmpUIObjArray.put(gIDX_initYoungMod, new Object[] {new double[]{1000.0f, 100000.0f, 100.0f}, 1.0*myMaterial.base_initYoungMod, "Young's Modulus", new boolean[]{false, false, true}});//gIDX_initYoungMod init 4.8e4f, 
		tmpUIObjArray.put(gIDX_poissonRatio, new Object[] {new double[]{.01f, 1.0f, .01f}, 1.0*myMaterial.base_poissonRatio, "Poisson Ratio", new boolean[]{false, false, true}});//gIDX_poissonRatio init 0.2f,  
		tmpUIObjArray.put(gIDX_hardeningCoeff , new Object[] {new double[]{1.0f, 100.0f, 1.0f}, 1.0*myMaterial.base_hardeningCoeff, "Hardening Coefficient", new boolean[]{false, false, true}});//gIDX_hardeningCoeff init 15.0f, 
		tmpUIObjArray.put(gIDX_criticalCompression, new Object[] {new double[]{0.001f, 0.1f, 0.001f}, 1.0*myMaterial.base_criticalCompression, "Critical Compression", new boolean[]{false, false, true}});//gIDX_criticalCompression  init .019f, 
		tmpUIObjArray.put(gIDX_criticalStretch , new Object[] {new double[]{0.0005f, 0.01f, 0.0005f}, 1.0*myMaterial.base_criticalStretch, "Critical Stretch",  new boolean[]{false, false, true}});//gIDX_criticalStretch init .0075f, 
		tmpUIObjArray.put(gIDX_alphaPicFlip, new Object[] {new double[]{0.0f, 1.0f, 0.01f}, 1.0*myMaterial.base_alphaPicFlip, "Particle PIC/FLIP Vel Ratio", new boolean[]{false, false, true}});//gIDX_alphaPicFlip init 0.95f, 
		tmpUIObjArray.put(gIDX_wallFricCoeff, new Object[] {new double[]{0.01f, 1.0f, 0.01f}, 1.0*MPM_Abs_CUDASim.wallFric, "Wall Friction Coefficient",  new boolean[]{false, false, true}});//gIDX_wallfricCoeffinit 1.0f  
		tmpUIObjArray.put(gIDX_collFricCoeff, new Object[] {new double[]{0.01f, 1.0f, 0.01f}, 1.0*MPM_Abs_CUDASim.collFric, "Collider Friction Coefficient", new boolean[]{false, false, true}});//gIDX_collfricCoeffinit 1.0f  
	}//setupGUIObjsAras
	
	@Override
	protected void setUIWinVals(int UIidx) {
		float val = (float)guiObjs[UIidx].getVal();
		if(val != uiVals[UIidx]){//if value has changed...
			uiVals[UIidx] = val;
			switch(UIidx){		
			case gIDX_timeStep 			:{
				currSim.setDeltaT(val);
				break;} 	
			case gIDX_numParticles				:{
				numParts = (int)val;
				currSim.setGridValsAndInit(numGridCells,  cellSize, numParts);
				break;}
			case gIDX_partMass 					:{
				partMass = val;
				//pa.outStr2Scr("Changing simulation parameter (particle mass) may require simulation reset; Not doing so may result in instability");
				break;}
			case gIDX_gridCellSize				:{
				cellSize = val;
				currSim.setGridValsAndInit(numGridCells,  cellSize, numParts);
				break;}
			case gIDX_gridCount					:{
				numGridCells = (int)val;
				currSim.setGridValsAndInit(numGridCells,  cellSize, numParts);				
				break;}
			case gIDX_initYoungMod 				:{
				currSim.mat.setYoungModulus(val);
				reinitSim();
				break;} 		
			case gIDX_poissonRatio 				:{
				currSim.mat.setPoissonRatio(val);
				reinitSim();
				break;} 		
			case gIDX_hardeningCoeff 			:{
				currSim.mat.setHC(val);
				reinitSim();
				break;} 	
			case gIDX_criticalCompression 		:{
				currSim.mat.setCC(val);
				reinitSim();
				break;}
			case gIDX_criticalStretch 			:{
				currSim.mat.setCS(val);
				reinitSim();
				break;} 	
			case gIDX_alphaPicFlip 				:{
				currSim.mat.setAlpha(val);
				break;} 		
			case gIDX_wallFricCoeff 				:{
				currSim.wallFric = val;
				reinitSim();
				break;} 
			case gIDX_collFricCoeff 				:{
				currSim.collFric = val;
				reinitSim();
				break;} 
			case gIDX_simStepsPerFrame 			:{
				MPM_Abs_CUDASim.simStepsPerFrame = (int)val;
				break;}
			default : {break;}
			}
		}
	}

	
	@Override
	public void initDrwnTrajIndiv(){}
	
//	public void setLights(){
//		pa.ambientLight(102, 102, 102);
//		pa.lightSpecular(204, 204, 204);
//		pa.directionalLight(180, 180, 180, 0, 1, -1);	
//	}	
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
		//pa.outStr2Scr("took : " + (pa.millis() - stVal) + " millis to simulate");
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
//		curMseLookVec = pa.c.getMse2DtoMse3DinWorld(pa.sceneCtrVals[pa.sceneIDX]);			//need to be here
//		curMseLoc3D = pa.c.getMseLoc(pa.sceneCtrVals[pa.sceneIDX]);
		
		//TODO draw simulation results here
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
	protected void stopMe() {	msgObj.dispInfoMessage("MPM_SimWindow","stopMe","Simulation Finished");	}	
	
	@Override
	public void processTrajIndiv(myDrawnSmplTraj drawnNoteTraj){	}
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


	@Override
	protected final void launchMenuBtnHndlr(int funcRow, int btn) {	}

	@Override
	public void handleSideMenuMseOvrDispSel(int btn, boolean val) {}

	@Override
	public void handleSideMenuDebugSel(int btn, int val) {	}
	
	
	@Override
	protected String[] getSaveFileDirNamesPriv() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected myPoint getMsePtAs3DPt(myPoint mseLoc) {
		// TODO Auto-generated method stub
		return new myPoint(mseLoc.x,mseLoc.y,0);
	}

	@Override
	protected void setVisScreenDimsPriv() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void setCustMenuBtnNames() {
		AppMgr.setAllMenuBtnNames(menuBtnNames);	
	}

	@Override
	public void hndlFileLoad(File file, String[] vals, int[] stIdx) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public ArrayList<String> hndlFileSave(File file) {		return null;	}

	@Override
	protected void drawOnScreenStuffPriv(float modAmtMillis) {	}

	@Override
	protected void drawRightSideInfoBarPriv(float modAmtMillis) {	}

	@Override
	protected base_UpdateFromUIData buildUIDataUpdateObject() {		return null;	}

	@Override
	protected void buildUIUpdateStruct_Indiv(TreeMap<Integer, Integer> intValues, TreeMap<Integer, Float> floatValues,TreeMap<Integer, Boolean> boolValues) {
	}


}//DESSimWindow

