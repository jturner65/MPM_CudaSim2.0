package MPM_CudaSim;

import java.util.*;

public class MPM_SimWindow extends myDispWindow {

	//simulator world within which simulation executes
	//TODO support multiple sim worlds possibly, with different configurations
	public MPM_ABS_Sim currSim;
	//grid variables
	private int numGridCells = 200;
	private float cellSize = .10f;
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
	//initial values - need one per object
	public float[] uiVals = new float[]{
		MPM_ABS_Sim.getDeltaT(),
		MPM_ABS_Sim.simStepsPerFrame,
		MPM_ABS_Sim.numPartsUI_Init,
		partMass,
		cellSize,
		numGridCells,
		myMaterial.base_initYoungMod,					//4.8e4f,
		myMaterial.base_poissonRatio,					//0.2f,
		myMaterial.base_hardeningCoeff,					//15.0f, 
		myMaterial.base_criticalCompression,			//.040f, 
		myMaterial.base_criticalStretch,				//.0075f, 
		myMaterial.base_alphaPicFlip,					//0.95f, 
		MPM_ABS_Sim.wallFric,							//init value for wall friction
		MPM_ABS_Sim.collFric,							//init value for collider friction
			//
	};			//values of 8 ui-controlled quantities
	public final int numGUIObjs = uiVals.length;											//# of gui objects for ui	
	
	//display variables
	//private float[] UIrectBox;	//box holding x,y,w,h values of black rectangle to hold UI sim display values
	
	/////////
	//custom debug/function ui button names -empty will do nothing
	public String[] menuDbgBtnNames = new String[] {};//must have literals for every button or this is ignored
	public String[] menuFuncBtnNames = new String[] {"Func 1", "Func 2", "Func 3", "Func 4", "Func 5"};//must have literals for every button or ignored
	
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
		
	public MPM_SimWindow(MPM_SimMain _p, String _n, int _flagIdx, int[] fc, int[] sc, float[] rd, float[] rdClosed, String _winTxt, boolean _canDrawTraj) {
		super(_p, _n, _flagIdx, fc, sc, rd, rdClosed, _winTxt, _canDrawTraj);
		float stY = rectDim[1]+rectDim[3]-4*yOff,stYFlags = stY + 2*yOff;
		trajFillClrCnst = MPM_SimMain.gui_DarkCyan;		
		trajStrkClrCnst = MPM_SimMain.gui_Cyan;
		super.initThisWin(_canDrawTraj, true, false);
	}//DancingBallWin
	
	@Override
	//initialize all private-flag based UI buttons here - called by base class
	public void initAllPrivBtns(){
		truePrivFlagNames = new String[]{								//needs to be in order of privModFlgIdxs
				"Visualization Debug","Resetting Simulation", "Showing Collider", "Showing Particle Vel", "Showing Grid", 
				"Showing Grid Vel", "Showing Grid Accel", "Showing Grid Mass", "Showing Active Nodes",
				"Showing Execution Time","Running Multi-Thd Sim"
		};
		falsePrivFlagNames = new String[]{			//needs to be in order of flags
				"Enable Debug","Reset Simulation", "Show Collider", "Show Particle Vel", "Show Grid", 
				"Show Grid Vel", "Show Grid Accel", "Show Grid Mass", "Show Active Nodes",
				"Show Execution Time","Run Multi-Thd Sim"
		};
		privModFlgIdxs = new int[]{
				debugAnimIDX, resetSimIDX, showCollider, showParticleVelArrows, showGrid, 
				showGridVelArrows, showGridAccelArrows, showGridMass, showActiveNodes,
				showExecTime, runMultiThdSim
		};
		numClickBools = privModFlgIdxs.length;	
		initPrivBtnRects(0,numClickBools);
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
		initPrivFlags(numPrivFlags);
		
		//init simulation construct here
		pa.outStr2Scr("Start building simulation : " + pa.millis());
		currSim = new MPM_BaseSim(numGridCells, cellSize);		
		setPrivFlags(showCollider, true);
		//initialize simulation here to simple world sim
		custMenuOffset = uiClkCoords[3];	//495	
	}//initMe	
	
	//call this to initialize or reinitialize simulation (on reset)
	protected void reinitSim() {
		pa.setFlags(pa.runSim, false);		
		currSim.resetSim();
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
				currSim.setSimFlags(MPM_ABS_Sim.showCollider, val);
				break;}
			case showParticleVelArrows	: {
				currSim.setSimFlags(MPM_ABS_Sim.showParticleVelArrows, val);
				break;} 
			case showGrid				: {
				currSim.setSimFlags(MPM_ABS_Sim.showGrid, val);
				break;} 				
			case showGridVelArrows 		: {
				currSim.setSimFlags(MPM_ABS_Sim.showGridVelArrows, val);
				break;} 		
			case showGridAccelArrows	: {
				currSim.setSimFlags(MPM_ABS_Sim.showGridAccelArrows, val);
				break;} 	
			case showGridMass  			: {
				currSim.setSimFlags(MPM_ABS_Sim.showGridMass, val);
				break;} 		
			case showActiveNodes 	: {
				currSim.setSimFlags(MPM_ABS_Sim.showActiveNodes, val);
				break;} 	
			case showExecTime			: { 
				break;}
			default:					{}
		}		
	}//setPrivFlags	
		
	//initialize structure to hold modifiable menu regions
	@Override
	protected void setupGUIObjsAras(){	
		//pa.outStr2Scr("setupGUIObjsAras start");
		guiMinMaxModVals = new double [][]{
			{.00005f,.0010f,.00005f},						//delta T for simulation    init  MPM_ABS_Sim.deltaT = 1e-3f;
			{1,20,1},									//gIDX_simStepsPerFrame  init 5
			{100,1000000,100},							//number of particles
			{.00005, 5.00, .00005},						//particle mass
			{.001, .1, .001},							//grid cell size
			{50,300,1},									//# of grid cells per side
			{1000.0f, 100000.0f, 100.0f},				//gIDX_initYoungMod 		init 4.8e4f, 
			{.01f,1.0f,.01f},							//gIDX_poissonRatio 		init 0.2f,  
			{1.0f,100.0f,1.0f},							//gIDX_hardeningCoeff 	    init 15.0f, 
			{0.001f, 0.1f, 0.001f},						//gIDX_criticalCompression  init .019f, 
			{0.0005f,0.01f,0.0005f},					//gIDX_criticalStretch 	    init .0075f,
			{0.0f,1.0f,0.01f},							//gIDX_alphaPicFlip 		init 0.95f, 
			{0.01f,1.0f,0.01f},							//gIDX_wallfricCoeff 			init 1.0f  
			{0.01f,1.0f,0.01f},							//gIDX_collfricCoeff 			init 1.0f  
		};		//min max mod values for each modifiable UI comp	

		guiStVals = new double[]{
			uiVals[gIDX_timeStep],
			uiVals[gIDX_simStepsPerFrame],
			uiVals[gIDX_numParticles],
			uiVals[gIDX_partMass],
			uiVals[gIDX_gridCellSize],
			uiVals[gIDX_gridCount],
			uiVals[gIDX_initYoungMod], 		
			uiVals[gIDX_poissonRatio],		
			uiVals[gIDX_hardeningCoeff], 	  
			uiVals[gIDX_criticalCompression],
			uiVals[gIDX_criticalStretch], 	  
			uiVals[gIDX_alphaPicFlip], 		
			uiVals[gIDX_wallFricCoeff], 	
			uiVals[gIDX_collFricCoeff], 	
		};								//starting value
		
		guiObjNames = new String[]{
				"Sim Time Step",
				"Sim Steps per Drawn Frame",
				"# of Particles",
				"Particle Mass",
				"Grid Cell Size",
				"Grid Cell Count Per Side",
				"Young's Modulus",
				"Poisson Ratio",
				"Hardening Coefficient",
				"Critical Compression",
				"Critical Stretch",
				"Particle PIC/FLIP Vel Ratio",
				"Wall Friction Coefficient",	
				"Collider Friction Coefficient",
		};								//name/label of component	
		
		//idx 0 is treat as int, idx 1 is obj has list vals, idx 2 is object gets sent to windows
		guiBoolVals = new boolean [][]{
			{false, false, true},	                     //delta T for simulation    
			{true, false, true},                         //gIDX_simStepsPerFrame  ini
			{true, false, true},	                     //number of particles       
			{false, false, true},	                     //particle mass             
			{false, false, true},	                     //grid cell size
			{true, false, true},	                     //# of grid cells per side
			{false, false, true},	                     //gIDX_initYoungMod 		
			{false, false, true},	                     //gIDX_poissonRatio 		
			{false, false, true},	                     //gIDX_hardeningCoeff 	    
			{false, false, true},	                     //gIDX_criticalCompression  
			{false, false, true},	                     //gIDX_criticalStretch 	    
			{false, false, true},	                     //gIDX_alphaPicFlip 		
			{false, false, true},                        //gIDX_wallFricCoeff 			
			{false, false, true},                        //gIDX_collFricCoeff 			
		};						//per-object  list of boolean flags
		
		//since horizontal row of UI comps, uiClkCoords[2] will be set in buildGUIObjs		
		guiObjs = new myGUIObj[numGUIObjs];			//list of modifiable gui objects
		if(numGUIObjs > 0){
			buildGUIObjs(guiObjNames,guiStVals,guiMinMaxModVals,guiBoolVals,new double[]{xOff,yOff});			//builds a horizontal list of UI comps
		}
		
//		setupGUI_XtraObjs();
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
				currSim.setNumPartsAndReset((int)val);
				break;}
			case gIDX_partMass 					:{
				partMass = val;
				//pa.outStr2Scr("Changing simulation parameter (particle mass) may require simulation reset; Not doing so may result in instability");
				break;}
			case gIDX_gridCellSize				:{
				cellSize = val;
				currSim.setGridValsAndInit(numGridCells,  cellSize);
				break;}
			case gIDX_gridCount					:{
				numGridCells = (int)val;
				currSim.setGridValsAndInit(numGridCells,  cellSize);				
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
				MPM_ABS_Sim.simStepsPerFrame = (int)val;
				break;}
			default : {break;}
			}
		}
	}
	//if any ui values have a string behind them for display
	@Override
	protected String getUIListValStr(int UIidx, int validx) {			
		switch(UIidx){
			//case gIDX_UAVTeamSize : {return uavTeamSizeVals[(validx % uavTeamSizeVals.length)];}
			default : {break;}
		}
		return "";
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
		pa.camera(camVals[0],camVals[1],camVals[2],camVals[3],camVals[4],camVals[5],camVals[6],camVals[7],camVals[8]);      
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
		} else if (getPrivFlags(runMultiThdSim)) {
			//done = currSim.simMeMThd(modAmtMillis);		
			done = currSim.simMeCuda(modAmtMillis);	
		} else {//base sim calculation		
			//done = currSim.simMe(modAmtMillis);
			done = currSim.simMeCuda(modAmtMillis);	
		}
		
		return done;	
	}//simMe
	
	
	@Override
	//draw 2d constructs over 3d area on screen - draws behind menu section
	//modAmtMillis is in milliseconds
	protected void drawRightSideInfoBar(float modAmtMillis) {
		pa.pushMatrix();pa.pushStyle();
		//display current simulation variables
		
		pa.popStyle();pa.popMatrix();				
	}//drawOnScreenStuff
	
	@Override
	//animTimeMod is in seconds.
	protected void drawMe(float animTimeMod) {
//		curMseLookVec = pa.c.getMse2DtoMse3DinWorld(pa.sceneCtrVals[pa.sceneIDX]);			//need to be here
//		curMseLoc3D = pa.c.getMseLoc(pa.sceneCtrVals[pa.sceneIDX]);
		
		//TODO draw simulation results here
		currSim.drawMe(pa, animTimeMod);
	}//drawMe	
	
	
	//draw custom 2d constructs below interactive component of menu
	@Override
	public void drawCustMenuObjs(){
		pa.pushMatrix();				pa.pushStyle();		
		//all sub menu drawing within push mat call
		pa.translate(0,custMenuOffset+yOff);		
		//draw any custom menu stuff here
		pa.popStyle();					pa.popMatrix();		
	}//drawCustMenuObjs

	
	@Override
	protected void closeMe() {
		//things to do when swapping this window out for another window - release objects that take up a lot of memory, for example.
	}	
	
	@Override
	protected void showMe() {
		//things to do when swapping into this window - reinstance released objects, for example.
		pa.setMenuDbgBtnNames(menuDbgBtnNames);
		pa.setMenuFuncBtnNames(menuFuncBtnNames);
	}
	
	@Override
	//stopping simulation
	protected void stopMe() {
		System.out.println("Simulation Finished");	
	}
	
	//custom functions launched by UI input
	//if launching threads for custom functions, need to remove clearFuncBtnState call in function below and call clearFuncBtnState when thread ends
	private void custFunc0(){
		clearFuncBtnState(0,false);
	}			
	private void custFunc1(){
		clearFuncBtnState(1,false);
	}	
	
	private void custFunc2(){	
		//custom function code here
		clearFuncBtnState(2,false);
	}			
	private void custFunc3(){	
		//custom function code here
		clearFuncBtnState(3,false);
	}			
	private void custFunc4(){	
		//custom function code here
		clearFuncBtnState(4,false);
	}		
	@Override
	public void clickFunction(int btnNum) {
		//pa.outStr2Scr("click cust function in "+name+" : btn : " + btnNum);
		switch(btnNum){
			case 0 : {	custFunc0();	break;}
			case 1 : {	custFunc1();	break;}
			case 2 : {	custFunc2();	break;}
			case 3 : {	custFunc3();	break;}
			case 4 : {	custFunc4();	break;}
			default : {break;}
		}	
	}		//only for display windows
	
	private void toggleDbgBtn(int idx, boolean val) {
		setPrivFlags(idx, !getPrivFlags(idx));
	}
	
	//debug function
	//if launching threads for debugging, need to remove clearDBGState call in function below and call clearDBGState when thread ends
	private void dbgFunc0() {			
		clearDBGBtnState(0,false);
	}	
	private void dbgFunc1(){
		clearDBGBtnState(1,false);
	}	
	private void dbgFunc2(){
		clearDBGBtnState(2,false);
	}	
	private void dbgFunc3(){	
		clearDBGBtnState(3,false);
	}	
	private void dbgFunc4(){	
		clearDBGBtnState(4,false);
	}	
	private void dbgFunc5(){	
		clearDBGBtnState(5,false);
	}	

	@Override
	public void clickDebug(int btnNum){
		pa.outStr2Scr("click debug in "+name+" : btn : " + btnNum);
		switch(btnNum){
			case 0 : {	dbgFunc0();	break;}//verify priority queue functionality
			case 1 : {	dbgFunc1();	break;}//verify FEL pq integrity
			case 2 : {	dbgFunc2();	break;}
			case 3 : {	dbgFunc3();	break;}
			case 4 : {	dbgFunc4();	break;}
			case 5 : {	dbgFunc5();	break;}
			default : {break;}
		}		
	}//clickDebug
	
	@Override
	public void hndlFileLoadIndiv(String[] vals, int[] stIdx) {}
	@Override
	public List<String> hndlFileSaveIndiv() {List<String> res = new ArrayList<String>();return res;}
	@Override
	protected void processTrajIndiv(myDrawnSmplTraj drawnNoteTraj){	}
	@Override
	protected myPoint getMsePtAs3DPt(int mouseX, int mouseY){return pa.P(mouseX,mouseY,0);}
	@Override
	protected boolean hndlMouseMoveIndiv(int mouseX, int mouseY, myPoint mseClckInWorld){
		return false;
	}
	
	//alt key pressed handles trajectory
	
	//cntl key pressed handles unfocus of spherey
	@Override
	protected boolean hndlMouseClickIndiv(int mouseX, int mouseY, myPoint mseClckInWorld, int mseBtn) {	
		boolean res = checkUIButtons(mouseX, mouseY);	
		return res;}//hndlMouseClickIndiv
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
}//DESSimWindow
