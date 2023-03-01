package MPM_SimMain;

import java.util.HashMap;

import MPM_CPUSim.ui.MPM_CPUSimWindow;
import MPM_CudaSim.ui.MPM_CudaSimWindow;
import base_UI_Objects.GUI_AppManager;
import base_UI_Objects.windowUI.sidebar.SidebarMenu;
import base_Utils_Objects.io.messaging.MsgCodes;

/**
* MPM Snow Simulation in CUDA
* @author john turner
*/

public class MPM_SimMain extends GUI_AppManager {
	
	public String prjNmShrt = "MPM_SnowSim_cpu_cuda";
	public String prjNmLong = "MPM Snow Simulation Multi-Threaded CPU and CUDA 9.1"; 
	public String projDesc = "Simulate numerous snow balls colliding using Material Point Method solved via MT CPU solver and CUDA 9.1 kernel.";
	
	//use sphere background for this program
	private boolean useSphereBKGnd = true;
	
	private String bkSkyBox = "winter.jpg";
	//bground color
	private final int[] bground = new int[]{244,244,244,255};				
	
	private final int gridDim = 1500;

	/**
	 * idx's in dispWinFrames for each window - 0 is always left side menu window
	 * Side menu is dispMenuIDX == 0
	 */
	private static final int
		dispMPMCudaWinIDX = 1,
		dispMPMCPUWinIDX = 2;
	
	/**
	 * # of visible windows including side menu (always at least 1 for side menu)
	 */
	private static final int numVisWins = 3;
	
///////////////
//CODE STARTS
///////////////	
//////////////////////////////////////////////// code
	//needs main to run project - do not modify this code in any way
	public static void main(String[] passedArgs) {		
		MPM_SimMain me = new MPM_SimMain();
	    MPM_SimMain.invokeProcessingMain(me, passedArgs);
	}//main	

	/**
	 * Set various relevant runtime arguments in argsMap
	 * @param _passedArgs command-line arguments
	 */
	@Override
	protected HashMap<String,Object> setRuntimeArgsVals(HashMap<String, Object> _passedArgsMap) {
		return  _passedArgsMap;
	}
	
	/**
	 * Called in pre-draw initial setup, before first init
	 * potentially override setup variables on per-project basis.
	 * Do not use for setting background color or Skybox anymore.
	 *  	(Current settings in my_procApplet) 	
	 *  	strokeCap(PROJECT);
	 *  	textSize(txtSz);
	 *  	textureMode(NORMAL);			
	 *  	rectMode(CORNER);	
	 *  	sphereDetail(4);	 * 
	 */
	@Override
	protected void setupAppDims_Indiv() {
		//Set grid dimensions
		setDesired3DGridDims(gridDim);		
	}
	@Override
	protected boolean getUseSkyboxBKGnd(int winIdx) {	return useSphereBKGnd;}
	@Override
	protected String getSkyboxFilename(int winIdx) {	return bkSkyBox;}
	@Override
	protected int[] getBackgroundColor(int winIdx) {return bground;}
	@Override
	protected int getNumDispWindows() {	return numVisWins;	}
	
	/**
	 * whether or not we want to restrict window size on widescreen monitors
	 * 
	 * @return 0 - use monitor size regardless
	 * 			1 - use smaller dim to be determine window 
	 * 			2+ - TBD
	 */
	@Override
	protected int setAppWindowDimRestrictions() {	return 1;}	
	
	@Override
	public String getPrjNmShrt() {		return prjNmShrt;}
	@Override
	public String getPrjNmLong() {		return prjNmLong;}
	@Override
	public String getPrjDescr() {		return projDesc;}	
	/**
	 * Set minimum level of message object console messages to display for this application. If null then all messages displayed
	 * @return
	 */
	@Override
	protected final MsgCodes getMinConsoleMsgCodes() {return MsgCodes.info1;}
	/**
	 * Set minimum level of message object log messages to save to log for this application. If null then all messages saved to log.
	 * @return
	 */
	@Override
	protected final MsgCodes getMinLogMsgCodes() {return MsgCodes.info1;}

	
	@Override
	protected void initBaseFlags_Indiv() {
		setBaseFlagToShow_debugMode(true);
		setBaseFlagToShow_saveAnim(true); 
		setBaseFlagToShow_runSim(true);
		setBaseFlagToShow_singleStep(true);
		setBaseFlagToShow_showRtSideMenu(true);
		setBaseFlagToShow_showDrawableCanvas(false);
	}
	@Override
	protected void initAllDispWindows() {
		showInfo = true;
		//titles and descs, need to be set before sidebar menu is defined
		String[] _winTitles = new String[]{"","Snow Balls!","Snow Ball"},
				_winDescr = new String[] {"", "Display Colliding Snowballs Simulated via MPM CUDA Solver", "Display Falling Snowball Simulated via CPU/MT Solver"};

		//instanced window dims when open and closed - only showing 1 open at a time - and init cam vals
		float[][] _floatDims  = new float[][] {getDefaultWinDimOpen(), getDefaultWinDimClosed(), getInitCameraValues()};	

		//Builds sidebar menu button config - application-wide menu button bar titles and button names
		String[] menuBtnTitles =  new String[]{"Functions 1","Functions 2","Functions 3","Functions 4"};
		String[][] menuBtnNames = new String[][] {	//each must have literals for every button defined in side bar menu, or ignored
			{"Restore Init Vals", "Func 01", "Func 02"},				//row 1
			{"Func 10", "Func 11", "Func 12", "Func 13"},	//row 2
			{"Func 10", "Func 11", "Func 12", "Func 13"},	//row 3
			{"Func 20", "Func 21", "Func 22", "Func 23"},	//row 4		
		};
		String[] dbgBtnNames = new String[] {"Debug 0","Debug 1","Debug 2","Debug 3","Debug 4"};
		buildSideBarMenu(_winTitles, menuBtnTitles, menuBtnNames, dbgBtnNames, true, false);

		initXORWins(new int[]{dispMPMCudaWinIDX, dispMPMCPUWinIDX},new int[]{dispMPMCudaWinIDX, dispMPMCPUWinIDX});
		//define windows
		/**
		 *  _winIdx The index in the various window-descriptor arrays for the dispWindow being set
		 *  _title string title of this window
		 *  _descr string description of this window
		 *  _dispFlags Essential flags describing the nature of the dispWindow for idxs : 
		 * 		0 : dispWinIs3d, 
		 * 		1 : canDrawInWin; 
		 * 		2 : canShow3dbox (only supported for 3D); 
		 * 		3 : canMoveView
		 *  _floatVals an array holding float arrays for 
		 * 				rectDimOpen(idx 0),
		 * 				rectDimClosed(idx 1),
		 * 				initCameraVals(idx 2)
		 *  _intClrVals and array holding int arrays for
		 * 				winFillClr (idx 0),
		 * 				winStrkClr (idx 1),
		 * 				winTrajFillClr(idx 2),
		 * 				winTrajStrkClr(idx 3),
		 * 				rtSideFillClr(idx 4),
		 * 				rtSideStrkClr(idx 5)
		 *  _sceneCenterVal center of scene, for drawing objects (optional)
		 *  _initSceneFocusVal initial focus target for camera (optional)
		 */
		
		int wIdx = dispMPMCudaWinIDX;
		//setInitDispWinVals(wIdx, _dimOpen, _dimClosed,new boolean[]{false,true,true,true}, new int[]{255,245,255,255},new int[]{0,0,0,255},new int[]{180,180,180,255},new int[]{100,100,100,255}); 		
		setInitDispWinVals(wIdx, _winTitles[wIdx], _winDescr[wIdx], getDfltBoolAra(true), _floatDims,		
				new int[][] {new int[]{255,245,255,255},new int[]{0,0,0,255},
					new int[]{180,180,180,255},new int[]{100,100,100,255},
					new int[]{0,0,0,200},new int[]{255,255,255,255}});

		dispWinFrames[wIdx] = new MPM_CudaSimWindow(ri, this, wIdx);
		
		wIdx = dispMPMCPUWinIDX;
		//setInitDispWinVals(wIdx, _dimOpen, _dimClosed,new boolean[]{false,true,true,true}, new int[]{255,245,255,255},new int[]{0,0,0,255},new int[]{180,180,180,255},new int[]{100,100,100,255}); 		
		setInitDispWinVals(wIdx, _winTitles[wIdx], _winDescr[wIdx], getDfltBoolAra(true), _floatDims,		
				new int[][] {new int[]{255,255,255,255},new int[]{0,0,0,255},
					new int[]{180,180,180,255},new int[]{100,100,100,255},
					new int[]{0,0,0,200},new int[]{255,255,255,255}});
		
		dispWinFrames[wIdx] = new MPM_CPUSimWindow(ri, this, wIdx);
			
		
	}//initVisOnce_Priv
	
	@Override
	protected void initOnce_Indiv() {
		//which objects to initially show
		setVisFlag(dispMPMCudaWinIDX, true);		
	}
	
	@Override
	protected void initProgram_Indiv() {}

	@Override
	public String[] getMouseOverSelBtnLabels() {
		return new String[0];
	}

	
	@Override
	protected void handleKeyPress(char key, int keyCode) {
		switch (key){
			case ' ' : {toggleSimIsRunning(); break;}							//run sim
			case 'f' : {dispWinFrames[curFocusWin].setInitCamView();break;}					//reset camera
			case 'a' :
			case 'A' : {toggleSaveAnim();break;}						//start/stop saving every frame for making into animation
			case 's' :
			case 'S' : {break;}//save(getScreenShotSaveName(prjNmShrt));break;}//save picture of current image			
			default : {	}
		}//switch	
	}
	/**
	 * gives multiplier based on whether shift, alt or cntl (or any combo) is pressed
	 */
	@Override
	public double clickValModMult(){return ((altIsPressed() ? .1 : 1.0) * (shiftIsPressed() ? 10.0 : 1.0));}	
	@Override
	public boolean isClickModUIVal() {
		//TODO change this to manage other key settings for situations where multiple simultaneous key presses are not optimal or conventient
		return altIsPressed() || shiftIsPressed();		
	}
	
	/**
	 * get the ui rect values of the "master" ui region (another window).
	 * this is so ui objects of one window can be made, clicked, and shown displaced from those of the parent window
	 */
	@Override
	public float[] getUIRectVals_Indiv(int idx, float[] menuClickDim){
		switch(idx){
			case dispMPMCudaWinIDX 	: {return menuClickDim;}
			case dispMPMCPUWinIDX	: {return menuClickDim;}
			default :  return menuClickDim;
			}
	}//getUIRectVals

	/**
	 * these tie using the UI buttons to modify the window in with using the boolean tags - PITA but currently necessary
	 */
	@Override
	public void handleShowWin(int btn, int val, boolean callFlags){//display specific windows - multi-select/ always on if sel
		if(!callFlags){//called from setflags - only sets button state in UI to avoid infinite loop
			setMenuBtnState(SidebarMenu.btnShowWinIdx,btn, val);
		} else {//called from clicking on buttons in UI
			//val is btn state before transition 
			boolean bVal = (val == 1?  false : true);
			//each entry in this array should correspond to a clickable window
			setVisFlag(winFlagsXOR[btn], bVal);
		}
	}//handleShowWin

	/**
	 * Individual extending Application Manager post-drawMe functions
	 * @param modAmtMillis
	 * @param is3DDraw
	 */
	@Override
	protected void drawMePost_Indiv(float modAmtMillis, boolean is3DDraw) {}

	//////////////////////////////////////////
	/// graphics and base functionality utilities and variables
	//////////////////////////////////////////
	
	/**
	 * return the number of visible window flags for this application
	 * @return
	 */
	@Override
	public int getNumVisFlags() {return numVisWins;}
	/**
	 * address all flag-setting here, so that if any special cases need to be addressed they can be easily
	 */
	@Override
	protected void setVisFlag_Indiv(int idx, boolean val ){
		switch (idx){
			case dispMPMCudaWinIDX		: {setWinFlagsXOR(dispMPMCudaWinIDX, val); break;}
			case dispMPMCPUWinIDX		: {setWinFlagsXOR(dispMPMCPUWinIDX, val); break;}
			default : {break;}
		}
	}//setFlags  
	
	@Override
	public int[] getClr_Custom(int colorVal, int alpha) {	return new int[] {255,255,255,alpha};}

	@Override
	protected void setSmoothing() {
		ri.setSmoothing(0);		
	}

}//class MPM_SimMain
