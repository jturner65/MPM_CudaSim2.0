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

	private final int
		showUIMenu = 0,
		showMPMCudawin = 1,
		showMPMCPUwin = 2;
	public final int numVisFlags = 3;
	
	//idx's in dispWinFrames for each window - 0 is always left side menu window
	private static final int
		dispMPMCudaWinIDX = 1,
		dispMPMCPUWinIDX = 2;
	
	private final int[] bground = new int[]{244,244,244,255};		//bground color		
	
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
	
	/**
	 * instance-specific setup code
	 */
	protected void setup_Indiv() {
		//modify default grid dims to be 1500x1500x1500
		setDesired3DGridDims(1500);
		//TODO move to window to set up specific background for each different "scene" type
		//PImage bgrndTex = loadImage("bkgrndTex.jpg"); 
		//PImage bgrndTex = loadImage("sky_1.jpg");
		if(useSphereBKGnd) {			pa.loadBkgndSphere("winter.jpg");	} else {		setBkgrnd();	}
	}	
	@Override
	public void setBkgrnd(){
		//TODO move to Base_DispWindow	
		if(useSphereBKGnd) { pa.setBkgndSphere();	} else {pa.setRenderBackground(bground[0],bground[1],bground[2],bground[3]);		}
	}//setBkgrnd

	@Override
	protected void initBaseFlags_Indiv() {
		setBaseFlagToShow_debugMode(false);
		setBaseFlagToShow_saveAnim(true); 
		setBaseFlagToShow_runSim(true);
		setBaseFlagToShow_singleStep(true);
		setBaseFlagToShow_showRtSideMenu(true);
	}
	@Override
	protected void initAllDispWindows() {
		showInfo = true;
		
		//includes 1 for menu window (never < 1) - always have same # of visFlags as Base_DispWindows
		int numWins = numVisFlags;		
		//titles and descs, need to be set before sidebar menu is defined
		String[] _winTitles = new String[]{"","Snow Balls!","Snow Ball"},
				_winDescr = new String[] {"", "Display Colliding Snowballs Simulated via MPM CUDA Solver", "Display Falling Snowball Simulated via CPU/MT Solver"};
		initWins(numWins,_winTitles, _winDescr);
		//call for menu window
		buildInitMenuWin();
		//instanced window dimensions when open and closed - only showing 1 open at a time
		float[] _dimOpen  = getDefaultWinDimOpen(), 
				_dimClosed  = getDefaultWinDimClosed();	
		//menu bar init
		int wIdx = dispMenuIDX,fIdx=showUIMenu;
		String[] menuBtnTitles =  new String[]{"Functions 1","Functions 2","Functions 3","Functions 4"};
		String[][] menuBtnNames = new String[][] {	//each must have literals for every button defined in side bar menu, or ignored
			{"Restore Init Vals", "Func 01", "Func 02"},				//row 1
			{"Func 10", "Func 11", "Func 12", "Func 13"},	//row 2
			{"Func 10", "Func 11", "Func 12", "Func 13"},	//row 3
			{"Func 20", "Func 21", "Func 22", "Func 23"},	//row 4		
		};
		String[] dbgBtnNames = new String[] {"Debug 0","Debug 1","Debug 2","Debug 3","Debug 4"};
		dispWinFrames[wIdx] = buildSideBarMenu(wIdx, fIdx,menuBtnTitles, menuBtnNames, dbgBtnNames, true, false);

		//setInitDispWinVals : use this to define the values of a display window
		//int _winIDX, 
		//float[] _dimOpen, float[] _dimClosed  : dimensions opened or closed
		//String _ttl, String _desc 			: window title and description
		//boolean[] _dispFlags 					: 
		//   flags controlling display of window :  idxs : 0 : canDrawInWin; 1 : canShow3dbox; 2 : canMoveView; 3 : dispWinIs3d
		//int[] _fill, int[] _strk, 			: window fill and stroke colors
		//int _trajFill, int _trajStrk)			: trajectory fill and stroke colors, if these objects can be drawn in window (used as alt color otherwise)
		//specify windows that cannot be shown simultaneously here
		initXORWins(new int[]{showMPMCudawin, showMPMCPUwin},new int[]{dispMPMCudaWinIDX, dispMPMCPUWinIDX});

		wIdx = dispMPMCudaWinIDX; fIdx= showMPMCudawin;
		setInitDispWinVals(wIdx, _dimOpen, _dimClosed,new boolean[]{false,true,true,true}, new int[]{255,245,255,255},new int[]{0,0,0,255},new int[]{180,180,180,255},new int[]{100,100,100,255}); 		
		dispWinFrames[wIdx] = new MPM_CudaSimWindow(pa, this, wIdx, fIdx);
		
		wIdx = dispMPMCPUWinIDX; fIdx= showMPMCPUwin;
		setInitDispWinVals(wIdx, _dimOpen, _dimClosed,new boolean[]{false,true,true,true}, new int[]{255,245,255,255},new int[]{0,0,0,255},new int[]{180,180,180,255},new int[]{100,100,100,255}); 		
		dispWinFrames[wIdx] = new MPM_CPUSimWindow(pa, this, wIdx, fIdx);
			
		
	}//initVisOnce_Priv
	
	@Override
	protected void initOnce_Indiv() {
		//which objects to initially show
		setVisFlag(showUIMenu, true);					//show input UI menu	
		setVisFlag(showMPMCudawin, true);		
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
	public float[] getUIRectVals(int idx){
		switch(idx){
			case dispMenuIDX 		: {return new float[0];}			//idx 0 is parent menu sidebar
			case dispMPMCudaWinIDX 	: {return dispWinFrames[dispMenuIDX].uiClkCoords;}
			case dispMPMCPUWinIDX	: {return dispWinFrames[dispMenuIDX].uiClkCoords;}
			default :  return dispWinFrames[dispMenuIDX].uiClkCoords;
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
	public int getNumVisFlags() {return numVisFlags;}
	/**
	 * address all flag-setting here, so that if any special cases need to be addressed they can be easily
	 */
	@Override
	protected void setVisFlag_Indiv(int idx, boolean val ){
		switch (idx){
			case showUIMenu 		: { dispWinFrames[dispMenuIDX].dispFlags.setShowWin(val);    break;}											//whether or not to show the main ui window (sidebar)			
			case showMPMCudawin		: {setWinFlagsXOR(dispMPMCudaWinIDX, val); break;}
			case showMPMCPUwin		: {setWinFlagsXOR(dispMPMCPUWinIDX, val); break;}
			default : {break;}
		}
	}//setFlags  
	
	@Override
	public int[] getClr_Custom(int colorVal, int alpha) {	return new int[] {255,255,255,alpha};}

	@Override
	protected void setSmoothing() {
		pa.setSmoothing(0);		
	}

}//class MPM_SimMain
