package MPM_CudaSim;


import base_UI_Objects.GUI_AppManager;
import base_UI_Objects.windowUI.base.myDispWindow;

/**
* MPM Snow Simulation in CUDA 2.0
* @author john turner
*/

public class MPM_SimMain extends GUI_AppManager {
	
	public String prjNmLong = "MPM Simulation CUDA 2.0", prjNmShrt = "MPM_SnowSim_cuda_2.0";
	public String authorString = "John Turner";

	private final int
		showUIMenu = 0,
		showMPMwin = 1;
	public final int numVisFlags = 2;
	
	//idx's in dispWinFrames for each window - 0 is always left side menu window
	private static final int
		dispMPMWinIDX = 1;
	
	//private boolean cyclModCmp;										//comparison every draw of cycleModDraw			
	private final int[] bground = new int[]{244,244,244,255};		//bground color	
	
	
///////////////
//CODE STARTS
///////////////	
//////////////////////////////////////////////// code																		//set array of vector values (sceneFcsVals) based on application
	//needs main to run project - do not modify this code in any way
	public static void main(String[] passedArgs) {		
		MPM_SimMain me = new MPM_SimMain();
	    MPM_SimMain.invokeProcessingMain(me, passedArgs);
	}//main	
	
	/**
	 * whether or not we want to restrict window size on widescreen monitors
	 * 
	 * @return 0 - use monitor size regardless
	 * 			1 - use smaller dim to be determine window 
	 * 			2+ - TBD
	 */
	@Override
	protected int setAppWindowDimRestrictions() {	return 1;}	
	
	//instance-specific setup code
	protected void setup_Indiv() {
		//modify default grid dims to be 1500x1500x1500
		setDesired3DGridDims(1500);

		setBkgrnd();
	}	
	@Override
	public void setBkgrnd(){pa.setRenderBackground(bground[0],bground[1],bground[2],bground[3]);}//setBkgrnd

	@Override
	protected void initMainFlags_Indiv() {
		setMainFlagToShow_debugMode(false);
		setMainFlagToShow_saveAnim(true); 
		setMainFlagToShow_runSim(true);
		setMainFlagToShow_singleStep(true);
		setMainFlagToShow_showRtSideMenu(true);
	}
	@Override
	protected void initAllDispWindows() {
		showInfo = true;
		
		//includes 1 for menu window (never < 1) - always have same # of visFlags as myDispWindows
		int numWins = numVisFlags;		
		//titles and descs, need to be set before sidebar menu is defined
		String[] _winTitles = new String[]{"","Snow Balls"},
				_winDescr = new String[] {"", "Display Two Colliding Snowballs Simulated via MPM"};
		initWins(numWins,_winTitles, _winDescr);
		//call for menu window
		buildInitMenuWin(showUIMenu);
		//menu bar init
		int wIdx = dispMenuIDX,fIdx=showUIMenu;
		dispWinFrames[wIdx] = this.buildSideBarMenu(wIdx, fIdx, new String[]{"Functions 1","Functions 2","Functions 3","Functions 4"}, new int[] {3,4,4,4}, 5, false, false);
				
		//new mySideBarMenu(this, winTitles[wIdx], fIdx, winFillClrs[wIdx], winStrkClrs[wIdx], winRectDimOpen[wIdx], winRectDimClose[wIdx], winDescr[wIdx]);	
		//instanced window dimensions when open and closed - only showing 1 open at a time
		float[] _dimOpen  =  new float[]{menuWidth, 0, pa.getWidth()-menuWidth, pa.getHeight()}, 
				_dimClosed  =  new float[]{menuWidth, 0, hideWinWidth, pa.getHeight()};	
		//setInitDispWinVals : use this to define the values of a display window
		//int _winIDX, 
		//float[] _dimOpen, float[] _dimClosed  : dimensions opened or closed
		//String _ttl, String _desc 			: window title and description
		//boolean[] _dispFlags 					: 
		//   flags controlling display of window :  idxs : 0 : canDrawInWin; 1 : canShow3dbox; 2 : canMoveView; 3 : dispWinIs3d
		//int[] _fill, int[] _strk, 			: window fill and stroke colors
		//int _trajFill, int _trajStrk)			: trajectory fill and stroke colors, if these objects can be drawn in window (used as alt color otherwise)

		wIdx = dispMPMWinIDX; fIdx= showMPMwin;
		setInitDispWinVals(wIdx, _dimOpen, _dimClosed,new boolean[]{false,true,true,true}, new int[]{255,245,255,255},new int[]{0,0,0,255},new int[]{180,180,180,255},new int[]{100,100,100,255}); 		
		dispWinFrames[wIdx] = new MPM_SimWindow(pa, this, winTitles[wIdx], fIdx, winFillClrs[wIdx], winStrkClrs[wIdx], winRectDimOpen[wIdx], winRectDimClose[wIdx], winDescr[wIdx]);
		
		//specify windows that cannot be shown simultaneously here
		initXORWins(new int[]{showMPMwin},new int[]{dispMPMWinIDX});
		
	}//initVisOnce_Priv
	
	@Override
	protected void initOnce_Indiv() {
		//which objects to initially show
		setVisFlag(showUIMenu, true);					//show input UI menu	
		setVisFlag(showMPMwin, true);		
	}
	
	@Override
	protected void initVisProg_Indiv() {}
	@Override
	protected void initProgram_Indiv() {}

	@Override
	public String[] getMouseOverSelBtnNames() {
		// TODO Auto-generated method stub
		return new String[0];
	}
	

	@Override
	protected String getPrjNmLong() {		return prjNmLong;}
	@Override
	protected String getPrjNmShrt() {		return prjNmShrt;}

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
	@Override
	//gives multiplier based on whether shift, alt or cntl (or any combo) is pressed
	public double clickValModMult(){return ((altIsPressed() ? .1 : 1.0) * (shiftIsPressed() ? 10.0 : 1.0));}	
	@Override
	public boolean isClickModUIVal() {
		//TODO change this to manage other key settings for situations where multiple simultaneous key presses are not optimal or conventient
		return altIsPressed() || shiftIsPressed();		
	}
	
	//get the ui rect values of the "master" ui region (another window) -> this is so ui objects of one window can be made, clicked, and shown displaced from those of the parent windwo
	@Override
	public float[] getUIRectVals(int idx){
			//this.pr("In getUIRectVals for idx : " + idx);
		switch(idx){
			case dispMenuIDX 		: {return new float[0];}			//idx 0 is parent menu sidebar
			case dispMPMWinIDX 	: {return dispWinFrames[dispMenuIDX].uiClkCoords;}
			//case disp2ndWinIDX 	: {	return dispWinFrames[dispMenuIDX].uiClkCoords;}
			default :  return dispWinFrames[dispMenuIDX].uiClkCoords;
			}
	}//getUIRectVals

	@Override
	//these tie using the UI buttons to modify the window in with using the boolean tags - PITA but currently necessary
	public void handleShowWin(int btn, int val, boolean callFlags){//display specific windows - multi-select/ always on if sel
		if(!callFlags){//called from setflags - only sets button state in UI to avoid infinite loop
			//setMenuBtnState(mySideBarMenu.btnShowWinIdx,btn, val);
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
	@Override
	//address all flag-setting here, so that if any special cases need to be addressed they can be
	protected void setVisFlag_Indiv(int idx, boolean val ){
		switch (idx){
			case showUIMenu 	    : { dispWinFrames[dispMenuIDX].setFlags(myDispWindow.showIDX,val);    break;}											//whether or not to show the main ui window (sidebar)			
			case showMPMwin			: {setWinFlagsXOR(dispMPMWinIDX, val); break;}
			default : {break;}
		}
	}//setFlags  
	
	@Override
	public int[] getClr_Custom(int colorVal, int alpha) {	return new int[] {255,255,255,alpha};}

	@Override
	protected void setSmoothing() {		pa.setSmoothing(4);		}

}//class MPM_SimMain
