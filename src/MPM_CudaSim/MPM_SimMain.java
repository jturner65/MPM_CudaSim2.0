package MPM_CudaSim;

import java.io.File;
import java.util.*;
import java.util.concurrent.*;


import processing.core.*;
import processing.event.*;
import processing.opengl.*;


/**
 * MPM Snow Simulation in CUDA 2.0
 * @author john turner
 */

public class MPM_SimMain extends PApplet{

	public String prjNmLong = "MPM Simulation CUDA 2.0", prjNmShrt = "MPM_SnowSim_cuda_2.0";
	public String authorString = "John Turner";
	
	//don't use sphere background for this program
	private boolean useSphereBKGnd = false;
	
//	//epsilon value for calculations
	public final float epsValCalc = .00000001f;

//	//how much force is exerted at click location -decays with distance sqrd
	public final float msClickForce = 100000000;
	
	public final int drawnTrajEditWidth = 10; //TODO make config file var?			//width in cntl points of the amount of the drawn trajectory deformed by dragging
	public final float
				//PopUpWinOpenFraction = .40f,				// No Popup in this proj //fraction of screen not covered by popwindow
				wScale = frameRate/5.0f,					//velocity drag scaling	
				trajDragScaleAmt = 100.0f;					//amt of displacement when dragging drawn trajectory to edit
			
	public String msClkStr = "";
	
///////////////
//CODE STARTS
///////////////	
	//////////////////////////////////////////////// code
	public static void main(String[] passedArgs) {
	    String[] appletArgs = new String[] { "MPM_CudaSim.MPM_SimMain" };
	    if (passedArgs != null) {	    	PApplet.main(PApplet.concat(appletArgs, passedArgs));  } else {	    	PApplet.main(appletArgs);	    }
	}
	public void settings(){
		size((int)(displayWidth*.95f), (int)(displayHeight*.92f),P3D);
		noSmooth();
	}		

	public void setup() {
		colorMode(RGB, 255, 255, 255, 255);
		frameRate(frate);
		if(useSphereBKGnd) {
			setBkgndSphere();
		} else {
			setBkgrnd();
		}
		initVisOnce();
		//call this in first draw loop?
		initOnce();
	}// setup
	
	private void setBkgndSphere() {
		sphereDetail(100);
		//TODO move to window to set up specific background for each different "scene" type
		PImage bgrndTex = loadImage("bkgrndTex.jpg");
		bgrndSphere = createShape(SPHERE, 10000);
		bgrndSphere.setTexture(bgrndTex);
		bgrndSphere.rotate(HALF_PI,-1,0,0);
		bgrndSphere.setStroke(false);	
		//TODO move to myDispWindow
		background(bground[0],bground[1],bground[2],bground[3]);		
		shape(bgrndSphere);	
	}
	
	public void setBkgrnd(){
		background(bground[0],bground[1],bground[2],bground[3]);		
	}//setBkgrnd

	private void initOnce() {
		numThreadsAvail = Runtime.getRuntime().availableProcessors();
		pr("# threads Total on this machine : "+ numThreadsAvail);
		th_exec = Executors.newCachedThreadPool();		
		initDispWins();
		setFlags(showUIMenu, true);					//show input UI menu	
		setFlags(showDESwin, true);					//show first window
		initProgram();
		setFlags(finalInitDone, true);
	}//	initOnce
	
	//called multiple times, whenever re-initing
	public void initProgram(){
		initVisProg();				//always first
		drawCount = 0;
	}//initProgram
	
	//get difference between frames and set both glbl times
	private float getModAmtMillis() {
		glblStartSimFrameTime = millis();
		float modAmtMillis = (glblStartSimFrameTime - glblLastSimFrameTime);
		glblLastSimFrameTime = millis();
		return modAmtMillis;
	}
	
	public void draw(){	
		if(!flags[finalInitDone]) {initOnce(); return;}	
		float modAmtMillis = getModAmtMillis();
		//simulation section
		if(flags[runSim] ){
			//run simulation
			drawCount++;									//needed to stop draw update so that pausing sim retains animation positions	
			for(int i =1; i<numDispWins; ++i){if((isShowingWindow(i)) && (dispWinFrames[i].getFlags(myDispWindow.isRunnable))){dispWinFrames[i].simulate(modAmtMillis);}}
			if(flags[singleStep]){setFlags(runSim,false);}
			simCycles++;
		}		//play in current window

		//drawing section
		pushMatrix();pushStyle();
		drawSetup();																//initialize camera, lights and scene orientation and set up eye movement		
		if((curFocusWin == -1) || (dispWinIs3D[curFocusWin])){	//allow for single window to have focus, but display multiple windows	
			//if refreshing screen, this clears screen, sets background
			setBkgrnd();				
			draw3D_solve3D(modAmtMillis);
			c.buildCanvas();			
			if(canShow3DBox[curFocusWin]){drawBoxBnds();}
			if(dispWinFrames[curFocusWin].chkDrawMseRet()){			c.drawMseEdge();	}			
			popStyle();popMatrix(); 
		} else {	//either/or 2d window
			//2d windows paint window box so background is always cleared
			c.buildCanvas();
			c.drawMseEdge();
			popStyle();popMatrix(); 
			for(int i =1; i<numDispWins; ++i){if (isShowingWindow(i) && !(dispWinFrames[i].getFlags(myDispWindow.is3DWin))){dispWinFrames[i].draw2D(modAmtMillis);}}
		}
		drawUI(modAmtMillis);																	//draw UI overlay on top of rendered results			
		if (flags[saveAnim]) {	savePic();}
		updateConsoleStrs();
		surface.setTitle(prjNmLong + " : " + (int)(frameRate) + " fps|cyc curFocusWin : " + curFocusWin);
	}//draw
	
	private void updateConsoleStrs(){
		++drawCount;
		if(drawCount % cnslStrDecay == 0){drawCount = 0;	consoleStrings.poll();}			
	}//updateConsoleStrs
	
	public void draw3D_solve3D(float modAmtMillis){
		//System.out.println("drawSolve");
		pushMatrix();pushStyle();
		for(int i =1; i<numDispWins; ++i){
			if((isShowingWindow(i)) && (dispWinFrames[i].getFlags(myDispWindow.is3DWin))){
				dispWinFrames[i].draw3D(modAmtMillis);
			}
		}
		popStyle();popMatrix();
		//fixed xyz rgb axes for visualisation purposes and to show movement and location in otherwise empty scene
		drawAxes(100,3, new myPoint(-c.viewDimW/2.0f+40,0.0f,0.0f), 200, false); 		
	}//draw3D_solve3D
	
	//if should show problem # i
	public boolean isShowingWindow(int i){return flags[(i+this.showUIMenu)];}//showUIMenu is first flag of window showing flags
	public void drawUI(float modAmtMillis){					
		//for(int i =1; i<numDispWins; ++i){if ( !(dispWinFrames[i].dispFlags[myDispWindow.is3DWin])){dispWinFrames[i].draw(sceneCtrVals[sceneIDX]);}}
		//dispWinFrames[0].draw(sceneCtrVals[sceneIDX]);
		for(int i =1; i<numDispWins; ++i){dispWinFrames[i].drawHeader(modAmtMillis);}
		//menu always idx 0
		dispWinFrames[0].draw2D(modAmtMillis);
		dispWinFrames[0].drawHeader(modAmtMillis);
		drawOnScreenData();				//debug and on-screen data
	}//drawUI	
	
	//handle pressing keys 0-9
	//keyVal is actual value of key (screen character as int)
	//keyPressed is actual key pressed (shift-1 gives keyVal 33 ('!') but keyPressed 49 ('1')) 
	//need to subtract 48 from keyVal or keyPressed to get actual number
	private void handleNumberKeyPress(int keyVal, int keyPressed) {
		//use key if want character 
		if (key == '0') {
			setFlags(showUIMenu,true);
		} 
		
	}//handleNumberKeyPress

//////////////////////////////////////////////////////
/// user interaction
//////////////////////////////////////////////////////	
	//key is key pressed
	//keycode is actual physical key pressed == key if shift/alt/cntl not pressed.,so shift-1 gives key 33 ('!') but keycode 49 ('1')
	public void keyPressed(){
		if(key==CODED) {
			if(!flags[shiftKeyPressed]){setFlags(shiftKeyPressed,(keyCode  == 16));} //16 == KeyEvent.VK_SHIFT
			if(!flags[cntlKeyPressed]){setFlags(cntlKeyPressed,(keyCode  == 17));}//17 == KeyEvent.VK_CONTROL			
			if(!flags[altKeyPressed]){setFlags(altKeyPressed,(keyCode  == 18));}//18 == KeyEvent.VK_ALT
		} else {	
			//handle pressing keys 0-9 (with or without shift,alt, cntl)
			if ((keyCode>=48) && (keyCode <=57)) { handleNumberKeyPress(((int)key),keyCode);}
			else {					//handle all other (non-numeric) keys
				switch (key){
					case ' ' : {setFlags(runSim,!flags[runSim]); break;}							//run sim
					case 'f' : {dispWinFrames[curFocusWin].setInitCamView();break;}					//reset camera
					case 'a' :
					case 'A' : {setFlags(saveAnim,!flags[saveAnim]);break;}						//start/stop saving every frame for making into animation
					case 's' :
					case 'S' : {save(sketchPath() +File.separatorChar+prjNmShrt+"_"+dateStr+File.separatorChar+prjNmShrt+"_img"+timeStr + ".jpg");break;}//save picture of current image			
					default : {	}
				}//switch	
			}
		}
	}
	public void keyReleased(){
		if(key==CODED) {
			if((flags[shiftKeyPressed]) && (keyCode == 16)){endShiftKey();}
			if((flags[cntlKeyPressed]) && (keyCode == 17)){endCntlKey();}
			if((flags[altKeyPressed]) && (keyCode == 18)){endAltKey();}
		}
	}		
	public void endShiftKey(){
		clearFlags(new int []{shiftKeyPressed, modView});
		for(int i =0; i<numDispWins; ++i){dispWinFrames[i].endShiftKey();}
	}
	public void endAltKey(){
		clearFlags(new int []{altKeyPressed});
		for(int i =0; i<numDispWins; ++i){dispWinFrames[i].endAltKey();}			
	}
	public void endCntlKey(){
		clearFlags(new int []{cntlKeyPressed});
		for(int i =0; i<numDispWins; ++i){dispWinFrames[i].endCntlKey();}			
	}
	
	//2d range checking of point
	public boolean ptInRange(double x, double y, double minX, double minY, double maxX, double maxY){return ((x > minX)&&(x <= maxX)&&(y > minY)&&(y <= maxY));}	
	//gives multiplier based on whether shift, alt or cntl (or any combo) is pressed
	public double clickValModMult(){
		return ((flags[altKeyPressed] ? .1 : 1.0) * (flags[cntlKeyPressed] ? 10.0 : 1.0));			
	}
	public void mouseMoved(){for(int i =0; i<numDispWins; ++i){if (dispWinFrames[i].handleMouseMove(mouseX, mouseY)){return;}}}
	public void mousePressed() {
		//verify left button if(mouseButton == LEFT)
		setFlags(mouseClicked, true);
		if(mouseButton == LEFT){			mouseClicked(0);} 
		else if (mouseButton == RIGHT) {	mouseClicked(1);}
		//for(int i =0; i<numDispWins; ++i){	if (dispWinFrames[i].handleMouseClick(mouseX, mouseY,c.getMseLoc(sceneCtrVals[sceneIDX]))){	return;}}
	}// mousepressed	
	
	private void mouseClicked(int mseBtn){ for(int i =0; i<numDispWins; ++i){if (dispWinFrames[i].handleMouseClick(mouseX, mouseY,mseBtn)){return;}}}
	public void mouseDragged(){
		if(mouseButton == LEFT){			mouseDragged(0);}
		else if (mouseButton == RIGHT) {	mouseDragged(1);}
	}//mouseDragged()
	//only for zooming
	public void mouseWheel(MouseEvent event) {
		if (dispWinFrames[curFocusWin].getFlags(myDispWindow.canChgView)) {// (canMoveView[curFocusWin]){	
			float mult = (flags[shiftKeyPressed]) ? 5.0f * mouseWhlSens : mouseWhlSens;
			dispWinFrames[curFocusWin].handleViewChange(true,(mult * event.getCount()),0);
		}
	}
	private void mouseDragged(int mseBtn){
		for(int i =0; i<numDispWins; ++i){if (dispWinFrames[i].handleMouseDrag(mouseX, mouseY, pmouseX, pmouseY,c.getMseDragVec(),mseBtn)) {return;}}		
	}
	
	public void mouseReleased(){
		clearFlags(new int[]{mouseClicked, modView});
		msClkStr = "";
		for(int i =0; i<numDispWins; ++i){dispWinFrames[i].handleMouseRelease();}
		flags[drawing] = false;
		//c.clearMsDepth();
	}//mouseReleased
	
	private void setMenuBtnState(int row, int col, int val) {
		((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiBtnSt[row][col] = val;	
		if (val == 1) {
			outStr2Scr("turning on button row : " + row + "  col " + col);
			((mySideBarMenu)dispWinFrames[dispMenuIDX]).setWaitForProc(row,col);}//if programmatically (not through UI) setting button on, then set wait for proc value true 
	}//setMenuBtnState
	
	//these tie using the UI buttons to modify the window in with using the boolean tags - PITA but currently necessary
	public void handleShowWin(int btn, int val){handleShowWin(btn, val, true);}					//display specific windows - multi-select/ always on if sel
	public void handleShowWin(int btn, int val, boolean callFlags){//display specific windows - multi-select/ always on if sel
	if(!callFlags){//called from setflags - only sets button state in UI to avoid infinite loop
		//setMenuBtnState(mySideBarMenu.btnShowWinIdx,btn, val);
	} else {//called from clicking on buttons in UI
		//val is btn state before transition 
		boolean bVal = (val == 1?  false : true);
		switch(btn){
			case 0 : {setFlags(showDESwin, bVal);break;}
			//case 1 : {setFlags(show2ndWinIDX, bVal);break;}
			}
		}
	}//handleShowWin
	
	//call menu from instance of dispwindow to update primary custom function button names with window-relevant entries
	public void setMenuFuncBtnNames(String[] btnNames) {
		((mySideBarMenu)dispWinFrames[dispMenuIDX]).setBtnNames(((mySideBarMenu)dispWinFrames[dispMenuIDX]).btnAuxFuncIdx,btnNames);
	}
	
	//isSlowProc means original calling process lasted longer than mouse click release and so button state should be forced to be off
	private void clearBtnState(int row, int col, boolean isSlowProc) {
		((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiBtnWaitForProc[row][col] = false;
		if(isSlowProc) {((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiBtnSt[row][col] = 0;}		
	}
	
	//turn off specific function button that might have been kept on during processing - btn must be in range of size of guiBtnSt[mySideBarMenu.btnAuxFuncIdx]
	//isSlowProc means function this was waiting on is a slow process and escaped the click release in the window (i.e. if isSlowProc then we must force button to be off)
	public void clearFuncBtnSt(int btn, boolean isSlowProc) {clearBtnState(mySideBarMenu.btnAuxFuncIdx,btn, isSlowProc);}
	//process to delete an existing component
	public void handleFuncSelCmp(int btn, int val){handleFuncSelCmp(btn, val, true);}					//display specific windows - multi-select/ always on if sel
	public void handleFuncSelCmp(int btn, int val, boolean callFlags){
		if(!callFlags){
			setMenuBtnState(mySideBarMenu.btnAuxFuncIdx,btn, val);
		} else {
			dispWinFrames[curFocusWin].clickFunction(btn) ;
		}
	}//handleAddDelSelCmp	
	
	//call menu from instance of dispwindow to update primary debug button names with window-relevant entries
	public void setMenuDbgBtnNames(String[] btnNames) {
		((mySideBarMenu)dispWinFrames[dispMenuIDX]).setBtnNames(((mySideBarMenu)dispWinFrames[dispMenuIDX]).btnDBGSelCmpIdx,btnNames);
	}
	//turn off specific debug button that might have been kept on during processing - btn must be in range of size of guiBtnSt[mySideBarMenu.btnDBGSelCmpIdx]
	//isSlowProc means function this was waiting on is a slow process and escaped the click release in the window (i.e. if isSlowProc then we must force button to be off)
	public void clearDBGBtnSt(int btn, boolean isSlowProc)  {clearBtnState(mySideBarMenu.btnDBGSelCmpIdx,btn, isSlowProc);}
	//process to delete an existing component
	public void handleDBGSelCmp(int btn, int val){handleDBGSelCmp(btn, val, true);}					//display specific windows - multi-select/ always on if sel
	public void handleDBGSelCmp(int btn, int val, boolean callFlags){
		if(!callFlags){
			setMenuBtnState(mySideBarMenu.btnDBGSelCmpIdx,btn, val);
		} else {
			dispWinFrames[curFocusWin].clickDebug(btn) ;
		}
	}//handleAddDelSelCmp	
	
	//process to handle file io	- TODO	
	public void handleFileCmd(int btn, int val){handleFileCmd(btn, val, true);}					//display specific windows - multi-select/ always on if sel
	public void handleFileCmd(int btn, int val, boolean callFlags){//{"Load","Save"},							//load an existing score, save an existing score - momentary	
		if(!callFlags){
			setMenuBtnState(mySideBarMenu.btnFileCmdIdx,btn, val);
		} else {
			switch(btn){
				case 0 : {selectInput("Select a txt file to load parameters from : ", "loadFromFile");break;}
				case 1 : {selectOutput("Select a txt file to save parameters to : ", "saveToFile");break;}
			}
			((mySideBarMenu)dispWinFrames[dispMenuIDX]).hndlMouseRelIndiv();
		}
	}//handleFileCmd
	
	//load strings of data from a text File named file
	public void loadFromFile(File file){
		if (file == null) {
		    outStr2Scr("Load was cancelled.");
		    return;
		} 
		String[] res = loadStrings(file.getAbsolutePath());
		//stop any simulations while file is loaded
		//load - iterate through for each window
		int[] stIdx = {0};//start index for a particular window - make an array so it can be passed by ref and changed by windows
		for(int i =0; i<numDispWins; ++i){
			while(!res[stIdx[0]].contains(dispWinFrames[i].name)){++stIdx[0];}			
			dispWinFrames[i].hndlFileLoad(res,stIdx);
		}//accumulate array of params to save
		//resume simulations		
	}//loadFromFile
	
	public void saveToFile(File file){
		if (file == null) {
		    outStr2Scr("Save was cancelled.");
		    return;
		} 
		ArrayList<String> res = new ArrayList<String>();
		//save - iterate through for each window
		for(int i =0; i<numDispWins; ++i){
			res.addAll(dispWinFrames[i].hndlFileSave());	
		}//accumulate array of params to save
		saveStrings(file.getAbsolutePath(), res.toArray(new String[0]));  
	}//saveToFile

	public void initDispWins(){
		//float popUpWinHeight = PopUpWinOpenFraction * height;		//how high is the InstEdit window when shown
		//instanced window dimensions when open and closed - only showing 1 open at a time
		winRectDimOpen[dispDESWinIDX] =  new float[]{menuWidth, 0,width-menuWidth,height};			
		//winRectDimOpen[disp2ndWinIDX] =   new float[]{menuWidth, 0,width-menuWidth,height};				
		//hidden
		winRectDimClose[dispDESWinIDX] =  new float[]{menuWidth, 0, hideWinWidth, height};				
		//winRectDimClose[disp2ndWinIDX] =  new float[]{menuWidth, 0, hideWinWidth, height};				
		
		winTrajFillClrs = new int []{gui_Black,gui_LightGray,gui_LightGray};		//set to color constants for each window
		winTrajStrkClrs = new int []{gui_Black,gui_DarkGray,gui_DarkGray};				//set to color constants for each window			
		//			//display window initialization	
		int wIdx = dispDESWinIDX , fIdx = showDESwin;
		dispWinFrames[wIdx] = new MPM_SimWindow(this, winTitles[wIdx], fIdx, winFillClrs[wIdx], winStrkClrs[wIdx], winRectDimOpen[wIdx], winRectDimClose[wIdx], winDescr[wIdx],canDrawInWin[wIdx]);
		//wIdx = disp2ndWinIDX; fIdx = show2ndWinIDX;
		//dispWinFrames[wIdx] = new altWindow(this, winTitles[wIdx], fIdx,winFillClrs[wIdx], winStrkClrs[wIdx], winRectDimOpen[wIdx], winRectDimClose[wIdx], winDescr[wIdx],canDrawInWin[wIdx]);

		for(int i =0; i < numDispWins; ++i){
			//dispWinFrames[i].initDrwnTrajs();		//drawn trajectories not used in this application (so far)
			//outStr2Scr("i:"+i);
			//init focus and scene center variables for each window

			int scIdx = dispWinIs3D[i] ? 1 : 0;
			dispWinFrames[i].finalInit(dispWinIs3D[i], canMoveView[i], sceneCtrValsBase[scIdx], sceneFcsValsBase[scIdx]);
			dispWinFrames[i].setTrajColors(winTrajFillClrs[i], winTrajStrkClrs[i]);
		}	
		//set initial state to be true - show info window
		setFlags(showRtSideMenu, true);

	}//initDispWins
	
	//get the ui rect values of the "master" ui region (another window) -> this is so ui objects of one window can be made, clicked, and shown displaced from those of the parent windwo
	public float[] getUIRectVals(int idx){
			//this.pr("In getUIRectVals for idx : " + idx);
		switch(idx){
			case dispMenuIDX 		: {return new float[0];}			//idx 0 is parent menu sidebar
			case dispDESWinIDX 	: {return dispWinFrames[dispMenuIDX].uiClkCoords;}
			//case disp2ndWinIDX 	: {	return dispWinFrames[dispMenuIDX].uiClkCoords;}
			default :  return dispWinFrames[dispMenuIDX].uiClkCoords;
			}
	}//getUIRectVals
	
	//find mouse "force" exerted upon a particular location - distance from mouse to passed location
	public myVectorf mouseForceAtLoc(myPointf _loc, boolean attractMode){
		myPointf mouseFrcLoc = c.getTransMseLoc(new myPointf(gridDimX/2.0f, gridDimY/2.0f,gridDimZ/2.0f));// new myPointf(c.dfCtr.x+gridDimX/2.0f,c.dfCtr.y+gridDimY/2.0f,c.dfCtr.z+gridDimZ/2.0f);// new myVector(lstClkX,0,lstClkY);//translate click location to where the space where the boids are	
		myVectorf resFrc = new myVectorf(_loc, mouseFrcLoc);		
		float sqDist = resFrc.sqMagn;
		if(sqDist<epsValCalc){sqDist=epsValCalc;}
		float mag = (attractMode? 1 : -1) * msClickForce / sqDist;
		resFrc._scale(mag);
		return resFrc;	
	}//mouseForceAtLoc

//////////////////////////////////////////
/// graphics and base functionality utilities and variables
//////////////////////////////////////////
	
	public int glblStartSimFrameTime,			//begin of draw
	glblLastSimFrameTime,					//begin of last draw
	glblStartProgTime;					//start of program
	
	//size of printed text (default is 12)
	public static final int txtSz = 11;
	//mouse wheel sensitivity
	public static final float mouseWhlSens = 1.0f;
	
	//display-related size variables
	public final int grid2D_X=800, grid2D_Y=800;	
	public final int gridDimX = 1500, gridDimY = 1500, gridDimZ = 1500;				//dimensions of 3d region
	
	public final myVectorf gridHalfDim = new myVectorf(gridDimX*.5f,gridDimY*.5f,gridDimZ*.5f );
	//boundary regions for enclosing cube - given as min and difference of min and max
	public float[][] cubeBnds = new float[][]{//idx 0 is min, 1 is diffs
	new float[]{-gridDimX/2.0f,-gridDimY/2.0f,-gridDimZ/2.0f},//mins
	new float[]{gridDimX,gridDimY,gridDimZ}};					//diffs
	
	public final float fsqrt2 = (float)(Math.sqrt(2.0));
	
	private final int cnslStrDecay = 3;							//how long a message should last before it is popped from the console strings deque
	
	public int scrWidth, scrHeight;								//set to be applet.width and applet.height unless otherwise specified below
	public final int scrWidthMod = 200, 
			scrHeightMod = 0;
	public final float frate = 120.0f;								//frame rate - # of playback updates per second
	
	//2D, 3D
	private myVector[] sceneFcsValsBase = new myVector[]{						//set these values to be different targets of focus
			new myVector(-grid2D_X/2,-grid2D_Y/1.75f,0),
			new myVector(0,0,0)
	};
	//2D, 3D
	private myPoint[] sceneCtrValsBase = new myPoint[]{				//set these values to be different display center translations -
		new myPoint(0,0,0),										// to be used to calculate mouse offset in world for pick
		new myPoint(-gridDimX/2.0,-gridDimY/2.0,-gridDimZ/2.0)
	};
	
	//static variables - put obj constructor counters here
	public static int GUIObjID = 0;										//counter variable for gui objs
	
	//visualization variables
	// boolean flags used to control various elements of the program 
	public boolean[] flags;
	//dev/debug flags
	public final int debugMode 			= 0;			//whether we are in debug mode or not	
	public final int finalInitDone		= 1;			//used only to call final init in first draw loop, to avoid stupid timeout error processing 3.x's setup introduced
	public final int saveAnim 			= 2;			//whether we are saving or not
	//interface flags	
	public final int shiftKeyPressed 	= 3;			//shift pressed
	public final int altKeyPressed  	= 4;			//alt pressed
	public final int cntlKeyPressed  	= 5;			//cntrl pressed
	public final int mouseClicked 		= 6;			//mouse left button is held down	
	public final int drawing			= 7; 			//currently drawing
	public final int modView	 		= 8;			//shift+mouse click+mouse move being used to modify the view
	//simulation
	public final int runSim				= 9;			//run simulation
	public final int singleStep			= 10;			//run single sim step
	public final int showRtSideMenu		= 11;			//display the right side info menu for the current window, if it supports that display
	public final int flipDrawnTraj  	= 12;			//whether or not to flip the direction of the drawn trajectory
	//window control
	public final int showUIMenu 		= 13;			//whether or not to show sidebar menu
	public final int showDESwin			= 14;			//whether to show 1st window
	//public final int show2ndWinIDX		= 14;			//whether to show 2nd window
	
	public final int numFlags = 15;
	
//	public boolean[] debugVisFlags;
//	public final int showParticleVelArrows 	= 0;			//plot velocity arrows for each particle	
//	public final int showGrid				= 1;			//plot the computational grid
//	public final int showGridVelArrows 		= 2;			//plot velocity arrows for each gridNode
//	public final int showGridAccelArrows 	= 3;			//plot acceleration arrows for each gridNode
//	public final int showGridMass  			= 4;			//plot variable sized spheres proportional to gridnode mass
//	public final int showParticleNeighbors  = 5;			//show the grid nodes influenced by each particle
	
	public final int numDebugVisFlags = 6;
	//flags to actually display in menu as clickable text labels - order does matter
	public List<Integer> flagsToShow = Arrays.asList( 
		debugMode, 			
		saveAnim,
		runSim,
		singleStep,
		showRtSideMenu
		);
	
	public final int numFlagsToShow = flagsToShow.size();
	
	public List<Integer> stateFlagsToShow = Arrays.asList( 
		shiftKeyPressed,			//shift pressed
		altKeyPressed,				//alt pressed
		cntlKeyPressed,				//cntrl pressed
		mouseClicked,				//mouse left button is held down	
		drawing, 					//currently drawing
		modView	 					//shift+mouse click+mouse move being used to modify the view					
			);
	public final int numStFlagsToShow = stateFlagsToShow.size();	
	
	///////////////////////////////////
	// window related variables	
	public String[] winTitles = new String[]{"","Snowball on Ball","Alternate Window"},
			winDescr = new String[] {"", "Snow falls on sphere collider","Alternate Window"};
	//individual display/HUD windows for gui/user interaction
	public myDispWindow[] dispWinFrames;
	//idx's in dispWinFrames for each window
	public static final int dispMenuIDX = 0,
							dispDESWinIDX = 1;
							//disp2ndWinIDX = 2;
	
	public static final int numDispWins = 2;	
			
	public int curFocusWin;				//which myDispWindow currently has focus 
	
	//whether or not the display windows will accept a drawn trajectory
	public boolean[] canDrawInWin = new boolean[]{false,false,false};		
	public boolean[] canShow3DBox = new boolean[]{false,true,true};		
	public boolean[] canMoveView = new boolean[]{false,true,true};		
	public static final boolean[] dispWinIs3D = new boolean[]{false,true,true};
	
	public static final int[][] winFillClrs = new int[][]{          
		new int[]{255,255,255,255},                                 	// dispMenuIDX = 0,
		new int[]{255,255,255,255},                                        	// dispJTWinIDX = 1;
		new int[]{255,255,255,255}                                        	// dispYYResIDX = 2
	};
	public static final int[][] winStrkClrs = new int[][]{
		new int[]{0,0,0,255},                                    		//dispMenuIDX = 0,
		new int[]{0,0,0,255},                               		//dispJTWinIDX = 1
		new int[]{255,255,255,255}                               		//dispYYResIDX = 2
	};
	
	public static int[] winTrajFillClrs = new int []{0,0,0};		//set to color constants for each window
	public static int[] winTrajStrkClrs = new int []{0,0,0};		//set to color constants for each window	
	
	//unblocked window dimensions - location and dim of window if window is one\
	public float[][] winRectDimOpen;// = new float[][]{new float[]{0,0,0,0},new float[]{0,0,0,0},new float[]{0,0,0,0},new float[]{0,0,0,0}};
	//window dimensions if closed -location and dim of all windows if this window is closed
	public float[][] winRectDimClose;// = new float[][]{new float[]{0,0,0,0},new float[]{0,0,0,0},new float[]{0,0,0,0},new float[]{0,0,0,0}};
	
	public boolean showInfo;										//whether or not to show start up instructions for code		
	//private myVector focusTar;										//target of focus - used in translate to set where the camera is looking - 
																//set array of vector values (sceneFcsVals) based on application
	//private boolean cyclModCmp;										//comparison every draw of cycleModDraw			
	public final int[] bground = new int[]{244,244,255,255};		//bground color
	private PShape bgrndSphere;										//giant sphere encapsulating entire scene
			
	public myPoint mseCurLoc2D;
	//how many frames to wait to actually refresh/draw
	//public int cycleModDraw = 1;
	public final int maxCycModDraw = 20;	//max val for cyc mod draw		
	
	// path and filename to save pictures for animation
	public String screenShotPath;
	public int animCounter;	
	public final int scrMsgTime = 50;									//5 seconds to delay a message 60 fps (used against draw count)
	public ArrayDeque<String> consoleStrings;							//data being printed to console - show on screen
	
	public int drawCount,simCycles;												// counter for draw cycles		
	public float menuWidth,menuWidthMult = .15f, hideWinWidth, hideWinWidthMult = .03f, //side menu is 15% of screen grid2D_X, 
			hidWinHeight, hideWinHeightMult = .05f;								//leave a bound at bottom of screen for hidden pop-up window				
	
	public ArrayList<String> DebugInfoAra;										//enable drawing dbug info onto screen
	public String debugInfoString;
	
	//animation control variables	
	//public float animCntr = 0, animModMult = 1.0f;
	public final float maxAnimCntr = PI*1000.0f, baseAnimSpd = 1.0f;	
	public float msSclX, msSclY;											//scaling factors for mouse movement		
	my3DCanvas c;															//3d interaction stuff and mouse tracking	
	public float[] camVals;		
	public String dateStr, timeStr;											//used to build directory and file names for screencaps
	
	public PGraphicsOpenGL pg; 
	public PGL pgl;
	//public GL2 gl;
	
	public double eps = .000000001, msClkEps = 40;							//calc epsilon, distance within which to check if clicked from a point
	public float feps = .000001f;
	public float SQRT2 = sqrt(2.0f);
	
	public int[] rgbClrs = new int[]{gui_Red,gui_Green,gui_Blue};
	//3dbox stuff
	public myVector[] boxNorms = new myVector[] {new myVector(1,0,0),new myVector(-1,0,0),new myVector(0,1,0),new myVector(0,-1,0),new myVector(0,0,1),new myVector(0,0,-1)};//normals to 3 d bounding boxes
	private final float hGDimX = gridDimX/2.0f, hGDimY = gridDimY/2.0f, hGDimZ = gridDimZ/2.0f;
	private final float tGDimX = gridDimX*10, tGDimY = gridDimY*10, tGDimZ = gridDimZ*20;
	public myPoint[][] boxWallPts = new myPoint[][] {//pts to check if intersection with 3D bounding box happens
			new myPoint[] {new myPoint(hGDimX,tGDimY,tGDimZ), new myPoint(hGDimX,-tGDimY,tGDimZ), new myPoint(hGDimX,tGDimY,-tGDimZ)  },
			new myPoint[] {new myPoint(-hGDimX,tGDimY,tGDimZ), new myPoint(-hGDimX,-tGDimY,tGDimZ), new myPoint(-hGDimX,tGDimY,-tGDimZ) },
			new myPoint[] {new myPoint(tGDimX,hGDimY,tGDimZ), new myPoint(-tGDimX,hGDimY,tGDimZ), new myPoint(tGDimX,hGDimY,-tGDimZ) },
			new myPoint[] {new myPoint(tGDimX,-hGDimY,tGDimZ),new myPoint(-tGDimX,-hGDimY,tGDimZ),new myPoint(tGDimX,-hGDimY,-tGDimZ) },
			new myPoint[] {new myPoint(tGDimX,tGDimY,hGDimZ), new myPoint(-tGDimX,tGDimY,hGDimZ), new myPoint(tGDimX,-tGDimY,hGDimZ)  },
			new myPoint[] {new myPoint(tGDimX,tGDimY,-hGDimZ),new myPoint(-tGDimX,tGDimY,-hGDimZ),new myPoint(tGDimX,-tGDimY,-hGDimZ)  }
	};
	//for multithreading 
	public ExecutorService th_exec;
	public int numThreadsAvail;	
	
	public Calendar now;
	
	///////////////////////////////////
	/// generic graphics functions and classes
	///////////////////////////////////
	//1 time initialization of things that won't change
	public void initVisOnce(){	
		//date and time of program launch
		dateStr = "_"+day() + "-"+ month()+ "-"+year();
		timeStr = "_"+hour()+"-"+minute()+"-"+second();
		now = Calendar.getInstance();
		scrWidth = width + scrWidthMod;
		scrHeight = height + scrHeightMod;		//set to be applet.width and applet.height unless otherwise specified below
		
		msSclX = PI/width;
		msSclY = PI/height;

		consoleStrings = new ArrayDeque<String>();				//data being printed to console		
		menuWidth = width * menuWidthMult;						//grid2D_X of menu region	
		hideWinWidth = width * hideWinWidthMult;				//dims for hidden windows
		hidWinHeight = height * hideWinHeightMult;
		c = new my3DCanvas(this);			
		winRectDimOpen = new float[numDispWins][];
		winRectDimClose = new float[numDispWins][];
		winRectDimOpen[0] =  new float[]{0,0, menuWidth, height};
		winRectDimClose[0] =  new float[]{0,0, hideWinWidth, height};
		
		strokeCap(SQUARE);//makes the ends of stroke lines squared off
		
		//display window initialization
		dispWinFrames = new myDispWindow[numDispWins];		
		//menu bar init
		dispWinFrames[dispMenuIDX] = new mySideBarMenu(this, "UI Window", showUIMenu,  winFillClrs[dispMenuIDX], winStrkClrs[dispMenuIDX], winRectDimOpen[dispMenuIDX],winRectDimClose[dispMenuIDX], "User Controls",canDrawInWin[dispMenuIDX]);			
		
		colorMode(RGB, 255, 255, 255, 255);
		mseCurLoc2D = new myPoint(0,0,0);	
		initBoolFlags();
		//initDebugVisFlags();
		//camVals = new float[]{width/2.0f, height/2.0f, (height/2.0f) / tan(PI/6.0f), width/2.0f, height/2.0f, 0, 0, 1, 0};
		camVals = new float[]{0, 0, (height/2.0f) / tan(PI/6.0f), 0, 0, 0, 0,1,0};
		showInfo = true;
		textSize(txtSz);
		outStr2Scr("Current sketchPath " + sketchPath());
		textureMode(NORMAL);			
		rectMode(CORNER);	
		sphereDetail(4);
		simCycles = 0;
		screenShotPath = sketchPath() +File.separatorChar+prjNmShrt+"_" + (int) random(1000)+File.separatorChar;		
		glblStartProgTime = millis();
		glblStartSimFrameTime = glblStartProgTime;
		glblLastSimFrameTime =  glblStartProgTime;	
	}//initVisOnce
	
		//init boolean state machine flags for program
	public void initBoolFlags(){
		flags = new boolean[numFlags];
		for (int i = 0; i < numFlags; ++i) { flags[i] = false;}	
		((mySideBarMenu)dispWinFrames[dispMenuIDX]).initPFlagColors();			//init sidebar window flags
	}	

	
	//address all flag-setting here, so that if any special cases need to be addressed they can be
	public void setFlags(int idx, boolean val ){
		flags[idx] = val;
		switch (idx){
			case debugMode 			: { break;}//anything special for debugMode 	
			case finalInitDone		: {//flag to handle long setup - processing seems to time out if setup takes too long, so this will continue setup in the first draw loop
				break;}
			case saveAnim 			: { break;}//anything special for saveAnim 			
			case altKeyPressed 		: { break;}//anything special for altKeyPressed 	
			case shiftKeyPressed 	: { break;}//anything special for shiftKeyPressed 	
			case mouseClicked 		: { break;}//anything special for mouseClicked 		
			case modView	 		: { break;}//anything special for modView	 	
			case drawing			: { break;}
			case runSim				: {
				//send out command to stop
				if(val==false) {for(int i=1; i<dispWinFrames.length;++i){if(isShowingWindow(i)) {dispWinFrames[i].stopMe();}}}
				break;}
			case singleStep			: {break;}////anything special for single step	
			case showRtSideMenu		: {	for(int i =1; i<dispWinFrames.length;++i){dispWinFrames[i].setRtSideInfoWinSt(val);}break;}	//set value for every window - to show or not to show info window
			//case flipDrawnTraj		: { dispWinFrames[dispPianoRollIDX].rebuildDrawnTraj();break;}						//whether or not to flip the drawn trajectory, width-wise
			case flipDrawnTraj		: { for(int i =1; i<dispWinFrames.length;++i){dispWinFrames[i].rebuildAllDrawnTrajs();}break;}						//whether or not to flip the drawn melody trajectory, width-wise
			case showUIMenu 	    : { dispWinFrames[dispMenuIDX].setFlags(myDispWindow.showIDX,val);    break;}											//whether or not to show the main ui window (sidebar)			

			case showDESwin			: {setWinFlagsXOR(dispDESWinIDX, val); break;}
			//case show2ndWinIDX		: {setWinFlagsXOR(disp2ndWinIDX, val); break;}
			//following for popup window
			//case showPopUpWinUI 		: {dispWinFrames[dispYYWinIDX].setShow(val);handleShowWin(dispYYWinIDX-1 ,(val ? 1 : 0),false); setWinsHeight(); break;}	//show InstEdit window
			//case useDrawnVels 		: {for(int i =1; i<dispWinFrames.length;++i){dispWinFrames[i].rebuildAllDrawnTrajs();}break;}
			default : {break;}
		}
	}//setFlags  	

	//specify mutually exclusive flags here
	public int[] winFlagsXOR = new int[]{showDESwin};//,show2ndWinIDX};//showSequence,showSphereUI};
	//specify windows that cannot be shown simultaneously here
	public int[] winDispIdxXOR = new int[]{dispDESWinIDX};//, disp2ndWinIDX};//dispPianoRollIDX,dispSphereUIIDX};
	public void setWinFlagsXOR(int idx, boolean val){
		//outStr2Scr("SetWinFlagsXOR : idx " + idx + " val : " + val);
		if(val){//turning one on
			//turn off not shown, turn on shown				
			for(int i =0;i<winDispIdxXOR.length;++i){//check windows that should be mutually exclusive during display
				if(winDispIdxXOR[i]!= idx){dispWinFrames[winDispIdxXOR[i]].setShow(false);handleShowWin(i ,0,false); flags[winFlagsXOR[i]] = false;}
				else {
					dispWinFrames[idx].setShow(true);
					handleShowWin(i ,1,false); 
					flags[winFlagsXOR[i]] = true;
					curFocusWin = winDispIdxXOR[i];
					//setCamView();	//camera now handled by individual windows
					dispWinFrames[idx].setInitCamView();
				}
			}
		} else {//if turning off a window - need a default uncloseable window - for now just turn on next window
			//idx is dispXXXIDX idx of allowable windows (1+ since idx 0 is sidebar menu), so use idx-1 for mod function
			//add 1 to (idx-1) to get next window index, modulo for range adherence, and then add 1 to move back to 1+ from 0+ result from mod		
			//setWinFlagsXOR((((idx-1) + 1) % winFlagsXOR.length)+1, true);
			setWinFlagsXOR((idx % winFlagsXOR.length)+1, true);
		}			
	}//setWinFlagsXOR
	
	//set flags appropriately when only 1 can be true 
	public void setFlagsXOR(int tIdx, int[] fIdx){for(int i =0;i<fIdx.length;++i){if(tIdx != fIdx[i]){flags[fIdx[i]] =false;}}}				
	public void clearFlags(int[] idxs){		for(int idx : idxs){flags[idx]=false;}	}
		//called every time re-initialized
	public void initVisProg(){	drawCount = 0;		debugInfoString = "";		reInitInfoStr();}
	public void reInitInfoStr(){		DebugInfoAra = new ArrayList<String>();		DebugInfoAra.add("");	}	
	public int addInfoStr(String str){return addInfoStr(DebugInfoAra.size(), str);}
	public int addInfoStr(int idx, String str){	
		int lstIdx = DebugInfoAra.size();
		if(idx >= lstIdx){		for(int i = lstIdx; i <= idx; ++i){	DebugInfoAra.add(i,"");	}}
		setInfoStr(idx,str);	return idx;
	}
	public void setInfoStr(int idx, String str){DebugInfoAra.set(idx,str);	}
	
	public void drawInfoStr(float sc, int[] fillClr){//draw text on main part of screen
		pushMatrix();		pushStyle();
			fill(fillClr[0],fillClr[1],fillClr[2],fillClr[3]);
			translate((menuWidth),0);
			scale(sc,sc);
			for(int i = 0; i < DebugInfoAra.size(); ++i){		text((flags[debugMode]?(i<10?"0":"")+i+":     " : "") +"     "+DebugInfoAra.get(i)+"\n\n",0,(10+(12*i)));	}
		popStyle();	popMatrix();
	}		
	//vector and point functions to be compatible with earlier code from jarek's class or previous projects	
	//draw bounding box for 3d
	public void drawBoxBnds(){
		pushMatrix();	pushStyle();
		strokeWeight(3f);
		noFill();
		setColorValStroke(gui_TransGray);
		
		box(gridDimX,gridDimY,gridDimZ);
		popStyle();	popMatrix();
	}	
	//drawsInitial setup for each draw
	public void drawSetup(){
		perspective(PI/3.0f, (1.0f*width)/(1.0f*height), .5f, camVals[2]*100.0f);
	    turnOnLights();
	    dispWinFrames[curFocusWin].drawSetupWin(camVals);
	}//drawSetup	
	//turn on lights for this sketch
	public void turnOnLights(){
	    lights(); 		
	}
	//used to handle camera location/motion - handled by individual windows because each maintains their own camera
	public void setCamOrient(){dispWinFrames[curFocusWin].setCamOrient();		}//sets the rx, ry, pi/2 orientation of the camera eye	
	//used to draw text on screen without changing mode - reverses camera orientation setting - handled by individual windows because each maintains their own camera
	public void unSetCamOrient(){dispWinFrames[curFocusWin].unSetCamOrient(); }//reverses the rx,ry,pi/2 orientation of the camera eye - paints on screen and is unaffected by camera movement
	//public void unSetCamOrient(){rotateY(-ry);   rotateX(-rx); }//reverses the rx,ry,pi/2 orientation of the camera eye - paints on screen and is unaffected by camera movement
	public void drawAxes(double len, float stW, myPoint ctr, int alpha, boolean centered){//axes using current global orientation
		pushMatrix();pushStyle();
			strokeWeight(stW);
			stroke(255,0,0,alpha);
			if(centered){
				double off = len*.5f;
				line(ctr.x-off,ctr.y,ctr.z,ctr.x+off,ctr.y,ctr.z);stroke(0,255,0,alpha);line(ctr.x,ctr.y-off,ctr.z,ctr.x,ctr.y+off,ctr.z);stroke(0,0,255,alpha);line(ctr.x,ctr.y,ctr.z-off,ctr.x,ctr.y,ctr.z+off);} 
			else {		line(ctr.x,ctr.y,ctr.z,ctr.x+len,ctr.y,ctr.z);stroke(0,255,0,alpha);line(ctr.x,ctr.y,ctr.z,ctr.x,ctr.y+len,ctr.z);stroke(0,0,255,alpha);line(ctr.x,ctr.y,ctr.z,ctr.x,ctr.y,ctr.z+len);}
		popStyle();	popMatrix();	
	}//	drawAxes
	public void drawAxes(double len, double stW, myPoint ctr, myVectorf[] _axis, int alpha){
		pushMatrix();pushStyle();
			strokeWeight((float)stW);stroke(255,0,0,alpha);line(ctr.x,ctr.y,ctr.z,ctr.x+(_axis[0].x)*len,ctr.y+(_axis[0].y)*len,ctr.z+(_axis[0].z)*len);stroke(0,255,0,alpha);line(ctr.x,ctr.y,ctr.z,ctr.x+(_axis[1].x)*len,ctr.y+(_axis[1].y)*len,ctr.z+(_axis[1].z)*len);	stroke(0,0,255,alpha);	line(ctr.x,ctr.y,ctr.z,ctr.x+(_axis[2].x)*len,ctr.y+(_axis[2].y)*len,ctr.z+(_axis[2].z)*len);
		popStyle();	popMatrix();	
	}//	drawAxes

	public void drawAxes(double len, float stW, myPoint ctr, myVector[] _axis, int alpha, boolean drawVerts){//RGB -> XYZ axes
		pushMatrix();pushStyle();
		if(drawVerts){
			show(ctr,3,gui_Black,gui_Black, false);
			for(int i=0;i<_axis.length;++i){show(myPoint._add(ctr, myVector._mult(_axis[i],len)),3,rgbClrs[i],rgbClrs[i], false);}
		}
		strokeWeight(stW);
		for(int i =0; i<3;++i){	setColorValStroke(rgbClrs[i]);	showVec(ctr,len, _axis[i]);	}
		popStyle();	popMatrix();	
	}//	drawAxes
	public void drawAxes(double len, float stW, myPoint ctr, myVector[] _axis, int[] clr, boolean drawVerts){//all axes same color
		pushMatrix();pushStyle();
			if(drawVerts){
				show(ctr,2,gui_Black,gui_Black, false);
				for(int i=0;i<_axis.length;++i){show(myPoint._add(ctr, myVector._mult(_axis[i],len)),2,rgbClrs[i],rgbClrs[i], false);}
			}
			strokeWeight(stW);stroke(clr[0],clr[1],clr[2],clr[3]);
			for(int i =0; i<3;++i){	showVec(ctr,len, _axis[i]);	}
		popStyle();	popMatrix();	
	}//	drawAxes
	
	public void drawText(String str, double x, double y, double z, int clr){
		int[] c = getClr(clr);
		pushMatrix();	pushStyle();
			fill(c[0],c[1],c[2],c[3]);
			unSetCamOrient();
			translate((float)x,(float)y,(float)z);
			text(str,0,0,0);		
		popStyle();	popMatrix();	
	}//drawText	
	
	//save screenshot
	public void savePic(){		save(screenShotPath + prjNmShrt + ((animCounter < 10) ? "000" : ((animCounter < 100) ? "00" : ((animCounter < 1000) ? "0" : ""))) + animCounter + ".jpg");		animCounter++;		}
	public void line(double x1, double y1, double z1, double x2, double y2, double z2){line((float)x1,(float)y1,(float)z1,(float)x2,(float)y2,(float)z2 );}
	public void line(myPoint p1, myPoint p2){line((float)p1.x,(float)p1.y,(float)p1.z,(float)p2.x,(float)p2.y,(float)p2.z);}
	public void line(myPointf p1, myPointf p2){line(p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);}
	public void line(myPointf a, myPointf b, int stClr, int endClr){
		beginShape();
		this.strokeWeight(1.0f);
		this.setColorValStroke(stClr);
		this.vertex((float)a.x,(float)a.y,(float)a.z);
		this.setColorValStroke(endClr);
		this.vertex((float)b.x,(float)b.y,(float)b.z);
		endShape();
	}
	
	public void drawOnScreenData(){
		if(flags[debugMode]){
			pushMatrix();pushStyle();			
			reInitInfoStr();
			addInfoStr(0,"mse loc on screen : " + new myPoint(mouseX, mouseY,0) + " mse loc in world :"+c.mseLoc +"  Eye loc in world :"+ c.eyeInWorld+ dispWinFrames[curFocusWin].getCamDisp());//" camera rx :  " + rx + " ry : " + ry + " dz : " + dz);
			String[] res = ((mySideBarMenu)dispWinFrames[dispMenuIDX]).getDebugData();		//get debug data for each UI object
			int numToPrint = min(res.length,80);
			for(int s=0;s<numToPrint;++s) {	addInfoStr(res[s]);}				//add info to string to be displayed for debug
			drawInfoStr(1.0f, new int[] {255,255,255,100}); 	
			popStyle();	popMatrix();		
		}
		else if(showInfo){
			pushMatrix();pushStyle();			
			reInitInfoStr();	
			String[] res = consoleStrings.toArray(new String[0]);
			int dispNum = min(res.length, 80);
			for(int i=0;i<dispNum;++i){addInfoStr(res[i]);}
		    drawInfoStr(1.1f, new int[] {255,255,255,100}); 
			popStyle();	popMatrix();	
		}
	}
	//print out multiple-line text to screen
	public void ml_text(String str, float x, float y){
		String[] res = str.split("\\r?\\n");
		float disp = 0;
		for(int i =0; i<res.length; ++i){
			text(res[i],x, y+disp);		//add console string output to screen display- decays over time
				disp += 12;
			}
	}
	
	//print out string in display window
	public void outStr2Scr(String str){outStr2Scr(str,true);}
	//print informational string data to console, and to screen
	public void outStr2Scr(String str, boolean showDraw){
		if(trim(str) != ""){	System.out.println(str);}
		String[] res = str.split("\\r?\\n");
		if(showDraw){
			for(int i =0; i<res.length; ++i){
				consoleStrings.add(res[i]);		//add console string output to screen display- decays over time
			}
		}
	}
	//build a date with each component separated by token
	public String getDateTimeString(){return getDateTimeString(true, false,".");}
	public String getDateTimeString(boolean useYear, boolean toSecond, String token){
		String result = "";
		int val;
		if(useYear){val = now.get(Calendar.YEAR);		result += ""+val+token;}
		val = now.get(Calendar.MONTH)+1;				result += (val < 10 ? "0"+val : ""+val)+ token;
		val = now.get(Calendar.DAY_OF_MONTH);			result += (val < 10 ? "0"+val : ""+val)+ token;
		val = now.get(Calendar.HOUR_OF_DAY);					result += (val < 10 ? "0"+val : ""+val)+ token;
		val = now.get(Calendar.MINUTE);					result += (val < 10 ? "0"+val : ""+val);
		if(toSecond){val = now.get(Calendar.SECOND);	result += token + (val < 10 ? "0"+val : ""+val);}
		return result;
	}
	//utilities
	
	//handle user-driven file load or save - returns a filename + filepath string
	public String FileSelected(File selection){
		if (null==selection){return null;}
		return selection.getAbsolutePath();		
	}//FileSelected
	
	//		//s-cut to print to console
	public void pr(String str){outStr2Scr(str);}
	
	public String getFName(String fNameAndPath){
		String[] strs = fNameAndPath.split("/");
		return strs[strs.length-1];
	}
	
	//load a file as text strings
	public String[] loadFileIntoStringAra(String fileName, String dispYesStr, String dispNoStr){
		String[] strs = null;
		try{
			strs = loadStrings(fileName);
			System.out.println(dispYesStr+"\tLength : " + strs.length);
	} catch (Exception e){System.out.println("!!"+dispNoStr);return null;}
		return strs;		
	}//loadFileIntoStrings
	
	//public void scribeHeaderRight(String s) {scribeHeaderRight(s, 20);} // writes black on screen top, right-aligned
	//public void scribeHeaderRight(String s, float y) {fill(0); text(s,width-6*s.length(),y); noFill();} // writes black on screen top, right-aligned
	public void displayHeader(int[] clr) { // Displays title and authors face on screen
		float stVal = 17;
		int idx = 1;	
		pushMatrix();pushStyle();
		translate(0,10,0);
		setFill(clr);
		text("Shift-Click-Drag to change view.",width-190, stVal*idx++);
		text("Shift-RClick-Drag to zoom.",width-160, stVal*idx++); 
		text(authorString,width-75, stVal*idx++); 
		popStyle();popMatrix();
	}
	
	//project passed point onto box surface based on location - to help visualize the location in 3d
	public void drawProjOnBox(myPoint p){
		//myPoint[]  projOnPlanes = new myPoint[6];
		myPoint prjOnPlane;
		//public myPoint intersectPl(myPoint E, myVector T, myPoint A, myPoint B, myPoint C) { // if ray from E along T intersects triangle (A,B,C), return true and set proposal to the intersection point
		pushMatrix();
		translate(-p.x,-p.y,-p.z);
		for(int i  = 0; i< 6; ++i){				
			prjOnPlane = bndChkInCntrdBox3D(intersectPl(p, boxNorms[i], boxWallPts[i][0],boxWallPts[i][1],boxWallPts[i][2]));				
			show(prjOnPlane,5,rgbClrs[i/2],rgbClrs[i/2], false);				
		}
		popMatrix();
	}//drawProjOnBox
	
	public double min3(double x, double y, double z) { return (x < y ? (x < z ? x : z) : (y < z ? y : z)); }
	public double max3(double x, double y, double z) { return (x > y ? (x > z ? x : z) : (y > z ? y : z)); }

	private static final double lcl_third = 1.0/3.0;
	//return a random position within a sphere
	public myVectorf getRandPosInSphereSlow(double rad, myVectorf ctr){
		myVectorf pos = new myVectorf();
		do{
			double u = ThreadLocalRandom.current().nextDouble(0,1), r = rad * Math.pow(u, lcl_third),
					cosTheta = ThreadLocalRandom.current().nextDouble(-1,1), sinTheta =  Math.sin(Math.acos(cosTheta)),
					phi = ThreadLocalRandom.current().nextDouble(0,PConstants.TWO_PI);
			pos.set(sinTheta * Math.cos(phi), sinTheta * Math.sin(phi),cosTheta);
			pos._mult(r);
			pos._add(ctr);
		} while (pos.z < 0);
		return pos;
	}
	//return a random position within a sphere of radius rad centered at ctr
	public myVectorf getRandPosInSphere(double rad, myVectorf ctr){
		myVectorf pos = new myVectorf();
		double u = ThreadLocalRandom.current().nextDouble(0,1), r = rad * Math.pow(u, lcl_third);
		do{
			pos.set(ThreadLocalRandom.current().nextDouble(-1,1), ThreadLocalRandom.current().nextDouble(-1,1),ThreadLocalRandom.current().nextDouble(-1,1));
		} while (pos.sqMagn > 1.0f);
		pos._mult(r);
		pos._add(ctr);
		return pos;
	}
	
	/** 
	 * convert from spherical coords to cartesian
	 * @param rad
	 * @param thet
	 * @param phi
	 * @return ara : norm, surface point == x,y,z of coords passed
	 */
	public myVectorf[] getXYZFromRThetPhi(float rad, float thet, float phi, float scaleZ) {
		double sinThet = Math.sin(thet);	
		myVectorf[] res = new myVectorf[2];
		res[1] = new myVectorf(sinThet * Math.cos(phi) * rad, sinThet * Math.sin(phi) * rad,Math.cos(thet)*rad*scaleZ);
		res[0] = myVectorf._normalize(res[1]);
		return res;
	}//
	
	/**
	 * builds a list of N regularly placed vertices and normals for a sphere of radius rad centered at ctr
	 */	
	public myVectorf[][] getRegularSphereList(float rad, int N, float scaleZ) {
		ArrayList<myVectorf[]> res = new ArrayList<myVectorf[]>();
		//choose 1 point per dArea, where dArea is area of sphere parsed into N equal portions
		float lclA = 4*PI/N, lclD = sqrt(lclA);
		int Mthet = round(PI/lclD), Mphi;
		float dThet = PI/Mthet, dPhi = lclA/dThet, thet, phi, twoPiOvDPhi = TWO_PI/dPhi;
		for(int i=0;i<Mthet;++i) {
			thet = dThet * (i + 0.5f);
			Mphi = round(twoPiOvDPhi * sin(thet));
			for (int j=0;j<Mphi; ++j) { 
				phi = (TWO_PI*j)/Mphi;		
				res.add(getXYZFromRThetPhi(rad, thet, phi, scaleZ));
			}
		}
		return res.toArray(new myVectorf[0][]);
	}//getRegularSphereList	
	
	
	//find random point on a sphere of radius rad centered at ctr, and norm vector from center to point
	public myVectorf[] getRandPosOnSphere(double rad, float scale) {//, myVectorf ctr){
		myVectorf[] res = new myVectorf[2];//idx 0 norm from center, idx 1 location in space (norm * rad + ctr
		double 	cosTheta = ThreadLocalRandom.current().nextDouble(-1,1), sinTheta =  Math.sin(Math.acos(cosTheta)),
				phi = ThreadLocalRandom.current().nextDouble(0,PConstants.TWO_PI);
		res[1] = new myVectorf(sinTheta * Math.cos(phi)*rad, sinTheta * Math.sin(phi)*rad,cosTheta *rad* scale);
		res[0] = myVectorf._normalize(res[1]);
		//res[1]._add(ctr);
		return res;
	}

	
	//very fast mechanism for setting an array of doubles to a specific val - takes advantage of caching
	public void dAraFill(double[] ara, double val){
		int len = ara.length;
		if (len > 0){ara[0] = val; }
		for (int i = 1; i < len; i += i){  System.arraycopy(ara, 0, ara, i, ((len - i) < i) ? (len - i) : i);  }		
	}
	
	//convert a world location within the bounded cube region to be a 4-int color array
	public int[] getClrFromCubeLoc(float[] t){
		return new int[]{(int)(255*(t[0]-cubeBnds[0][0])/cubeBnds[1][0]),(int)(255*(t[1]-cubeBnds[0][1])/cubeBnds[1][1]),(int)(255*(t[2]-cubeBnds[0][2])/cubeBnds[1][2]),255};
	}
	
	//return a string holding the hex format of a passed value
	public String intToHexStr(int val){return  String.format("0x%08X", val);}

	
	//random location within coords[0] and coords[1] extremal corners of a cube - bnds is to give a margin of possible random values
	public myVectorf getRandPosInCube(float[][] coords, float bnds){
		return new myVectorf(
				ThreadLocalRandom.current().nextDouble(coords[0][0]+bnds,(coords[0][0] + coords[1][0] - bnds)),
				ThreadLocalRandom.current().nextDouble(coords[0][1]+bnds,(coords[0][1] + coords[1][1] - bnds)),
				ThreadLocalRandom.current().nextDouble(coords[0][2]+bnds,(coords[0][2] + coords[1][2] - bnds)));}		
	public myPoint getScrLocOf3dWrldPt(myPoint pt){	return new myPoint(screenX((float)pt.x,(float)pt.y,(float)pt.z),screenY((float)pt.x,(float)pt.y,(float)pt.z),screenZ((float)pt.x,(float)pt.y,(float)pt.z));}
	
	public myPoint bndChkInBox2D(myPoint p){p.set(Math.max(0,Math.min(p.x,grid2D_X)),Math.max(0,Math.min(p.y,grid2D_Y)),0);return p;}
	public myPoint bndChkInBox3D(myPoint p){p.set(Math.max(0,Math.min(p.x,gridDimX)), Math.max(0,Math.min(p.y,gridDimY)),Math.max(0,Math.min(p.z,gridDimZ)));return p;}	
	public myPoint bndChkInCntrdBox3D(myPoint p){
		p.set(Math.max(-hGDimX,Math.min(p.x,hGDimX)), 
				Math.max(-hGDimY,Math.min(p.y,hGDimY)),
				Math.max(-hGDimZ,Math.min(p.z,hGDimZ)));return p;}	
	 
	public void translate(myPoint p){translate((float)p.x,(float)p.y,(float)p.z);}
	public void translate(myPointf p){translate(p.x,p.y,p.z);}
	public void translate(myVector p){translate((float)p.x,(float)p.y,(float)p.z);}
	public void translate(double x, double y, double z){translate((float)x,(float)y,(float)z);}
	public void translate(double x, double y){translate((float)x,(float)y);}
	public void rotate(float thet, myPoint axis){rotate(thet, (float)axis.x,(float)axis.y,(float)axis.z);}
	public void rotate(float thet, double x, double y, double z){rotate(thet, (float)x,(float)y,(float)z);}
	//************************************************************************
	//**** SPIRAL
	//************************************************************************
	//3d rotation - rotate P by angle a around point G and axis normal to plane IJ
	public myPoint R(myPoint P, double a, myVector I, myVector J, myPoint G) {
		double x= myVector._dot(new myVector(G,P),U(I)), y=myVector._dot(new myVector(G,P),U(J)); 
		double c=Math.cos(a), s=Math.sin(a); 
		double iXVal = x*c-x-y*s, jYVal= x*s+y*c-y;			
		return myPoint._add(P,iXVal,I,jYVal,J); }; 
		
	public cntlPt R(cntlPt P, double a, myVector I, myVector J, myPoint G) {
		double x= myVector._dot(new myVector(G,P),U(I)), y=myVector._dot(new myVector(G,P),U(J)); 
		double c=Math.cos(a), s=Math.sin(a); 
		double iXVal = x*c-x-y*s, jYVal= x*s+y*c-y;		
		return new cntlPt(this, P(P,iXVal,I,jYVal,J), P.r, P.w); };
		
	public myPoint PtOnSpiral(myPoint A, myPoint B, myPoint C, double t) {
		//center is coplanar to A and B, and coplanar to B and C, but not necessarily coplanar to A, B and C
		//so center will be coplanar to mp(A,B) and mp(B,C) - use mpCA midpoint to determine plane mpAB-mpBC plane?
		myPoint mAB = new myPoint(A,.5f, B);
		myPoint mBC = new myPoint(B,.5f, C);
		myPoint mCA = new myPoint(C,.5f, A);
		myVector mI = U(mCA,mAB);
		myVector mTmp = myVector._cross(mI,U(mCA,mBC));
		myVector mJ = U(mTmp._cross(mI));	//I and J are orthonormal
		double a =spiralAngle(A,B,B,C); 
		double s =spiralScale(A,B,B,C);
		
		//myPoint G = spiralCenter(a, s, A, B, mI, mJ); 
		myPoint G = spiralCenter(A, mAB, B, mBC); 
		return new myPoint(G, Math.pow(s,t), R(A,t*a,mI,mJ,G));
	  }
	public double spiralAngle(myPoint A, myPoint B, myPoint C, myPoint D) {return myVector._angleBetween(new myVector(A,B),new myVector(C,D));}
	public double spiralScale(myPoint A, myPoint B, myPoint C, myPoint D) {return myPoint._dist(C,D)/ myPoint._dist(A,B);}
	
	public myPoint R(myPoint Q, myPoint C, myPoint P, myPoint R) { // returns rotated version of Q by angle(CP,CR) parallel to plane (C,P,R)
		myVector I0=U(C,P), I1=U(C,R), V=new myVector(C,Q); 
		double c=myPoint._dist(I0,I1), s=Math.sqrt(1.-(c*c)); 
		if(Math.abs(s)<0.00001) return Q;
		myVector J0=V(1./s,I1,-c/s,I0);  
		myVector J1=V(-s,I0,c,J0);  
		double x=V._dot(I0), y=V._dot(J0);  
		return P(Q,x,M(I1,I0),y,M(J1,J0)); 
	} 	
	// spiral given 4 points, AB and CD are edges corresponding through rotation
	public myPoint spiralCenter(myPoint A, myPoint B, myPoint C, myPoint D) {         // new spiral center
		myVector AB=V(A,B), CD=V(C,D), AC=V(A,C);
		double m=CD.magn/AB.magn, n=CD.magn*AB.magn;		
		myVector rotAxis = U(AB._cross(CD));		//expect ab and ac to be coplanar - this is the axis to rotate around to find f		
		myVector rAB = myVector._rotAroundAxis(AB, rotAxis, PConstants.HALF_PI);
		double c=AB._dot(CD)/n, 
				s=rAB._dot(CD)/n;
		double AB2 = AB._dot(AB), a=AB._dot(AC)/AB2, b=rAB._dot(AC)/AB2;
		double x=(a-m*( a*c+b*s)), y=(b-m*(-a*s+b*c));
		double d=1+m*(m-2*c);  if((c!=1)&&(m!=1)) { x/=d; y/=d; };
		return P(P(A,x,AB),y,rAB);
	  }
	
	private static float cyl_da = TWO_PI/36.0f;
	public void cylinder(myPoint A, myPoint B, float r, int c1, int c2) {
		myPoint P = A;
		myVector V = V(A,B);
		myVector I = myVector.UP;//c.drawSNorm;//U(Normal(V)); -- need vector that will not be singular on xprod below
		myVector Nvec = N(I,V);
		if(Math.abs(Nvec.magn) < this.epsValCalc) {//singular - cylinder wanting to go up
			I = myVector.RIGHT;Nvec = N(I,V);
		}
		myVector J = U(Nvec);
		beginShape(QUAD_STRIP);
		float rcA, rsA;
		setColorValStroke(c1);
		noFill();
		for(float a=0; a<=TWO_PI+cyl_da; a+=cyl_da) {rcA = r*cos(a); rsA = r*sin(a);setColorValStroke(c1); gl_vertex(P(P,rcA,I,rsA,J,0,V)); setColorValStroke(c2); gl_vertex(P(P,rcA,I,rsA,J,1,V));}
		endShape();
	}
	
	public void cylinder(myPointf A, myPointf B, float r, int c1, int c2) {
		cylinder(new myPoint(A.x, A.y, A.z), new myPoint(B.x, B.y, B.z), r, c1, c2);
	}
	
	public void drawCylinder(ArrayList<myPointf>[] pts, int c1, int c2) {
		beginShape(QUAD_STRIP);
		setColorValStroke(c1);
		noFill();
		for(int i=0;i<pts[0].size();++i) {
			setColorValStroke(c1); 
			gl_vertex(pts[0].get(i)); 
			setColorValStroke(c2); 
			gl_vertex(pts[1].get(i)); 
		}		
		endShape();
	}	
	
	//point functions
	public myPoint P() {return new myPoint(); };                                                                          // point (x,y,z)
	public myPoint P(double x, double y, double z) {return new myPoint(x,y,z); };                                            // point (x,y,z)
	public myPoint P(myPoint A) {return new myPoint(A.x,A.y,A.z); };                                                           // copy of point P
	public myPoint P(myPoint A, double s, myPoint B) {return new myPoint(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };        // A+sAB
	public myPoint L(myPoint A, double s, myPoint B) {return new myPoint(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };        // A+sAB
	public myPoint P(myPoint A, myPoint B) {return P((A.x+B.x)/2.0,(A.y+B.y)/2.0,(A.z+B.z)/2.0); }                             // (A+B)/2
	public myPoint P(myPoint A, myPoint B, myPoint C) {return new myPoint((A.x+B.x+C.x)/3.0,(A.y+B.y+C.y)/3.0,(A.z+B.z+C.z)/3.0); };     // (A+B+C)/3
	public myPoint P(myPoint A, myPoint B, myPoint C, myPoint D) {return P(P(A,B),P(C,D)); };                                            // (A+B+C+D)/4
	public myPoint P(double s, myPoint A) {return new myPoint(s*A.x,s*A.y,s*A.z); };                                            // sA
	public myPoint A(myPoint A, myPoint B) {return new myPoint(A.x+B.x,A.y+B.y,A.z+B.z); };                                         // A+B
	public myPoint P(double a, myPoint A, double b, myPoint B) {return A(P(a,A),P(b,B));}                                        // aA+bB 
	public myPoint P(double a, myPoint A, double b, myPoint B, double c, myPoint C) {return A(P(a,A),P(b,B,c,C));}                     // aA+bB+cC 
	public myPoint P(double a, myPoint A, double b, myPoint B, double c, myPoint C, double d, myPoint D){return A(P(a,A,b,B),P(c,C,d,D));}   // aA+bB+cC+dD
	public myPoint P(myPoint P, myVector V) {return new myPoint(P.x + V.x, P.y + V.y, P.z + V.z); }                                 // P+V
	public myPoint P(myPoint P, double s, myVector V) {return new myPoint(P.x+s*V.x,P.y+s*V.y,P.z+s*V.z);}                           // P+sV
	public myPoint P(myPoint O, double x, myVector I, double y, myVector J) {return P(O.x+x*I.x+y*J.x,O.y+x*I.y+y*J.y,O.z+x*I.z+y*J.z);}  // O+xI+yJ
	public myPoint P(myPoint O, double x, myVector I, double y, myVector J, double z, myVector K) {return new myPoint(O.x+x*I.x+y*J.x+z*K.x,O.y+x*I.y+y*J.y+z*K.y,O.z+x*I.z+y*J.z+z*K.z);}  // O+xI+yJ+kZ
	public myPointf Pf(myPointf O, float x, myVectorf I, double y, myVectorf J, double z, myVectorf K) {return new myPointf(O.x+x*I.x+y*J.x+z*K.x,O.y+x*I.y+y*J.y+z*K.y,O.z+x*I.z+y*J.z+z*K.z);}  // O+xI+yJ+zK
	void makePts(myPoint[] C) {for(int i=0; i<C.length; i++) C[i]=P();}
	
	//draw a circle - JT
	void circle(myPoint P, float r, myVector I, myVector J, int n) {myPoint[] pts = new myPoint[n];pts[0] = P(P,r,U(I));float a = (2*PI)/(1.0f*n);for(int i=1;i<n;++i){pts[i] = R(pts[i-1],a,J,I,P);}pushMatrix(); pushStyle();noFill(); show(pts);popStyle();popMatrix();}; // render sphere of radius r and center P
	
	void circle(myPoint p, float r){ellipse((float)p.x, (float)p.y, r, r);}
	void circle(float x, float y, float r1, float r2){ellipse(x,y, r1, r2);}
	
	void noteArc(float[] dims, int[] noteClr){
		noFill();
		setStroke(noteClr);
		strokeWeight(1.5f*dims[3]);
		arc(0,0, dims[2], dims[2], dims[0] - this.HALF_PI, dims[1] - this.HALF_PI);
	}
	//draw a ring segment from alphaSt in radians to alphaEnd in radians
	void noteArc(myPoint ctr, float alphaSt, float alphaEnd, float rad, float thickness, int[] noteClr){
		noFill();
		setStroke(noteClr);
		strokeWeight(thickness);
		arc((float)ctr.x, (float)ctr.y, rad, rad, alphaSt - this.HALF_PI, alphaEnd- this.HALF_PI);
	}
	
	
	void bezier(myPoint A, myPoint B, myPoint C, myPoint D) {bezier((float)A.x,(float)A.y,(float)A.z,(float)B.x,(float)B.y,(float)B.z,(float)C.x,(float)C.y,(float)C.z,(float)D.x,(float)D.y,(float)D.z);} // draws a cubic Bezier curve with control points A, B, C, D
	void bezier(myPoint [] C) {bezier(C[0],C[1],C[2],C[3]);} // draws a cubic Bezier curve with control points A, B, C, D
	myPoint bezierPoint(myPoint[] C, float t) {return P(bezierPoint((float)C[0].x,(float)C[1].x,(float)C[2].x,(float)C[3].x,(float)t),bezierPoint((float)C[0].y,(float)C[1].y,(float)C[2].y,(float)C[3].y,(float)t),bezierPoint((float)C[0].z,(float)C[1].z,(float)C[2].z,(float)C[3].z,(float)t)); }
	myVector bezierTangent(myPoint[] C, float t) {return V(bezierTangent((float)C[0].x,(float)C[1].x,(float)C[2].x,(float)C[3].x,(float)t),bezierTangent((float)C[0].y,(float)C[1].y,(float)C[2].y,(float)C[3].y,(float)t),bezierTangent((float)C[0].z,(float)C[1].z,(float)C[2].z,(float)C[3].z,(float)t)); }
	
	
	public myPoint Mouse() {return new myPoint(mouseX, mouseY,0);}                                          			// current mouse location
	public myVector MouseDrag() {return new myVector(mouseX-pmouseX,mouseY-pmouseY,0);};                     			// vector representing recent mouse displacement
	
	//public int color(myPoint p){return color((int)p.x,(int)p.z,(int)p.y);}	//needs to be x,z,y for some reason - to match orientation of color frames in z-up 3d geometry
	public int color(myPoint p){return color((int)p.x,(int)p.y,(int)p.z);}	
	
	// =====  vector functions
	public myVector V() {return new myVector(); };                                                                          // make vector (x,y,z)
	public myVector V(double x, double y, double z) {return new myVector(x,y,z); };                                            // make vector (x,y,z)
	public myVector V(myVector V) {return new myVector(V.x,V.y,V.z); };                                                          // make copy of vector V
	public myVector A(myVector A, myVector B) {return new myVector(A.x+B.x,A.y+B.y,A.z+B.z); };                                       // A+B
	public myVector A(myVector U, float s, myVector V) {return V(U.x+s*V.x,U.y+s*V.y,U.z+s*V.z);};                               // U+sV
	public myVector M(myVector U, myVector V) {return V(U.x-V.x,U.y-V.y,U.z-V.z);};                                              // U-V
	public myVector M(myVector V) {return V(-V.x,-V.y,-V.z);};                                                              // -V
	public myVector V(myVector A, myVector B) {return new myVector((A.x+B.x)/2.0,(A.y+B.y)/2.0,(A.z+B.z)/2.0); }                      // (A+B)/2
	public myVector V(myVector A, float s, myVector B) {return new myVector(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };      // (1-s)A+sB
	public myVector V(myVector A, myVector B, myVector C) {return new myVector((A.x+B.x+C.x)/3.0,(A.y+B.y+C.y)/3.0,(A.z+B.z+C.z)/3.0); };  // (A+B+C)/3
	public myVector V(myVector A, myVector B, myVector C, myVector D) {return V(V(A,B),V(C,D)); };                                         // (A+B+C+D)/4
	public myVector V(double s, myVector A) {return new myVector(s*A.x,s*A.y,s*A.z); };                                           // sA
	public myVector V(double a, myVector A, double b, myVector B) {return A(V(a,A),V(b,B));}                                       // aA+bB 
	public myVector V(double a, myVector A, double b, myVector B, double c, myVector C) {return A(V(a,A,b,B),V(c,C));}                   // aA+bB+cC
	public myVector V(myPoint P, myPoint Q) {return new myVector(P,Q);};                                          // PQ
	public myVector N(myVector U, myVector V) {return V( U.y*V.z-U.z*V.y, U.z*V.x-U.x*V.z, U.x*V.y-U.y*V.x); };                  // UxV cross product (normal to both)
	public myVectorf Nf(myVectorf U, myVectorf V) {return new myVectorf( U.y*V.z-U.z*V.y, U.z*V.x-U.x*V.z, U.x*V.y-U.y*V.x); };                  // UxV cross product (normal to both)
	public myVector N(myPoint A, myPoint B, myPoint C) {return N(V(A,B),V(A,C)); };                                                   // normal to triangle (A,B,C), not normalized (proportional to area)
	public myVector B(myVector U, myVector V) {return U(N(N(U,V),U)); }        

	//calculate the normal, tangent, binormal components of passed partVec compared to a passed normal (needs to be normalized)
	public myVectorf[] getVecFrame(myVectorf partVec, myVectorf norm) {
		myVectorf[] result = new myVectorf[3];//(2, myVector(0, 0, 0));
		result[0] = myVectorf._mult(norm,(norm._dot(partVec)));//norm dir
		result[1] = myVectorf._sub(partVec, result[0]);		//tan dir
		result[2] = myVectorf._cross(result[0], result[1]);
		return result;
	}
	
	public double d(myVector U, myVector V) {return U.x*V.x+U.y*V.y+U.z*V.z; };                                            //U*V dot product
	public double dot(myVector U, myVector V) {return U.x*V.x+U.y*V.y+U.z*V.z; };                                            //U*V dot product
	public double det2(myVector U, myVector V) {return -U.y*V.x+U.x*V.y; };                                       		// U|V det product
	public double det3(myVector U, myVector V) {double dist = d(U,V); return Math.sqrt(d(U,U)*d(V,V) - (dist*dist)); };                                // U|V det product
	public double m(myVector U, myVector V, myVector W) {return d(U,N(V,W)); };                                                 // (UxV)*W  mixed product, determinant - measures 6x the volume of the parallelapiped formed by myVectortors
	public double m(myPoint E, myPoint A, myPoint B, myPoint C) {return m(V(E,A),V(E,B),V(E,C));}                                    // det (EA EB EC) is >0 when E sees (A,B,C) clockwise
	public double n2(myVector V) {return (V.x*V.x)+(V.y*V.y)+(V.z*V.z);};                                                   // V*V    norm squared
	public double n(myVector V) {return  Math.sqrt(n2(V));};                                                                // ||V||  norm
	public double d(myPoint P, myPoint Q) {return  myPoint._dist(P, Q); };                            // ||AB|| distance
	public double area(myPoint A, myPoint B, myPoint C) {return n(N(A,B,C))/2; };                                               // area of triangle 
	public double volume(myPoint A, myPoint B, myPoint C, myPoint D) {return m(V(A,B),V(A,C),V(A,D))/6; };                           // volume of tet 
	public boolean parallel (myVector U, myVector V) {return n(N(U,V))<n(U)*n(V)*0.00001; }                              // true if U and V are almost parallel
	public double angle(myPoint A, myPoint B, myPoint C){return angle(V(A,B),V(A,C));}												//angle between AB and AC
	public double angle(myPoint A, myPoint B, myPoint C, myPoint D){return angle(U(A,B),U(C,D));}							//angle between AB and CD
	public double angle(myVector U, myVector V){double angle = Math.atan2(n(N(U,V)),d(U,V)),sign = m(U,V,V(0,0,1));if(sign<0){    angle=-angle;}	return angle;}
	public boolean cw(myVector U, myVector V, myVector W) {return m(U,V,W)>0; };                                               // (UxV)*W>0  U,V,W are clockwise
	public boolean cw(myPoint A, myPoint B, myPoint C, myPoint D) {return volume(A,B,C,D)>0; };                                     // tet is oriented so that A sees B, C, D clockwise 
	public boolean projectsBetween(myPoint P, myPoint A, myPoint B) {return dot(V(A,P),V(A,B))>0 && dot(V(B,P),V(B,A))>0 ; };
	public double distToLine(myPoint P, myPoint A, myPoint B) {double res = det3(U(A,B),V(A,P)); return Double.isNaN(res) ? 0 : res; };		//MAY RETURN NAN IF point P is on line
	public myPoint projectionOnLine(myPoint P, myPoint A, myPoint B) {return P(A,dot(V(A,B),V(A,P))/dot(V(A,B),V(A,B)),V(A,B));}
	public boolean isSame(myPoint A, myPoint B) {return (A.x==B.x)&&(A.y==B.y)&&(A.z==B.z) ;}                                         // A==B
	public boolean isSame(myPoint A, myPoint B, double e) {return ((Math.abs(A.x-B.x)<e)&&(Math.abs(A.y-B.y)<e)&&(Math.abs(A.z-B.z)<e));}                   // ||A-B||<e
	
	public myVector W(double s,myVector V) {return V(s*V.x,s*V.y,s*V.z);}                                                      // sV
	
	public myVector U(myVector v){myVector u = new myVector(v); return u._normalize(); }
	public myVector U(myVector v, float d, myVector u){myVector r = new myVector(v,d,u); return r._normalize(); }
	public myVector Upt(myPoint v){myVector u = new myVector(v); return u._normalize(); }
	public myVector U(myPoint a, myPoint b){myVector u = new myVector(a,b); return u._normalize(); }
	public myVectorf Uf(myPoint a, myPoint b){myVectorf u = new myVectorf(a,b); return u._normalize(); }
	public myVector U(double x, double y, double z) {myVector u = new myVector(x,y,z); return u._normalize();}
	
	public myVector normToPlane(myPoint A, myPoint B, myPoint C) {return myVector._cross(new myVector(A,B),new myVector(A,C)); };   // normal to triangle (A,B,C), not normalized (proportional to area)
	
	public void gl_normal(myVector V) {normal((float)V.x,(float)V.y,(float)V.z);}                                          // changes normal for smooth shading
	public void gl_vertex(myPoint P) {vertex((float)P.x,(float)P.y,(float)P.z);}                                           // vertex for shading or drawing
	public void gl_normal(myVectorf V) {normal(V.x,V.y,V.z);}                                          // changes normal for smooth shading
	public void gl_vertex(myPointf P) {vertex(P.x,P.y,P.z);}                                           // vertex for shading or drawing
	///show functions
	public void showVec( myPoint ctr, double len, myVector v){line(ctr.x,ctr.y,ctr.z,ctr.x+(v.x)*len,ctr.y+(v.y)*len,ctr.z+(v.z)*len);}
	public void show(myPoint P, double r, int clr, String txt) {
		pushMatrix(); pushStyle(); 
		if(clr!= -1){setColorValFill(clr); setColorValStroke(clr);}
		sphereDetail(5);
		translate((float)P.x,(float)P.y,(float)P.z); 
		sphere((float)r); 
		setColorValFill(gui_Black);setColorValStroke(gui_Black);
		double d = 1.1 * r;
		show(myPoint.ZEROPT, txt, new myVector(d,d,d));
		popStyle(); popMatrix();} // render sphere of radius r and center P)

	public void show(myPoint P, double r, int clr) {pushMatrix(); pushStyle(); if(clr!= -1){setColorValFill(clr); setColorValStroke(clr);}
			sphereDetail(5);translate((float)P.x,(float)P.y,(float)P.z); sphere((float)r); popStyle(); popMatrix();} // render sphere of radius r and center P)
	public void show(myPoint P, double r,int fclr, int sclr, boolean flat) {//TODO make flat circles for points if flat
		pushMatrix(); pushStyle(); 
		if((fclr!= -1) && (sclr!= -1)){setColorValFill(fclr); setColorValStroke(sclr);}
		if(!flat){
			translate((float)P.x,(float)P.y,(float)P.z); 
			sphereDetail(5);
			sphere((float)r);
		} else {
			translate((float)P.x,(float)P.y,0); 
			this.circle(0,0,(float)r,(float)r);				
		}
	popStyle(); popMatrix();} // render sphere of radius r and center P)
	
	public void show(myPoint P, double rad, int fclr, int sclr, int tclr, String txt) {
		pushMatrix(); pushStyle(); 
		if((fclr!= -1) && (sclr!= -1)){setColorValFill(fclr); setColorValStroke(sclr);}
		sphereDetail(5);
		translate((float)P.x,(float)P.y,(float)P.z); 
		setColorValFill(tclr);setColorValStroke(tclr);
		showOffsetText(1.2f * (float)rad,tclr, txt);
		popStyle(); popMatrix();
	} // render sphere of radius r and center P)
	
	public void show(myPoint P, double r, int fclr, int sclr) {
		pushMatrix(); pushStyle(); 
		if((fclr!= -1) && (sclr!= -1)){setColorValFill(fclr); setColorValStroke(sclr);}
		sphereDetail(5);
		translate((float)P.x,(float)P.y,(float)P.z); 
		sphere((float)r); 
		popStyle(); popMatrix();
	} // render sphere of radius r and center P)
	
	public void vTextured(myPointf P, float u, float v) {vertex((float)P.x,(float)P.y,(float)P.z,(float)u,(float)v);};                          // vertex with texture coordinates
	public void show(myPoint P, double r){show(P,r, gui_Black, gui_Black, false);}
	public void show(myPoint P, String s) {text(s, (float)P.x, (float)P.y, (float)P.z); } // prints string s in 3D at P
	public void show(myPoint P, String s, myVector D) {text(s, (float)(P.x+D.x), (float)(P.y+D.y), (float)(P.z+D.z));  } // prints string s in 3D at P+D
	public void show(myPoint P, double r, String s, myVector D){show(P,r, gui_Black, gui_Black, false);pushStyle();setColorValFill(gui_Black);show(P,s,D);popStyle();}
	public void show(myPoint P, double r, String s, myVector D, int clr, boolean flat){show(P,r, clr, clr, flat);pushStyle();setColorValFill(clr);show(P,s,D);popStyle();}
	public void show(myPoint[] ara) {beginShape(); for(int i=0;i<ara.length;++i){gl_vertex(ara[i]);} endShape(CLOSE);};                     
	public void show(myPoint[] ara, myVector norm) {beginShape();gl_normal(norm); for(int i=0;i<ara.length;++i){gl_vertex(ara[i]);} endShape(CLOSE);};                     
	
	public void showVec( myPointf ctr, float len, myVectorf v){line(ctr.x,ctr.y,ctr.z,ctr.x+(v.x)*len,ctr.y+(v.y)*len,ctr.z+(v.z)*len);}
	
	public void showOffsetText(float d, int tclr, String txt){
		setColorValFill(tclr);setColorValStroke(tclr);
		text(txt, d, d,d); 
	}
	public void showOffsetText(float[] c, int tclr, String txt){
		setColorValFill(tclr);setColorValStroke(tclr);
		text(txt, c[0], c[1], c[2]); 
	}

	public void showFlat(myPointf P, float r,int fclr, int sclr, int tclr, String txt) {
		pushMatrix(); pushStyle(); 
		if((fclr!= -1) && (sclr!= -1)){setColorValFill(fclr); setColorValStroke(sclr);}
		translate(P.x,P.y,0); 
		circle(0,0,r,r);	
		setColorValFill(tclr);setColorValStroke(tclr);
		text(txt, r, r,r); 
		popStyle(); popMatrix();
	} // render sphere of radius r and center P)
	

	public void show(myPointf P, float r,int fclr, int sclr, boolean flat) {
		pushMatrix(); pushStyle(); 
		if((fclr!= -1) && (sclr!= -1)){setColorValFill(fclr); setColorValStroke(sclr);}
		if(!flat){
			translate(P.x,P.y,P.z); 
			sphereDetail(5);
			sphere(r);
		} else {
			translate(P.x,P.y,0); 
			circle(0,0,r,r);				
		}
		popStyle(); popMatrix();
	} // render sphere of radius r and center P)
	
	
	public void show(myPointf P, float rad, int fclr, int sclr, int tclr, String txt) {
		pushMatrix(); pushStyle(); 
		if((fclr!= -1) && (sclr!= -1)){setColorValFill(fclr); setColorValStroke(sclr);}
		sphereDetail(5);
		translate(P.x,P.y,P.z);
		sphere(rad); 
		showOffsetText(1.2f * rad,tclr, txt);
		popStyle(); popMatrix();
	} // render sphere of radius r and center P)
	
	public void show(myPointf P, float rad, int[] fclr, int[] sclr, int tclr, String txt) {
		pushMatrix(); pushStyle(); 
		if((fclr!= null) && (sclr!= null)){setFill(fclr,255); setStroke(sclr,255);}
		sphereDetail(5);
		translate(P.x,P.y,P.z);
		sphere(rad); 
		showOffsetText(1.2f * rad,tclr, txt);
		popStyle(); popMatrix();
	} // render sphere of radius r and center P)
	
	public void show(myPointf P, float rad, int det, int[] fclr, int[] sclr, int tclr, String txt) {
		pushMatrix(); pushStyle(); 
		setFill(fclr); 
		setStroke(sclr);
		sphereDetail(det);
		translate(P.x,P.y,P.z); 
		sphere(rad); 
		showOffsetText(1.2f * rad,tclr, txt);
		popStyle(); popMatrix();
	} // render sphere of radius r and center P)
	
	//show sphere of certain radius
	public void show(myPointf P, float rad, int det, int[] fclr, int[] sclr) {
		pushMatrix(); pushStyle(); 
		if((fclr!= null) && (sclr!= null)){setFill(fclr,255); setStroke(sclr,255);}
		sphereDetail(det);
		translate(P.x,P.y,P.z); 
		sphere(rad); 
		popStyle(); popMatrix();
	} // render sphere of radius r and center P)
	public void show(myPointf P, float rad, int det, int fclr, int sclr) {//only call with set fclr and sclr
		pushMatrix(); pushStyle(); 
		setColorValFill(fclr,255); 
		setColorValStroke(sclr,255);
		sphereDetail(det);
		translate(P.x,P.y,P.z); 
		sphere(rad); 
		popStyle(); popMatrix();
	} // render sphere of radius r and center P)
	
	public void show(myPointf P, float rad, int det){			
		pushMatrix(); pushStyle(); 
		fill(0,0,0,255); 
		stroke(0,0,0,255);
		sphereDetail(det);
		translate(P.x,P.y,P.z); 
		sphere(rad); 
		popStyle(); popMatrix();
	}
	
	public void show(myPointf P, String s) {text(s, P.x, P.y, P.z); } // prints string s in 3D at P
	public void show(myPointf P, String s, myVectorf D) {text(s, (P.x+D.x), (P.y+D.y),(P.z+D.z));  } // prints string s in 3D at P+D
	public void show(myPointf P, float r, String s, myVectorf D){show(P,r, gui_Black, gui_Black, false);pushStyle();setColorValFill(gui_Black);show(P,s,D);popStyle();}
	public void show(myPointf P, float r, String s, myVectorf D, int clr, boolean flat){show(P,r, clr, clr, flat);pushStyle();setColorValFill(clr);show(P,s,D);popStyle();}
	public void show(myPointf[] ara) {beginShape(); for(int i=0;i<ara.length;++i){gl_vertex(ara[i]);} endShape(CLOSE);};                     
	public void show(myPointf[] ara, myVectorf norm) {beginShape();gl_normal(norm); for(int i=0;i<ara.length;++i){gl_vertex(ara[i]);} endShape(CLOSE);};                     
	
	public void showNoClose(myPoint[] ara) {beginShape(); for(int i=0;i<ara.length;++i){gl_vertex(ara[i]);} endShape();};                     
	public void showNoClose(myPointf[] ara) {beginShape(); for(int i=0;i<ara.length;++i){gl_vertex(ara[i]);} endShape();};                     
	///end show functions
	
	public void curveVertex(myPoint P) {curveVertex((float)P.x,(float)P.y);};                                           // curveVertex for shading or drawing
	public void curve(myPoint[] ara) {if(ara.length == 0){return;}beginShape(); curveVertex(ara[0]);for(int i=0;i<ara.length;++i){curveVertex(ara[i]);} curveVertex(ara[ara.length-1]);endShape();};                      // volume of tet 	
	
	public boolean intersectPl(myPoint E, myVector T, myPoint A, myPoint B, myPoint C, myPoint X) { // if ray from E along T intersects triangle (A,B,C), return true and set proposal to the intersection point
		myVector EA=new myVector(E,A), AB=new myVector(A,B), AC=new myVector(A,C); 		double t = (float)(myVector._mixProd(EA,AC,AB) / myVector._mixProd(T,AC,AB));		X.set(myPoint._add(E,t,T));		return true;
	}	
	public myPoint intersectPl(myPoint E, myVector T, myPoint A, myPoint B, myPoint C) { // if ray from E along T intersects triangle (A,B,C), return true and set proposal to the intersection point
		myVector EA=new myVector(E,A), AB=new myVector(A,B), AC=new myVector(A,C); 		
		double t = (float)(myVector._mixProd(EA,AC,AB) / myVector._mixProd(T,AC,AB));		
		return (myPoint._add(E,t,T));		
	}	
	// if ray from E along V intersects sphere at C with radius r, return t when intersection occurs
	public double intersectPt(myPoint E, myVector V, myPoint C, double r) { 
		myVector Vce = V(C,E);
		double CEdCE = Vce._dot(Vce), VdV = V._dot(V), VdVce = V._dot(Vce), b = 2 * VdVce, c = CEdCE - (r*r),
				radical = (b*b) - 4 *(VdV) * c;
		if(radical < 0) return -1;
		double t1 = (b + Math.sqrt(radical))/(2*VdV), t2 = (b - Math.sqrt(radical))/(2*VdV);			
		return ((t1 > 0) && (t2 > 0) ? Math.min(t1, t2) : ((t1 < 0 ) ? ((t2 < 0 ) ? -1 : t2) : t1) );
		
	}
	
	public void rect(float[] a){rect(a[0],a[1],a[2],a[3]);}				//rectangle from array of floats : x, y, w, h	
	public myPoint WrldToScreen(myPoint wPt){return new myPoint(screenX((float)wPt.x,(float)wPt.y,(float)wPt.z),screenY((float)wPt.x,(float)wPt.y,(float)wPt.z),screenZ((float)wPt.x,(float)wPt.y,(float)wPt.z));}
	
	/////////////////////		
	///color utils
	/////////////////////
	public final int  // set more colors using Menu >  Tools > Color Selector
	  black=0xff000000, 
	  white=0xffFFFFFF,
	  red=0xffFF0000, 
	  green=0xff00FF00, 
	  blue=0xff0000FF, 
	  yellow=0xffFFFF00, 
	  cyan=0xff00FFFF, 
	  magenta=0xffFF00FF,
	  grey=0xff818181, 
	  orange=0xffFFA600, 
	  brown=0xffB46005, 
	  metal=0xffB5CCDE, 
	  dgreen=0xff157901;
	//set color based on passed point r= x, g = z, b=y
	public void fillAndShowLineByRBGPt(myPoint p, float x,  float y, float w, float h){
		fill((int)p.x,(int)p.y,(int)p.z);
		stroke((int)p.x,(int)p.y,(int)p.z);
		rect(x,y,w,h);
		//show(p,r,-1);
	}
	
	public int[][] triColors = new int[][] {
		{gui_DarkMagenta,gui_DarkBlue,gui_DarkGreen,gui_DarkCyan}, 
		{gui_LightMagenta,gui_LightBlue,gui_LightGreen,gui_TransCyan}};
	
	public int getRndClrInt(){return ThreadLocalRandom.current().nextInt(0,23);}		//return a random color flag value from below
	public int[] getRndClr(int alpha){return new int[]{ThreadLocalRandom.current().nextInt(0,250),ThreadLocalRandom.current().nextInt(0,250),ThreadLocalRandom.current().nextInt(0,250),alpha};	}
	public int[] getRndClr(){return getRndClr(255);	}		
	public Integer[] getClrMorph(int a, int b, double t){return getClrMorph(getClr(a), getClr(b), t);}    
	public Integer[] getClrMorph(int[] a, int[] b, double t){
		if(t==0){return new Integer[]{a[0],a[1],a[2],a[3]};} else if(t==1){return new Integer[]{b[0],b[1],b[2],b[3]};}
		return new Integer[]{(int)(((1.0f-t)*a[0])+t*b[0]),(int)(((1.0f-t)*a[1])+t*b[1]),(int)(((1.0f-t)*a[2])+t*b[2]),(int)(((1.0f-t)*a[3])+t*b[3])};
	}
	public void setFill(int[] clr){setFill(clr,clr[3]);}
	public void setStroke(int[] clr){setStroke(clr,clr[3]);}		
	public void setFill(int[] clr, int alpha){fill(clr[0],clr[1],clr[2], alpha);}
	public void setStroke(int[] clr, int alpha){stroke(clr[0],clr[1],clr[2], alpha);}
	
	public void setColorValFillAmbSh(PShape sh, int colorVal){ setColorValFillAmbSh(sh,colorVal,255);}
	public void setColorValFillAmbSh(PShape sh, int colorVal, int alpha){
		switch (colorVal){
			case gui_rnd				: { sh.fill(ThreadLocalRandom.current().nextInt(255),ThreadLocalRandom.current().nextInt(255),ThreadLocalRandom.current().nextInt(255),alpha); sh.ambient(120,120,120);break;}
	    	case gui_White  			: { sh.fill(255,255,255,alpha); sh.ambient(255,255,alpha); break; }
	    	case gui_Gray   			: { sh.fill(120,120,120,alpha); sh.ambient(120,120,120); break;}
	    	case gui_Yellow 			: { sh.fill(255,255,0,alpha); sh.ambient(255,255,0); break; }
	    	case gui_Cyan   			: { sh.fill(0,255,255,alpha); sh.ambient(0,255,255); break; }
	    	case gui_Magenta			: { sh.fill(255,0,255,alpha); sh.ambient(255,0,255); break; }
	    	case gui_Red    			: { sh.fill(255,0,0,alpha); sh.ambient(255,0,0); break; }
	    	case gui_Blue				: { sh.fill(0,0,255,alpha); sh.ambient(0,0,255); break; }
	    	case gui_Green				: { sh.fill(0,255,0,alpha); sh.ambient(0,255,0); break; } 
	    	case gui_DarkGray   		: { sh.fill(80,80,80,alpha); sh.ambient(80,80,80); break;}
	    	case gui_DarkRed    		: { sh.fill(120,0,0,alpha); sh.ambient(120,0,0); break;}
	    	case gui_DarkBlue   		: { sh.fill(0,0,120,alpha); sh.ambient(0,0,120); break;}
	    	case gui_DarkGreen  		: { sh.fill(0,120,0,alpha); sh.ambient(0,120,0); break;}
	    	case gui_DarkYellow 		: { sh.fill(120,120,0,alpha); sh.ambient(120,120,0); break;}
	    	case gui_DarkMagenta		: { sh.fill(120,0,120,alpha); sh.ambient(120,0,120); break;}
	    	case gui_DarkCyan   		: { sh.fill(0,120,120,alpha); sh.ambient(0,120,120); break;}		   
	    	case gui_LightGray   		: { sh.fill(200,200,200,alpha); sh.ambient(200,200,200); break;}
	    	case gui_LightRed    		: { sh.fill(255,110,110,alpha); sh.ambient(255,110,110); break;}
	    	case gui_LightBlue   		: { sh.fill(110,110,255,alpha); sh.ambient(110,110,255); break;}
	    	case gui_LightGreen  		: { sh.fill(110,255,110,alpha); sh.ambient(110,255,110); break;}
	    	case gui_LightYellow 		: { sh.fill(255,255,110,alpha); sh.ambient(255,255,110); break;}
	    	case gui_LightMagenta		: { sh.fill(255,110,255,alpha); sh.ambient(255,110,255); break;}
	    	case gui_LightCyan   		: { sh.fill(110,255,255,alpha); sh.ambient(110,255,255); break;}	    	
	    	case gui_Black			 	: { sh.fill(0,0,0,alpha); sh.ambient(0,0,0); break;}//
	    	case gui_TransBlack  	 	: { sh.fill(0x00010100); sh.ambient(0,0,0); break;}//	have to use hex so that alpha val is not lost    	
	    	case gui_FaintGray 		 	: { sh.fill(77,77,77,77); sh.ambient(77,77,77); break;}//
	    	case gui_FaintRed 	 	 	: { sh.fill(110,0,0,alpha/2); sh.ambient(110,0,0); break;}//
	    	case gui_FaintBlue 	 	 	: { sh.fill(0,0,110,alpha/2); sh.ambient(0,0,110); break;}//
	    	case gui_FaintGreen 	 	: { sh.fill(0,110,0,alpha/2); sh.ambient(0,110,0); break;}//
	    	case gui_FaintYellow 	 	: { sh.fill(110,110,0,alpha/2); sh.ambient(110,110,0); break;}//
	    	case gui_FaintCyan  	 	: { sh.fill(0,110,110,alpha/2); sh.ambient(0,110,110); break;}//
	    	case gui_FaintMagenta  	 	: { sh.fill(110,0,110,alpha/2); sh.ambient(110,0,110); break;}//
	    	case gui_TransGray 	 	 	: { sh.fill(120,120,120,alpha/8); sh.ambient(120,120,120); break;}//
	    	case gui_TransRed 	 	 	: { sh.fill(255,0,0,alpha/2); sh.ambient(255,0,0); break;}//
	    	case gui_TransBlue 	 	 	: { sh.fill(0,0,255,alpha/2); sh.ambient(0,0,255); break;}//
	    	case gui_TransGreen 	 	: { sh.fill(0,255,0,alpha/2); sh.ambient(0,255,0); break;}//
	    	case gui_TransYellow 	 	: { sh.fill(255,255,0,alpha/2); sh.ambient(255,255,0); break;}//
	    	case gui_TransCyan  	 	: { sh.fill(0,255,255,alpha/2); sh.ambient(0,255,255); break;}//
	    	case gui_TransMagenta  	 	: { sh.fill(255,0,255,alpha/2); sh.ambient(255,0,255); break;}//
//	    	case gui_boatBody1 	  		: { sh.fill(110, 65, 30,alpha); sh.ambient(130, 75, 40); 	sh.specular(130, 75, 40); break;}
//	    	case gui_boatBody2 	  		: { sh.fill(0, 0, 0,alpha); 	sh.ambient(222, 222, 222);	sh.specular(255,255,255); break;}
//	    	case gui_boatBody3 	  		: { sh.fill(130, 22, 10,alpha); sh.ambient(255, 100, 100); sh.specular(255, 75, 40); break;}//fill(255, 0, 0,255); 	ambient(130, 75, 40); break;}
//	    	case gui_boatBody4 	  		: { sh.fill(22, 130, 10,alpha); sh.ambient(100, 255, 100); sh.specular(75, 255, 40); break;}//fill(255, 0, 0,255); 	ambient(130, 75, 40); break;}
//	    	case gui_boatBody5 	  		: { sh.fill(22, 10, 130,alpha); sh.ambient(100, 100, 255); sh.specular(75, 40, 255); break;}//fill(255, 0, 0,255); 	ambient(130, 75, 40); break;}
	    	case gui_boatBody1 	  		: { sh.fill(110, 65, 30,alpha); break;}
	    	case gui_boatBody2 	  		: { sh.fill(0, 0, 0, alpha); 	break;}
	    	case gui_boatBody3 	  		: { sh.fill(130, 22, 10,alpha); break;}//fill(255, 0, 0,255); 	ambient(130, 75, 40); break;}
	    	case gui_boatBody4 	  		: { sh.fill(22, 230, 10,alpha); break;}//fill(255, 0, 0,255); 	ambient(130, 75, 40); break;}
	    	case gui_boatBody5 	  		: { sh.fill(22, 10, 130,alpha); break;}//fill(255, 0, 0,255); 	ambient(130, 75, 40); break;}
	    	case gui_boatStrut 	  		: { sh.fill(80, 40, 25, alpha); sh.ambient(130, 75, 40);break;}
	    	
	    	default         			: { sh.fill(255,255,255,alpha); sh.ambient(255,255,255); break; }	    	    	
		}//switch	
	}//setShcolorValFill
	

	public void setColorValFill(int colorVal){ setColorValFill(colorVal,255);}
	public void setColorValFill(int colorVal, int alpha){
		switch (colorVal){
			case gui_rnd				: { fill(random(255),random(255),random(255),alpha);break;}
	    	case gui_White  			: { fill(255,255,255,alpha);break; }
	    	case gui_Gray   			: { fill(120,120,120,alpha); break;}
	    	case gui_Yellow 			: { fill(255,255,0,alpha);break; }
	    	case gui_Cyan   			: { fill(0,255,255,alpha);  break; }
	    	case gui_Magenta			: { fill(255,0,255,alpha);break; }
	    	case gui_Red    			: { fill(255,0,0,alpha); break; }
	    	case gui_Blue				: { fill(0,0,255,alpha); break; }
	    	case gui_Green				: { fill(0,255,0,alpha);  break; } 
	    	case gui_DarkGray   		: { fill(80,80,80,alpha); break;}
	    	case gui_DarkRed    		: { fill(120,0,0,alpha);break;}
	    	case gui_DarkBlue   		: { fill(0,0,120,alpha); break;}
	    	case gui_DarkGreen  		: { fill(0,120,0,alpha); break;}
	    	case gui_DarkYellow 		: { fill(120,120,0,alpha); break;}
	    	case gui_DarkMagenta		: { fill(120,0,120,alpha); break;}
	    	case gui_DarkCyan   		: { fill(0,120,120,alpha); break;}	   
	    	case gui_LightGray   		: { fill(200,200,200,alpha); break;}
	    	case gui_LightRed    		: { fill(255,110,110,alpha); break;}
	    	case gui_LightBlue   		: { fill(110,110,255,alpha); break;}
	    	case gui_LightGreen  		: { fill(110,255,110,alpha); break;}
	    	case gui_LightYellow 		: { fill(255,255,110,alpha); break;}
	    	case gui_LightMagenta		: { fill(255,110,255,alpha); break;}
	    	case gui_LightCyan   		: { fill(110,255,255,alpha); break;}    	
	    	case gui_Black			 	: { fill(0,0,0,alpha);break;}//
			case gui_TransBlack  	 	: { fill(0x00010100);  break;}//	have to use hex so that alpha val is not lost    	
			case gui_FaintGray 		 	: { fill(77,77,77,alpha/3); break;}
			case gui_FaintRed 	 	 	: { fill(110,0,0,alpha/2);  break;}
			case gui_FaintBlue 	 	 	: { fill(0,0,110,alpha/2);  break;}
			case gui_FaintGreen 	 	: { fill(0,110,0,alpha/2);  break;}
			case gui_FaintYellow 	 	: { fill(110,110,0,alpha/2); break;}
			case gui_FaintCyan  	 	: { fill(0,110,110,alpha/2); break;}
			case gui_FaintMagenta  	 	: { fill(110,0,110,alpha/2); break;}
			case gui_TransGray 	 	 	: { fill(120,120,120,alpha/8); break;}//
			case gui_TransRed 	 	 	: { fill(255,0,0,alpha/2);  break;}
			case gui_TransBlue 	 	 	: { fill(0,0,255,alpha/2);  break;}
			case gui_TransGreen 	 	: { fill(0,255,0,alpha/2);  break;}
			case gui_TransYellow 	 	: { fill(255,255,0,alpha/2);break;}
			case gui_TransCyan  	 	: { fill(0,255,255,alpha/2);break;}
			case gui_TransMagenta  	 	: { fill(255,0,255,alpha/2);break;}
			case gui_OffWhite			: { fill(248,248,255,alpha);break; }
	    	case gui_boatBody1 	  		: { fill(110, 65, 30,255); 	break;}
	    	case gui_boatBody2 	  		: { fill(0, 0, 0,255);  break;}
	    	case gui_boatBody3 	  		: { fill(130, 0, 0,255); break;}//fill(255, 0, 0,255); 	ambient(130, 75, 40); break;}
	    	case gui_boatBody4 	  		: { fill(0, 130, 0,255); break;}//fill(255, 0, 0,255); 	ambient(130, 75, 40); break;}
	    	case gui_boatBody5 	  		: { fill(50, 0, 130,255); break;}//fill(255, 0, 0,255); 	ambient(130, 75, 40); break;}
	    	case gui_boatStrut 	  		: { fill(80, 40, 25, 255); break;}	    	
	    	default         			: { fill(255,255,255,255);  break; }	    	    	
		}//switch	
	}//setcolorValFill
	
	public void setColorValStroke(int colorVal){ setColorValStroke(colorVal, 255);}
	public void setColorValStroke(int colorVal, int alpha){
		switch (colorVal){
	    	case gui_White  	 	    : { stroke(255,255,255,alpha); break; }
	    	case gui_Gray   	 	    : { stroke(120,120,120,alpha); break;}
	    	case gui_Yellow      	    : { stroke(255,255,0,alpha); break; }
	    	case gui_Cyan   	 	    : { stroke(0,255,255,alpha); break; }
	    	case gui_Magenta	 	    : { stroke(255,0,255,alpha);  break; }
	    	case gui_Red    	 	    : { stroke(255,120,120,alpha); break; }
	    	case gui_Blue		 	    : { stroke(120,120,255,alpha); break; }
	    	case gui_Green		 	    : { stroke(120,255,120,alpha); break; }
	    	case gui_DarkGray    	    : { stroke(80,80,80,alpha); break; }
	    	case gui_DarkRed     	    : { stroke(120,0,0,alpha); break; }
	    	case gui_DarkBlue    	    : { stroke(0,0,120,alpha); break; }
	    	case gui_DarkGreen   	    : { stroke(0,120,0,alpha); break; }
	    	case gui_DarkYellow  	    : { stroke(120,120,0,alpha); break; }
	    	case gui_DarkMagenta 	    : { stroke(120,0,120,alpha); break; }
	    	case gui_DarkCyan    	    : { stroke(0,120,120,alpha); break; }	   
	    	case gui_LightGray   	    : { stroke(200,200,200,alpha); break;}
	    	case gui_LightRed    	    : { stroke(255,110,110,alpha); break;}
	    	case gui_LightBlue   	    : { stroke(110,110,255,alpha); break;}
	    	case gui_LightGreen  	    : { stroke(110,255,110,alpha); break;}
	    	case gui_LightYellow 	    : { stroke(255,255,110,alpha); break;}
	    	case gui_LightMagenta	    : { stroke(255,110,255,alpha); break;}
	    	case gui_LightCyan   		: { stroke(110,255,255,alpha); break;}		   
	    	case gui_Black				: { stroke(0,0,0,alpha); break;}
	    	case gui_TransBlack  		: { stroke(1,1,1,1); break;}	    	
	    	case gui_FaintGray 			: { stroke(120,120,120,250); break;}
	    	case gui_FaintRed 	 		: { stroke(110,0,0,alpha); break;}
	    	case gui_FaintBlue 	 		: { stroke(0,0,110,alpha); break;}
	    	case gui_FaintGreen 		: { stroke(0,110,0,alpha); break;}
	    	case gui_FaintYellow 		: { stroke(110,110,0,alpha); break;}
	    	case gui_FaintCyan  		: { stroke(0,110,110,alpha); break;}
	    	case gui_FaintMagenta  		: { stroke(110,0,110,alpha); break;}
	    	case gui_TransGray 	 		: { stroke(150,150,150,alpha/4); break;}
	    	case gui_TransRed 	 		: { stroke(255,0,0,alpha/2); break;}
	    	case gui_TransBlue 	 		: { stroke(0,0,255,alpha/2); break;}
	    	case gui_TransGreen 		: { stroke(0,255,0,alpha/2); break;}
	    	case gui_TransYellow 		: { stroke(255,255,0,alpha/2); break;}
	    	case gui_TransCyan  		: { stroke(0,255,255,alpha/2); break;}
	    	case gui_TransMagenta  		: { stroke(255,0,255,alpha/2); break;}
	    	case gui_OffWhite			: { stroke(248,248,255,alpha);break; }
	    	case gui_boatBody1 	  		: { stroke(80, 40, 25,alpha); break;}
	    	case gui_boatBody2 	  		: { stroke(0, 0, 0,alpha); break;}
	    	case gui_boatBody3 	  		: { stroke(40, 0, 0,alpha); break;}//stroke(222, 0, 0); break;}
	    	case gui_boatBody4 	  		: { stroke(0, 40, 0,alpha); break;}//stroke(222, 0, 0); break;}
	    	case gui_boatBody5 	  		: { stroke(0, 0, 40,alpha); break;}//stroke(222, 0, 0); break;}
	    	case gui_boatStrut 	  		: { stroke(80, 40, 25,alpha); break;}
	    	default         			: { stroke(55,55,255,alpha); break; }
		}//switch	
	}//setcolorValStroke		

	
	//returns one of 30 predefined colors as an array (to support alpha)
	//processing's hex clr format is 0xAARRGGBB
	public int getHexClr(int colorVal){
		int[] clrAra = getClr(colorVal, 255);
		int res = (byte)255 << 24;
		for(int i =2;i>=0;--i){//RGB
			res+= ((byte)clrAra[i] << (i*8));
		}
		return res;
	}
	public int[] getClr(int colorVal){		return getClr(colorVal, 255);	}//getClr
	public int[] getClr(int colorVal, int alpha){
		switch (colorVal){
		case gui_Gray   		         : { return new int[] {120,120,120,alpha}; }
		case gui_White  		         : { return new int[] {255,255,255,alpha}; }
		case gui_Yellow 		         : { return new int[] {255,255,0,alpha}; }
		case gui_Cyan   		         : { return new int[] {0,255,255,alpha};} 
		case gui_Magenta		         : { return new int[] {255,0,255,alpha};}  
		case gui_Red    		         : { return new int[] {255,0,0,alpha};} 
		case gui_Blue			         : { return new int[] {0,0,255,alpha};}
		case gui_Green			         : { return new int[] {0,255,0,alpha};}  
		case gui_DarkGray   	         : { return new int[] {80,80,80,alpha};}
		case gui_DarkRed    	         : { return new int[] {120,0,0,alpha};}
		case gui_DarkBlue  	 	         : { return new int[] {0,0,120,alpha};}
		case gui_DarkGreen  	         : { return new int[] {0,120,0,alpha};}
		case gui_DarkYellow 	         : { return new int[] {120,120,0,alpha};}
		case gui_DarkMagenta	         : { return new int[] {120,0,120,alpha};}
		case gui_DarkCyan   	         : { return new int[] {0,120,120,alpha};}	   
		case gui_LightGray   	         : { return new int[] {200,200,200,alpha};}
		case gui_LightRed    	         : { return new int[] {255,110,110,alpha};}
		case gui_LightBlue   	         : { return new int[] {110,110,255,alpha};}
		case gui_LightGreen  	         : { return new int[] {110,255,110,alpha};}
		case gui_LightYellow 	         : { return new int[] {255,255,110,alpha};}
		case gui_LightMagenta	         : { return new int[] {255,110,255,alpha};}
		case gui_LightCyan   	         : { return new int[] {110,255,255,alpha};}
		case gui_Black			         : { return new int[] {0,0,0,alpha};}
		case gui_FaintGray 		         : { return new int[] {110,110,110,alpha};}
		case gui_FaintRed 	 	         : { return new int[] {110,0,0,alpha};}
		case gui_FaintBlue 	 	         : { return new int[] {0,0,110,alpha};}
		case gui_FaintGreen 	         : { return new int[] {0,110,0,alpha};}
		case gui_FaintYellow 	         : { return new int[] {110,110,0,alpha};}
		case gui_FaintCyan  	         : { return new int[] {0,110,110,alpha};}
		case gui_FaintMagenta  	         : { return new int[] {110,0,110,alpha};}    	
		case gui_TransBlack  	         : { return new int[] {1,1,1,alpha/2};}  	
		case gui_TransGray  	         : { return new int[] {110,110,110,alpha/2};}
		case gui_TransLtGray  	         : { return new int[] {180,180,180,alpha/2};}
		case gui_TransRed  	         	 : { return new int[] {110,0,0,alpha/2};}
		case gui_TransBlue  	         : { return new int[] {0,0,110,alpha/2};}
		case gui_TransGreen  	         : { return new int[] {0,110,0,alpha/2};}
		case gui_TransYellow  	         : { return new int[] {110,110,0,alpha/2};}
		case gui_TransCyan  	         : { return new int[] {0,110,110,alpha/2};}
		case gui_TransMagenta  	         : { return new int[] {110,0,110,alpha/2};}	
		case gui_TransWhite  	         : { return new int[] {220,220,220,alpha/2};}	
		case gui_OffWhite				 : { return new int[] {255,255,235,alpha};}		
    	case gui_boatBody1 	  			: {return new int[] {80, 40, 25,alpha};}
    	case gui_boatBody2 	  			: {return new int[] {0, 0, 0,alpha};}
    	case gui_boatBody3 	  			: {return new int[] {40, 0, 0,alpha};}
    	case gui_boatBody4 	  			: {return new int[] {0, 80, 0,alpha};}
    	case gui_boatBody5 	  			: {return new int[] {40, 0, 80,alpha};}
    	case gui_boatStrut 	  			: {return new int[] {80, 40, 25,alpha};}
    	default         		         : { return new int[] {255,255,255,alpha};}    
		}//switch
	}//getClr	
	
	//used to generate random color
	public static final int gui_rnd = -1;
	//color indexes
	public static final int gui_Black 	= 0;
	public static final int gui_White 	= 1;	
	public static final int gui_Gray 	= 2;
	
	public static final int gui_Red 	= 3;
	public static final int gui_Blue 	= 4;
	public static final int gui_Green 	= 5;
	public static final int gui_Yellow 	= 6;
	public static final int gui_Cyan 	= 7;
	public static final int gui_Magenta = 8;
	
	public static final int gui_LightRed = 9;
	public static final int gui_LightBlue = 10;
	public static final int gui_LightGreen = 11;
	public static final int gui_LightYellow = 12;
	public static final int gui_LightCyan = 13;
	public static final int gui_LightMagenta = 14;
	public static final int gui_LightGray = 15;

	public static final int gui_DarkCyan = 16;
	public static final int gui_DarkYellow = 17;
	public static final int gui_DarkGreen = 18;
	public static final int gui_DarkBlue = 19;
	public static final int gui_DarkRed = 20;
	public static final int gui_DarkGray = 21;
	public static final int gui_DarkMagenta = 22;
	
	public static final int gui_FaintGray = 23;
	public static final int gui_FaintRed = 24;
	public static final int gui_FaintBlue = 25;
	public static final int gui_FaintGreen = 26;
	public static final int gui_FaintYellow = 27;
	public static final int gui_FaintCyan = 28;
	public static final int gui_FaintMagenta = 29;
	
	public static final int gui_TransBlack = 30;
	public static final int gui_TransGray = 31;
	public static final int gui_TransMagenta = 32;	
	public static final int gui_TransLtGray = 33;
	public static final int gui_TransRed = 34;
	public static final int gui_TransBlue = 35;
	public static final int gui_TransGreen = 36;
	public static final int gui_TransYellow = 37;
	public static final int gui_TransCyan = 38;	
	public static final int gui_TransWhite = 39;	
	public static final int gui_OffWhite = 40;
	
	public static final int gui_boatBody1 = 41;
	public static final int gui_boatBody2 = 42;
	public static final int gui_boatBody3 = 43;	
	public static final int gui_boatBody4 = 44;	
	public static final int gui_boatBody5 = 45;	

	public static final int gui_boatStrut = 46;
	
}//papplet class





