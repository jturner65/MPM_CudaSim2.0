package MPM_Simulation_PKG;


import java.util.HashMap;

import MPM.CPUSim.ui.MPM_CPUSimWindow;
import MPM.CudaSim.ui.MPM_CudaSimWindow;
import base_Render_Interface.IRenderInterface;
import base_UI_Objects.GUI_AppManager;
import base_Utils_Objects.io.messaging.MsgCodes;

/**
* MPM Snow Simulation in CUDA
* @author john turner
*/

public class MPM_SimMain extends GUI_AppManager {
    
    public final String prjNmShrt = "MPM_SnowSim_cpu_cuda";
    public final String prjNmLong = "MPM Snow Simulation Multi-Threaded CPU and CUDA/JCUDA 9.2"; 
    public final String projDesc = "Simulate numerous snow balls colliding using Material Point Method solved via MT CPU solver or a CUDA kernel.";
    
    private String bkSkyBox = "winter.jpg";
    //bground color
    private final int[] bground = new int[]{0,0,0,255};                
    
    private final int gridDim = 1500;

    /**
     * idx's in dispWinFrames for each window - 0 is always left side menu window
     * Side menu is dispMenuIDX == 0
     */
    private static final int
        dispMPMCudaWinIDX = 1,
        dispMPMCudaWin2IDX = 2,
        dispMPMCPUWinIDX = 3;
    
    /**
     * # of visible windows including side menu (always at least 1 for side menu)
     */
    private static final int numVisWins = 4;
    
///////////////
//CODE STARTS
///////////////    
//////////////////////////////////////////////// code
    //needs main to run project - do not modify this code in any way
    public static void main(String[] passedArgs) {        
        MPM_SimMain me = new MPM_SimMain();
        MPM_SimMain.invokeProcessingMain(me, passedArgs);
    }//main    
    
    protected MPM_SimMain(){super();}

    @Override
    protected boolean showMachineData() {return true;}
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
     *      (Current settings in ProcessingRenderer)     
     *      strokeCap(PROJECT);
     *      textSize(txtSz);
     *      textureMode(NORMAL);            
     *      rectMode(CORNER);    
     *      sphereDetail(4);     * 
     */
    @Override
    protected void setupAppDims_Indiv() {
        //Set grid dimensions
        setDesired3DGridDims(gridDim);        
    }
    @Override
    protected boolean getUseSkyboxBKGnd(int winIdx) {    return winIdx == dispMPMCudaWin2IDX;}
    @Override
    protected String getSkyboxFilename(int winIdx) {    return bkSkyBox;}
    @Override
    protected int[] getBackgroundColor(int winIdx) {return bground;}
    @Override
    protected int getNumDispWindows() {    return numVisWins;    }
    
    /**
     * whether or not we want to restrict window size on widescreen monitors
     * 
     * @return 0 - use monitor size regardless
     *             1 - use smaller dim to be determine window 
     *             2+ - TBD
     */
    @Override
    protected int setAppWindowDimRestrictions() {    return 1;}    
    
    @Override
    public String getPrjNmShrt() {        return prjNmShrt;}
    @Override
    public String getPrjNmLong() {        return prjNmLong;}
    @Override
    public String getPrjDescr() {        return projDesc;}    
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
        setBaseFlagToShow_showStatusBar(true);
        setBaseFlagToShow_showDrawableCanvas(false);
    }
    @Override
    protected void initAllDispWindows() {
        showInfo = true;
        //titles and descs, need to be set before sidebar menu is defined
        String[] _winTitles = new String[]{"","Snow Balls!","Snow Balls In The Alps!","CPU Snow Ball"},
                _winDescr = new String[]{"", "Display Colliding Snowballs Simulated via MPM CUDA Solver", "Display Colliding Snowballs Simulated via MPM CUDA Solver", "Display Falling Snowball Simulated via CPU/MT Solver"};

        //instanced window dims when open and closed - only showing 1 open at a time - and init cam vals
        float[][] _floatDims  = getDefaultWinAndCameraDims();    

        //Builds sidebar menu button config - application-wide menu button bar titles and button names
        String[] menuBtnTitles =  new String[]{"Functions 1","Functions 2","Functions 3","Functions 4"};
        String[][] menuBtnNames = new String[][] {    //each must have literals for every button defined in side bar menu, or ignored
            {"Restore Init Vals", "Func 01", "Func 02"},                //row 1
            {"Func 10", "Func 11", "Func 12", "Func 13"},    //row 2
            {"Func 10", "Func 11", "Func 12", "Func 13"},    //row 3
            {"Func 20", "Func 21", "Func 22", "Func 23"},    //row 4        
        };
        String[] dbgBtnNames = new String[]{"Debug 0","Debug 1","Debug 2","Debug 3","Debug 4"};
        buildSideBarMenu(_winTitles, menuBtnTitles, menuBtnNames, dbgBtnNames, true, false);

        initXORWins(new int[]{dispMPMCudaWinIDX, dispMPMCudaWin2IDX, dispMPMCPUWinIDX},new int[]{dispMPMCudaWinIDX, dispMPMCudaWin2IDX, dispMPMCPUWinIDX});
        //define windows
        /**
         *  _winIdx The index in the various window-descriptor arrays for the dispWindow being set
         *  _title string title of this window
         *  _descr string description of this window
         *  _dispFlags Essential flags describing the nature of the dispWindow for idxs : 
         *         0 : dispWinIs3d, 
         *         1 : canDrawInWin; 
         *         2 : canShow3dbox (only supported for 3D); 
         *         3 : canMoveView
         *  _floatVals an array holding float arrays for 
         *                 rectDimOpen(idx 0),
         *                 rectDimClosed(idx 1),
         *                 initCameraVals(idx 2)
         *  _intClrVals and array holding int arrays for
         *                 winFillClr (idx 0),
         *                 winStrkClr (idx 1),
         *                 winTrajFillClr(idx 2),
         *                 winTrajStrkClr(idx 3),
         *                 rtSideFillClr(idx 4),
         *                 rtSideStrkClr(idx 5)
         *  _sceneCenterVal center of scene, for drawing objects (optional)
         *  _initSceneFocusVal initial focus target for camera (optional)
         */        
        int wIdx = dispMPMCudaWinIDX;
        //setInitDispWinVals(wIdx, _dimOpen, _dimClosed,new boolean[]{false,true,true,true}, new int[]{255,245,255,255},new int[]{0,0,0,255},new int[]{180,180,180,255},new int[]{100,100,100,255});         
        setInitDispWinVals(wIdx, _winTitles[wIdx], _winDescr[wIdx], getDfltBoolAra(true), _floatDims,        
                new int[][] {ri.getClr(IRenderInterface.gui_White, 255),ri.getClr(IRenderInterface.gui_White, 255),
                    ri.getClr(IRenderInterface.gui_LightGray, 255),ri.getClr(IRenderInterface.gui_FaintGray, 255),
                    ri.getClr(IRenderInterface.gui_Black, 200),ri.getClr(IRenderInterface.gui_White, 255)});

        setDispWindow(wIdx, new MPM_CudaSimWindow(ri, this, wIdx));
        
        wIdx = dispMPMCudaWin2IDX;
        //setInitDispWinVals(wIdx, _dimOpen, _dimClosed,new boolean[]{false,true,true,true}, new int[]{255,245,255,255},new int[]{0,0,0,255},new int[]{180,180,180,255},new int[]{100,100,100,255});         
        setInitDispWinVals(wIdx, _winTitles[wIdx], _winDescr[wIdx], getDfltBoolAra(true), _floatDims,        
                new int[][] {ri.getClr(IRenderInterface.gui_White, 255),ri.getClr(IRenderInterface.gui_Black, 255),
                    ri.getClr(IRenderInterface.gui_LightGray, 255),ri.getClr(IRenderInterface.gui_FaintGray, 255),
                    ri.getClr(IRenderInterface.gui_Black, 200),ri.getClr(IRenderInterface.gui_White, 255)});

        setDispWindow(wIdx, new MPM_CudaSimWindow(ri, this, wIdx));    
        
        wIdx = dispMPMCPUWinIDX;
        //setInitDispWinVals(wIdx, _dimOpen, _dimClosed,new boolean[]{false,true,true,true}, new int[]{255,245,255,255},new int[]{0,0,0,255},new int[]{180,180,180,255},new int[]{100,100,100,255});         
        setInitDispWinVals(wIdx, _winTitles[wIdx], _winDescr[wIdx], getDfltBoolAra(true), _floatDims,        
                new int[][] {ri.getClr(IRenderInterface.gui_White, 255),ri.getClr(IRenderInterface.gui_White, 255),
                    ri.getClr(IRenderInterface.gui_LightGray, 255),ri.getClr(IRenderInterface.gui_FaintGray, 255),
                    ri.getClr(IRenderInterface.gui_Black, 200),ri.getClr(IRenderInterface.gui_White, 255)});
        
        setDispWindow(wIdx, new MPM_CPUSimWindow(ri, this, wIdx));
    }//initAllDispWindows
    
    /**
     * Map indexed by window ID, holding an array of the titles (idx 0) and descriptions (idx 1) for every sub window
     * return null if none exist, and only put an entry in the map if one exists for that window
     * @return
     */
    @Override
    protected final HashMap<Integer, String[]> getSubWindowTitles(){ return null;}
    
    @Override
    protected void initOnce_Indiv() {
        //which objects to initially show
        setWinVisFlag(dispMPMCudaWinIDX, true);    
        setShowStatusBar(true);    
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
            case ' ' : {toggleSimIsRunning(); break;}                            //run sim
            case 'f' : {getCurFocusDispWindow().setInitCamView();break;}                    //reset camera
            case 'a' :
            case 'A' : {toggleSaveAnim();break;}                        //start/stop saving every frame for making into animation
            case 's' :
            case 'S' : {break;}//save(getScreenShotSaveName(prjNmShrt));break;}//save picture of current image            
            default : {    }
        }//switch    
    }
    
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
            case dispMPMCudaWinIDX     : {return menuClickDim;}
            case dispMPMCudaWin2IDX : {return menuClickDim;}
            case dispMPMCPUWinIDX    : {return menuClickDim;}
            default :  return menuClickDim;
            }
    }//getUIRectVals

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
            case dispMPMCudaWinIDX        : {setWinFlagsXOR(dispMPMCudaWinIDX, val); break;}
            case dispMPMCudaWin2IDX        : {setWinFlagsXOR(dispMPMCudaWin2IDX, val); break;}
            case dispMPMCPUWinIDX        : {setWinFlagsXOR(dispMPMCPUWinIDX, val); break;}
            default : {break;}
        }
    }//setFlags  
    
    @Override
    public void setSmoothing() {        ri.setSmoothing(0);}

}//class MPM.BaseSim
