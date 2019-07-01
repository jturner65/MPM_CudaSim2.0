package MPM_CudaSim;

import base_UI_Objects.my_procApplet;
import base_UI_Objects.windowUI.BaseBarMenu;

/**
 * class to manage buttons used by sidebar window - overrides base setup, allows for custom config
 * @author john
 *
 */
public class mySideBarMenu extends BaseBarMenu {
	
	public static final int 
		btnAuxFunc1Idx = 0,			//aux functionality 1
		btnAuxFunc2Idx = 1,			//aux functionality 2
		btnAuxFunc3Idx = 2,			//aux functionality 3
		btnAuxFunc4Idx = 3,			//load prebuilt maps
		btnDBGSelCmpIdx = 4;			//debug

	public mySideBarMenu(my_procApplet _p, String _n, int _flagIdx, int[] fc, int[] sc, float[] rd, float[] rdClosed, String _winTxt, boolean _canDrawTraj) {
		super(_p, _n, _flagIdx, fc, sc, rd, rdClosed, _winTxt, _canDrawTraj);
	}

	@Override
	protected void initSideBarMenuBtns_Priv() {
		guiBtnRowNames = new String[]{"Raw Data Conversion/Processing","Load Post Proc Data","Console Exec Testing","Load Prebuilt Maps","DEBUG"};//,"File"};

		//names for each row of buttons - idx 1 is name of row
		guiBtnNames = new String[][]{
			//new String[]{pa.winTitles[1], pa.winTitles[2], pa.winTitles[3], pa.winTitles[4]},							//display specific windows - multi-select/ always 
			new String[]{"Func 1","Func 2","Func 3"},						//per-window user functions - momentary
			new String[]{"Func 1","Func 2","Func 3","Func 4"},						//per-window user functions - momentary
			new String[]{"Func 1","Func 2","Func 3","Func 4"},						//2nd row per-window user functions - momentary
			new String[]{"Func 1","Func 2","Func 3","Func 4"},						//2nd row per-window user functions - momentary
			new String[]{"Dbg 1","Dbg 2","Dbg 3","Dbg 4","Dbg 5"},						//DEBUG - momentary
//			new String[]{"Load Txt File","Save Txt File"}							//load an existing score, save an existing score - momentary		
		};
		//default names, to return to if not specified by user
		defaultUIBtnNames = new String[][]{
			//new String[]{pa.winTitles[1], pa.winTitles[2], pa.winTitles[3], pa.winTitles[4]},							//display specific windows - multi-select/ always 
			new String[]{"Func 1","Func 2","Func 3"},					//per-window user functions - momentary
			new String[]{"Func 1","Func 2","Func 3","Func 4"},			//per-window user functions - momentary
			new String[]{"Func 1","Func 2","Func 3","Func 4"},			//per-window user functions - momentary
			new String[]{"Func 1","Func 2","Func 3","Func 4"},			//per-window user functions - momentary
			new String[]{"Dbg 1","Dbg 2","Dbg 3","Dbg 4","Dbg 5"},						//DEBUG - momentary
//			new String[]{"Load Txt File","Save Txt File"}							//load an existing score, save an existing score - momentary		
		};
		//whether buttons are momentary or not (on only while being clicked)
		guiBtnInst = new boolean[][]{
			//new boolean[]{false,false,false,false},         						//display specific windows - multi-select/ always on if sel
			new boolean[]{false,false,false},                   //functionality - momentary
			new boolean[]{false,false,false,false},                   //functionality - momentary
			new boolean[]{false,false,false,false},                   //functionality - momentary
			new boolean[]{false,false,false,false},                   //functionality - momentary
			new boolean[]{false,false,false,false,false},                   		//debug - momentary
//			new boolean[]{true,true},			              			//load an existing score, save an existing score - momentary	
		};		
		//whether buttons are waiting for processing to complete (for non-momentary buttons)
		guiBtnWaitForProc = new boolean[][]{
			//new boolean[]{false,false,false,false},         						//display specific windows - multi-select/ always on if sel
			new boolean[]{false,false,false},                   //functionality - momentary
			new boolean[]{false,false,false,false},                   //functionality - momentary
			new boolean[]{false,false,false,false},                   //functionality - momentary
			new boolean[]{false,false,false,false},                   //functionality - momentary
			new boolean[]{false,false,false,false,false},                   		//debug - momentary
//			new boolean[]{false,false},			              			//load an existing score, save an existing score - momentary	
		};			
		
		//whether buttons are disabled(-1), enabled but not clicked/on (0), or enabled and on/clicked(1)
		guiBtnSt = new int[][]{
			//new int[]{0,0,1,0},                    					//display specific windows - multi-select/ always on if sel
			new int[]{0,0,0},                   					//debug - momentary
			new int[]{0,0,0,0},                   					//debug - momentary
			new int[]{0,0,0,0},                   					//debug - momentary
			new int[]{0,0,0,0},                   					//debug - momentary
			new int[]{0,0,0,0,0},                   					//debug - momentary
//			new int[]{0,0}			              					//load an existing score, save an existing score - momentary	
		};	}

	@Override
	public void handleButtonClick(int row, int col){
		int val = guiBtnSt[row][col];//initial state, before being changed
		guiBtnSt[row][col] = (guiBtnSt[row][col] + 1)%2;//change state
		//if not momentary buttons, set wait for proc to true
		setWaitForProc(row,col);
		switch(row){
			//case btnShowWinIdx 		: {pa.handleShowWin(col, val);break;}
			case btnAuxFunc1Idx 		: //{pa.handleMenuBtnSelCmp(btnAuxFunc1Idx,col, val);break;}
			case btnAuxFunc2Idx 		: //{pa.handleMenuBtnSelCmp(btnAuxFunc2Idx,col, val);break;}
			case btnAuxFunc3Idx 		: //{pa.handleMenuBtnSelCmp(btnAuxFunc2Idx,col, val);break;}
			case btnAuxFunc4Idx			:
			case btnDBGSelCmpIdx  		: {pa.handleMenuBtnSelCmp(row, col, val);break;}//{pa.handleMenuBtnSelCmp(btnDBGSelCmpIdx,col, val);break;}
//			case btnFileCmdIdx 			: {pa.handleFileCmd(btnFileCmdIdx, col, val);break;}
		}				
	}	

	@Override
	protected void launchMenuBtnHndlr() {
		switch(curCustBtnType) {
		case btnAuxFunc1Idx : {break;}//row 1 of menu side bar buttons
		case btnAuxFunc2Idx : {break;}//row 2 of menu side bar buttons
		case btnAuxFunc3Idx : {break;}//row 2 of menu side bar buttons
		case btnAuxFunc4Idx : {break;}
		case btnDBGSelCmpIdx : {break;}//row 3 of menu side bar buttons (debug)			
		}		
	}

}//mySideBarMenu