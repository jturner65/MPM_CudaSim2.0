package MPM_CudaSim;

import base_UI_Objects.my_procApplet;
import base_UI_Objects.windowUI.BaseBarMenu;

/**
 * class to manage buttons used by sidebar window - overrides base setup, allows for custom config
 * @author john
 *
 */
public class mySideBarMenu extends BaseBarMenu {

	public mySideBarMenu(my_procApplet _p, String _n, int _flagIdx, int[] fc, int[] sc, float[] rd, float[] rdClosed, String _winTxt) {
		super(_p, _n, _flagIdx, fc, sc, rd, rdClosed, _winTxt);
	}

	@Override
	protected void initSideBarMenuBtns_Priv() {
		/**
		 * set row names for each row of ui action buttons getMouseOverSelBtnNames()
		 * @param _funcRowNames array of names for each row of functional buttons 
		 * @param _numBtnsPerFuncRow array of # of buttons per row of functional buttons
		 * @param _numDbgBtns # of debug buttons
		 * @param _inclWinNames include the names of all the instanced windows
		 * @param _inclMseOvValues include a row for possible mouse over values
		 */
		//protected void setBtnData(String[] _funcRowNames, int[] _numBtnsPerFuncRow, int _numDbgBtns, boolean _inclWinNames, boolean _inclMseOvValues) {

		setBtnData(new String[]{"Functions 1","Functions 2","Functions 3","Functions 4"}, new int[] {3,4,4,4}, 5, false, false);
		
//		guiBtnRowNames = new String[]{"Raw Data Conversion/Processing","Load Post Proc Data","Console Exec Testing","Load Prebuilt Maps","DEBUG"};//,"File"};
//
//		//names for each row of buttons - idx 1 is name of row
//		guiBtnNames = new String[][]{
//			//new String[]{pa.winTitles[1], pa.winTitles[2], pa.winTitles[3], pa.winTitles[4]},							//display specific windows - multi-select/ always 
//			new String[]{"Func 1","Func 2","Func 3"},						//per-window user functions - momentary
//			new String[]{"Func 1","Func 2","Func 3","Func 4"},						//per-window user functions - momentary
//			new String[]{"Func 1","Func 2","Func 3","Func 4"},						//2nd row per-window user functions - momentary
//			new String[]{"Func 1","Func 2","Func 3","Func 4"},						//2nd row per-window user functions - momentary
//			new String[]{"Dbg 1","Dbg 2","Dbg 3","Dbg 4","Dbg 5"},						//DEBUG - momentary
////			new String[]{"Load Txt File","Save Txt File"}							//load an existing score, save an existing score - momentary		
//		};
//		//default names, to return to if not specified by user
//		defaultUIBtnNames = new String[][]{
//			//new String[]{pa.winTitles[1], pa.winTitles[2], pa.winTitles[3], pa.winTitles[4]},							//display specific windows - multi-select/ always 
//			new String[]{"Func 1","Func 2","Func 3"},					//per-window user functions - momentary
//			new String[]{"Func 1","Func 2","Func 3","Func 4"},			//per-window user functions - momentary
//			new String[]{"Func 1","Func 2","Func 3","Func 4"},			//per-window user functions - momentary
//			new String[]{"Func 1","Func 2","Func 3","Func 4"},			//per-window user functions - momentary
//			new String[]{"Dbg 1","Dbg 2","Dbg 3","Dbg 4","Dbg 5"},						//DEBUG - momentary
////			new String[]{"Load Txt File","Save Txt File"}							//load an existing score, save an existing score - momentary		
//		};
//		//whether buttons are momentary or not (on only while being clicked)
//		guiBtnInst = new boolean[][]{
//			//new boolean[]{false,false,false,false},         						//display specific windows - multi-select/ always on if sel
//			new boolean[]{false,false,false},                   //functionality - momentary
//			new boolean[]{false,false,false,false},                   //functionality - momentary
//			new boolean[]{false,false,false,false},                   //functionality - momentary
//			new boolean[]{false,false,false,false},                   //functionality - momentary
//			new boolean[]{false,false,false,false,false},                   		//debug - momentary
////			new boolean[]{true,true},			              			//load an existing score, save an existing score - momentary	
//		};		
//		//whether buttons are waiting for processing to complete (for non-momentary buttons)
//		guiBtnWaitForProc = new boolean[][]{
//			//new boolean[]{false,false,false,false},         						//display specific windows - multi-select/ always on if sel
//			new boolean[]{false,false,false},                   //functionality - momentary
//			new boolean[]{false,false,false,false},                   //functionality - momentary
//			new boolean[]{false,false,false,false},                   //functionality - momentary
//			new boolean[]{false,false,false,false},                   //functionality - momentary
//			new boolean[]{false,false,false,false,false},                   		//debug - momentary
////			new boolean[]{false,false},			              			//load an existing score, save an existing score - momentary	
//		};			
//		
//		//whether buttons are disabled(-1), enabled but not clicked/on (0), or enabled and on/clicked(1)
//		guiBtnSt = new int[][]{
//			//new int[]{0,0,1,0},                    					//display specific windows - multi-select/ always on if sel
//			new int[]{0,0,0},                   					//debug - momentary
//			new int[]{0,0,0,0},                   					//debug - momentary
//			new int[]{0,0,0,0},                   					//debug - momentary
//			new int[]{0,0,0,0},                   					//debug - momentary
//			new int[]{0,0,0,0,0},                   					//debug - momentary
////			new int[]{0,0}			              					//load an existing score, save an existing score - momentary	
//		};	
		}


}//mySideBarMenu