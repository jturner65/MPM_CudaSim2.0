package MPM_CudaSim.utils;

import java.util.TreeMap;

import base_UI_Objects.windowUI.base.base_UpdateFromUIData;
import base_UI_Objects.windowUI.base.myDispWindow;

public class MPM_SimUpdateFromUIData extends base_UpdateFromUIData {

	public MPM_SimUpdateFromUIData(myDispWindow _win) {		super(_win);	}
	public MPM_SimUpdateFromUIData(myDispWindow _win, TreeMap<Integer, Integer> _iVals, TreeMap<Integer, Float> _fVals,
			TreeMap<Integer, Boolean> _bVals) {
		super(_win, _iVals, _fVals, _bVals);
	}
	public MPM_SimUpdateFromUIData(base_UpdateFromUIData _otr) {		super(_otr);	}


	/**
	 * access app-specific ints
	 */
	
	
	/**
	 * access app-specific floats
	 */
	
	
	/**
	 * access app-specific booleans
	 */
		

	
}//class MPM_SimUpdateFromUIData
