package MPM_CudaSim.sim;

import java.util.HashMap;
import java.util.Map;


/**
 * This enum specifies the various simulation initialization that should be performed, 
 * depending on what UI values have changed. The higher values have higher priority
 */
public enum SimResetProcess{
	DoNothing(0),RemakeKernel(1),ResetSim(2),RebuildSim(3);
	private int value; 
	private String[] _valExplanation = new String[] {
			"No need to modify simulation in any way due to UI Input",
			"Remap Kernel functions but do not modify simulation environment",
			"Reset existing simulation points using some new values",
			"Rebuild simulation environment entirely"};
	private static String[] _valName = new String[] {"No Modification","Remap Kernel","Reset Existing Sim","Rebuild Simulation"};
	public static String[] getListOfValNames() {return _valName;}
	private static Map<Integer, SimResetProcess> map = new HashMap<Integer, SimResetProcess>(); 
	static { for (SimResetProcess enumV : SimResetProcess.values()) { map.put(enumV.value, enumV);}}
	private SimResetProcess(int _val){value = _val;} 
	public int getVal(){return value;}
	public static SimResetProcess getVal(int idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
	public String getName() {return _valName[value];}
	@Override
    public String toString() { return ""+value + ":"+_valExplanation[value]; }	
}//SimResetProcess