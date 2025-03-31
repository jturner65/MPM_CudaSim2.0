package MPM_SimMain.sim;

import java.util.HashMap;
import java.util.Map;


/**
 * This enum specifies the various simulation initialization that should be performed, 
 * depending on what UI values have changed. The higher values have higher priority
 */
public enum SimResetProcess{
	DoNothing,RemakeKernel,ResetSim,RebuildSim;
	private String[] _valExplanation = new String[] {
			"No need to modify simulation in any way due to UI Input",
			"Remap Solver/Kernel functions but do not modify simulation environment",
			"Reset existing simulation points using some new values",
			"Rebuild simulation environment entirely"};
	private static String[] _valName = new String[] {"No Modification","Remap Kernel","Reset Existing Sim","Rebuild Simulation"};
	public static String[] getListOfValNames() {return _valName;}
	private static Map<Integer, SimResetProcess> map = new HashMap<Integer, SimResetProcess>(); 
	static { for (SimResetProcess enumV : SimResetProcess.values()) { map.put(enumV.ordinal(), enumV);}}
	public int getVal(){return ordinal();}
	public static SimResetProcess getEnumByIndex(int idx){return map.get(idx);}
	public static SimResetProcess getEnumFromValue(int idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
	public String getName() {return _valName[ordinal()];}
	@Override
   public String toString() { return ""+name() + ":"+_valExplanation[ordinal()]; }	
}//SimResetProcess