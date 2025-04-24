package MPM.CPUSim.sim.threads;

import MPM.CPUSim.sim.base.Base_MPMCPUSim;
import MPM.CPUSim.sim.threads.base.Base_MPMCPUSimThdBldr;

/**
 * class to build a sphere of particles of specified size, dispersed amongst specified # of threads
 */

public class MPM_CPUPartBuilder extends Base_MPMCPUSimThdBldr {
	
	public MPM_CPUPartBuilder(Base_MPMCPUSim _sim, int _numThds) {super(_sim, _numThds);}
	
	@Override
	public void initSimBuilder() {	
		initWrkrIDXs(sim.parts.length, "PartBuilder");
	}//initSimBuilder	
	
	//build arraylist of thread worker callables
	@Override
	protected void initWrkrs() {
		//build thread list structure - this won't change unless simulation structure will change, so on reinitialization should reuse these
		for(int i =0;i<numThds;++i) {
			//System.out.println("i:"+i+" | start idx : " + stIDXs[i] + " | end : " + endIDXs[i]);
			// TODO Finish this conversion
			callExecs.add(new MPM_CPUPartBuildWrkr(sim, stIDXs[i],endIDXs[i], i));
		}	
	}//initWrkrs

	//after all threads are finished this method is called, with specific code to execute depending on which sim step was run
	@Override
	protected void endRun() {
		switch(simStep) {//execute multiple different steps within this environment
			case 0 : {//initialization of all particles		
//				sim.setSimPartsAreBuilt();
//				sim.resetSimEnd();
				break;}
			case 1 : {//project to grid				
				//once done, need to aggregate mass and vels in every grid node 
		      	sim.gridThdMgr.setSimStep(2);//aggregate mass and vels
		      	sim.gridThdMgr.synchRun();					
				break;}
			case 2 : {//compute grid forces from particles
				//once done, need to aggregate forces in every grid node 
				sim.gridThdMgr.setSimStep(3);//aggregate forces sim step
				sim.gridThdMgr.synchRun();				      	
				break;}
			case 3 :{
				break;}
			case 4 : {//NOP				
				break;}
			default : {}
		}
	}//endRun

}//class myPartBuilder