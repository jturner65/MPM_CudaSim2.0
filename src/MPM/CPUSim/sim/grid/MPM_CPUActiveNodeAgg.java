package MPM.CPUSim.sim.grid;

import base_Math_Objects.vectorObjs.floats.myVectorf;
/**
 * Class to hold an active node's aggregation construct - there should be 1 of these per node.
 * This will hold data that is calculated per-thread for the node, which will then be aggregated 
 * once the calculation is finished.  this is cheaper than putting these data structures in the 
 * grid node class, particularly since the grid will be so sparse 
 * this object will be made and put in the active nodes conccurenthashmap array
 */
public class MPM_CPUActiveNodeAgg {
	public MPM_CPUGridNode node;
	//aggregate force array(per thd for each _particle thread_) for a particular thread's worth of weighted particle forces
	private myVectorf[] aggForces;
	//aggregate force(per thd for each _particle thread_) for a particular thread's worth of weighted particle masses
	private float[] aggMass;
	//aggregate force(per thd for each _particle thread_) for a particular thread's worth of weighted particle velocity contributions (really momentum)
	private myVectorf[] aggVel;
	//size of arrays == # of threads available
	public static int araSz;
	
	//DO NOT NEED TO CLEAR THESE VALUES BECAUSE THIS OBJECT IS TOSSED AFTER 1 USE
	
	//# threads here is making space for each -particle- thread. should always be the same as # of grid threads
	public MPM_CPUActiveNodeAgg(MPM_CPUGridNode _node, int numThds) {
		node = _node;araSz = numThds;
		aggForces = new myVectorf[araSz];
		aggMass = new float[araSz];
		aggVel = new myVectorf[araSz];
		for(int i=0;i<araSz;++i) {
			aggForces[i] = new myVectorf();
			aggVel[i] = new myVectorf();
		}
	}
	
	public void addPartFrc(myVectorf _frc, int _thd) {aggForces[_thd]._add(_frc);}
	//thd is particle thread id
	public void addPartMassAndWtdVel(float _m, myVectorf _vel, int _thd) {
		aggMass[_thd] += _m;
		aggVel[_thd]._add(_vel);
	}//addMassAndVel
	
	//aggregate the collected forces for this node
	public void aggregateForces() {
		myVectorf tmpFrces = new myVectorf();
		for(int i=0;i<aggForces.length;++i) {
			tmpFrces._add(aggForces[i]);
		}
		//we subtract the aggregated forces from the existing force in the node.
		node.forces._sub(tmpFrces);
	}
	
	public void aggregateMassAndVel() {
		float vx = 0, vy = 0, vz = 0, tmpM = 0;
		for(int i=0;i<aggVel.length;++i) {
			tmpM += aggMass[i];
			vx += aggVel[i].x;
			vy += aggVel[i].y;
			vz += aggVel[i].z;
		}
		//we add the mass and velocity to the existing mass and velocity in the node (TODO check if they are always 0, we can just set them)
		node.addMass(tmpM);
		node.addWtdVelocity(new myVectorf(vx,vy,vz));
	}
	
}//activeNodeAgg