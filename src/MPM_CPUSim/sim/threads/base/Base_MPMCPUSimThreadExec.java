package MPM_CPUSim.sim.threads.base;

import java.util.concurrent.*;

import MPM_CPUSim.sim.base.Base_MPMCPUSim;
import base_Math_Objects.vectorObjs.floats.myVectorf;

//1 thread's worth of execution for executing work
public abstract class Base_MPMCPUSimThreadExec implements Callable<Boolean> {
	//ref to owning sim
	protected Base_MPMCPUSim sim;

	//sim step now running
	protected int simStep;
	//idxs in component array that this thread will map to
	protected int stIDX, endIDX, thIDX;
	/**
	 * Wall friction for collision against sim bounds
	 */
	protected float wallFriction;
	/**
	 * Collider friction
	 */
	protected float colFriction;
	
	/**
	 * Min sim bounds in x,y,z
	 */
	protected myVectorf minSimBnds;
	/**
	 * Max sim bounds in x,y,z
	 */
	protected myVectorf maxSimBnds;
	
	public Base_MPMCPUSimThreadExec(Base_MPMCPUSim _sim, int _stIDX, int _endIDX, int _thIDX) {
		sim=_sim;stIDX = _stIDX; endIDX = _endIDX;thIDX = _thIDX;		
		simStep=0;
		wallFriction = sim.getWallFric();
		colFriction = sim.getCollFric();
		minSimBnds = sim.getMinSimBnds();
		maxSimBnds = sim.getMaxSimBnds();
	}
	//set what step of sim is being executed
	public void setSimStep(int _s) {simStep = _s;}	
			
	//methods to execute different simulations steps - not all thread exec processes will use all step functs
	protected abstract void execSimStep0();
	protected abstract void execSimStep1();
	protected abstract void execSimStep2();
	protected abstract void execSimStep3();
	protected abstract void execSimStep4();
	protected abstract void execSimStep5();
	
	@Override
	public Boolean call() throws Exception { 
		switch(simStep) {
			case 0 : {	 execSimStep0();return true;		}
			case 1 : {	 execSimStep1();return true;		}
			case 2 : {	 execSimStep2();return true;		}
			case 3 : {	 execSimStep3();return true;		}
			case 4 : {	 execSimStep4();return true;		}
			case 5 : {	 execSimStep5();return true;		}
			default : { return true;}
		}
	}
	
	
	//cube wall normals
	protected final myVectorf[] wallNorms = new myVectorf[] {
			new myVectorf( 1.0, 0.0, 0.0),new myVectorf(-1.0, 0.0, 0.0),
			new myVectorf( 0.0, 1.0, 0.0),new myVectorf( 0.0,-1.0, 0.0),
			new myVectorf( 0.0, 0.0, 1.0),new myVectorf( 0.0, 0.0,-1.0)
	};
	
	/**
	 * calc collision velocity if collision is going to occur - norm is collider normal
	 * @param fricCoeff
	 * @param velocity
	 * @param colNorm
	 * @return
	 */
	protected myVectorf calcCollVel(float fricCoeff, myVectorf velocity, myVectorf colNorm) {		
		float velNormDir=myVectorf._dot(velocity, colNorm);
		if(velNormDir >=0){return velocity;}
		//velocity in opposite direction of normal if dot prod is <0
		myVectorf velTanDir = myVectorf._sub(velocity, myVectorf._mult(colNorm,velNormDir));
		float fricNormVel = fricCoeff*velNormDir;
		if(velTanDir.magn <= -fricNormVel) {//coulomb friction law
			//no bounce, so just stops
			return new myVectorf(0.0,0.0,0.0);
		}		
		return myVectorf._add(velTanDir,(myVectorf._normalize(velTanDir))._mult(fricNormVel));
		 
	}//calcCollVel
	
	//collision detection for wall collisions, returns idx in wallNorms of collision
	protected int checkWallCollision(myVectorf pos) {
		if(pos.x<=minSimBnds.x) {		return 0;	} else if(pos.x>=maxSimBnds.x) {	return 1;	}			
		if(pos.y<=minSimBnds.y) {		return 2;	} else if(pos.y>=maxSimBnds.y) {	return 3;	}
		if(pos.z<=minSimBnds.z) {		return 4;	} else if(pos.z>=maxSimBnds.z) {	return 5;	}
		return -1;
	}//checkWallCollision
	
	//looking ahead at position, to return new velocity that responds to collision if position has collided
	protected myVectorf applyCollisions(myVectorf pos, myVectorf velocity) {
		//check collider collisions
		myVectorf sphereColNorm = sim.checkColliderCollision(pos);
		if(sphereColNorm.magn != 0) {//colliding with central collider			
			return calcCollVel(colFriction, velocity, sphereColNorm);
			
		} else {//not colliding with sphere collider, check walls	
			//check wall collisions
			int colType = checkWallCollision(pos);
			if(-1==colType) {return velocity;}
			//calc collision velocity if collision is going to occur
			return calcCollVel(wallFriction, velocity, wallNorms[colType]);
		}
	}//applyCollisions

	
}//class simThreadExec

//1 thread's worth of execution for building particles

