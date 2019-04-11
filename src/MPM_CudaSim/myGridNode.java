package MPM_CudaSim;

import org.jblas.FloatMatrix;

//class holding a single node of the grid structure encapsulating the material.  this is projected to with each particle's quantities
public class myGridNode {
	public final int ID;
	public static int IDgen = 0;
	
	private float mass;
	private myVectorf velocity;
	public myVectorf newvelocity, accels;
	public FloatMatrix forces;	
	//indices
	public int[] idxs;
	//position in world
	myVectorf pos;
	
	public myGridNode(int[] _idxs, myVectorf _pos) {
		ID = IDgen++;
		idxs=_idxs;
		pos =_pos;
		velocity = new myVectorf();
		newvelocity = new myVectorf();
		forces = FloatMatrix.zeros(3);
		accels = new myVectorf();
		reset();
	}//myGridNode
	
	//call before each sim cycle
	public void reset() {
		mass = 0;
		accels.set(0,0,0);
		forces = FloatMatrix.zeros(3);//new float[] {0,0,0};
		velocity.set(0,0,0);
		newvelocity.set(0,0,0);
	}//reset	

	public boolean hasMass() {return mass>0;}
	
	public void addMass(float _m) {
		mass += _m;
	}
	
	public void addWtdVelocity(myVectorf _newVel) {
		velocity._add(_newVel);
	}
	
	public myVectorf getWtDiffVel(float wt) {
		return myVectorf._sub(newvelocity, velocity)._mult(wt);
	}
	
	
	//get particle's contribution to weight
	public float getWtdMass(float _pWt) {
		return mass * _pWt;
	}
	
	//calculate grid velocity
	public void calcVel(float deltaT, FloatMatrix gravBase) {
		//FloatMatrix gravity = gravBase.mul(mass);
		forces.addi(gravBase.mul(mass));
		FloatMatrix toAdd = forces.mul(deltaT / mass);
		myVectorf toAddV = new myVectorf(toAdd.get(0), toAdd.get(1), toAdd.get(2));
		accels=myVectorf._div(toAddV,deltaT);
		velocity._div(mass);
		newvelocity = myVectorf._add(velocity, toAddV);		
	}//calcVel
		
	
	public void drawMe(MPM_SimMain pa, boolean showVel, boolean showAccel, boolean showMass) {
		pa.pushMatrix();pa.pushStyle();
		pa.translate(pos);
		if(showVel && velocity.magn>0) {
			pa.stroke(255,100,60,255);
			pa.line(0,0,0,velocity.x, velocity.y, velocity.z);
		}
		if(showAccel && accels.magn>0) {
			pa.stroke(0,255,0,255);
			pa.line(0,0,0,accels.x, accels.y, accels.z);
		}
		if(showMass && mass>0) {
			pa.stroke(0,255,255,255);
			pa.line(0,0,0,0,0,mass*.1f);
		}
		pa.popStyle();pa.popMatrix();
		
	}//drawMe

}//myGridNode
