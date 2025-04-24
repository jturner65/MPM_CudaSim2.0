package MPM.CPUSim.sim.grid;

import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import base_Render_Interface.IRenderInterface;

//class holding a single node of the grid structure encapsulating the material.  this is projected to with each particle's quantities
public class MPM_CPUGridNode {
	public final int ID;
	public static int IDgen = 0;
	
	public float mass;
	public myVectorf velocity;
	public myVectorf newvelocity, accels;
	public myVectorf forces;	
	//indices
	public int[] idxs;
	//position in world
	public myPointf pos;
	
	public MPM_CPUGridNode(int[] _idxs, myPointf _pos) {
		ID = IDgen++;
		idxs=_idxs;
		pos =_pos;
		velocity = new myVectorf();
		newvelocity = new myVectorf();
		forces = new myVectorf();
		accels = new myVectorf();
		reset();
	}//myGridNode
	
	//call before each sim cycle
	public void reset() {
		mass = 0;
		accels.set(0,0,0);
		forces.set(0,0,0);
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
	public void calcVel(float deltaT, myVectorf gravBase) {
		forces._add(myVectorf._mult(gravBase, mass));
		myVectorf toAddV = myVectorf._mult(forces, deltaT / mass);
		accels = myVectorf._div(toAddV,deltaT);
		velocity._div(mass);
		newvelocity = myVectorf._add(velocity, toAddV);		
	}//calcVel
		
	public void drawMe(IRenderInterface pa) {
		pa.pushMatState();
		pa.translate(pos);
		pa.setStrokeWt(1.0f);
		pa.drawSphere(2.0f);
		pa.popMatState();
		
	}//drawMe
	
	public void drawVel(IRenderInterface pa) {
		if(velocity.magn>0) {
			pa.pushMatState();
			pa.translate(pos);
			pa.setStroke(255,100,60,255);
			pa.drawLine(0,0,0,velocity.x, velocity.y, velocity.z);
			pa.popMatState();
		}
	}

	public void drawAccel(IRenderInterface pa) {
		if(accels.magn>0) {
			pa.pushMatState();
			pa.translate(pos);
			pa.setStroke(0,255,0,255);
			pa.drawLine(0,0,0,accels.x, accels.y, accels.z);	
			pa.popMatState();
		}
	}
	
	public void drawMass(IRenderInterface pa) {
		if(mass>0) {
			pa.pushMatState();
			pa.translate(pos);
			pa.setStroke(0,255,255,255);
			pa.drawLine(0,0,0,0,0,mass*.1f);		
			pa.popMatState();
		}
	}
	
}//myGridNode
