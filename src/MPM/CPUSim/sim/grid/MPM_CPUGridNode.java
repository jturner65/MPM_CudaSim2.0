package MPM.CPUSim.sim.grid;

import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import base_Render_Interface.IGraphicsAppInterface;

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
        // add force due to gravity
        forces._add(myVectorf._mult(gravBase, mass));
        // find velocity by multiplying accel by time
        myVectorf toAddV = myVectorf._mult(forces, deltaT / mass);
        accels = myVectorf._div(toAddV,deltaT);
        velocity._div(mass);
        newvelocity = myVectorf._add(velocity, toAddV);        
    }//calcVel
        
    public void drawMe(IGraphicsAppInterface ri) {
        ri.pushMatState();
        ri.translate(pos);
        ri.setStrokeWt(1.0f);
        ri.drawSphere(2.0f);
        ri.popMatState();
        
    }//drawMe
    
    public void drawVel(IGraphicsAppInterface ri) {
        if(velocity.magn>0) {
            ri.pushMatState();
            ri.translate(pos);
            ri.setStroke(255,100,60,255);
            ri.drawLine(0,0,0,velocity.x, velocity.y, velocity.z);
            ri.popMatState();
        }
    }

    public void drawAccel(IGraphicsAppInterface ri) {
        if(accels.magn>0) {
            ri.pushMatState();
            ri.translate(pos);
            ri.setStroke(0,255,0,255);
            ri.drawLine(0,0,0,accels.x, accels.y, accels.z);    
            ri.popMatState();
        }
    }
    
    public void drawMass(IGraphicsAppInterface ri) {
        if(mass>0) {
            ri.pushMatState();
            ri.translate(pos);
            ri.setStroke(0,255,255,255);
            ri.drawLine(0,0,0,0,0,mass*.1f);        
            ri.popMatState();
        }
    }
    
}//myGridNode
