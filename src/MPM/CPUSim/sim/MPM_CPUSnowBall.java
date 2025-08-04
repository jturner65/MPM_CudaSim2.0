package MPM.CPUSim.sim;

import MPM.BaseSim.sim.Base_MPMSimFlags;
import MPM.BaseSim.ui.Base_MPMSimWindow;
import MPM.BaseSim.utils.MPM_SimUpdateFromUIData;
import MPM.CPUSim.sim.base.Base_MPMCPUSim;
import MPM.CPUSim.sim.base.MPM_CPUSimFlags;
import MPM.CPUSim.sim.particles.MPM_CPURndrdPart;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;
import base_Render_Interface.IGraphicsAppInterface;

/**
 * @author John Turner
 *
 */
public class MPM_CPUSnowBall extends Base_MPMCPUSim {
    
    //collider to demonstrate behavior
    /**
     * location of spherical collider
     */
    protected myVectorf colLocation;
    /**
     * radius of sphere collider
     */
    protected float colRad, colSqRad;
    
    /**
     * @param _pa
     * @param _win
     * @param _simName
     * @param _gravity
     * @param _currUIVals
     */
    public MPM_CPUSnowBall(IGraphicsAppInterface _pa, Base_MPMSimWindow _win, MPM_SimUpdateFromUIData _currUIVals) {
        super(_pa, _win, "Snowball Fall", _currUIVals);
        colLocation = new myVectorf();
    }
    
    @Override
    protected void updateCPUSimVals_FromUI_Indiv(MPM_SimUpdateFromUIData upd) {
    }//updateCPUSimVals_FromUI_Indiv


    @Override
    protected void buildPartLayoutArray(MPM_CPURndrdPart[] parts) {
        // TODO Auto-generated method stub
        
    }

    @Override
    protected void reinitSimObjects(MPM_CPURndrdPart[] parts) {
        // TODO Auto-generated method stub
        
    }
    
    
    @Override
    protected void drawColliders_Indiv(float animTimeMod) {
        ri.pushMatState();    
        ri.setSphereDetail(20);
        ri.translate(colLocation);
        ri.drawSphere(this.colRad);
        ri.popMatState();    
    }

    /**
     * sim method to show execution time and debug information for each sim step
     * @param modAmtMillis
     * @return
     */
    @Override
    public final boolean simMeDebug(float modAmtMillis) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    protected Base_MPMSimFlags buildSimFlags() {return new MPM_CPUSimFlags(this);    }

    @Override
    public void handleSimFlagsDebug_Indiv(boolean val) {
        // TODO Auto-generated method stub
        
    }

    //check central floating collider - return 0 vec if no collision, otherwise return normal of collision
    @Override
    public myVectorf checkColliderCollision(myPointf pos) {
        if(myPointf._SqrDist(pos, colLocation) <= colSqRad) {
            myVectorf tmp = new myVectorf(pos);
            tmp._sub(colLocation);
            return tmp._normalize();
        }
        return new myVectorf(0,0,0);
    }


}//class SnowBallSim
