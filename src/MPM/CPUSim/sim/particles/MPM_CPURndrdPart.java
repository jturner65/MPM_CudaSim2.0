package MPM.CPUSim.sim.particles;

import MPM.CPUSim.sim.particles.base.Base_MPMCPUParticle;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Render_Interface.IRenderInterface;

public class MPM_CPURndrdPart extends Base_MPMCPUParticle{
    public int[] color;
    //so we can reset color if we change it
    public static final int[] origColor = new int[] {255,255,255,255};
    
    public MPM_CPURndrdPart(myPointf _iPos) {
        super(_iPos);
        color = new int[] {255,255,255,255};
    }
    
    public void drawMe(IRenderInterface ri, float partRad) {
        ri.pushMatState();
        ri.translate(pos);
        ri.drawSphere(partRad);
        ri.popMatState();

    }//drawMe
    
    public void drawMeVel(IRenderInterface ri, float partRad) {
        ri.pushMatState();
        ri.translate(pos);
        ri.drawSphere(partRad);
        
        ri.setStroke(190,20,70,155);
        ri.drawLine(0,0,0,vel.x,vel.y,vel.z);
        ri.popMatState();
    }//drawMeVel
}//class myRndrdPart
