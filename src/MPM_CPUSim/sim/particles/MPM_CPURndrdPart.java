package MPM_CPUSim.sim.particles;

import MPM_CPUSim.sim.particles.base.Base_MPMCPUParticle;
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
	
	public void drawMe(IRenderInterface pa, float partRad) {
		pa.pushMatState();
		pa.translate(pos);
		pa.drawSphere(partRad);
		pa.popMatState();

	}//drawMe
	
	public void drawMeVel(IRenderInterface pa, float partRad) {
		pa.pushMatState();
		pa.translate(pos);
		pa.drawSphere(partRad);
		
		pa.setStroke(190,20,70,155);
		pa.drawLine(0,0,0,vel.x,vel.y,vel.z);
		pa.popMatState();
	}//drawMeVel
}//class myRndrdPart
