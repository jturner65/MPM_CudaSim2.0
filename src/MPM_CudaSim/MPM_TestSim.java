package MPM_CudaSim;

//class to manage testing functionality.
public class MPM_TestSim extends MPM_ABS_Sim {

	public MPM_TestSim(int _gridCount, float _h) {
		super(_gridCount, _h);
	}

	@Override
	//all code to run before sim loop, at end of init/reset of sim world
	protected void resetSimPriv() {
		// TODO Auto-generated method stub
		
	}
	@Override
	//check central floating collider - return 0 vec if no collision, otherwise return normal of collision
	public synchronized myVectorf checkColliderCollision(myVectorf pos) {
		return new myVectorf(0,0,0);
		
	}

	@Override
	public boolean simMe(float modAmtMillis) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean simMeDebug(float modAmtMillis) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean simMeMThd(float modAmtMillis) {
		// TODO Auto-generated method stub
		return false;
	}
	
	@Override
	public void drawMe(MPM_SimMain pa, float animTimeMod) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean simMeCuda(float modAmtMillis) {
		// TODO Auto-generated method stub
		return false;
	}
	
}//MPM_TestSim