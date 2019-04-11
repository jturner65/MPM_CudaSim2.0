package MPM_CudaSim;

//class for the material that the simulation consists of
public class myMaterial {
	private float initYoungMod, poissonRatio, lambda0, mu0;
	public float hardeningCoeff, criticalCompression, criticalStretch, alphaPicFlip;
	

	//initial values of material quantities
	public static float	
		base_initYoungMod 			 = 1.4e5f,
		base_poissonRatio 			 = 0.2f,
		base_hardeningCoeff 		 = 15.0f, 
		base_criticalCompression 	 = 0.010f, 
		base_criticalStretch 		 = 0.0025f, 
		base_alphaPicFlip 			 = 0.95f;
	
	//TODO only using this so that initial material quantities will always match what is displayed in UI
	private myMaterial(float _ym, float _pr, float _hc, float _cc, float _cs, float _aPF) {
		initYoungMod = _ym; 
		poissonRatio = _pr;   
		hardeningCoeff = _hc;     
		criticalCompression = _cc;
		criticalStretch = _cs;   
		alphaPicFlip = _aPF;
		recalcParams();
	}
		
	public myMaterial() {
		this(base_initYoungMod, base_poissonRatio, base_hardeningCoeff, base_criticalCompression,base_criticalStretch, base_alphaPicFlip );
	}
	
	public void setYoungModulus(float _ym) {
		initYoungMod = _ym;
		recalcParams();
	}
	
	//return all current material properties, to populate UI on reset
	public float[] getMatProps() {
		return new float[] {initYoungMod, poissonRatio,hardeningCoeff, criticalCompression, criticalStretch, alphaPicFlip};
	}
	
	public void setPoissonRatio(float _pr) {
		poissonRatio = _pr;
		recalcParams();
	}
	public void setHC(float _hc) {hardeningCoeff = _hc;}
	public void setCC(float _cc) {criticalCompression = _cc;}
	public void setCS(float _cs) {criticalStretch = _cs;}
	public void setAlpha(float _ap) {alphaPicFlip = _ap;}
	
	
	//recalc lambda0 and mu0 when young mod or poisson ratio changes
	private void recalcParams() {	
		lambda0 = (initYoungMod * poissonRatio) / ((1.0f+poissonRatio) * (1.0f-2.0f*poissonRatio));
		mu0 = initYoungMod / (2.0f * (1.0f+poissonRatio));
	}
	//recalc mu0 when ym or poisson ratio changes
	
	//private to govern recalculation when youngmod and poisson ratio change via UI input
	public float getLambda0() {return lambda0;}
	public float getMu0() {return mu0;}	
	public float getYoungMod() {return initYoungMod;}
	public float getPoissonRatio() {return poissonRatio;}
	
	
	
}//myMaterial
