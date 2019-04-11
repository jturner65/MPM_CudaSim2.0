package MPM_CudaSim;

public class myPoint {
	public double x,y,z;
	public static final myPoint ZEROPT = new myPoint(0,0,0);

	myPoint(double _x, double _y, double _z){this.x = _x; this.y = _y; this.z = _z;}         //constructor 3 args  
	myPoint(myPoint p){ this(p.x, p.y, p.z); }                                                                                                           	//constructor 1 arg  
	myPoint(myPoint A, myVector B) {this(A.x+B.x,A.y+B.y,A.z+B.z); };
	
	myPoint(myPoint A, double s, myPoint B) {this(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };		//builds a point somewhere in between a and b
	myPoint(){ this(0,0,0);}                                                                                                                               //constructor 0 args
	public void clear() {this.x = 0; this.y = 0; this.z = 0;}	
	public void set(double _x, double _y, double _z){ this.x = _x;  this.y = _y;  this.z = _z; }                                               //set 3 args 
	public myPoint set(myPoint p){ this.x = p.x; this.y = p.y; this.z = p.z; return this;}                                                                   //set 1 args
	public void set(double _x, double _y, double _z, double _sqMagn){ this.x = _x;  this.y = _y;  this.z = _z; }                                                                     //set 3 args 
	
	public myPoint _mult(double n){ this.x *= n; this.y *= n; this.z *= n; return this; }                                                     //_mult 3 args  
	public static myPoint _mult(myPoint p, double n){ myPoint result = new myPoint(p.x * n, p.y * n, p.z * n); return result;}                          //1 pt, 1 double
	public static myPoint _mult(myPoint p, myPoint q){ myPoint result = new myPoint(p.x *q.x, p.y * q.y, p.z * q.z); return result;}           //return elementwise product
	public static void _mult(myPoint p, myPoint q, myPoint r){ myPoint result = new myPoint(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}           //2 pt src, 1 pt dest  

	public void _div(double q){this.x /= q; this.y /= q; this.z /= q; }  
	public static myPoint _div(myPoint p, double n){ if(n==0) return p; myPoint result = new myPoint(p.x / n, p.y / n, p.z / n); return result;}                          //1 pt, 1 double
	
	public void _add(double _x, double _y, double _z){ this.x += _x; this.y += _y; this.z += _z;   }                                            //_add 3 args
	public void _add(myPoint v){ this.x += v.x; this.y += v.y; this.z += v.z;   }                                                 //_add 1 arg  
	public static myPoint _add(myPoint O, myVector I){														return new myPoint(O.x+I.x,O.y+I.y,O.z+I.z);}  
	public static myPoint _add(myPoint O, double a, myVector I){												return new myPoint(O.x+a*I.x,O.y+a*I.y,O.z+a*I.z);}                						//2 vec
	public static myPoint _add(myPoint O, double a, myVector I, double b, myVector J) {						return new myPoint(O.x+a*I.x+b*J.x,O.y+a*I.y+b*J.y,O.z+a*I.z+b*J.z);}  					// O+xI+yJ
	public static myPoint _add(myPoint O, double a, myVector I, double b, myVector J, double c, myVector K) {	return new myPoint(O.x+a*I.x+b*J.x+c*K.x,O.y+b*I.y+b*J.y+c*K.y,O.z+b*I.z+b*J.z+c*K.z);} // O+xI+yJ+kZ
	
	public static void _add(myPoint p, myPoint q, myPoint r){ myPoint result = new myPoint(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}       	//2 pt src, 1 pt dest  
	public static myPoint _add(myPoint p, myPoint q){ myPoint result = new myPoint(p.x + q.x, p.y + q.y, p.z + q.z); return result;}
	
	public void _sub(double _x, double _y, double _z){ this.x -= _x; this.y -= _y; this.z -= _z;  }                                                                   //_sub 3 args
	public void _sub(myPoint v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  }                                                                           //_sub 1 arg 
	public static myPoint _sub(myPoint p, myPoint q){ myPoint result = new myPoint(p.x - q.x, p.y - q.y, p.z - q.z); return result; }
	public static void _sub(myPoint p, myPoint q, myPoint r){ myPoint result = new myPoint(p.x - q.x, p.y - q.y, p.z - q.z); r.set(result);}       //2 pt src, 1 pt dest  	

	public myPoint cloneMe(){myPoint retVal = new myPoint(this.x, this.y, this.z); return retVal;}  
	
	public double _L1Dist(myPoint q){return Math.abs((this.x - q.x)) + Math.abs((this.y - q.y)) + Math.abs((this.z - q.z)); }
	public static double _L1Dist(myPoint q, myPoint r){ return q._L1Dist(r);}
	
	public double _SqrDist(myPoint q){ return (((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z))); }
	public double _SqrDist(myPointf q){ return (((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z))); }
	public static double _SqrDist(myPoint q, myPoint r){  return (((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z)));}
	
	public double _dist(myPoint q){ return (double)Math.sqrt( ((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z)) ); }
	public static double _dist(myPoint q, myPoint r){  return (double)Math.sqrt(((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z)));}
	
	public double _dist(double qx, double qy, double qz){ return (double)Math.sqrt( ((this.x - qx)*(this.x - qx)) + ((this.y - qy)*(this.y - qy)) + ((this.z - qz)*(this.z - qz)) ); }
	public static double _dist(myPoint r, double qx, double qy, double qz){  return(double) Math.sqrt(((r.x - qx) *(r.x - qx)) + ((r.y - qy) *(r.y - qy)) + ((r.z - qz) *(r.z - qz)));}	
	
	public double[] asArray(){return new double[]{x,y,z};}
	
	public float[] asFltArray(){return new float[]{(float)x,(float)y,(float)z};}
	
	public boolean clickIn(myPoint p, double eps) { return(_dist(p) < eps);}
	/**
	 * returns if this pttor is equal to passed pttor
	 * @param b myPoint to check
	 * @return whether they are equal
	 */
	public boolean equals(Object b){
		if (this == b) return true;
		if (!(b instanceof myPoint)) return false;
		myPoint v = (myPoint)b;
		return ((this.x == v.x) && (this.y == v.y) && (this.z == v.z));		
	}				
	public String toStrBrf(){return "|(" + String.format("%.4f",this.x) + ", " + String.format("%.4f",this.y) + ", " + String.format("%.4f",this.z)+")";}	
	public String toString(){return "|(" + String.format("%.4f",this.x) + ", " + String.format("%.4f",this.y) + ", " + String.format("%.4f",this.z)+")";}
}


class myPointf {
	public float x,y,z;
	public static final myPointf ZEROPT = new myPointf(0,0,0);

	myPointf(float _x, float _y, float _z){this.x = _x; this.y = _y; this.z = _z;}         //constructor 3 args  
	myPointf(double _x, double _y, double _z){this((float)_x, (float)_y,(float)_z);}         //constructor 3 args  
	myPointf(myPointf p){ this(p.x, p.y, p.z); }                                                                                                           	//constructor 1 arg  
	myPointf(myPointf A, myVectorf B) {this(A.x+B.x,A.y+B.y,A.z+B.z); };
	
	myPointf(myPointf A, float s, myPointf B) {this(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };		//builds a point somewhere in between a and b
	myPointf(){ this(0,0,0);}                                                                                                                               //constructor 0 args
	public void clear() {this.x = 0; this.y = 0; this.z = 0;}
	public void set(float _x, float _y, float _z){ this.x = _x;  this.y = _y;  this.z = _z; }                                               //set 3 args 
	public myPointf set(myPointf p){ this.x = p.x; this.y = p.y; this.z = p.z; return this;}                                                                   //set 1 args
	public void set(float _x, float _y, float _z, float _sqMagn){ this.x = _x;  this.y = _y;  this.z = _z; }                                                                     //set 3 args 
	
	public myPointf _mult(float n){ this.x *= n; this.y *= n; this.z *= n; return this; }                                                     //_mult 3 args  
	public static myPointf _mult(myPointf p, float n){ myPointf result = new myPointf(p.x * n, p.y * n, p.z * n); return result;}                          //1 pt, 1 float
	public static myPointf _mult(myPointf p, myPointf q){ myPointf result = new myPointf(p.x *q.x, p.y * q.y, p.z * q.z); return result;}           //return elementwise product
	public static void _mult(myPointf p, myPointf q, myPointf r){ myPointf result = new myPointf(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}           //2 pt src, 1 pt dest  

	public void _div(float q){this.x /= q; this.y /= q; this.z /= q; }  
	public static myPointf _div(myPointf p, float n){ if(n==0) return p; myPointf result = new myPointf(p.x / n, p.y / n, p.z / n); return result;}                          //1 pt, 1 float
	
	public void _add(float _x, float _y, float _z){ this.x += _x; this.y += _y; this.z += _z;   }                                            //_add 3 args
	public void _add(myPointf v){ this.x += v.x; this.y += v.y; this.z += v.z;   }                                                 //_add 1 arg  
	public static myPointf _add(myPointf O, myVectorf I){														return new myPointf(O.x+I.x,O.y+I.y,O.z+I.z);}  
	public static myPointf _add(myPointf O, float a, myVectorf I){												return new myPointf(O.x+a*I.x,O.y+a*I.y,O.z+a*I.z);}                						//2 vec
	public static myPointf _add(myPointf O, float a, myVectorf I, float b, myVectorf J) {						return new myPointf(O.x+a*I.x+b*J.x,O.y+a*I.y+b*J.y,O.z+a*I.z+b*J.z);}  					// O+xI+yJ
	public static myPointf _add(myPointf O, float a, myVectorf I, float b, myVectorf J, float c, myVectorf K) {	return new myPointf(O.x+a*I.x+b*J.x+c*K.x,O.y+b*I.y+b*J.y+c*K.y,O.z+b*I.z+b*J.z+c*K.z);} // O+xI+yJ+kZ
	
	public static void _add(myPointf p, myPointf q, myPointf r){ myPointf result = new myPointf(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}       	//2 pt src, 1 pt dest  
	public static myPointf _add(myPoint p, myPointf q){ myPointf result = new myPointf(p.x + q.x, p.y + q.y, p.z + q.z); return result;}
	public static myPointf _add(myPointf p, myPointf q){ myPointf result = new myPointf(p.x + q.x, p.y + q.y, p.z + q.z); return result;}
	
	public void _sub(float _x, float _y, float _z){ this.x -= _x; this.y -= _y; this.z -= _z;  }                                                                   //_sub 3 args
	public void _sub(myPointf v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  }                                                                           //_sub 1 arg 
	public static myPointf _sub(myPointf p, myPointf q){ myPointf result = new myPointf(p.x - q.x, p.y - q.y, p.z - q.z); return result; }
	public static void _sub(myPointf p, myPointf q, myPointf r){ myPointf result = new myPointf(p.x - q.x, p.y - q.y, p.z - q.z); r.set(result);}       //2 pt src, 1 pt dest  	

	public myPointf cloneMe(){myPointf retVal = new myPointf(this.x, this.y, this.z); return retVal;}  
	
	public float _L1Dist(myPointf q){return Math.abs((this.x - q.x)) + Math.abs((this.y - q.y)) + Math.abs((this.z - q.z)); }
	public static float _L1Dist(myPointf q, myPointf r){ return q._L1Dist(r);}
	
	public float _SqrDist(myPointf q){ return (((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z))); }
	public static float _SqrDist(myPointf q, myPointf r){  return (((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z)));}
	
	public float _dist(myPointf q){ return (float)Math.sqrt( ((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z)) ); }
	public static float _dist(myPointf q, myPointf r){  return (float)Math.sqrt(((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z)));}
	
	public float _dist(float qx, float qy, float qz){ return (float)Math.sqrt( ((this.x - qx)*(this.x - qx)) + ((this.y - qy)*(this.y - qy)) + ((this.z - qz)*(this.z - qz)) ); }
	public static float _dist(myPointf r, float qx, float qy, float qz){  return(float) Math.sqrt(((r.x - qx) *(r.x - qx)) + ((r.y - qy) *(r.y - qy)) + ((r.z - qz) *(r.z - qz)));}	
	
	public float[] asArray(){return new float[]{x,y,z};}
	public double[] asDblArray(){return new double[]{x,y,z};}
	
	public boolean clickIn(myPointf p, float eps) { return(_dist(p) < eps);}
	/**
	 * returns if this pttor is equal to passed pttor
	 * @param b myPointf to check
	 * @return whether they are equal
	 */
	public boolean equals(Object b){
		if (this == b) return true;
		if (!(b instanceof myPointf)) return false;
		myPointf v = (myPointf)b;
		return ((this.x == v.x) && (this.y == v.y) && (this.z == v.z));		
	}				
	public String toStrBrf(){return "|(" + String.format("%.4f",this.x) + ", " + String.format("%.4f",this.y) + ", " + String.format("%.4f",this.z)+")";}	
	public String toString(){return "|(" + String.format("%.4f",this.x) + ", " + String.format("%.4f",this.y) + ", " + String.format("%.4f",this.z)+")";}
}

class myVector extends myPoint{
	public double sqMagn, magn;
	
		//vector constants available to all consumers of myVector
	public static final myVector ZEROVEC = new myVector(0,0,0);
//	public static final myVector UP	=	new myVector(0,1,0);
//	public static final myVector RIGHT = new myVector(1,0,0);
//	public static final myVector FORWARD = new myVector(0,0,1);
	public static final myVector UP	=	new myVector(0,0,1);
	public static final myVector RIGHT = new myVector(0,1,0);
	public static final myVector FORWARD = new myVector(1,0,0);
	
	myVector(double _x, double _y, double _z){super(_x,_y,_z); this._mag();}         //constructor 3 args  
	myVector(myVector p){ this(p.x, p.y, p.z); }                                                                                                           	//constructor 1 arg  
	myVector(myVectorf p){ this(p.x, p.y, p.z); }                                                                                                           	//constructor 1 arg  
	myVector(){ }//this(0,0,0);}                                                                                                                               //constructor 0 args
	myVector(myPoint a, myPoint b){this(b.x-a.x,b.y-a.y,b.z-a.z);}			//vector from a->b
	myVector(myPoint a){this(a.x,a.y,a.z);}			//vector from 0->a
	
	myVector(myVector a, double _y, myVector b) {super(a,_y,b);this._mag();	}//interp cnstrctr
	public void clear() {super.clear();this.magn = 0; this.sqMagn=0;}
	public void set(double _x, double _y, double _z){ super.set(_x, _y, _z); this._mag(); }                                               //set 3 args 
	public void set(myVector p){ this.x = p.x; this.y = p.y; this.z = p.z;  this._mag();}                                                                   //set 1 args
	public void set(myPoint p, myPoint q){ this.x = q.x - p.x; this.y = q.y - p.y; this.z = q.z - p.z;  this._mag();}                                                                   //set 1 args
	public void set(double _x, double _y, double _z, double _sqMagn){ super.set(_x, _y, _z); this.sqMagn = _sqMagn; }                                                                     //set 3 args 
	
	public myVector _mult(double n){ super._mult(n); this._mag(); return this; }                                                     //_mult 3 args  
	public static myVector _mult(myVector p, double n){ myVector result = new myVector(p.x * n, p.y * n, p.z * n); return result;}                          //1 vec, 1 double
	public static myVector _mult(myVectorf p, double n){ myVector result = new myVector(p.x * n, p.y * n, p.z * n); return result;}                          //1 vec, 1 double
	public static myVector _ewise_mult(myVector p, myVector q){ myVector result = new myVector(p.x *q.x, p.y * q.y, p.z * q.z); return result;}                   //2 vec - point-wise
	public static void _mult(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}           //2 vec src, 1 vec dest  

	public void _div(double q){super._div(q); this._mag();}  
	public static myVector _div(myVector p, double n){ if(n==0) return p; myVector result = new myVector(p.x / n, p.y / n, p.z / n); return result;}                          //1 pt, 1 double
	
	public void _add(double _x, double _y, double _z){ super._add(_x, _y, _z); this._mag(); }                                            //_add 3 args
	public void _add(myVector v){ this.x += v.x; this.y += v.y; this.z += v.z;  this._mag();  }                                                 //_add 1 arg  
	public void _add(myVectorf v){ this.x += v.x; this.y += v.y; this.z += v.z;  this._mag();  }                                                 //_add 1 arg  
	public static myVector _add(myVector p, myVector q){ myVector result = new myVector(p.x + q.x, p.y + q.y, p.z + q.z); return result;}                	//2 vec
	public static myVector _add(myPoint p, myVector q){ myVector result = new myVector(p.x + q.x, p.y + q.y, p.z + q.z); return result;}                	//2 vec
	public static void _add(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}       	//2 vec src, 1 vec dest  
	
	public void _sub(double _x, double _y, double _z){ super._sub(_x, _y, _z);  this._mag(); }                                                                   //_sub 3 args
	public void _sub(myVector v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  this._mag(); }                                                                           //_sub 1 arg 
	public static myVector _sub(myPoint p, myPoint q){ return new myVector(p.x - q.x, p.y - q.y, p.z - q.z);}					             //2 pts or 2 vectors or any combo
	public static void _sub(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x - q.x, p.y - q.y, p.z - q.z); r.set(result);}       //2 vec src, 1 vec dest  
	
	public double _mag(){ this.magn = Math.sqrt(this._SqMag()); return magn; }  
	public double _SqMag(){ this.sqMagn =  ((this.x*this.x) + (this.y*this.y) + (this.z*this.z)); return this.sqMagn; }  							//squared magnitude
	
	public void _scale(double _newMag){this._normalize()._mult(_newMag);}
	
	public myVector _normalize(){this._mag();if(magn==0){return this;}this.x = this.x/magn; this.y = this.y/magn; this.z = this.z/magn; _mag();return this;}
	//returns normal but doesn't modify v
	public static myVector _normalize(myVector v){double magn = v._mag(); if(magn==0){return v;} myVector newVec = new myVector( v.x / magn, v.y / magn, v.z / magn); return newVec;}// newVec._mag(); return newVec;}
	
	public double _norm(){return _mag();}
	public static double _L2Norm(myVector v){return Math.sqrt(v._SqMag());}
	public static double _L2SqNorm(myVector v){return v._SqMag();}
	
	public myVector _normalized(){double magn = this._mag(); myVector newVec = (magn == 0) ? (new myVector(0,0,0)) : (new myVector( this.x /= magn, this.y /= magn, this.z /= magn)); newVec._mag(); return newVec;}

	public myVector cloneMe(){myVector retVal = new myVector(this.x, this.y, this.z); return retVal;}  
		
	public double _SqrDist(myVector q){ return (((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z))); }
	public double _SqrDist(myVectorf q){ return (((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z))); }
	public static double _SqrDist(myVector q, myVector r){  return ((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z));}
	
	public double _dist(myVector q){ return Math.sqrt( ((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z)) ); }
	public static double _dist(myVector q, myVector r){  return Math.sqrt(((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z)));}
	
	public static double _dist(myVector r, double qx, double qy, double qz){  return Math.sqrt(((r.x - qx) *(r.x - qx)) + ((r.y - qy) *(r.y - qy)) + ((r.z - qz) *(r.z - qz)));}	
	
	public myVector _cross(myVector b){ return new myVector((this.y * b.z) - (this.z*b.y), (this.z * b.x) - (this.x*b.z), (this.x * b.y) - (this.y*b.x));}		//cross product 
	public static myVector _cross(myVector a, myVector b){		return a._cross(b);}
	public static myVector _cross(double ax, double ay, double az, double bx, double by, double bz){		return new myVector((ay*bz)-(az*by), (az*bx)-(ax*bz), (ax*by)-(ay*bx));}
	
	
	public static myVector _elemMult(myVector a, myVector b){return new myVector(a.x*b.x, a.y*b.y, a.z*b.z);}
	public static myVector _elemDiv(myVector a, myVector b){return new myVector(a.x/b.x, a.y/b.y, a.z/b.z);}
	
	public static double _det3(myVector U, myVector V) {double udv = U._dot(V); return (Math.sqrt(U._dot(U)*V._dot(V) - (udv*udv))); };                                // U|V det product
	public static double _mixProd(myVector U, myVector V, myVector W) {return U._dot(myVector._cross(V,W)); };                                                 // (UxV)*W  mixed product, determinant - measures 6x the volume of the parallelapiped formed by myVectortors
	public double _dot(myVector b){return ((this.x * b.x) + (this.y * b.y) + (this.z * b.z));}																	//dot product
	public double _dot(myVectorf b){return ((this.x * b.x) + (this.y * b.y) + (this.z * b.z));}																	//dot product
	public static double _dot(myVector a, myVector b){		return a._dot(b);}
	
	public static double _angleBetween(myVector v1, myVector v2) {
		double 	_v1Mag = v1._mag(), 
				_v2Mag = v2._mag(), 
				dotProd = v1._dot(v2),
				cosAngle = dotProd/(_v1Mag * _v2Mag),
				angle = Math.acos(cosAngle);
		return angle;
	}//_angleBetween
	public static myVector _rotAroundAxis(myVector v1, myVector u){return _rotAroundAxis(v1, u, Math.PI*.5);}
	//rotate v1 around axis unit vector u, by give angle thet, around origin
	public static myVector _rotAroundAxis(myVector v1, myVector u, double thet){		
		double cThet = Math.cos(thet), sThet = Math.sin(thet), oneMC = 1-cThet,
				ux2 = u.x*u.x, uy2 = u.y*u.y, uz2 = u.z*u.z,
				uxy = u.x * u.y, uxz = u.x * u.z, uyz = u.y*u.z,
				uzS = u.z*sThet, uyS = u.y*sThet, uxS = u.x*sThet,
				uxzC1 = uxz *oneMC, uxyC1 = uxy*oneMC, uyzC1 = uyz*oneMC;
		//build rot matrix in vector def
		myVector res = new myVector(
				(ux2*oneMC+cThet) * v1.x + (uxyC1-uzS) 		* v1.y + (uxzC1+uyS) *v1.z,
				(uxyC1+uzS) 	  * v1.x + (uy2*oneMC+cThet)* v1.y + (uyzC1-uxS) *v1.z,
				(uxzC1-uyS) 	  * v1.x + (uyzC1+uxS)		* v1.y + (uz2*oneMC + cThet) * v1.z);
		
		return res;		
	}
	
	/**
	 * returns if this vector is equal to passed vector
	 * @param b vector to check
	 * @return whether they are equal
	 */
	public boolean equals(Object b){
		if (this == b) return true;
		if (!(b instanceof myVector)) return false;
		myVector v = (myVector)b;
		return ((this.x == v.x) && (this.y == v.y) && (this.z == v.z));		
	}				
	public String toString(){return super.toString()+ " | Mag:" + String.format("%.4f",this.magn)+ " | sqMag:" + String.format("%.4f",this.sqMagn);}
}//myVector

class myVectorf extends myPointf{
	public float sqMagn, magn;
	
		//vector constants available to all consumers of myVectorf
	public static final myVectorf ZEROVEC = new myVectorf(0,0,0);
//	public static final myVectorf UP	=	new myVectorf(0,1,0);
//	public static final myVectorf RIGHT = new myVectorf(1,0,0);
//	public static final myVectorf FORWARD = new myVectorf(0,0,1);
	public static final myVectorf UP	=	new myVectorf(0,0,1);
	public static final myVectorf RIGHT = new myVectorf(0,1,0);
	public static final myVectorf FORWARD = new myVectorf(1,0,0);
	
	myVectorf(float _x, float _y, float _z){super(_x,_y,_z); this._mag();}         //constructor 3 args  
	myVectorf(double _x, double _y, double _z){this((float)_x,(float)_y,(float)_z); }         //constructor with doubles
	myVectorf(myVectorf p){ this(p.x, p.y, p.z); }                                                                                                           	//constructor 1 arg  
	myVectorf(myVector p){ this((float)p.x, (float)p.y, (float)p.z); }                                                                                                           	//constructor 1 arg  
	myVectorf(){ }//this(0,0,0);}                                                                                                                               //constructor 0 args
	myVectorf(myPointf a, myPointf b){this(b.x-a.x,b.y-a.y,b.z-a.z);}			//vector from a->b
	myVectorf(myPoint a, myPointf b){this(b.x-a.x,b.y-a.y,b.z-a.z);}			//vector from a->b
	myVectorf(myPoint a, myPoint b){this(b.x-a.x,b.y-a.y,b.z-a.z);}			//vector from a->b
	myVectorf(myPointf a){this(a.x,a.y,a.z);}			//vector from 0->a
	
	myVectorf(myVectorf a, float _y, myVectorf b) {super(a,_y,b);this._mag();	}//interp cnstrctr
	public void clear() {super.clear();this.magn = 0; this.sqMagn=0;}
	public void set(float _x, float _y, float _z){ super.set(_x, _y, _z); this._mag(); }                                               //set 3 args 
	public void set(double _x, double _y, double _z){ this.set((float)_x,(float)_y,(float)_z); }                                               //set 3 args 
	public void set(myVectorf p){ this.x = p.x; this.y = p.y; this.z = p.z;  this._mag();}                                                                   //set 1 args
	public void set(myPointf p, myPointf q){ this.x = q.x - p.x; this.y = q.y - p.y; this.z = q.z - p.z;  this._mag();}                                                                   //set 1 args
	public void set(float _x, float _y, float _z, float _sqMagn){ super.set(_x, _y, _z); this.sqMagn = _sqMagn; }                                                                     //set 3 args 
	
	public myVectorf _mult(float n){ super._mult(n); this._mag(); return this; }                                                     //_mult 3 args  
	public myVectorf _mult(double n){ super._mult((float)n); this._mag(); return this; }                                                     //_mult 3 args  
	public static myVectorf _mult(myVectorf p, float n){ myVectorf result = new myVectorf(p.x * n, p.y * n, p.z * n); return result;}                          //1 vec, 1 float
	public static myVectorf _mult(myVectorf p, double n){ myVectorf result = new myVectorf(p.x * n, p.y * n, p.z * n); return result;}                          //1 vec, 1 float
	public static myVectorf _mult(myVectorf p, myVectorf q){ myVectorf result = new myVectorf(p.x *q.x, p.y * q.y, p.z * q.z); return result;}                   //2 vec
	public static void _mult(myVectorf p, myVectorf q, myVectorf r){ myVectorf result = new myVectorf(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}           //2 vec src, 1 vec dest  

	public void _div(float q){super._div(q); this._mag();}  
	public static myVectorf _div(myVectorf p, float n){ if(n==0) return p; myVectorf result = new myVectorf(p.x / n, p.y / n, p.z / n); return result;}                          //1 pt, 1 float
	
	public void _add(float _x, float _y, float _z){ super._add(_x, _y, _z); this._mag(); }                                            //_add 3 args
	public void _add(myVectorf v){ this.x += v.x; this.y += v.y; this.z += v.z;  this._mag();  }                                                 //_add 1 arg  
	public static myVectorf _add(myVectorf p, myVectorf q){ myVectorf result = new myVectorf(p.x + q.x, p.y + q.y, p.z + q.z); return result;}                	//2 vec
	public static myVectorf _add(myPointf p, myVectorf q){ myVectorf result = new myVectorf(p.x + q.x, p.y + q.y, p.z + q.z); return result;}                	//2 vec
	public static void _add(myVectorf p, myVectorf q, myVectorf r){ myVectorf result = new myVectorf(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}       	//2 vec src, 1 vec dest  
	
	public void _sub(float _x, float _y, float _z){ super._sub(_x, _y, _z);  this._mag(); }                                                                   //_sub 3 args
	public void _sub(myVectorf v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  this._mag(); }                                                                           //_sub 1 arg 
	public static myVectorf _sub(myPointf p, myPointf q){ return new myVectorf(p.x - q.x, p.y - q.y, p.z - q.z);}					             //2 pts or 2 vectors or any combo
	public static void _sub(myVectorf p, myVectorf q, myVectorf r){ myVectorf result = new myVectorf(p.x - q.x, p.y - q.y, p.z - q.z); r.set(result);}       //2 vec src, 1 vec dest  
	
	public float _mag(){ this.magn = (float)Math.sqrt(this._SqMag()); return magn; }  
	public float _SqMag(){ this.sqMagn =  ((this.x*this.x) + (this.y*this.y) + (this.z*this.z)); return this.sqMagn; }  							//squared magnitude
	
	public void _scale(float _newMag){this._normalize()._mult(_newMag);}
	
	public myVectorf _normalize(){this._mag();if(magn==0){return this;}this.x = this.x/magn; this.y = this.y/magn; this.z = this.z/magn; _mag();return this;}
	//public static myVectorf _normalize(myVectorf v){float magn = v._mag(); if(magn==0){return v;} myVectorf newVec = new myVectorf( v.x /= magn, v.y /= magn, v.z /= magn); newVec._mag(); return newVec;}
	public static myVectorf _normalize(myVectorf v){double magn = v._mag(); if(magn==0){return v;} myVectorf newVec = new myVectorf( v.x / magn, v.y / magn, v.z / magn);return newVec;}// newVec._mag(); return newVec;}

	public float _norm(){return _mag();}
	public static float _L2Norm(myVectorf v){return (float)Math.sqrt(v._SqMag());}
	public static float _L2SqNorm(myVectorf v){return v._SqMag();}
	
	public myVectorf _normalized(){float magn = this._mag(); myVectorf newVec = (magn == 0) ? (new myVectorf(0,0,0)) : (new myVectorf( this.x /= magn, this.y /= magn, this.z /= magn)); newVec._mag(); return newVec;}

	public myVectorf cloneMe(){myVectorf retVal = new myVectorf(this.x, this.y, this.z); return retVal;}  
		
	public float _SqrDist(myVectorf q){ return (((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z))); }
	public static float _SqrDist(myVectorf q, myVectorf r){  return ((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z));}
	
	public float _dist(myVectorf q){ return (float)Math.sqrt( ((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z)) ); }
	public static float _dist(myVectorf q, myVectorf r){  return (float)Math.sqrt(((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z)));}
	
	public static float _dist(myVectorf r, float qx, float qy, float qz){  return (float)Math.sqrt(((r.x - qx) *(r.x - qx)) + ((r.y - qy) *(r.y - qy)) + ((r.z - qz) *(r.z - qz)));}	
	
	public myVectorf _cross(myVectorf b){ return new myVectorf((this.y * b.z) - (this.z*b.y), (this.z * b.x) - (this.x*b.z), (this.x * b.y) - (this.y*b.x));}		//cross product 
	public static myVectorf _cross(myVectorf a, myVectorf b){		return a._cross(b);}
	public static myVectorf _cross(float ax, float ay, float az, float bx, float by, float bz){		return new myVectorf((ay*bz)-(az*by), (az*bx)-(ax*bz), (ax*by)-(ay*bx));}
	
	
	public static myVectorf _elemMult(myVectorf a, myVectorf b){return new myVectorf(a.x*b.x, a.y*b.y, a.z*b.z);}
	public static myVectorf _elemDiv(myVectorf a, myVectorf b){return new myVectorf(a.x/b.x, a.y/b.y, a.z/b.z);}
	
	public static float _det3(myVectorf U, myVectorf V) {float udv = U._dot(V); return (float)(Math.sqrt(U._dot(U)*V._dot(V) - (udv*udv))); };                                // U|V det product
	public static float _mixProd(myVectorf U, myVectorf V, myVectorf W) {return U._dot(myVectorf._cross(V,W)); };
	public float _dot(myVectorf b){return ((this.x * b.x) + (this.y * b.y) + (this.z * b.z));}
	public static float _dot(myVectorf a, myVectorf b){		return a._dot(b);}
	
	public static float _angleBetween(myVectorf v1, myVectorf v2) {
		float 	_v1Mag = v1._mag(), 
				_v2Mag = v2._mag(), 
				dotProd = v1._dot(v2),
				cosAngle = dotProd/(_v1Mag * _v2Mag),
				angle = (float)(Math.acos(cosAngle));
		return angle;
	}//_angleBetween
	public static myVectorf _rotAroundAxis(myVectorf v1, myVectorf u){return _rotAroundAxis(v1, u, (float)Math.PI*.5f);}
	//rotate v1 around axis unit vector u, by give angle thet, around origin
	public static myVectorf _rotAroundAxis(myVectorf v1, myVectorf u, float thet){		
		float cThet = (float)(Math.cos(thet)), sThet = (float)(Math.sin(thet)), oneMC = 1-cThet,
				ux2 = u.x*u.x, uy2 = u.y*u.y, uz2 = u.z*u.z,
				uxy = u.x * u.y, uxz = u.x * u.z, uyz = u.y*u.z,
				uzS = u.z*sThet, uyS = u.y*sThet, uxS = u.x*sThet,
				uxzC1 = uxz *oneMC, uxyC1 = uxy*oneMC, uyzC1 = uyz*oneMC;
		//build rot matrix in vector def
		myVectorf res = new myVectorf(
				(ux2*oneMC+cThet) * v1.x + (uxyC1-uzS) 		* v1.y + (uxzC1+uyS) *v1.z,
				(uxyC1+uzS) 	  * v1.x + (uy2*oneMC+cThet)* v1.y + (uyzC1-uxS) *v1.z,
				(uxzC1-uyS) 	  * v1.x + (uyzC1+uxS)		* v1.y + (uz2*oneMC + cThet) * v1.z);
		
		return res;		
	}
	
	/**
	 * returns if this vector is equal to passed vector
	 * @param b vector to check
	 * @return whether they are equal
	 */
	public boolean equals(Object b){
		if (this == b) return true;
		if (!(b instanceof myVectorf)) return false;
		myVectorf v = (myVectorf)b;
		return ((this.x == v.x) && (this.y == v.y) && (this.z == v.z));		
	}				
	public String toString(){return super.toString()+ " | Mag:" + String.format("%.4f",this.magn)+ " | sqMag:" + String.format("%.4f",this.sqMagn);}
}//myVectorf



