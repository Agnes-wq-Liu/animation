package comp559.particle;

import javax.vecmath.Vector2d;

public class SymplecticEuler implements Integrator {

    @Override
    public String getName() {
        return "symplectic Euler";
    }

    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
        // TODO: Objective 7, complete the symplectic Euler integration method.
    	// note you'll need to know how p is packed to properly implement this, so go
    	// look at ParticleSystem.getPhaseSpace()
    	Vector2d v1 = new Vector2d();
    	Vector2d x1 = new Vector2d();
    	double[] dpdt = new double [n];
		derivs.derivs(t,p,dpdt);
    	for (int i=0;i<n-3;i+=4) {
    		//v1 = v0+h(v0-x0)
    		v1.set(new Vector2d(p[i+2],p[i+3]));
    		v1.scaleAdd(h, new Vector2d(p[i+2]-p[i],p[i+3]-p[i+1]), v1);
    		//x1 = x0+h*v1
    		x1.set(new Vector2d(p[i],p[i+1]));
    		x1.scaleAdd(h, v1, x1);
    		//update velocity and positions
    		pout[i] = x1.x;
    		pout[i+1] = x1.y;
    		pout[i+2] = v1.x;
    		pout[i+3] = v1.y;
    	}
    }

}
