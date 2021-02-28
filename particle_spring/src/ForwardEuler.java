package comp559.particle;
import javax.vecmath.Vector2d;
/**
 * TODO: finish this class!
 */
public class ForwardEuler implements Integrator {
    
    @Override
    public String getName() {
        return "Forward Euler";
    }
    
    /** 
     * Advances the system at t by h 
     * @param p The state at time h
     * @param n The dimension of the state (i.e., p.length)
     * @param t The current time (in case the derivs function is time dependent)
     * @param h The step size
     * @param pout  The state of the system at time t+h
     * @param derivs The object which computes the derivative of the system state
     */
    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
        // TODO: Objective 3, implement the forward Euler method
        //state: p.x,p.y,v.x,v.y
    	Vector2d v1 = new Vector2d();
    	Vector2d x1 = new Vector2d();
    	Vector2d a = new Vector2d();
    	double[] dpdt = new double [n];
		derivs.derivs(t,p,dpdt);
    	for (int i=0;i<n-3;i+=4) {
    		//x1 = x0+h*v0
    		x1.set(new Vector2d(p[i],p[i+1]));
    		x1.scaleAdd(h, new Vector2d(p[i+2],p[i+3]), x1);
    		//v1 = v0+h*a
    		a.set(dpdt[i+2],dpdt[i+3]);
    		v1.set(new Vector2d(p[i+2],p[i+3]));
    		v1.scaleAdd(h, a, v1);
    		//update velocity and positions
    		pout[i] = x1.x;
    		pout[i+1] = x1.y;
    		pout[i+2] = v1.x;
    		pout[i+3] = v1.y;
    	}
    	    	
    }

}
