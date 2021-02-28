package comp559.particle;

import javax.vecmath.Vector2d;

public class Midpoint implements Integrator {

    @Override
    public String getName() {
        return "midpoint";
    }

    private double[] tmp;
    
    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
        // TODO: Objective 4, implement midpoint method

    	// You will probably want a temporary array in this method and perhaps
    	// multiple temporary arrays in other higher order explicit integrators.
    	// Avoid thrashing memory by reallocating only when necessary.
    	if ( tmp == null || tmp.length != n ) {
            tmp = new double[n];
    	}
    	//state: p.x,p.y,v.x,v.y
    	Vector2d v1 = new Vector2d();
    	Vector2d x1 = new Vector2d();
    	Vector2d a = new Vector2d();
		derivs.derivs(t,p,tmp);
		//gives me f(x0) 
    	for (int i=0;i<p.length-3;i+=4) {
    		//x1 = x0+h(x0+1/2h(x0)')'
    		double[] tuples = new double[] {tmp[i]*1/2*h+p[i],tmp[i+1]*1/2*h+p[i+1],
    				tmp[i+2]*1/2*h+p[i+2],tmp[i+3]*1/2*h+p[i+3]};
    		double[] newSet = new double[4];
    		derivs.derivs(t,tuples,newSet);
    		x1.set(new Vector2d(p[i],p[i+1]));
    		x1.scaleAdd(h, new Vector2d(newSet[i],newSet[i+1]), x1);
    		//v1 = v0+h*a
    		a.set(newSet[i+2],newSet[i+3]);
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
