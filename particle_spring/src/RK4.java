package comp559.particle;

import javax.vecmath.Vector2d;

public class RK4 implements Integrator {
    
    @Override
    public String getName() {
        return "RK4";
    }
    private double[] tmp;
    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
        // TODO: Objective 6, implement the RK4 integration method
    	// see also efficient memory management suggestion in provided code for the Midpoint method.
    	if ( tmp == null || tmp.length != n ) {
            tmp = new double[n];
    	}
    	//state: p.x,p.y,v.x,v.y
    	Vector2d v1 = new Vector2d();
    	Vector2d x1 = new Vector2d();
    	Vector2d a = new Vector2d();
		derivs.derivs(t,p,tmp);
		//gives me k1
    	for (int i=0;i<p.length-3;i+=4) {
    		double[] k2 = new double[4];
    		double[] k3 = new double[4];
    		double[] k4 = new double[4];
    		//x1 = x0+h(x0+1/2h(x0)')'
    		double[] tuples = new double[] {tmp[i]*1/2*h+p[i],tmp[i+1]*1/2*h+p[i+1],
    				tmp[i+2]*1/2*h+p[i+2],tmp[i+3]*1/2*h+p[i+3]};
    		derivs.derivs(t,tuples,k2);
    		tuples = new double[] {k2[i]*1/2*h+p[i],k2[i+1]*1/2*h+p[i+1],
    				k2[i+2]*1/2*h+p[i+2],k2[i+3]*1/2*h+p[i+3]};
    		derivs.derivs(t,tuples,k3);
    		tuples = new double[] {k3[i]*h+p[i],k3[i+1]*h+p[i+1],
    				k3[i+2]*h+p[i+2],k3[i+3]*h+p[i+3]};
    		derivs.derivs(t,tuples,k4);
    		x1.set(new Vector2d(p[i],p[i+1]));
    		x1.scaleAdd(h/6, new Vector2d(tmp[i],tmp[i+1]), x1);
    		x1.scaleAdd(h/3, new Vector2d(k2[i],k2[i+1]), x1);
    		x1.scaleAdd(h/3, new Vector2d(k3[i],k3[i+1]), x1);
    		x1.scaleAdd(h/6, new Vector2d(k4[i],k4[i+1]), x1);
    		//v1 = v0+h*a
    		a.set(tmp[i+2],tmp[i+3]);
    		v1.set(new Vector2d(p[i+2],p[i+3]));
    		v1.scaleAdd(h/6, a, v1);
    		v1.scaleAdd(h/3, new Vector2d(k2[i+2],k2[i+3]), v1);
    		v1.scaleAdd(h/3, new Vector2d(k3[i+2],k3[i+3]), v1);
    		v1.scaleAdd(h/6, new Vector2d(k4[i+2],k4[i+3]), v1);
    		//update velocity and positions
    		pout[i] = x1.x;
    		pout[i+1] = x1.y;
    		pout[i+2] = v1.x;
    		pout[i+3] = v1.y;
    	}
    }
}
