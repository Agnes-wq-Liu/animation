package comp559.particle;

/**
 * Interface for a class that computes an unknown function's derivative
 * and checks that a provided state is valid.
 * @author kry
 */
public interface Function {
    
    /**
     * Evaluates derivatives for ODE integration.  
     * The forces could be time varying, which is why t is provided, but
     * in the main objectives of the assignment you will note that there 
     * is no time dependence for the forces.
     * 
     * @param t time 
     * @param p phase space state
     * @param dpdt to be filled with the derivative
     */
    public void derivs( double t, double p[], double dpdt[] );

}




