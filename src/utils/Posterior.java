import cern.jet.stat.*;

public class Posterior {
	public Posterior(){
	}
	public Posterior(double ft1, double ft2, double nt){
	 	m1 = ft1;
		m2 = ft2;
		n  = nt;
	}
	double m1;
	double m2;
	double n;
	
        public static void main(String[] args) {
		double b = Gamma.incompleteBeta(500.,1000.,1./3.);
		System.out.println("\n"+b);
        }
}
