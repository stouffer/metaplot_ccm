#include <R.h>
#include <Rmath.h>
void getorder(int neighbors[], double distances[], int E, int from, int to, int i, int LibUse[]);
void getrho(double A[], double Aest[], double rho[], int from, int to, int l, int LibLength, double *varrho);
void CCM_bootstrap(double *A, double *Aest, double *B, double *rho, int *pE, int *ptau, int *plengtht, int *pLibLength, int *DesiredL, int *plengthDesL, int *piterations, double *varrho, int *acceptablelib, int *plengthacceptablelib) {
    int i, j, k, l, from, to, slide, count, slidetrip=0, lindex;
    double distsv, sumu, sumaest, sumw;
    int E = *pE;
	int nneigh=E+1;
    int tau = *ptau;
    int lengtht= *plengtht;
    int lengthDesL= *plengthDesL;
    int LibLength= *pLibLength;
    int neighbors[E+1];
    double u[E+1], w[E+1], distances[LibLength];
    int LibUse[LibLength];
    int iterations = *piterations;
    int lengthacceptablelib = *plengthacceptablelib;
    int integerpos;
    GetRNGstate(); // Load seed for random number generation, R base

    
    /* Code to implement Sugihara&al 2012 CCM algorithm to determing causality */
    /* Checks to see that *A causes *B */
    /* Note that *A and *B must be the same length, and standardized to same timestep */
    /* *distances is length(*A), and *Aest are length (length(*A or *B)- tau -(E+1)) */
    /* *E, *tau, and *lengtht are all single integers */
    /* *neighbors, *u, and *w are all (E+1) long */
    
    from=(tau*(E-1)); // starting point for most vectors
    for(lindex=0; lindex<lengthDesL; lindex++) { /* try out all desired library sizes, within ((E+1) to LibLengh-tau*(E-1)) */
        l=DesiredL[lindex];
        if(l<(tau*(E-1)+(E+1))) {
            l=(tau*(E-1)+(E+1));
        } /* Catch values that fall outside of feasible library */
        if(l>=lengtht) {
            l=(lengtht-1);
        } /* Catch values that fall outside of feasible library */
        to=l; // Set "end" of library for each iteration
        if(slidetrip==0) { // Trigger to end function when we reach the end of the library
            for(slide=from; slide<iterations; slide++) { // Run a step for each desired iteration
		if(slidetrip==0) {
			integerpos=floor(runif(0,1)*(lengthacceptablelib));
      			count=acceptablelib[integerpos]; //Random number generator - populates Count with an entry that has enough of a history for give E.
			for(j=from; j<=to; j++){ /* create first round of L assignments */
				// Use "skipvector" to identify regions with holes, and leave them out of the library.


				if(count<LibLength) { // If we have not yet "wrapped around" to the beginning of the data
					LibUse[j]=count;
				}else{ // Otherwise, account for "jump" over the end of the library
					LibUse[j]=count-(LibLength-(tau*(E-1)));
				} // LibUse is now a vector of positions, based on acceptable starting points (don't "jump" over gaps), including any wrapping around the Library that was required.

				integerpos=floor(runif(0,1)*(lengthacceptablelib));
      				count=acceptablelib[integerpos]; //Random number generator
			}
                    
                    for(int ii=0; ii<lengthacceptablelib; ii++) { //Predict all points in A using information from the minimized library
                    	i=acceptablelib[ii]; // Make sure we only predict points that have suitable time lag information
                        for(j=from; j<=to; j++) { // scroll across elements in minimized L (based on lengthDesL, including wrapping)
                            distances[LibUse[j]]=0;
                            for(k=0; k<E; k++) {
                                distances[LibUse[j]]=distances[LibUse[j]]+pow((B[i-tau*k]-B[LibUse[j]-tau*k]),2); //calculate distances between points on shadow manifold for all E lagged dimensions
                            }
                            distances[LibUse[j]]=sqrt(distances[LibUse[j]]);
                        }
                        
                        getorder(neighbors, distances, E, from, to, i, LibUse); /*find "tme" position of (E+1) closest points to B[i] on the shadow manifold*/
                        
                        distsv=distances[neighbors[0]];
                        
                        sumaest=0.; /* find w, and weighted Aest variables */
                        if(distsv!=0) { /* check whether minimum distance is zero */
                            sumu=0.; /* find u for all neighbors, and sumu */
                            for(j=0; j<(nneigh); j++) {
                                u[j]=exp(-distances[neighbors[j]]/distsv);
                                sumu=sumu+u[j];
                            }
                            sumw=0.;
                            for(j=0; j<(nneigh); j++) {
                                w[j]=u[j]/sumu;
                                if(w[j]<0.000001) { //minimum cap on weights, taken from Hao Ye code (Sugihara CCM paper)
                                    w[j]=0.000001;
                                }
                                sumw=sumw+w[j];
                            }
                            for(j=0; j<(nneigh); j++) {
                                w[j]=w[j]/sumw;
                                sumaest=sumaest+(A[neighbors[j]])*(w[j]);
                            }
                        } else {
                            sumw=0.;
                            for(j=0; j<(nneigh); j++) {
                                if(distances[neighbors[j]]==0) {
                                    w[j]=1;
                                } else {
                                    w[j]=0.000001;
                                }
                                sumw=sumw+w[j];
                            }
                            for(j=0; j<(nneigh); j++) {
                                w[j]=w[j]/sumw;
                                sumaest=sumaest+(A[neighbors[j]])*(w[j]);
                            }
                        }
                        Aest[i]=sumaest; /* calculate Aest */
                    }
                    getrho(A, Aest, rho, from, to, l, LibLength, varrho);
                    if(to==LibLength-1) {
                        slidetrip=1;
                    }
		}
            }
        }
        if(slidetrip==0) {
            rho[l]=rho[l]/iterations; // Calcualte mean of rho values
            varrho[l]=varrho[l]/iterations; // Calculate mean of rho squared
            varrho[l]=varrho[l]-pow(rho[l], 2); // Calculate variance of rho
        }
    }
    PutRNGstate(); // Free up state of R random number generator
}

void getorder(int neighbors[], double distances[], int E, int from, int to, int i, int LibUse[]) {
    /*find "tme" position of (Ego+1) closest points to B[i] on the shadow manifold*/
    /* include output only for distancesgo >0 */
    int trip, n, ii, j, k, nneigh;
    nneigh=1;
    n=0;
    if(LibUse[from]==i) { /* if first element is a self-reference */
        n=n+1; /* #### scroll across elements in L (including wrapping) #### */
    }
	neighbors[0]=LibUse[from+n];
    for(ii=(from+n); ii<=to; ii++) { /* #### scroll across elements in L (including wrapping) #### */
        trip=0;
        for(j=0; j<nneigh; j++) {
            if((distances[LibUse[ii]]<distances[neighbors[j]])&(LibUse[ii]!=i)) {
                for(k=(nneigh); k>(j); k--) {
                    if(k<(E+1)) {
                        neighbors[k]=neighbors[k-1];
                    }
                }
                neighbors[j]=LibUse[ii];
                if(nneigh<(E+1)) {
                    nneigh=nneigh+1;
                }
                trip=1;
                break;
            }
        }
        if((trip==0)&(nneigh<(E+1))&(LibUse[ii]!=i)&(neighbors[nneigh-1]!=LibUse[ii])) { /* Add element if vector is not yet full, but distance is greater than all already recorded */
            neighbors[nneigh]=LibUse[ii];
            if(nneigh<(E+1)) {
				nneigh=nneigh+1;
            }
        }
    }
}

void getrho(double A[], double Aest[], double rho[], int from, int to, int l, int LibLength, double *varrho) {
    /*Calculate Pearson correlation coefficient between A and Aest*/
    int j, n=0;
    double xbar=0, ybar=0, rhocalc=0, xyxybar=0, xxbarsq=0, yybarsq=0;
    for(j=from; j<LibLength; j++) { /* #### scroll across elements in L (including wrapping) #### */
        xbar=xbar+A[j];
        ybar=ybar+Aest[j];
        n=n+1;
    }
    xbar=xbar/n;
    ybar=ybar/n;
    for(j=from; j<LibLength; j++) { /* #### scroll across elements in L (including wrapping) #### */
        xyxybar=xyxybar+((A[j]-xbar)*(Aest[j]-ybar));
        xxbarsq=xxbarsq+pow(A[j]-xbar,2);
        yybarsq=yybarsq+pow(Aest[j]-ybar,2);
    }
    rhocalc=xyxybar/(sqrt(xxbarsq)*sqrt(yybarsq));
    rho[l]=rho[l]+rhocalc;
    varrho[l]=varrho[l]+pow(rhocalc, 2); //Save rho squared, for calculating variance
}

