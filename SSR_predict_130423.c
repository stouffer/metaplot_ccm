#include <R.h>
#include <Rmath.h>
void getorder(int neighbors[], double distances[], int E, int from, int i, int Blength, int predstep);
void SSR_predict_130423(double *A, double *Aest, double *B, int *pE, int *ptau, int *pBlength, int *pAlength, int *ppredstep, int *prepvec, int *pmatchSugi) {
    int i, j, k, from;
    double distsv, sumu, sumaest, sumw;
    int E = *pE;
	int nneigh=E+1;
    int tau = *ptau;
    int Blength= *pBlength;
    int Alength= *pAlength;
    int repvec= *prepvec;
    int neighbors[E+1];
    int predstep= *ppredstep;
    int matchSugi=*pmatchSugi;
    double u[E+1], w[E+1], distances[Blength];
    /* Code to implement Sugihara&al SSR algorithm */
    /* Uses information in B to predict next predstep steps of A */
    /* Note that *A and *B need not be the same size - but if they are the same vector, then repvec must = "1" */
    
    from=(tau*(E-1)); /* starting point in most vectors */
    for(i=from; i<Alength; i++) {	/*  include all indices  */
        if(repvec==1) { /* A and B are same vector - use "leave one out cross validation" */
            for(j=from; j<Blength-predstep; j++) { /* scroll across elements in L */
                distances[j]=0;
                
                if(matchSugi==1) {
                if(i!=j) {
                    for(k=0; k<E; k++) {
                        distances[j]=distances[j]+pow((A[i-tau*k]-B[j-tau*k]),2); /*calculate distances between points on A and B for all E lagged dimensions*/
                    }
                    distances[j]=sqrt(distances[j]);
                } else {
                    distances[j]=999;
                }
                } else{
                    if((j>i+predstep)|(j<=(i-E))) {
                        for(k=0; k<E; k++) {
                            distances[j]=distances[j]+pow((A[i-tau*k]-B[j-tau*k]),2); /*calculate distances between points on A and B for all E lagged dimensions*/
                        }
                        distances[j]=sqrt(distances[j]);
                    } else {
                        distances[j]=999;
                    }
                }
            }
        } else {
            for(j=from; j<Blength; j++) { /* scroll across elements in L */
                distances[j]=0;
                for(k=0; k<E; k++) {
                    distances[j]=distances[j]+pow((A[i-tau*k]-B[j-tau*k]),2); /*calculate distances between points on A and B for all E lagged dimensions*/
                }
                distances[j]=sqrt(distances[j]);
            }
        }
        
        getorder(neighbors, distances, E, from, i, Blength, predstep); /*find "tme" position of (E+1) closest points to A[i] on B manifold*/
        
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
                if(w[j]<0.000001) {
                    w[j]=0.000001;
                }
                sumw=sumw+w[j];
            }
            for(j=0; j<(nneigh); j++) {
                w[j]=w[j]/sumw;
                sumaest=sumaest+(B[neighbors[j]+predstep])*(w[j]);
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
}


void getorder(int neighbors[], double distances[], int E, int from, int i, int Blength, int predstep) {
    int trip, n, ii, j, k, nneigh;
    nneigh=1;
    n=0;
    if(from==i) { /* if first element is a self-reference */
        n=n+1;
    }
	neighbors[0]=from+n;
    
    for(ii=(from+n); ii<Blength-predstep; ii++) { /* scroll across elements in L */
        trip=0;
        for(j=0; j<nneigh; j++) {
            if((distances[ii]<distances[neighbors[j]])&(ii!=i)&(distances[ii]>0)) {
                for(k=(nneigh); k>(j); k--) {
                    if(k<(E+1)) {
                        neighbors[k]=neighbors[k-1];
                    }
                }
                neighbors[j]=ii;
                if(nneigh<(E+1)) {
                    nneigh=nneigh+1;
                }
                trip=1;
                break;
            }
        }
        if((trip==0)&(nneigh<(E+1))&(ii!=i)&(neighbors[nneigh-1]!=ii)&(distances[ii]>0)) { /* Add element if vector is not yet full, but distance is greater than all already recorded */
            neighbors[nneigh]=ii;
            if(nneigh<(E+1)) {
				nneigh=nneigh+1;
            }
        }
    }
}
