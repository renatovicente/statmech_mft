//**********************************************************************
//metropolis.cpp
//This code implements a  Metropolis simulation for the 
//hamiltonean V(h,b)=-hb[1-(1-\delta)\Theta(hb)] in a BA network.
//The main goal is obtaining histograms for 
//m_k=|J_k.Z|/(|J_k||Z|) as a function of delta and alpha 
//July 25, 2011
//R. Vicente
//**********************************************************************

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <malloc.h>
#include <time.h> 
#include <unistd.h>

#include </usr/include/gsl/gsl_math.h>
#include </usr/include/gsl/gsl_blas.h>
#include </usr/include/gsl/gsl_rng.h>
#include </usr/include/gsl/gsl_randist.h>
#include </usr/include/gsl/gsl_vector.h>
#include </usr/include/gsl/gsl_matrix.h>
#include </usr/include/gsl/gsl_histogram.h>
#include </usr/include/gsl/gsl_sf_erf.h>
#include </usr/include/igraph/igraph.h>

#include "moraldyn.h"

#define SIZE N*N	   
#define DIMENSION 2        //lattice dimension
#define TORUS 1            //periodic boundary conditions (=1)  
#define NBINS 200 	       //number of histogram bins 
#define EPS 5
#define EPS2 0.00
#define STEPTEMP 0.2
#define NTEMP 21  
#define NUMCONECT 4

void help();

main (int argc,char *argv[])
{
    int ia,ib,ic,id,it,inow,ineigh,icont,ix,iy,itemp,in=0;
    int ia2,ia3,irun,icurrent,iflag;
    int P,L,N,TMAX,NRUNS,RP;
    double u,pm,pr,c;
    double alfa,sigma,aux,Q1,Q2,QZ,RZQ,rho,R;
    double hm_sigma,hm_mean,hr_sigma,hr_mean;	
    double DELTA,BETA,TCOLECT;
    double H,Delta,E,DeltaE,mE,mE2;
    double ORTOGONALFLAG; 
    long seed;

    FILE * mhist;
    FILE * rhist;
    FILE * out;
    FILE * pajek;

//**********************************************************************
// Help
//**********************************************************************
   if (argc<13){
	   help();
	   return(1);
   }
   else{
	BETA = atof(argv[1]);
    	DELTA = atof(argv[2]);
   	P = atoi(argv[3]);
        RP = atoi(argv[4]);
	L = atof(argv[5]); 
	N = atof(argv[6]);
	TMAX = atoi(argv[7]);
	TCOLECT = atoi(argv[8]);	
	ORTOGONALFLAG = atoi (argv[9]);
	mhist=fopen(argv[10],"w");
	rhist=fopen(argv[11],"w");
	out=fopen(argv[12],"w");	
   }  	  
   
//**********************************************************************
// Alocate matrices                                              
//**********************************************************************
    gsl_matrix * sociedade = gsl_matrix_alloc(SIZE,L);
    gsl_matrix * issue = gsl_matrix_alloc(P,L);
    gsl_vector * current_issue = gsl_vector_alloc(L);
    gsl_vector * v0 = gsl_vector_alloc(L);
    gsl_vector * v_trial = gsl_vector_alloc(L);
    gsl_vector * rnd_vector = gsl_vector_alloc(L);
    gsl_vector * Z = gsl_vector_alloc(L);  
    gsl_vector * opinions = gsl_vector_alloc(SIZE);
   	  	
//**********************************************************************
// Inicialization                                               
//**********************************************************************
    const gsl_rng_type * T;
    gsl_rng * r; 
  
    gsl_rng_env_setup();  
    T = gsl_rng_default;
    r=gsl_rng_alloc (T);

    seed = time (NULL) * getpid();
    //seed = 13188839657852;
    gsl_rng_set(r,seed);
   	
    igraph_t graph;
    igraph_vector_t neighbors;
    igraph_vector_t result;
    igraph_vector_t dim_vector; 
    igraph_real_t res;
    igraph_bool_t C;


    igraph_vector_init(&result,0);
    igraph_vector_init(&dim_vector,DIMENSION);
    igraph_vector_fill(&dim_vector,N); 

  
    gsl_histogram * histm = gsl_histogram_alloc (NBINS);
    gsl_histogram_set_ranges_uniform (histm,-1.1,1.1);
    
    gsl_histogram * histr = gsl_histogram_alloc (NBINS);
    gsl_histogram_set_ranges_uniform (histr,0.0,1.1);

   //*******************************************************************
   // Social Graph                                             
   //*******************************************************************
   //Barabasi-Alberts network
    igraph_barabasi_game(&graph,SIZE,RP,&result,1,0);
    
    /* for (inow=0;inow<SIZE;inow++){
	 igraph_neighbors(&graph,&neighbors,inow,IGRAPH_OUT);
	 printf("%d ",inow);
	 for(ic=0;ic<igraph_vector_size(&neighbors);ic++)
	 {	
	    	ineigh=(int)VECTOR(neighbors)[ic];
		printf("%d ",ineigh);
	 }
	  printf("\n");	
     }*/

     pajek=fopen("graph.xml","w");
     igraph_write_graph_graphml(&graph,pajek);

     //igraph_write_graph_pajek(&graph, pajek);
     fclose(pajek);
  
   //*******************************************************************
   // Quenched issues set and Zeitgeist
   //*******************************************************************	 
   gsl_vector_set_zero(Z);  
   gera_config(Z,issue,P,L,r,1.0);  
   if (ORTOGONALFLAG==1) gsl_matrix_set_identity(issue);
   for (ib=0;ib<P;ib++){
	gsl_matrix_get_row(current_issue,issue,ib);
        gsl_blas_ddot(current_issue,current_issue,&Q1);
	gsl_vector_scale(current_issue,1/sqrt(Q1));	
	gsl_vector_add(Z,current_issue);
    }
    gsl_blas_ddot(Z,Z,&QZ);
    gsl_vector_scale(Z,1/sqrt(QZ));		

    //******************************************************************
    // Initial ordered configuration
    //******************************************************************
    gera_config(Z,sociedade,SIZE,L,r,0);
    E = hamiltoneana(sociedade,issue,SIZE,L,P,DELTA,graph);

    //printf("GROUND STATE=%f\n",E);
	 
    //******************************************************************		      	
    // MC sweeps
    //******************************************************************	
    it=0;
    iflag=0;
    gsl_histogram_reset(histm);	
    gsl_histogram_reset(histr);
    mE=0;
    mE2=0;

//------>SIMULATION    	 				        SIMULATION<-------******	
    do{	
//-------->SWEEP                                 SWEEP<-------------****       			
	for(inow=0;inow<SIZE;inow++){
		//**************************************************************
		// current site and neighbors
		//**************************************************************
                igraph_vector_init(&neighbors,SIZE);
		inow=gsl_rng_uniform_int(r,SIZE);
		gsl_matrix_get_row(v0,sociedade,inow);
                igraph_neighbors(&graph,&neighbors,inow,IGRAPH_ALL);
		
                /*printf("A. Os vizinhos de INOW %d são: ",inow);
		for(ic=0;ic<igraph_vector_size(&neighbors);ic++)
	 	{	
	    		ineigh=(int)VECTOR(neighbors)[ic];
			printf("%d ",ineigh);
		 }
                printf("\n");*/
		
		//**************************************************************		
	 	// trial
		//**************************************************************
		gera_vetor(v_trial,L,r);
                //printf("B. Vetor aleatorio gerado\n");
		
		//**************************************************************
		// Energy change
		//**************************************************************
		DeltaE=variacaoE(v0,v_trial,inow,sociedade,issue,
					 	   N,L,P,DELTA,graph,neighbors);	      		
                //printf("C. DeltaE=%f\n",DeltaE);
		
		//**************************************************************
		// test new configuration
		//**************************************************************
		u=gsl_rng_uniform(r);		
		if(DeltaE<0){ 
			gsl_matrix_set_row(sociedade,inow,v_trial);	
                        //printf("D. Vetor aceito pois DeltaE<0\n\n");
		}			
		else if (u<exp(-BETA*DeltaE)){   	
			gsl_matrix_set_row(sociedade,inow,v_trial);
                        //printf("D. Vetor aceito u=%f\n\n",u);
		     }
                     //else  printf("D. Vetor rejeitado \n\n");
                 igraph_vector_destroy(&neighbors);
	}
//-------->SWEEP                                 SWEEP<-------------****
        E=hamiltoneana(sociedade,issue,SIZE,L,P,DELTA,graph);
	//printf("****** ------>   SWEEP%d E=%f\n",it,E);			
	it++; 
		
	//******************************************************************	
    //Gather statistics
    //******************************************************************
	if(it>TMAX){		       	 
		for(ia=0;ia<SIZE;ia++){
			gsl_matrix_get_row(v0,sociedade,ia);
			gsl_blas_ddot(v0,v0,&Q1);
			gsl_blas_ddot(v0,Z,&RZQ);
			pm = RZQ/sqrt(Q1*QZ);
			pr = fabs(RZQ)/sqrt(Q1*QZ);
			gsl_histogram_increment(histm,pm);
			gsl_histogram_increment(histr,pr);
       		}	
               mE +=E;
               mE2+=E*E;
               in++;
	} 
	if(it>TMAX+TCOLECT) iflag=1;
			    
    }while(iflag==0);   
//-------->SIMULATION                           SIMULATION<-----********

    mE=mE/(double)in;
    mE2=mE2/(double)in;

    c= ((BETA*BETA)/(double)SIZE)*(mE2-mE*mE);

   //******************************************************************* 
   //Statistical summary
   //*******************************************************************
   hm_sigma=gsl_histogram_sigma(histm);
   hm_mean=gsl_histogram_mean(histm); 
   
   hr_sigma=gsl_histogram_sigma(histr);
   hr_mean=gsl_histogram_mean(histr); 
   
   //*******************************************************************
   //Print results                                                  
   //*******************************************************************	
   fprintf(out,"%f %f %f %f %f %f\n",
							 BETA,hm_mean,hm_sigma,hr_mean,hr_sigma, c);
   gsl_histogram_fprintf(mhist,histm,"%g","%g"); 	  	   
   gsl_histogram_fprintf(rhist,histr,"%g","%g");	  	
	  	 
   //*******************************************************************
   // Finalizes                                                  
   //*******************************************************************
   igraph_destroy(&graph);
   igraph_vector_destroy(&result);
    
   gsl_matrix_free(issue);
   gsl_vector_free(current_issue);
   gsl_vector_free(opinions);
   gsl_vector_free(v0);
   gsl_vector_free(v_trial);
   gsl_matrix_free(sociedade);   	
   gsl_rng_free (r);
   
   fclose(mhist);
   fclose(rhist);
   fclose(out); 	
   
   return(0);
}

void help()
{
  printf("\n");
  printf("USAGE: metropolis BETA DELTA ISSUES M DIM SIZE ");
  printf("TMAX TCOLECT ORTFLAG FILEMHIST FILERHIST FILEOUT\n\n");
  printf("DELTA: confirmation bias\n");
  printf("ISSUES: number of simultaneously debated issues\n");
  printf("M: Barabasi-Alberts branching\n");
  printf("DIM: moral matrix dimension\n");
  printf("SIZE: lattice length (N_AGENTS=SIZE*SIZE) \n");
  printf("TMAX: number of sweeps TCOLECT: nb of samples \n");
  printf("ORTFLAG=1 for orthogonal issues \n\n");
  printf("OUTPUTs: BETA,hm_sigma,hm_mean,hr_sigma,hr_mean, c\n");
  printf("FILEs: mhist rhist out\n\n\n");
  exit(1);
}




