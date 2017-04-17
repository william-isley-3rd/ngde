// NGDE- Aerosol Dynamics solver.
// Authors: Anand Prakash, Ameya P. Bapat and Michael R. Zachariah, Mar 27, 2003.

// This C code is used to solve for nucleation, coagulation and surface growth problems. This example problem solves
// for characteristics of an Aluminum aerosol. However, it can be used for any other material whose properties are
// known. There are four different sections in this code:
//1) Coagulation, 
//2) Nucleation (classical) + coagulation, 
//3) Surface growth 
//4) Unified GDE with all the three phenomena combined. 

// The above phenomena constitute four different sections of the code. Each of the four sections are independent
// of each other. They have especially been written in such manner so that, one can identify the contribution of
// each phenomena to the GDE. Also, if one requires to use nucleation, coagulation or surface growth alone,
// it can be easily done.

// The solution algorithm involves a new approach to solve the Aerosol GDE, where the particle volume space is
// divided into nodes of zero width. Using nodes allows us to cover the entire size range of aerosol
// (1 nm to 10 micrometer) using just 40 nodes.A more detailed description of the theory and algorithm
// can be found in the reference:
// A. Prakash, A. P. Bapat and M.R. Zachariah,  Aerosol Science and Technology

// The main theoretical constraints of this implementation are:
// 1. Use of the Self Consistent Classical nucleation Model
// 2. Free molecule collision kernel
// 3. Free molecule surface growth


// The main variables in the code are v[MAX] and N [MAX]. These arrays are respectively the size (volume),
// and the number concentration of particles corresponding to that size. v1 is the monomer volume, which we assume,
// is equal to the volume of a molecule. We have tried to comment the code sufficiently to enable the user to make
// changes as needed. The main variables used are listed in the table below.

//**************************Note on compiling and running the code**************************
//To compile the code, use command "gcc gde.c -lm". To run the code the user should edit the input file as required
// and then the command "a.out > outputfile" should be used to direct the output to a file.
//******************************************************************************************

// *************************************NOTE ON UNITS***************************************
//                            SI Units are used throughout the code.
//******************************************************************************************
/*                                   ------------                                            */
/* 				                     NOMENCLATURE                                            */
/* 				                     ------------                                            */
/* __________________________________________________________________________________________*/
/* Variables            	            Meaning (Units)                                      */
/* __________________________________________________________________________________________*/
/* v[i]         	Volume of particles at ith node (m3/m3 of aerosol)                       */
/* N[i]         	Number of particles at ith node (#/m3 of aerosol)                        */
/* X[i][j][k]   	Splitting operator for coagulation (dimensionless)                       */
/* K[i][j]      	Collision frequency function(m3/s)                                       */
/* N1s[I]       	Saturation monomer concentration over particle of size I(#/m3)           */
/* sum1,sum2    	Addition and subtraction terms in Nk due to collision(#/(m3s))           */
/* Vav          	Average volume of particles (m3/m3 of aerosol)                           */
/* Vtotal       	Total volume of particles and monomer (m3/m3 of aerosol)                 */
/* T                Temperature	(K)                                                          */
/* coolrate     	Cooling rate	(K/s)                                                    */
/* t            	Time	(s)                                                              */
/* step	            Timestep for integration  (s)                                            */
/* Ninf	            Total number of particles as got by solving GDE	(#/m3 of aerosol)        */
/* Ninf_eq       	Total number of particles by eq 7.77, Friedlander (#/m3 of aerosol)      */
/* fi           	Total volume of particles (m3/m3 of aerosol)                             */
/* Vtot	            Total volume of particles (excluding monomers)	(m3/m3 of aerosol)       */
/* addterm,subterm	Addition and subtraction terms in Nk due to surface growth(#/(m3s))      */
/* Ntot         	Total number of particles (#/m3 of aerosol)                              */
/* dpav         	Mean Diameter	(m)                                                      */
/* v1           	Monomer volume	(m3/m3 of aerosol)                                       */
/* sigma        	Surface Tension	(N/m)                                                    */
/* Ps           	Saturation Pressure (Pa)                                                 */
/* ns           	Monomer concentration at saturation(#/m3 of aerosol)                     */
/* Ksp              Solublity product for solid in solution                                  */
/* Q                Initial concentration (M)                                                */
/* kstar        	Critical node size (dimensionless)                                       */
/* dpstar       	Critical particle size (m)                                               */
/* vstar        	Volume of particles in the critical bin	(m3/m3 of aerosol)               */
/* s1           	Surface area of monomer	(m2)                                             */
/* S            	Saturation ratio (dimensionless)                                         */
/* m1           	Monomer mass(Kg)                                                         */
/* zeta[i]      	Splitting operator for nucleation at the ith node (dimensionless)        */
/* t_SPD            Time required to reach a self-preserving distribution (s)                */

//*****************************EXAMPLE PROBLEM*********************************
// Aluminum vapor at a high temperature (1773 K) is being cooled at a rate ~1000 K/s. As the system cools down,
// supersaturation of vapor causes nucleation of particles, which then coagulate to form larger particles. Also, due to
// high concentration of monomer units (molecules), supersaturation causes surface growth of existing particles.
// The monomers either evaporate from the surface of existing particles, or condense on them.
//******************************************************************************

// Authors: Anand Prakash, Ameya P. Bapat and Michael R. Zachariah, Nov 27, 2002.

// Beginning of the code

#include<stdio.h>	
#include<math.h>
#include<time.h>
#define	pi	    3.1415926	    //pi
#define kb  	1.38e-23	    //Boltzmann constant
#define R       8.314           //Universal gas constant Joule/(K.mol)
#define Na      6.0225e23       //Avogadro's Number
#define MAX1    50              // Number of bins

int coagulation()
{
    //**********************************************
    //           PURE COAGULATION
    //**********************************************
    // This part of the code solves for coagulation problems. The code verifies that we obtain a self preserving
    // distribution. It verifies the problem stated in Friedlander, 2000 - page 218, line 4 - and shows that SPD
    // is obtained. And note that in this case node 1 is not the monomer. Its particle of size 1 nm. We do not have
    // to worry about monomers, since we are concerned with coagulation only.

	// Declaration of variables
	int i,j,k,choice,counter,counter2,opt,printopt,flag,MAX,beta_option;
	double v[MAX1+1],N[MAX1],X[MAX1][MAX1][MAX1],K[MAX1][MAX1],N1s[MAX1],sum1[MAX1],sum2[MAX1],temp,temp1,temp2,Vav,Vtotal;
	double t,step,temp3,Ninf,n,Vtot,addterm[MAX1],subterm[MAX1],Ntot,dpav,v1,fi,t_SPD,step_coag,q,P;
	double sigma,Ps,ns,theta,a,b,Jk,kstar,dpstar,vstar,s1,S,m1,c,zeta[MAX1],dp[MAX1],t1,t2,t3,t4,Kmin;
	double A,B,C,D,E,F,MW,rho,Ninf_eq,T,coolrate,d;
	double lambda,kn1,kn2,mu,D1,D2,m[MAX1],c1,c2,g1,g2,l1,l2,A_mu,B_mu;
    FILE *fptr,*fptr1,*fptr2;
    if((fptr=fopen("main.inp","r+"))==NULL)//Opening the input file that contains property data for the aerosol material.
    {
        printf("Warning: Could not open file");
        return(-1);
    }
    //Scanning the property data from the input file "main.inp"
    fscanf(fptr,"Number of nodes: %d Starting Temperature(Kelvin): %lf Cooling Rate(Kelvin/s): %lf Pressure(Atm): %lf Molecular Weight(Kg/mol): %lf Density(Kg/m3): %lf Surface Tension in the form A-BT(dyne/cm): A=%lf B=%lf Saturation Vapor Pressure of material in form Ps=exp(C-D/T)P: C=%lf D=%lf Enter your choice - 1)Coagulation,2)nucleation+coagulation,3)nucleation + coagulation + surface growth,4)surface growth:%d Choice of collision frequency function:1)Free Molecular Regime,2)Fuchs form for free molecular and transition regime:%d Sutherland's constants for carrier gas in the form A*T^1.5/(B+T): A=%lf B=%lf",&MAX,&T,&coolrate,&P,&MW,&rho,&A,&B,&C,&D,&choice,&beta_option,&A_mu,&B_mu);
    fclose(fptr);//Closing the main property data file.

    q=pow(10.0,(12.0/(MAX-2)));//Calculation of geometric spacing factor which depends on the number of nodes.
    v1 = MW/(rho*Na);// Volume of a monomer unit (a molecule); Na is the Avogadro's Number.

    if((fptr1=fopen("coag.inp","r+"))==NULL)//opening the input file for coagulation probelm.
    {
      printf("Warning: Could not open file");
      return(-1);
    }
      for(i=1;i<=MAX+1;i++)
	N[i]=0.0;
      //Scanning the initial monomer size and number concentration of the monodisperse distribution for coagulation.
      fscanf(fptr1,"\n Initial diameter of monodisperse particles(nm): %lf  Initial number concentration(#/m3): %lf Temperature(Kelvin): %lf",&d,&N[1],&T);
      v[1] = 1e-27*pi*pow(d,3)/6;// Setting diameter of particle that the user enters as the smallest node( coagulation would never lead to decrease in particle size).
      dp[1]=pow((6*v[1]/pi),0.3333333);
       m[1]=v[1]*rho;
      for(i=2;i<=MAX+1;i++) // Calculating volumes of larger nodes using a geometric spacing factor of q.
	{
	  v[i]=v[1]*pow(q,(i-1));
	  dp[i]=pow((6*v[i]/pi),0.3333333);
	  m[i]=v[i]*rho;
	}
      //Calculation of size splitting operators.
      for(k=1;k<=MAX;k++)
	{
	  for(i=1;i<=MAX;i++)
	    {
	      for(j=1;j<=MAX;j++)
		{
		  //Conditions in parentheses check if the combined volume of colliding particles is between k and k+1.
		  if (v[i]+v[j]>=v[MAX])
		    X[i][j][MAX]= (v[i]+v[j])/v[MAX];

		  else
		    {

		      if (v[k]<= (v[i]+v[j]) && (v[i]+v[j])<v[k+1])
			{
			  X[i][j][k]=(v[k+1]-v[i]-v[j])/(v[k+1]-v[k]);
			}
		      else
			{
			  if(v[k-1]<= (v[i]+v[j]) && (v[i]+v[j])<v[k])
			    {
			      X[i][j][k]=(v[i]+v[j]-v[k-1])/(v[k]-v[k-1]);
			    }
			  else
			    X[i][j][k]=0.0;
			}
		    }
		}
	    }
	}
      t=0.0; // Initializing time.
      step=1e-8;//Initializing timestep for integration.
      t_SPD=5/(pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(v[1],0.1666666)*N[1]);//Estimation of time required to reach SPD.
      printf("\n\nEstimated time required to reach SPD=%es\n",t_SPD);
      Ninf_eq=1e24;// A temporary variable that calculates N_total according to Friedlander, 2000, equation 7.77.
      printf("\nt\t\tNinf\t\tN[1]\t\tN[5]\t\tN[10]\t\tN[15]");
      while(t<1.e-2)//Running the coagulation mode of the code until SPD is reached. The dimensionless number distribution remains same after SPD has reached.
	{
	  //Calculation of collision frequency function ( same as in " nucleation + coagulation" section of the code).
	  Kmin=1e-9;//setting an arbitrary value for Kmin to start with.
	  if(beta_option==1)
	    {
	      for(i=1;i<=MAX;i++)
		{
		  for(j=1;j<=MAX;j++)
		    {

		      temp1=1/v[i] + 1/v[j];
		      temp2=pow(v[i],0.333333) + pow(v[j],0.333333);
		      K[i][j]= pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(temp1,0.5)*pow(temp2,2.0);
		      if (K[i][j]<Kmin)
			Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time.

		    }
		}
	    }
	  if(beta_option==2)
	    {
	      mu=A_mu*pow(T,1.5)/(B_mu+T);
	      lambda = (mu/(P*101325.0))* sqrt(pi*R*T/(2.0*0.04));
	      for(i=1;i<=MAX;i++)
		{
		  for(j=1;j<=MAX;j++)
		    {

		      kn1=(2.0*lambda)/dp[i];
		      kn2=(2.0*lambda)/dp[j];
		      D1=(kb*T)/(3.0*pi*mu*dp[i])*((5.0+4.0*kn1+6.0*kn1*kn1+18.0*kn1*kn1*kn1)/(5.0-kn1+(8.0+pi)*kn1*kn1));
		      D2=(kb*T)/(3.0*pi*mu*dp[j])*((5.0+4.0*kn2+6.0*kn2*kn2+18.0*kn2*kn2*kn2)/(5.0-kn2+(8.0+pi)*kn2*kn2));
		      c1=sqrt((8.0*kb*T)/(pi*m[i]));
		      c2=sqrt((8.0*kb*T)/(pi*m[j]));
		      l1=(8.0*D1)/(pi*c1);
		      l2=(8.0*D2)/(pi*c2);
		      g1=(pow((dp[i]+l1),3)-pow((dp[i]*dp[i]+l1*l1),1.5))/(3.0*dp[i]*l1)-dp[i];
		      g2=(pow((dp[j]+l2),3)-pow((dp[j]*dp[j]+l2*l2),1.5))/(3.0*dp[j]*l2)-dp[j];
		      K[i][j]=2.0*pi*(D1+D2)*(dp[i]+dp[j])/((dp[i]+dp[j])/(dp[i]+dp[j]+2.0*sqrt(g1*g1+g2*g2))+(8.0*(D1+D2))/(sqrt(c1*c1+c2*c2)*(dp[i]+dp[j])));
		      //  printf("kn %e\n",kn2);
		      if (K[i][j]<Kmin)
			Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time.

		    }
		}
	    }

	  //Calculating the gain and loss terms due to coagulation.
	  for(k=1;k<=MAX;k++)
	    {
	      sum1[k]=0;//Addition term when i and j collide to form a k sized particle.
	      sum2[k]=0;//Subtraction term when k collides with any other particle.
	      for(i=1;i<=MAX;i++)
		{
		  sum2[k]=sum2[k]+K[k][i]*N[i];
 		  for(j=1;j<=k;j++)
		    {
		      sum1[k]=sum1[k]+X[i][j][k]*K[i][j]*N[i]*N[j];
		    }

		}

	    }
	  Ninf=0;// Total number of particles, which is zero initially.
	  for(i=1;i<=MAX;i++)
	    {
	      N[i]= N[i] + step*(0.5*sum1[i] - N[i]*sum2[i]);
	      Ninf = Ninf + N[i];//N_infinity as we get by solving the coagulation part of GDE.
	    }
	  step = 1.e-3/(Kmin*Ninf);//Adaptive timestep for integration, based on characteristic coagulation time.
	  //step = 1.e-8;
	  fi=0;
	  for(i=1;i<=MAX;i++)
	    fi = fi + N[i]*v[i];
	  //N_infinity according to equation 7.77 in Friedlander.
	  Ninf_eq = Ninf_eq - step*(3.33*pow((3/(4*pi)),0.1666666)*pow(6*kb*T/rho,0.5)*pow(fi,0.1666666)*pow(Ninf_eq,(11.0/6.0)));
	  Vtot=0;
	  printf("\n%e\t%e\t%e\t%e\t%e\t%e",t,Ninf,N[1],N[5],N[10],fi);
	  t=t+step;
	}
      // Printing the final number distribution after SPD is reached.
      printf("\n\n***Size Distribution after reaching SPD****\n");
      printf("\nvolume\t\tnumber");
      for(i=1;i<=MAX-1;i++)
	printf("\n%e\t%e",v[i],N[i]);
      fclose(fptr1);
    return(0);
}
int coagulation_nucleation()
{
    // Declaration of variables
    int i,j,k,choice,counter,counter2,opt,printopt,flag,MAX,beta_option;
    double v[MAX1+1],N[MAX1],X[MAX1][MAX1][MAX1],K[MAX1][MAX1],N1s[MAX1],sum1[MAX1],sum2[MAX1],temp,temp1,temp2,Vav,Vtotal;
    double t,step,temp3,Ninf,n,Vtot,addterm[MAX1],subterm[MAX1],Ntot,dpav,v1,fi,t_SPD,step_coag,q,P;
    double sigma,Ps,ns,theta,a,b,Jk,kstar,dpstar,vstar,s1,S,m1,c,zeta[MAX1],dp[MAX1],t1,t2,t3,t4,Kmin;
    double A,B,C,D,E,F,MW,rho,Ninf_eq,T,coolrate,d;
    double lambda,kn1,kn2,mu,D1,D2,m[MAX1],c1,c2,g1,g2,l1,l2,A_mu,B_mu;
    FILE *fptr,*fptr1,*fptr2;
    if((fptr=fopen("main.inp","r+"))==NULL)//Opening the input file that contains property data for the aerosol material.
    {
        printf("Warning: Could not open file");
        return(-1);
    }
    //Scanning the property data from the input file "main.inp"
    fscanf(fptr,"Number of nodes: %d Starting Temperature(Kelvin): %lf Cooling Rate(Kelvin/s): %lf Pressure(Atm): %lf Molecular Weight(Kg/mol): %lf Density(Kg/m3): %lf Surface Tension in the form A-BT(dyne/cm): A=%lf B=%lf Saturation Vapor Pressure of material in form Ps=exp(C-D/T)P: C=%lf D=%lf Enter your choice - 1)Coagulation,2)nucleation+coagulation,3)nucleation + coagulation + surface growth,4)surface growth:%d Choice of collision frequency function:1)Free Molecular Regime,2)Fuchs form for free molecular and transition regime:%d Sutherland's constants for carrier gas in the form A*T^1.5/(B+T): A=%lf B=%lf",&MAX,&T,&coolrate,&P,&MW,&rho,&A,&B,&C,&D,&choice,&beta_option,&A_mu,&B_mu);
    fclose(fptr);//Closing the main property data file.
    q=pow(10.0,(12.0/(MAX-2)));//Calculation of geometric spacing factor which depends on the number of nodes.
    v1 = MW/(rho*Na);// Volume of a monomer unit (a molecule); Na is the Avogadro's Number.
    /*   //************************************************************ */
    /*   //             Nucleation + Coagulation                        */
    /*   //************************************************************ */

        v[1]=v1;//Assigning first node to monomer. Note that these are not considered as particles.
        for(i=1;i<=MAX+1;i++)
            N[i]=0.0;//Setting initial number of particles = 0 in all the nodes.
        for(i=2;i<=MAX+1;i++)
            //Defining volume spacing of bins. Using a geometric factor of 2 to cover 1nm to 10 micrometer size particles.
        {
            v[i]=v[1]*pow(q,(i-1));
            m[i]=rho*v[i];
            dp[i]=pow(6.0*v[i]/pi,0.3333333);
        }
        dp[1]=pow(6.0*v[1]/pi,0.3333333);
        m[1]=rho*v[1];
        //Calculation of splitting operator for coagulation, X[i][j][k].
        for(k=1;k<=MAX;k++)
        {
            for(i=1;i<=MAX;i++)
            {
                for(j=1;j<=MAX;j++)
                {
                    //Conditions in parentheses check if the combined volume of colliding particles is between k and k+1.
                    if (v[i]+v[j]>=v[MAX])
                        X[i][j][MAX]= (v[i]+v[j])/v[MAX];
                    else
                    {
                        if (v[k]<= (v[i]+v[j]) && (v[i]+v[j])<v[k+1])
                        {
                            X[i][j][k]=(v[k+1]-v[i]-v[j])/(v[k+1]-v[k]);
                        }
                        else
                        {
                            if(v[k-1]<= (v[i]+v[j]) && (v[i]+v[j])<v[k])
                            {
                                X[i][j][k]=(v[i]+v[j]-v[k-1])/(v[k]-v[k-1]);
                            }
                            else
                                X[i][j][k]=0.0;
                        }
                    }
                }
            }
        }

        t=0.0;//Initializing time.
        step=5e-4;//Timestep for integration.
        S=1.001;// Setting initial saturation ratio, a little larger than 1. Hereafter the saturation ratio is determined by cooling or heating rate of the aerosol. A heating rate would require a negative value.

        // Calculating surface area of the monomer, given the volume of monomer v1.
        temp1 = 6*v1/pi;
        s1 = pi*pow(temp1,0.66666666667);
        //Calculating mass of the monomer unit.
        m1=rho*v1;
        //Calculating saturation vapor pressure for the Aerosol.
        Ps=exp(13.07-36373.0/T)*101325.0*P;
        //Calculating number concentration of monomers at saturation.
        ns=Ps/(kb*T);
        // Calculating monomer concentration ( multiplying saturation ratio by saturation concentration of monomers).
        N[1]=S*ns;
        counter=0; // Temporary variable for printing after a large number of iterations.
        printf("\ntime\t\tMonomer\t\tJk\t\tS\t\tDp_mean\t\tk*\t\tNtot");
        while(T>300) // Calculating aerosol properties until the system cools down to 27 C.
        {
            sigma=(A-B*T)*1e-3; // Surface tension.
            Ps=exp(C-D/T)*101325.0*P; // Saturation pressure.
            ns=Ps/(kb*T); // Saturation concentration of monomers using ideal gas law.
            S=N[1]/ns; // Saturation ratio.
            theta = (s1*sigma)/(kb*T);// Dimensionless surface tension theta.
            a=(2*sigma)/(pi*m1); //Temporary variable.
            b= theta - (4*pow(theta,3))/(27*pow(log(S),2));//Temporary variable.
            Jk = pow(ns,2)*S*v1*pow(a,0.5)*exp(b); //Nucleation rate using the classical SCC model.
            c = 0.6666667*theta/log(S);// Temporary variable.
            kstar = pow(c,3.0);//Calculating critical cluster size that will determine the node at which nucleation occurs.
            dpstar= 4*sigma*v1/(kb*T*log(S));//Size of particles corresponding to the critical node.
            vstar=pi*pow(dpstar,3)/6; // Volume of particle corresponding to the critical node.
            Ntot=0.0; // Initializing total number of particles to zero, as there are no particles prior to nucleation.

            for(i=2;i<=MAX;i++)
                Ntot = Ntot + N[i];//Total number of particles (does not include monomers).
            Vtot=0.0;
            for(i=2;i<=MAX;i++)
                Vtot = Vtot + N[i]*v[i];//Total volume of particles. Note that loop runs from i=2 because i=1 corresponds to monomers, which we do not count as particles.
            Vtotal = Vtot + N[1]*v[1];//Total volume for mass conservation check (includes monomers)
            Vav = Vtot/Ntot; // Average volume of particles ( number average)
            dpav = pow((6*Vav/pi),0.3333333);//Volume based mean diameter of particles printing after every 50 times the loop runs.
            if(t>0.14)
                step=5e-8;
            counter=counter + 1;
            if (counter==50)
            {
                printf("\n%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e",t,N[1],Jk,S,dpav,kstar,Vtotal,Ntot);
                counter=0;
            }
            //Calculation of collision frequency function beta(i,j)=K(i,j).

            if(beta_option==1)
            {
                for(i=1;i<=MAX;i++)
                {
                    for(j=1;j<=MAX;j++)
                    {

                        temp1=1/v[i] + 1/v[j];
                        temp2=pow(v[i],0.333333) + pow(v[j],0.333333);
                        K[i][j]= pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(temp1,0.5)*pow(temp2,2.0);
                        if (K[i][j]<Kmin)
                            Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time.

                    }
                }
            }
            if(beta_option==2)
            {
                mu=A_mu*pow(T,1.5)/(B_mu+T);
                lambda = (mu/(P*101325.0))* sqrt(pi*R*T/(2.0*0.04));
                for(i=1;i<=MAX;i++)
                {
                    for(j=1;j<=MAX;j++)
                    {
                        kn1=(2.0*lambda)/dp[i];
                        kn2=(2.0*lambda)/dp[j];
                        D1=(kb*T)/(3.0*pi*mu*dp[i])*((5.0+4.0*kn1+6.0*kn1*kn1+18.0*kn1*kn1*kn1)/(5.0-kn1+(8.0+pi)*kn1*kn1));
                        D2=(kb*T)/(3.0*pi*mu*dp[j])*((5.0+4.0*kn2+6.0*kn2*kn2+18.0*kn2*kn2*kn2)/(5.0-kn2+(8.0+pi)*kn2*kn2));
                        c1=sqrt((8.0*kb*T)/(pi*m[i]));
                        c2=sqrt((8.0*kb*T)/(pi*m[j]));
                        l1=(8.0*D1)/(pi*c1);
                        l2=(8.0*D2)/(pi*c2);
                        g1=(pow((dp[i]+l1),3)-pow((dp[i]*dp[i]+l1*l1),1.5))/(3.0*dp[i]*l1)-dp[i];
                        g2=(pow((dp[j]+l2),3)-pow((dp[j]*dp[j]+l2*l2),1.5))/(3.0*dp[j]*l2)-dp[j];
                        K[i][j]=2.0*pi*(D1+D2)*(dp[i]+dp[j])/((dp[i]+dp[j])/(dp[i]+dp[j]+2.0*sqrt(g1*g1+g2*g2))+(8.0*(D1+D2))/(sqrt(c1*c1+c2*c2)*(dp[i]+dp[j])));
                        if (K[i][j]<Kmin)
                            Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time.
                    }
                }
            }

            //Operator to put nucleated particles in the bin just higher than k*.
            for(k=2;k<=MAX;k++)
            {
                if (vstar<v1)
                    zeta[2]=vstar/v[2];//Putting particles formed smaller than monomers in the smallest particle node (node 2). This situation arises when k* falls below monomer size.
                else if (v[k-1] <=vstar && vstar < v[k])
                    zeta[k]=vstar/v[k];// Putting particles in node just larger than k*.
                else
                    zeta[k]=0.0;
            }

            //Calculation of gain and loss terms for Nk due to coagulation.
            for(k=2;k<=MAX;k++)
            {
                //Initializing gain and loss terms for coagulation.
                sum1[k]=0;// Production term due to collision of two smaller particles to form a k size particle.
                sum2[k]=0;// Loss term due to collision of k size particle with any other particle.
                for(i=2;i<=MAX;i++)
                {
                    sum2[k]=sum2[k]+K[k][i]*N[i]; // k collides with any other particle to get out of node k, thus loss term.
                    for(j=2;j<=k;j++)
                        sum1[k]=sum1[k]+X[i][j][k]*K[i][j]*N[i]*N[j];// i and j collide to form particle of size k which is then multiplied by the size splitting operator to put the particles in the adjacent nodes after adjusting the volume.
                }
            }

            // Change in PSD due to nucleation + coagulation.
            for(k=2;k<=MAX;k++)
                N[k]= N[k] + step*(0.5*sum1[k] - N[k]*sum2[k] + Jk*zeta[k]);

            //Monomer balance. accounting for the monomers loss due to nucleation.
            N[1] = N[1] - Jk*kstar*step;
            T=T-step*coolrate;// Temperature would decrease as time progresses, due to cooling.
            t=t+step;//Time increment.
        }
        printf("\n\nFinal Particle Size Distribution\n\n");
        printf("\nvolume\t\tnumber");
        for(i=2;i<=MAX;i++)
            printf("\n%e\t%e",v[i],N[i]);// Printing the final number distribution obtained due to nucleation and coagulation. (number concentration or particles against nodes).
    return(0);
}
int coagulation_nucleation_surface()
{
    /*   //************************************************************ */
    /*   //             Nucleation + Coagulation + Surface              */
    /*   //************************************************************ */
    // Declaration of variables
    int i,j,k,choice,counter,counter2,MAX,beta_option,restart_option;
    // int opt,printopt,flag;
    double v[MAX1+1],N[MAX1],K[MAX1][MAX1],XxK[MAX1][MAX1][MAX1],
            N1s[MAX1],sum1[MAX1],sum2[MAX1], temp, temp1, temp2, Vav, Vtotal,
            t,step,Vtot,addterm[MAX1],subterm[MAX1],Ntot,dpav,v1,q,
            sigma,ns,theta,a,b,Jk,kstar,dpstar,vstar,s1,S,m1,c,zeta[MAX1],dp[MAX1],t1,t2,Kmin,
            A,B,MW,rho,T,coolrate,
            QM, Qsol, Ksp, viscosity, alpha, wetangle, temp_t1_k, temp_t2_k, time_taken,
            m[MAX1];
    // double kn1,kn2,c1,c2,D1,D2,g1,g2,l1,l2,A_mu,B_mu,lambda,mu,C,D,E,F,t3,t4,Ninf,Ps,P,n,fi,t_SPD,Ninf_eq,dstep_coag;
    int XxK_loop_list[1514],XxK_loop_counter=0,XxK_loop_val,XxK_looper;
    int XxK_loop_size = sizeof(XxK_loop_list) / sizeof(int);
    clock_t real_time;
    FILE *fptr,*fdout,*fbinsout,*frestart;
    if((fptr=fopen("main.inp","r+"))==NULL)//Opening the input file that contains property data for the aerosol material.
    {
        printf("Warning: Could not open file");
        return(-1);
    }
    //Scanning the property data from the input file "main.inp"
    printf("File Read");
    fscanf(fptr,"Number of nodes: %d Starting Temperature(Kelvin): %lf Cooling Rate(Kelvin/s): %lf Concentration(M): %lf "
                   "Molecular Weight(Kg/mol): %lf Density(Kg/m3): %lf Viscosity(Pa*s): %lf "
                   "Surface Tension in the form A-BT(dyne/cm): A=%lf B=%lf Solubility Product Ksp: %lf "
                   "Enter your choice - 1)Coagulation,2)nucleation+coagulation,3)nucleation + coagulation + surface growth,4)surface growth:%d "
                   "Choice of collision frequency function:1)Free Molecular Regime,2)Fuchs form for free molecular and transition regime,3)Stokes-Einstein:%d "
                   "Heterogeneous Wetting Angle: %lf Restart: %d",
           &MAX,&T,&coolrate,&QM,&MW,&rho,&viscosity,&A,&B,&Ksp,&choice,&beta_option,&wetangle,&restart_option);
    fclose(fptr);//Closing the main property data file.
    q=pow(10.0,(8.0/(MAX-2)));//Calculation of geometric spacing factor which depends on the number of nodes.
    v1 = MW/(rho*Na);// Volume of a monomer unit (a molecule); Na is the Avogadro's Number.

    if (restart_option == 1) {
        // need to read the bins to get initial # of particles and time
        frestart=fopen("nucl-coag-surf-restart.dat","r+");
        if(frestart == NULL)//Opening the input file that contains restart info
        {
            printf("Warning: Could not open restart file");
            return(-1);
        }
        double restart_array[MAX+1];
        for (i = 0; i <=MAX;i++) {
            fscanf(frestart, "%lf", &restart_array[i]);
        }
        t=restart_array[0];  // reading in time
        if(t>0.14)
            step=1e-9;
        else
            step=1e-5;
        for (i = 1; i<=MAX; i++){
            N[i] = restart_array[i];  // setting # of particles per bin
            N1s[i] = 0.0; // initializing values, will be compute below before they are used
        }
        fclose(frestart);
        fdout = fopen("nucl-coag-surf-out.dat", "w");
        if (fdout == NULL){
            printf("Warning: Could not open outfile");
            return(1);
        }
        fbinsout = fopen("nucl-coag-surf-bins.dat", "w");
        if (fbinsout == NULL){
            printf("Warning: Could not open bin file");
            return(1);
        }
    }
    else {
        for(i=1;i<=MAX+1;i++)
        {
            N[i]=0.0; // Setting number of particles = 0 at all the nodes initially.
            N1s[i]=0.0; // Setting number concentration of monomers over different size particles = 0 initially.
        }
        N[1]=QM*Na/1000;  // Initial #of metal monomers from input file (converted from M to molecules / m3)
        t=0.0;// Initializing time.
        step=1e-5; // Initial time step for integration.
        fdout = fopen("nucl-coag-surf-out.dat", "w");
        if (fdout == NULL){
            printf("Warning: Could not open outfile");
            return(1);
        }
        fbinsout = fopen("nucl-coag-surf-bins.dat", "w");
        if (fbinsout == NULL){
            printf("Warning: Could not open binfile");
            return(1);
        }
    }

    // Setting volume, diameters, masses for all Nodes
    v[1]=v1; // Setting volume of the first node equal to the monomer volume, i.e. volume of a molecule, which is v1.
    for(i=2;i<=MAX+1;i++)
    {
        v[i]=v[1]*pow(q,(i-1)); // Calculating volume of all larger nodes using a geometric factor of 2. Note that volume of node 1 is the volume of a molecule.
        dp[i]=pow((6.0*v[i]/pi),0.333333);
        m[i]=rho*v[i];
    }
    dp[1]=pow((6.0*v[1]/pi),0.333333); //diameter
    m[1]=rho*v[1];  //mass

    //for(i=2;i<MAX;i++){
    //    vfrac[i] = v[1]/(v[i+1]-v[i]);
    //}
    //Size splitting operator for coagulation ( Same as in "coagulation + nucleation" section of the code).
    //Now includes the collision frequency (as temperature isn't changing during simulation)
    for(k=1;k<=MAX;k++)
    {
        //printf('%f\n',XxK_loop_counter);
        for(i=1;i<=MAX;i++)
        {
            for(j=1;j<=MAX;j++)
            {
                //Collision frequency function for coagulation.
                // Removed beta options 1 and 2 to end of code for optimization...
                // Computing the collision frequency here
                //if(beta_option==3)
                temp1=1/pow(v[i],0.333333) + 1/pow(v[j],0.333333);
                temp2=pow(v[i],0.333333) + pow(v[j],0.333333);
                K[i][j]= ((2*kb*T)/(3*viscosity))*temp1*temp2;
                //Calculating the smallest collision frequency function to decide the characteristic coagulation time.
                if (K[i][j]<Kmin)
                    Kmin = K[i][j];

                if (v[i]+v[j]>=v[MAX]) {
                    XxK[i][j][MAX] = (v[i] + v[j]) / v[MAX] * K[i][j];
                    //XxK_loop_val = MAX * (MAX * MAX) + MAX * j + i;
                    //XxK_loop_list[XxK_loop_counter] = XxK_loop_val;
                    //XxK_loop_counter = XxK_loop_counter + 1;
                   // printf('loop index is %d \n', XxK_loop_counter);
                }
                else
                {
                    if(v[k]<= (v[i]+v[j]) && (v[i]+v[j])<v[k+1])
                    {
                        XxK[i][j][k]=(v[k+1]-v[i]-v[j])/(v[k+1]-v[k]) * K[i][j];
                        //XxK_loop_val = (int) k * (MAX * MAX) + MAX * j + i;
                        //XxK_loop_list[XxK_loop_counter] = XxK_loop_val;
                        //XxK_loop_counter = XxK_loop_counter + 1;
                    }
                    else
                    {
                        if(v[k-1]<= (v[i]+v[j]) && (v[i]+v[j])<v[k])
                        {
                            XxK[i][j][k]=(v[i]+v[j]-v[k-1])/(v[k]-v[k-1]) * K[i][j];
                            //XxK_loop_val = k * (MAX * MAX) + MAX * j + i;
                            //XxK_loop_list[XxK_loop_counter] = XxK_loop_val;
                            //XxK_loop_counter = XxK_loop_counter + 1;
                        }
                        else {
                            XxK[i][j][k] = 0;
                        }
                    }
                }// end else for XxK
            }// end for j
        }// end for i
    }// end for k


    Qsol=2*pow(Ksp/(4*27),0.2)*Na/1000;  //Compute number of soluble metal particles (converted from M to molecules /m3)
    ns=Qsol;                             // computing saturation concentration from equilibrium solubility product

    // Calculating surface area of the monomer, given the volume of monomer v1.
    temp1 = 6*v1/pi;                     // Temporary variable.
    s1 = pi*pow(temp1,0.66666666667);    // Surface area of monomer.

    m1=rho*v1;                           // Mass of a monomer.
    //S=N[1]/ns;                           // Saturation ratio

    alpha = (2+cos(wetangle*(pi/180)))*
            pow((1-cos(wetangle*(pi/180))),2)/4;  //heterogeneous modification to nucleation
    Kmin=1e-9;//setting an arbitrary value for Kmin to start with.
    sigma= A*1e-3;                                  // Surface tension.
    theta = (s1*sigma)/(kb*T);                      // Dimensionless surface tension theta.
    a=(2*sigma)/(pi*m1);                            // Temporary variable. (T dependence, surface tension)
    // Printing section for outfiles
    printf("\nTime\t\tN[1]\t\tJk\t\tS\t\tdiameter\tParticle Volume\tTotal Number\tReal t");// Printing result heading.
    fprintf(fdout, "\nTime\t\tN[1]\t\tJk\t\tS\t\tKstar\t\tdiameter\tParticle Volume\tTotal Number");// Printing result heading.
    fprintf(fbinsout, "\n%-16s", "Time"); // Printing bin header
    for (i=1;i<=MAX;i++) {
        //Calculating the saturation monomer concentration over an i sized particle that has a diameter dp[i].
        dp[i]=pow((6*v[i]/pi),0.333333); // Diameter of particle at the ith node.
        N1s[i]=ns*exp(4*sigma*MW/(R*T*rho*dp[i]));// Calculating saturation monomer concentration over a i size particle including the kelvin effect.
        fprintf(fbinsout, "N%-15d", i);  //Printing bin header
    }
    real_time = clock();
    // Counters below are temporary variables that determine number of time the loop runs before every print step.
    // Generally time steps are small and printing results is done after every 4000 steps.
    counter=1;
    counter2=1;
    while(t<240.0)
    {
        // Update terms with new Temperature & monomer pop Here
        S=N[1]/ns;                                      // Saturation ratio.
        b= alpha*(theta -
                  (4*pow(theta,3))/(27*pow(log(S),2))); // Temporary variable. (T dependence, surface tension, saturation)
        Jk = pow(ns,2)*S*v1*pow(a,0.5)*exp(b);          // Nucleation rate using the classical SCC model.
        c = 0.6666667*theta/log(S);                     // Temporary variable.
        kstar = pow(c,3.0);                             // Calculating critical cluster size that will determine the node at which nucleation occurs.
        dpstar= 4*sigma*v1/(kb*T*log(S));               // Size of particles corresponding to the critical node.
        vstar=pi*pow(dpstar,3)/6;                       // Volume of particle corresponding to the critical node.
        Ntot=0.0;                                       // Initializing total number of particles to zero, as there are no particles prior to nucleation.

        Vtot=0.0;
        for(i=2;i<=MAX;i++) {
            Ntot = Ntot + N[i];//Total number of particles.
            Vtot = Vtot + N[i] * v[i];//Total volume of particles . This volume is used to evaluate the mean volume and thus mean diameter of particles.
        }
        Vtotal = Vtot + N[1]*v[1];//Total mass volume (including monomers) This volume is used to check for mass conservation.
        Vav = Vtot/Ntot; // Mean volume of particles.
        dpav = pow((6*Vav/pi),0.3333333);//Mean diameter of particles.

        //Reducing the timestep for integration when there are sufficient number of particles such that surface growth
        // takes over and the code can no longer run for large time-steps.
        if(t>0.14)
            step=1e-9;
        counter = counter +1;
        if(counter == 50000)// Printing mean diameter etc. after every 10,000 steps.
        {
            real_time = clock() - real_time;
            time_taken = ((double)real_time)/CLOCKS_PER_SEC;
            //printf("\n%i while loops took %f seconds to execute\n", counter, time_taken);
            real_time = clock();
            printf("\n%.9f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f",t,N[1],Jk,S,dpav,Vtotal,Ntot,time_taken);
            fprintf(fdout,"\n%.9f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e",t,N[1],Jk,S,kstar,dpav,Vtotal,Ntot);
            counter = 0;
            counter2 = counter2 +1;
            if (counter2 == 10)// Printing PSD after every 500,000 steps.
            {
                fprintf(fbinsout,"\n%.9f\t",t);
                for (i=1;i<=MAX;i++)
                    fprintf(fbinsout, "%.6e\t", N[i]);
                fflush(fbinsout);
                fflush(fdout);
                counter2=1;
            }

        }// end print loop
        for(i=2;i<=MAX;i++)
        {
            // Putting nucleated particles in bins just larger than k*
            // ( explained in section " coagulation + nucleation" above).
            // Putting particles formed smaller than monomers in the smallest particle node (node 2).
            // This situation arises when k* falls below monomer size.
            if (vstar<v1) {
                zeta[2] = vstar / v[2];
            }
            else if (v[i-1] <=vstar && vstar < v[i])
                zeta[i]=vstar/v[i]; // Putting particles in node just larger than k*.
            else
                zeta[i]=0.0;
        } // end for i, zeta computing

        // t1 through t2 are terms used for monomer balance. These are loss and gain terms due to condensation and
        // evaporation respectively. Each of these terms has been explained below as I use them.
        // Initially setting these terms = 0.
        t1=0;  // Loss of monomers that have condensed.
        t2=0;  // Gain of monomers that have evaporated.

        // Note that in the following part "addterm" and "subterm" are addition and subtraction term in particle size
        // distribution due to condensation or evaporation (surface growth), where as terms t1 through t4 are addition
        // and subtraction term for monomer balance.
        for(k=2;k<=MAX;k++)
        {
            sum1[k]=0.0;//Addition term due to coagulation of size i+j to form k.
            sum2[k]=0.0;//Loss term due to coagulation of k with any other particle.

            addterm[k] = 0.0;
            subterm[k] = 0.0;

            // Activating surface growth when sufficient number of particles have nucleated. At smaller times there are
            // very small number of particles and surface growth hardly changes PSD. Doing so makes the code run faster.
            if(step < 2e-5)
            {
                // The following four "if" statements check if the vapor pressure of monomer in the system is larger or
                // smaller than saturation vapor pressure over the k size particle. The difference in the two vapor
                // pressures is the driving force for evaporation or condensation.

                // i.e. if actual monomer concentration in the aerosol is greater than saturation concentration over a k-1 size particle.
                if(N[1]>N1s[k-1])  // term1: condensation on k-1 particle, gain in k, loss in monomer
                {
                    if (k==2)
                        addterm[k]=0.0;
                    else
                    {
                        temp_t1_k = K[1][k-1]*(N[1]-N1s[k-1])*N[k-1];
                        addterm[k] = addterm[k] + (v[1]/(v[k]-v[k-1]))*temp_t1_k;//Growth of k due to condensation of monomers on k-1.
                        t1 = t1 + -temp_t1_k;//Loss of monomers that have condensed.
                    }
                } // end N[1]>N1s[k-1]

                // Similarly, we can explain the following terms below.
                if( N[1]<N1s[k+1] )  // term2: evaporation of k+1 particle, gain in k, gain in monomer
                {
                    if(k==MAX)
                        subterm[k] = 0.;
                    else
                    {
                        temp_t2_k = K[1][k+1]*(N[1]-N1s[k+1])*N[k+1];
                        addterm[k] = addterm[k] -(v[1]/(v[k+1]-v[k]))*temp_t2_k;//Growth of k due to evaporation of monomers from k+1.
                        t2 = t2 + -temp_t2_k;//Gain of monomers that have evaporated.
                    }
                } // end N[1]<N1s[k+1]

                if(N[1]>N1s[k])  //term 3, condensation on k particle, loss in k, loss in monomer (balanced in term2)
                {
                    if(k==MAX)
                        subterm[k] = 0.;
                    else
                        subterm[k] = subterm[k] + (v[1]/(v[k+1]-v[k]))*K[1][k]*(N[1]-N1s[k])*N[k];
                    //Loss of k due to condensation of monomers on k. (monomers balanced above)
                } // end N[1]>N1s[k]

                if(N[1]<N1s[k])  //term 4: evaporation from k particle, loss in k, gain in monomer (balanced in term1)
                {
                    if (k==2)
                    {
                        double temp_t2_k = K[1][k]*(N[1]-N1s[k])*N[k];
                        subterm[k] = subterm[k] - (v[1]/(v[k]-v[k-1]))*temp_t2_k;
                        t2 = t2 + -temp_t2_k;//Gain of monomers that have evaporated.
                    }
                    else
                        subterm[k] = subterm[k] - (v[1]/(v[k]-v[k-1]))*K[1][k]*(N[1]-N1s[k])*N[k];
                    //Loss of k due to evaporation of monomers from k. (monomers balanced above)
                } // end N[1]<N1s[k]
            } // end if step < 2e-5

            // Terms due to coagulation.
            for(i=2;i<=MAX;i++)
            {
                sum2[k]=sum2[k]+K[k][i]*N[i]; // k collides with any other particle to get out of node k, thus loss term.
                for(j=2;j<=k;j++) {
                    sum1[k] = sum1[k] + XxK[i][j][k] * N[i] * N[j];
                    // sum1[k] = sum1[k] + X[i][j][k] * K[i][j] * N[i] * N[j];
                    // i and j collide to form particle of size k which is then multiplied by the size splitting
                    // operator to put the particles in the adjacent nodes after adjusting the volume.
                }
            }// end for i, for coagulation
        }  // end for k, over surface and coag

        /*for(XxK_looper=0; XxK_looper< XxK_loop_size; XxK_looper++){
            XxK_loop_val = XxK_loop_list[XxK_looper];
            i = XxK_loop_val % MAX;
            j = ((XxK_loop_val - i)/MAX) % MAX;
            k = (((XxK_loop_val - i)/MAX) - j) / MAX;
            //printf('%d\n', XxK_loop_val);
            sum1[k] = sum1[k] + XxK[i][j][k] * N[i] * N[j];
        }*/

        // Change in PSD due to Nucleation, Coagulation and Surface Growth.
        for(k=2;k<=MAX;k++) {
            N[k] = N[k] + step * (0.5 * sum1[k] - N[k] * sum2[k] + Jk * zeta[k] + addterm[k] - subterm[k]);
        }
        // Change in the concentration of monomers as they are used up in nucleation and surface growth.
        N[1] = N[1] - Jk*kstar*step - step*(t1 -t2);//Monomer balance equation.

        T=T-step*coolrate;//Reducing the temperature at the given cooling rate.
        t=t+step;// Time increment.
    } // end while t < 240
    fclose(fdout);
    fclose(fbinsout);
    return(0);
}
int surface_only()
{
    // Declaration of variables
    int i,j,k,choice,counter,counter2,opt,printopt,flag,MAX,beta_option;
    double v[MAX1+1],N[MAX1],X[MAX1][MAX1][MAX1],K[MAX1][MAX1],N1s[MAX1],sum1[MAX1],sum2[MAX1],temp,temp1,temp2,Vav,Vtotal;
    double t,step,temp3,Ninf,n,Vtot,addterm[MAX1],subterm[MAX1],Ntot,dpav,v1,fi,t_SPD,step_coag,q,P;
    double sigma,Ps,ns,theta,a,b,Jk,kstar,dpstar,vstar,s1,S,m1,c,zeta[MAX1],dp[MAX1],t1,t2,t3,t4,Kmin;
    double A,B,C,D,E,F,MW,rho,Ninf_eq,T,coolrate,d;
    double lambda,kn1,kn2,mu,D1,D2,m[MAX1],c1,c2,g1,g2,l1,l2,A_mu,B_mu;
    FILE *fptr,*fptr1,*fptr2;
    if((fptr=fopen("main.inp","r+"))==NULL)//Opening the input file that contains property data for the aerosol material.
    {
        printf("Warning: Could not open file");
        return(-1);
    }
    //Scanning the property data from the input file "main.inp"
    fscanf(fptr,"Number of nodes: %d Starting Temperature(Kelvin): %lf Cooling Rate(Kelvin/s): %lf "
            "Pressure(Atm): %lf Molecular Weight(Kg/mol): %lf Density(Kg/m3): %lf "
            "Surface Tension in the form A-BT(dyne/cm): A=%lf B=%lf "
            "Saturation Vapor Pressure of material in form Ps=exp(C-D/T)P: C=%lf D=%lf "
            "Enter your choice - 1)Coagulation,2)nucleation+coagulation,3)nucleation + coagulation + surface growth,4)surface growth:%d "
            "Choice of collision frequency function:1)Free Molecular Regime,2)Fuchs form for free molecular and transition regime:%d "
            "Sutherland's constants for carrier gas in the form A*T^1.5/(B+T): A=%lf B=%lf",
           &MAX,&T,&coolrate,&P,&MW,&rho,&A,&B,&C,&D,&choice,&beta_option,&A_mu,&B_mu);
    fclose(fptr);//Closing the main property data file.
    q=pow(10.0,(12.0/(MAX-2)));//Calculation of geometric spacing factor which depends on the number of nodes.
    v1 = MW/(rho*Na);// Volume of a monomer unit (a molecule); Na is the Avogadro's Number.

    for(i=1;i<=MAX+1;i++)
    {
        N[i]=0.0;
        N1s[i]=0.0;
    }
    v[1]=v1;
    for(i=2;i<=MAX+1;i++)
        v[i]=v[1]*pow(q,(i-1));

    if((fptr2=fopen("grow.inp","r+"))==NULL)//Opening the input file for surface growth problems
    {
        printf("Warning: Could not open file");
        return(-1);
    }

    //Scanning the initial size and number concentration of particles whose surface growth is to be studied.
    fscanf(fptr2,"Enter the node in which you want to put in particles:%d Enter the particle number concentration in #/m3 initially:%lf",&i,&temp);
    fclose(fptr2);
    N[i]=temp;
    t=0.0;
    step=1e-8;
    S=1.001;
    temp1 = 6*v1/pi;
    s1 = pi*pow(temp1,0.66666666667);
    m1=rho*v1;
    Ps=exp(C-D/T)*101325.0*P;
    ns=Ps/(kb*T);
    N[1]=S*ns;
    counter=0;

    while(T>300)
    {
        sigma=(A-B*T)*1e-3;
        Ps=exp(C-D/T)*101325.0*P;
        ns=Ps/(kb*T);
        S=N[1]/ns;
        theta = (s1*sigma)/(kb*T);
        a=(2*sigma)/(pi*m1);
        b= theta - (4*pow(theta,3))/(27*pow(log(S),2));
        Jk = pow(ns,2)*S*v1*pow(a,0.5)*exp(b);
        c = 0.6666667*theta/log(S);
        kstar = pow(c,3.0);
        if(counter==1000)//printing after every 1000 steps
        {
            printf("\n\ntime=%es\n",t);
            for(i=1;i<=MAX-2;i++)
                printf("\nN[%d]\t%e",i,N[i]);
            counter=0;
        }
        counter=counter +1;

        //collision frequency functions for the monomers
        for(j=1;j<=MAX;j++)
        {
            temp1=1/v[1] + 1/v[j];
            temp2=pow(v[1],0.333333) + pow(v[j],0.333333);
            K[1][j]= pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(temp1,0.5)*pow(temp2,2.0);
        }

        //calculation of saturation concentration of monomers over an i sized particle
        for(i=1;i<=MAX;i++)
        {
            dp[i]=pow((6*v[i]/pi),0.333333);
            N1s[i]=ns*exp(4*sigma*MW/(R*T*rho*dp[i]));
        }

        // t1 through t4 are terms used for monomer balance. These are loss and gain terms due to condensation and evaporation respectively. Each of these terms has been explained below as I use them. Initially setting these terms = 0.
        t1=0;
        t2=0;

        for(k=2;k<=MAX;k++)
        {
            addterm[k] = 0;
            subterm[k] = 0;

            if(step!=1e-4)
            {
                // The following four "if" statements check if the vapor pressure of monomer in the system is larger or smaller than saturation vapor pressure over the k size particle. The difference in the two vapor pressures is the driving force for evaporation or condensation.

                if(N[1]>N1s[k-1]) // i.e. if actual monomer concentration in the aerosol is greater than saturation concentration over a k-1 size particle.
                {
                    if (k==2)
                        addterm[k]=0.0;
                    else
                    {
                        addterm[k] = addterm[k] + (v[1]/(v[k]-v[k-1]))*K[1][k-1]*(N[1]-N1s[k-1])*N[k-1];//Growth of k due to condensation of monomers on k-1.//
                        t1 = t1 + K[1][k-1]*(N[1]-N1s[k-1])*N[k-1];//Loss of monomers that have condensed.
                    }
                }

                // Similarly, we can explain the following terms below.
                if( N[1]<N1s[k+1] )
                {
                    if(k==MAX)
                        subterm[k] = 0.;
                    else
                    {
                        addterm[k] = addterm[k] -(v[1]/(v[k+1]-v[k]))*K[1][k+1]*(N[1]-N1s[k+1])*N[k+1];
                        //Growth of k due to evaporation of monomers from k+1.
                        t2 = t2 + K[1][k+1]*(-N[1]+N1s[k+1])*N[k+1];
                        //Gain of monomers that have evaporated.
                    }
                }

                if(N[1]>N1s[k])
                {
                    if(k==MAX)
                        subterm[k] = 0.;
                    else
                        subterm[k] = subterm[k] + (v[1]/(v[k+1]-v[k]))*K[1][k]*(N[1]-N1s[k])*N[k];
                    //Loss of k due to condensation of monomers on k.
                }

                if(N[1]<N1s[k])
                {
                    if (k==2)
                    {
                        subterm[k] = subterm[k] - (v[1]/(v[k]-v[k-1]))*K[1][k]*(N[1]-N1s[k])*N[k];
                        t2 = t2 + K[1][k]*(-N[1]+N1s[k])*N[k];//Gain of monomers that have evaporated.
                    }
                    else
                        subterm[k] = subterm[k] - (v[1]/(v[k]-v[k-1]))*K[1][k]*(N[1]-N1s[k])*N[k];
                    //Loss of k due to evaporation of monomers from k.
                }
            }

        }

        for(k=2;k<=MAX;k++)
            N[k]= N[k] + step*(addterm[k] - subterm[k]);//changing PSD due to surface growth/evaporation only

        T=T-step*coolrate;
        t=t+step;
    }
    return(0);
}
int main()
{
  // Declaration of variables
  int choice,MAX,beta_option,restart_option;
  double A,B,MW,rho,T,coolrate;
  double Ksp,QM,viscosity,wetangle;
  FILE *fptr,*fptr1,*fptr2;
  if((fptr=fopen("main.inp","r+"))==NULL)//Opening the input file that contains property data for the aerosol material.
  {
    printf("Warning: Could not open file");
    return(-1);
  }
  //Scanning the property data from the input file "main.inp"
    fscanf(fptr,"Number of nodes: %d Starting Temperature(Kelvin): %lf Cooling Rate(Kelvin/s): %lf Concentration(M): %lf "
                   "Molecular Weight(Kg/mol): %lf Density(Kg/m3): %lf Viscosity(Pa*s): %lf "
                   "Surface Tension in the form A-BT(dyne/cm): A=%lf B=%lf Solubility Product Ksp: %lf "
                   "Enter your choice - 1)Coagulation,2)nucleation+coagulation,3)nucleation + coagulation + surface growth,4)surface growth:%d "
                   "Choice of collision frequency function:1)Free Molecular Regime,2)Fuchs form for free molecular and transition regime,3)Stokes-Einstein:%d "
                   "Heterogeneous Wetting Angle: %lf Restart: %d",
           &MAX,&T,&coolrate,&QM,&MW,&rho,&viscosity,&A,&B,&Ksp,&choice,&beta_option,&wetangle,&restart_option);
    fclose(fptr);//Closing the main property data file.

     //**********************************************
    //           PURE COAGULATION
    //**********************************************
    // This part of the code solves for coagulation problems. The code verifies that we obtain a self preserving distribution.
    // It verifies the problem stated in Friedlander, 2000 - page 218, line 4 - and shows that SPD is obtained.
    // And note that in this case node 1 is not the monomer. Its particle of size 1 nm.
    // We do not have to worry about monomers, since we are concerned with coagulation only.
  if (choice==1)
  {
    coagulation();
  }
  //************************************************************
  //             Nucleation + Coagulation
  //************************************************************
  
 else if (choice==2) //If user chooses 2, code runs for nucleation + coagulation.
    {
      coagulation_nucleation();
    }
  
 //*****************************************
 // SURFACE GROWTH + NUCLEATION + COAGULATION
 //*****************************************
  else if(choice==3) // THE COMPLETE GDE
    {
      coagulation_nucleation_surface();
    }

  //*****************************************
  //         PURE SURFACE GROWTH
  //*****************************************
  
  else if (choice ==4)
    { 
        surface_only();
    } 
  
}


/*if(beta_option==1)
        {
            for(i=1;i<=MAX;i++)
            {
                for(j=1;j<=MAX;j++)
                {

                    temp1=1/v[i] + 1/v[j];
                    temp2=pow(v[i],0.333333) + pow(v[j],0.333333);
                    K[i][j]= pow(3.0/(4.0*pi),0.1666666)*pow((6.0*kb*T/rho),0.5)*pow(temp1,0.5)*pow(temp2,2.0);
                    if (K[i][j]<Kmin)
                        Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time.

                }
            }
        }
        if(beta_option==2)
        {
            mu=A_mu*pow(T,1.5)/(B_mu+T);
            lambda = (mu/(P*101325.0))* sqrt(pi*R*T/(2.0*0.04));
            for(i=1;i<=MAX;i++)
            {
                for(j=1;j<=MAX;j++)
                {
                    kn1=(2.0*lambda)/dp[i];
                    kn2=(2.0*lambda)/dp[j];
                    D1=(kb*T)/(3.0*pi*mu*dp[i])*((5.0+4.0*kn1+6.0*kn1*kn1+18.0*kn1*kn1*kn1)/(5.0-kn1+(8.0+pi)*kn1*kn1));
                    D2=(kb*T)/(3.0*pi*mu*dp[j])*((5.0+4.0*kn2+6.0*kn2*kn2+18.0*kn2*kn2*kn2)/(5.0-kn2+(8.0+pi)*kn2*kn2));
                    c1=sqrt((8.0*kb*T)/(pi*m[i]));
                    c2=sqrt((8.0*kb*T)/(pi*m[j]));
                    l1=(8.0*D1)/(pi*c1);
                    l2=(8.0*D2)/(pi*c2);
                    g1=(pow((dp[i]+l1),3)-pow((dp[i]*dp[i]+l1*l1),1.5))/(3.0*dp[i]*l1)-dp[i];
                    g2=(pow((dp[j]+l2),3)-pow((dp[j]*dp[j]+l2*l2),1.5))/(3.0*dp[j]*l2)-dp[j];
                    K[i][j]=2.0*pi*(D1+D2)*(dp[i]+dp[j])/((dp[i]+dp[j])/(dp[i]+dp[j]+2.0*sqrt(g1*g1+g2*g2))+(8.0*(D1+D2))/(sqrt(c1*c1+c2*c2)*(dp[i]+dp[j])));
                    if (K[i][j]<Kmin)
                        Kmin = K[i][j];//Calculating the smallest collision frequency function to decide the characteristic coagulation time.
                }
            }
        }*/
