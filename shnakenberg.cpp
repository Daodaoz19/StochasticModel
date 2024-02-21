//**********************************************
//the stochastic simulation for shnakenberg model using multigrid method
//
// Here we consider a special 1D RD system with only two variables. 
// In the future, we should make it more general. 
//
// The file was originally created by Fei Li, and then modified by Young Cao
// Spatial pattern forms when the diffusion rates D2 > 40 D1 in theory.
// But in real simulation, the pattern formation depends on the scale factor, the ratio, as well as random sequences. 
//
//
//**********************************************
//Usage: command [XBINNUM] [YBINNUM] [NUMBER OF RUNS]
//
//**********************************************

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

#define DefaultXBin 100 //default number of bins in x
#define DefaultYBin 100 // default number of bins in y
#define DefaultNumofRun 1 // we can choose to run it multiple times
// #define printtraj //Define this to allow trajectory output
// #define DEBUG // For debug purpose


int main(int argc, char* argv[]){

  int XBINNUM, YBINNUM, NUMofRUNS;
  if(argc>1) XBINNUM=atoi(argv[1]);
  else XBINNUM = DefaultXBin;
   
  if(argc>2) YBINNUM=atoi(argv[2]);
  else YBINNUM = DefaultYBin;

  if (argc > 3) NUMofRUNS = atoi(argv[3]); 
  else NUMofRUNS = DefaultNumofRun; 
    
  int pm = XBINNUM/YBINNUM; // To make the test simple, we require d2 > d1 and the ratio of the discretization sizes be integer

  long *x = new long[XBINNUM ]; // PDE people tends to use U and V as their variables, because X and Y are spatial indices. 
  long *y = new long[YBINNUM ];

    long Xtot, Ytot;
    int i;


  //*******************************************
  //model & parameter
  //
  // Shnakenberg model:
  // A -> X, with k1
  // X -> C, with k2
//   B -> Y, with k3
//  2X + Y -> 3X, with k4
  //*******************************************
  //Parameters
  //*******************************************

  ////////////////////////////////////////////////////
  //Parameter set
  ////////////////////////////////////////////////////

    double L = 10; // 1D interval length

    double w = 100; // 1000; //scale factor, interpreted as total population scale. When scale increases, the behavios is close to PDE. 
    double D1 = 5e-4, D2 = D1*10000; // Diffusion rates, D1 for X and D2 for Y
//    double D1 = 0, D2 = 0; // test no diffusion case
    double k1 = 1, k2 = 2, k3 = 3; // reaction constants
    
    double h1 = L/(XBINNUM), h2 = L/(YBINNUM); // interval length for the discretization
    double k4 = 1/(w*w); // reaction rate converted to propability rate and scaled with the scale factor
    double d1 = D1/(h1*h1), d2 = D2/(h2*h2); // jump rate in the SSA 
    const int RULENUM = 8; // Four types of chemical reactions + four diffusion events, all called rules

  //*****************************************
  //variables for SSA process
  //*****************************************
  double a[RULENUM];
  double a0;

  double r1, r2;
  double tau;
  int ita;
  double sum, sum_a, sum_x;
 
srand(time(NULL));
//    srand(0); 
  
  //*************************************************
  //get the population at each bin for PopZp at
  //the end of simulaiton
  //*************************************************  
  ofstream population("population", ios::out);
  
  //*************************************************
  //get the total population of M and P at the 
  //end of simulation
  //*************************************************
  ofstream total("total", ios::out);

  //*********************************
  //Get the firings of each rule
  //*********************************
  long firings[RULENUM];
  long location[XBINNUM]; 

  ofstream firefile("firings", ios::out);
  ofstream locationfile("location", ios::out); 

  //************************************************
  //Begin simulation, and set up timer
  //************************************************
  cout<<"Begin the model simulation ..."<<endl;
  cout<<"XBINNUM = "<< XBINNUM<<", YBINNUM = "<< YBINNUM<<endl;
  cout<<"Get the spatial population result in file: population"<<endl;
  cout<<"Get the total population result in file: total"<<endl;
  cout<<"total iterations "<< NUMofRUNS <<endl;



  for(int real=0; real< NUMofRUNS; real++){

    //*****************************************
    //*****************************************
    //Keep the trajectory of each run
    //*****************************************
#ifdef printtraj 
      char seq[8], trajfile[16];
    sprintf(seq, "%d", real);
    strcpy(trajfile, "poptraj");
    strcat(trajfile, seq);
    ofstream poptraj(trajfile, ios::out);
#endif
     // ****************************************/

    //*****************************************
    //set initial condition;
    //*****************************************
      for (i=0; i<XBINNUM; i++)
          x[i]= 2*w / XBINNUM;
      for (i=0; i<YBINNUM; i++)
          y[i]= 1*w / YBINNUM;
      
      Xtot= 2*w; Ytot= w;

      for( i=0; i<RULENUM; i++)
          firings[i]=0;
	for (i=0; i< XBINNUM; i++)
	  location[i] = 0; 
    
      double timeTracker=0.0, SimulationTime = 10;

    while(timeTracker<SimulationTime){

        //*****************************************
        // Calculate propensities. The indices are arranged so that diffusion rules are before reaction rules. 
	// This order helps in the search process. 
        //*****************************************

        a[0] = d1* (Xtot - x[0]); // X jumps to left
        a[1] = d1* (Xtot - x[XBINNUM - 1]); // X jumps to right
        a[2] = d2* (Ytot - y[0]); // Y jumps to left
        a[3] = d2* (Ytot - y[YBINNUM - 1]); // Y jumps to right
        a[4] = k1*w;  // k1 * A this is the total propensity for Rule 1
        a[5] = k2*Xtot; // total propensity for Rule 2
        a[6] = k3*w; // k3 * B total propensity for Rule 3
        
        a[7] = 0; // The multi-molecule reaction propensity calculation is tricky. 
        for (i= 0; i < XBINNUM; i++) {
            a[7] += x[i]*(x[i]-1)*y[i/pm];
        }
	a[7] *= k4*XBINNUM*YBINNUM;         

        a0=0;
        for(i=0; i<RULENUM; i++)
            a0 += a[i];
        
        do{r1=1.0*rand()/RAND_MAX;} while(r1<=0 || r1>=1);
        tau = -1.0/a0*log(r1);
        
        timeTracker += tau;

      //***********************************

        r2 = 1.0*rand()/RAND_MAX;
        sum = r2*a0;

        sum_a = a[0];
        for (ita = 0; sum_a < sum; ) {
            ita++; 
            sum_a += a[ita];
        }
        sum_a -= a[ita]; 
        
#ifdef DEBUG
	cout << "Y total = " << Ytot << "\t X total = " << Xtot << endl;
        cout << "rand num = " << r2 << endl; 
        for (i = 0; i < RULENUM; i ++)
            cout << a[i] << '\t';
        cout << endl;
        cout << "ita = " << ita << endl; 
	cout << "Rule " << ita << " fires! " ; 
#endif

        firings[ita]++;

        switch(ita){
            case 0: // X diffuse to left
                r2 = (sum - sum_a)/a[ita];
                sum = (Xtot - x[0])*r2;
                sum_x = x[1];
                for ( i = 1; sum_x < sum; ){  // Find the index for which bin an X jumps to left
                    i++;
                    sum_x += x[i];
                }
                x[i] --;
                x[i-1]++;
                break;
            case 1: // X diffuse to right
                r2 = (sum - sum_a)/a[ita];
                sum = (Xtot - x[XBINNUM -1])*r2;
                sum_x = x[0];
                for ( i = 0; sum_x < sum; ){  // Find the index for which bin an X jumps to right
                    i++;
                    sum_x += x[i];
                }
                x[i] --;
                x[i+1]++;
                break;
            case 2: // Y diffuse to left
                r2 = (sum - sum_a)/a[ita];
                sum = (Ytot - y[0])*r2;
                sum_x = y[1];
                for ( i = 1; sum_x < sum; ){  // Find the index for which bin an Y jumps to left
                    i++;
                    sum_x += y[i];
                }
                y[i] --;
                y[i-1]++;
                break;
            case 3: // Y diffuse to right
                r2 = (sum - sum_a)/a[ita];
                sum = (Ytot - y[YBINNUM-1])*r2;
                sum_x = y[0];
                for ( i = 0; sum_x < sum; ){  // Find the index for which bin an Y jumps to right
                    i++;
                    sum_x += y[i];
                }
                y[i] --;
                y[i+1]++;
                break;
            case 4: // reaction 1: A -> X, with k1
                r2 = (sum - sum_a)/a[ita];
                i = r2*XBINNUM;  // Find the index of the bin for the new X to appear at
                x[i] ++;
                Xtot ++;
                break;
            case 5:  // reaction 2: X -> C, with k2
                r2 = (sum - sum_a)/a[ita];
                sum_x = x[0];
                sum = r2 * Xtot;
                for ( i =0; sum_x < sum;){  // Find the index for which bin an X degradates
		    i++; 
                    sum_x += x[i];
                }
                x[i] --;
                Xtot --;
                break;
            case 6: //   B -> Y, with k3
                r2 = (sum - sum_a)/a[ita];
                i = r2*YBINNUM;  // Find the index of the bin for the new Y to appear at
                y[i] ++;
                Ytot ++;
#ifdef DEBUG
	cout << "r2 " << r2 << " round = " << r2*YBINNUM << " location " << i << ": x [i] =" << x[i] << " y[i] = " << y[i] << endl; 
#endif
                break;
            case 7:    //  2X + Y -> 3X, with k4; this is the most complicated one, be careful
                
                //r2 = (sum - sum_a)/a[ita];
                sum = (sum - sum_a)/(k4*XBINNUM*YBINNUM); // == a[3]*r2/k4
                sum_x = x[0]*(x[0]-1)*y[0];
                for ( i = 0; sum_x < sum; ){  // Find the index for which bin an X degradates
                    i++;
                    sum_x += x[i]*(x[i]-1)*y[i/pm];
                }
                x[i]++; Xtot++;
                y[i/pm]--; Ytot--;
                break; 

            default:
                cerr<<"In finding the rule number ..."<<endl;
                cerr<<"bin_ita = "<<ita<<" not found"<<endl;
                exit(1);
                break;
      }
    	location[i]++; 
        
        #ifdef DEBUG
		cout << " at location " << i << endl; 
		cout << " x population " << x[i] << " with total = " << Xtot << endl; 
		cout << " corresponding Y " << y[i/pm] << " with total = " << Ytot << endl; 
		// cout << "Time proceed to " << timeTracker << "\t tau = " << tau << endl;
		// cout << "rule " << ita << endl; 
	#endif
	#ifdef printtraj
	        poptraj << timeTracker <<"\t";
        	for(int i=0; i < XBINNUM; i++) poptraj << x[i] <<"\t";
        	for(int i=0; i < YBINNUM; i++) poptraj << y[i] << "\t";
        	poptraj << endl;
	#endif

    }
    
    
    //******************************************************
    //end of simulation
    //
    //collect the statistics
    //******************************************************
    total << Xtot << "\t" << Ytot << endl;
    for(int i=0; i<XBINNUM; i++) population << x[i]<<"\t";
    for(int i=0; i<YBINNUM; i++) population << y[i]<<"\t";

    population << endl;
    for(int i=0; i<RULENUM; i++) firefile<<firings[i]<<"\t";
    for(int i=0; i<XBINNUM; i++) locationfile<<location[i]<<'\t'; // note: location info is only meaningful when u talk about equal discretizations
    firefile<<endl;  
    locationfile << endl;

    }

    cout<<"end of simulation ..."<<endl;
}

