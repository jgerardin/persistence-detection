/****************************************************************
 * simulate.cc
 * Jaline Gerardin 2010
 *
 * Modeling enzymatic circuits, 3 nodes, OR/AND, 1 pulse of input
 * 
 * simulate.cc is given a number of command line arguments, simulates
 * a circuit's response to 1 pulse of input, and outputs to stdout in
 * one of 2 modes: metrics (no RECORD mode) and full timecourse
 * (RECORD mode).
 *
 * Command line arguments:
 *
 * [input duration] [basal input] [change in input]
 *
 * [node A logic] [node B logic] [node C logic] where 0 is OR and 1 is
 * AND
 *
 * 26 kcat's and Kms:
 * [kcat of input acting on node A] [Km of input acting on node A]
 * [kcat of A acting on A] [Km of A acting on A] [kcat of A acting on
 * B] [Km of A acting on B] [kcat of A acting on C] [Km of A acting on
 * C] [kcat of constitutive regulator of A][Km of constitutive
 * regulator of A]
 * + similar for node B and node C
 *
 * Note that regulation type (activator, inhibitor, absent) is coded
 * in kcat: positive kcat = activator, negative kcat = inhibitor, and
 * 0 kcat = absent.
 *
 * Optional arguments:

 * Initial concentrations of nodes A, B, and C. If none are given,
 * simulate.cc will initialize the circuit from initial concentration
 * of 0.1 for each node, then either output the initialized values
 * (metrics mode) or continue onward to apply input (fill timecourse
 * mode). If the circuit cannot reach steady state during
 * initialization, an error of "-1" is outputed. If initial
 * concentrations are given, simulate.cc will apply input without
 * pre-initializing and simulate from there.
 *
 * ODEs are defined in the getplus() function and integrated using a
 * fifth-order Runge-Kutta method with adaptive stepsize control.
 *
 * In metrics mode, simulate.cc will calculate and output metrics
 * immediately after the input pulse as well as after the circuit has
 * reached steady state once input has turned off. Metrics include
 * initial output value, final output value, output amplitude, and
 * other metrics. Normally simulate.cc is set to output an error
 * message of "-1" if the output node decreases in value in response
 * to input.
 ****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <iostream>
using namespace std;

/*** CONSTANTS AND PARAMETERS ***/
// mode toggler: uncomment to output timecourse; comment to output timecourse metrics
//#define RECORD 

#define RECNODE 2 // output node
#define NUMNODES 3 // number of nodes
#define INIT 0.1 // initial concentrations
#define FCONC 0.1 // constitutive enzyme conc
#define TOTALNODE 1 // total concentration of each node

// parameters for determining if steady state has been reached
#define ZERO 0.000001 // below 1e-6 -> set to 0.
#define SSCHECKSTART 50
#define SSCHECKEVERY 200

// adaptive Runge-Kutta parameters (ODE integrator)
#define DELTAZ 0.000001
#define SAFETY 0.9
#define ERRCON 0.000189
#define HMIN 0.001
#define HMAX 0.1

// parameters for metric measurement
#define INTTIME 50000 // duration to integrate over for measuring integration of output

// global parameters set in command line
double base_input; // basal input level
int Aand, Band, Cand; // logic of each node
double Ainit = INIT, Binit = INIT, Cinit = INIT; // initial value of each node

// kcat's and Km's for 3-node circuit
double kIA, KIA;
double kAA, KAA, kBA, KBA, kCA, KCA, kFAA, KFAA;
double kAB, KAB, kBB, KBB, kCB, KCB, kFBB, KFBB;
double kAC, KAC, kBC, KBC, kCC, KCC, kFCC, KFCC;


/*** FUNCTION DECLARATIONS ***/
// initializations
void start(vector<double> &ks, vector<double> &Ks, vector<int> &andtrack);

// Runge-Kutta ODE solver
void step_rkck(double h1, vector<double> &innodes, vector<double> &outnodes, vector<double> &yerr,
	       vector<double> &Ks, vector<double> &ks, vector<int> &andtrack);
double step_rkas(double h1, double &h2, vector<double> &nodes, vector<double> &Ks,
		 vector<double> &ks, vector<int> &andtrack);
double mymax(vector<double> &x);

// ODEs defined
void getplus(vector<double> &plus, vector<double> &nodes, vector<double> &Ks,
	     vector<double> &ks, vector<int> &andtrack);

// calculate timecourse metrics
int get_metrics(vector<double> &timecourse, vector<double> &timevec, double &init, double &final, double &peakheight,
		double &peaktime, double &halfup, double &halfdown, double &sstime, double &integral, double &ten, 
		double &ninety,  double &inttopeak, double &inttohalf);


/*** MAIN ***
 *
 * 1. Read in parameters from command line
 * 2. Initialize nodes to pre-input steady state
 * 3. Apply input and simulate response
 * 4. Remove input and simulate response until steady state is reached
 *
 */
int main(int argc, char** argv)
{

  /*** READ IN AND SET SIMULATION PARAMETERS ***/

  double input_duration, input_change;
  int initialization_duration = 1000000; // maximum allowable time for initialization
  bool initialize = false;

  if(argc < 33) // if not enough arguments in command line
    {
      cout << "Syntax: ./FILENAME input_duration basal_input change_in_input [node logic x 3] [26 kcats and Kms]. You have " 
	   << argc << "\n";
      return 1;
    }

  int k = 1;
  input_duration = atof(argv[k++]);
  base_input = atof(argv[k++]);
  input_change = atof(argv[k++]);

  // node logic. 0 = OR, 1 = AND
  Aand = atoi(argv[k++]);
  Band = atoi(argv[k++]);
  Cand = atoi(argv[k++]);

  // kcat and Km of input acting on node A
  kIA = atof(argv[k++]);
  KIA = atof(argv[k++]);

  // regulations by node A
  kAA = atof(argv[k++]);
  KAA = atof(argv[k++]);
  kAB = atof(argv[k++]);
  KAB = atof(argv[k++]);
  kAC = atof(argv[k++]);
  KAC = atof(argv[k++]);
  // constitutive regulator of A
  kFAA = atof(argv[k++]);
  KFAA = atof(argv[k++]);

  // regulations by node B
  kBA = atof(argv[k++]);
  KBA = atof(argv[k++]);
  kBB = atof(argv[k++]);
  KBB = atof(argv[k++]);
  kBC = atof(argv[k++]);
  KBC = atof(argv[k++]);
  // constitutive regulator of B
  kFBB = atof(argv[k++]);
  KFBB = atof(argv[k++]);

  // regulations by node C
  kCA = atof(argv[k++]);
  KCA = atof(argv[k++]);
  kCB = atof(argv[k++]);
  KCB = atof(argv[k++]);
  kCC = atof(argv[k++]);
  KCC = atof(argv[k++]);
  // constitutive regulator of C
  kFCC = atof(argv[k++]);
  KFCC = atof(argv[k++]);

  // if initialized node concentrations are given, use those; otherwise, simulation will initialize
  if(argc == 36) 
    {
      Ainit = atof(argv[k++]);
      Binit = atof(argv[k++]);
      Cinit = atof(argv[k++]);
    }
  else
    initialize = true;

 // starting node concentrations, parameter values, node logic set
  vector<double> nodes(NUMNODES);
  nodes[0] = Ainit;
  nodes[1] = Binit; 
  nodes[2] = Cinit;
  vector<double> ks(NUMNODES*(NUMNODES+1));
  vector<double> Ks(NUMNODES*(NUMNODES+1));
  vector<int> andtrack(NUMNODES);
  start(ks, Ks, andtrack);

  // variables for tracking steady state
  vector<double> sstrack = nodes;
  bool foundss = false;

  // variables for managing ODE solver
  double h, hnext = 0.05; // adaptive stepsize starting size
  double time = 0;

  /*** IF NECESSARY, INITIALIZE NODE CONCENTRATIONS BY LETTING THE SYSTEM COME TO STEADY STATE ***/
  if(initialize)
    {
      for(int i = 0; i < initialization_duration; i++)
	{
	  // one step in ODE solver
	  h = hnext;
	  h = step_rkas(h, hnext, nodes, Ks, ks, andtrack);
	  if(h < 0)
	    break;
	  time += h;

	  // check whether steady state has been reached; if so, stop initializing
	  if(i%SSCHECKEVERY == 0 && time > SSCHECKSTART)
	    if(fabs(nodes[0] - sstrack[0]) < ZERO && fabs(nodes[1] - sstrack[1]) < ZERO && fabs(nodes[2] - sstrack[2]) < ZERO)
	      {
		foundss = true;
		break;
	      }
	    else
	      for(int j = 0; j < NUMNODES; j++)
		sstrack[j] = nodes[j];
	}

      // if timecourse is not being outputed, output initialized node concentrations or output "-1" for error message
#ifndef RECORD
      if(foundss == false)
	{
	  cout << -1 << "\t" << -1 << "\t" << -1 << endl;
	  return 1;
	}
      else
	{
	  cout << nodes[0] << " " << nodes[1] << " " << nodes[2] << endl;
	  return 0;	
	}
#endif
    }


  /*** APPLY INPUT AND SIMULATE RESPONSE ***/

  // variables for measuring timecourse metrics
  double initC, finalC, peakC, integral;
  double peaktime, halfup, halfdown, sstime;
  double ten, ninety;
  double inttopeak, inttohalf;
  vector<double> timecourse;
  vector<double> timevec;

  // input is changed
  base_input += input_change;

  // record starting output concentration and starting time
  time = 0;
  timecourse.push_back(nodes[RECNODE]);
  timevec.push_back(time);

  // ODE solver management
  h = 0.01;

  // simulate for duration input_duration
  while(time < input_duration)
    {
      // one step in ODE solver
      double myh = step_rkas(h, hnext, nodes, Ks, ks, andtrack);
      if(myh < 0) // solver error
	{
	  cout << -1 << "\t" << -1 << "\t" << -1 << endl;
	  return 1;
	}
      time += myh;
      h = hnext;

      // record output concentration and time
      timecourse.push_back(nodes[RECNODE]);
      timevec.push_back(time);

      // if timecourse is being outputed, do that
#ifdef RECORD
      cout << time << "\t" << nodes[0] << "\t" << nodes[1] << "\t" << nodes[2] << endl;
#endif

    }

  // calculate timecourse metrics
  int sign = get_metrics(timecourse, timevec, initC, finalC, peakC, peaktime, halfup, halfdown, sstime, integral, ten, ninety, inttopeak, inttohalf);
  // if metrics are being outputed, do that
#ifndef RECORD
  cout << initC << "\t" << finalC << "\t" << peakC << "\t" << peaktime << "\t" << halfup << "\t" << halfdown
       << "\t" << sstime << "\t" << integral << "\t" << ten << "\t" << ninety << "\t" << inttopeak << "\t" 
       << inttohalf << "\t";
#endif


  /*** REMOVE INPUT AND SIMULATE RESPONSE ***/

  // steady state tracking variables
  foundss = false;
  sstrack = nodes;

  // return input to initial level
  base_input -= input_change;

  // ODE solver
  h = 0.01;

  // simulate until steady state is reached or 86400 seconds, whichever comes first
  int i = 0;
  while(time < 86400)
    {
      i++;
      // one step in ODE solver
      double myh = step_rkas(h, hnext, nodes, Ks, ks, andtrack);
      if(myh < 0)
	{
	  cout << -1 << "\t" << -1 << "\t" << -1 << endl;
	  return 1;
	}
      time += myh;
      h = hnext;

      // record output concentration and time
      timecourse.push_back(nodes[RECNODE]);
      timevec.push_back(time);
      // if timecourse is being outputed, do that
#ifdef RECORD
      cout << time << "\t" << nodes[0] << "\t" << nodes[1] << "\t" << nodes[2] << endl;
#endif

      // check whether steady state has been reached
      if(i%SSCHECKEVERY == 0 && i > SSCHECKSTART)
	if(fabs(nodes[0] - sstrack[0]) < ZERO && fabs(nodes[1] - sstrack[1]) < ZERO && fabs(nodes[2] - sstrack[2]) < ZERO
	   && fabs(nodes[RECNODE] - timecourse[timecourse.size()-2]) < ZERO)
	  {
	    foundss = true;
	    break;
	  }
	else
	  for(int j = 0; j < NUMNODES; j++)
	    sstrack[j] = nodes[j];
    }

  // calculate timecourse metrics
  sign = get_metrics(timecourse, timevec, initC, finalC, peakC, peaktime, halfup, halfdown, sstime, integral, ten, ninety, inttopeak, inttohalf);

  // if metrics are being outputed, do that
#ifndef RECORD
  if(foundss == false || sign < 0) // if steady state was not reached, or if output turned off in response to input, output error
    cout << "-1\n";
  else
    cout << initC << "\t" << finalC << "\t" << peakC << "\t" << peaktime << "\t" << halfup << "\t" << halfdown
	 << "\t" << sstime << "\t" << integral << "\t" << ten << "\t" << ninety << "\t" << inttopeak << "\t" 
	 << inttohalf << endl;
#endif

  return 0;
}


/*** FUNCTION DEFINITIONS ***/
/* find maximum element in vector */
double mymax(vector<double> &x)
{
  double max = x[0];
  for(unsigned int i = 1; i < x.size(); i++)
    if(x[i] > max)
      max = x[i];
  return max;
}

/* adaptive stepsize control for Runge-Kutta ODE solver */
double step_rkas(double h1, double &h2, vector<double> &nodes, vector<double> &Ks,
		 vector<double> &ks, vector<int> &andtrack)
{
  double delta1, h = h1, htemp;
  vector<double> outnodes(NUMNODES);
  vector<double> yerr(NUMNODES);

  time_t rkstart = time(NULL);

  for(;;)
    {
      time_t rknow = time(NULL);
      if (rknow - rkstart > 100)
	return -1;
      step_rkck(h, nodes, outnodes, yerr, Ks, ks, andtrack);
      delta1 = mymax(yerr);
      if(h < HMIN || delta1 < DELTAZ)
	break;
      if(delta1 <= DELTAZ)
	break;
      htemp = h*SAFETY*pow(fabs(DELTAZ/delta1), 0.2);
      if(htemp < 0.1*h)
	h = 0.1*h;
      else
	h = htemp;
    }
  if(delta1 > ERRCON)
    h2 = h*SAFETY*pow(fabs(DELTAZ/delta1), 0.25);
  else
    h2 = 5*h;
  if(h2 < HMIN)
    h2 = HMIN;
  if(h2 > HMAX)
    h2 = HMAX;

  for(int i = 0; i < NUMNODES; i++)
    {
      nodes[i] = outnodes[i];
      if(nodes[i] < ZERO)
	nodes[i] = 0;
      else if(nodes[i] > TOTALNODE)
	nodes[i] = TOTALNODE;
    }

  return h;

}
/* Fifth-order Runge-Kutta step */
void step_rkck(double h1, vector<double> &innodes, vector<double> &outnodes, vector<double> &yerr,
	       vector<double> &Ks, vector<double> &ks, vector<int> &andtrack)
{
  vector<double> plus(NUMNODES);
  getplus(plus, innodes, Ks, ks, andtrack);
  vector<double> rkk1(NUMNODES);
  vector<double> rkk2(NUMNODES);
  vector<double> rkk3(NUMNODES);
  vector<double> rkk4(NUMNODES);
  vector<double> rkk5(NUMNODES);
  vector<double> rkk6(NUMNODES);
  vector<double> tempnodes(NUMNODES);

  for(int i = 0; i < NUMNODES; i++)
    {  
      rkk1[i] = h1*plus[i]; // k1
      tempnodes[i] = innodes[i] + (1./5)*rkk1[i];
    }
  getplus(plus, tempnodes, Ks, ks, andtrack);
  for(int i = 0; i < NUMNODES; i++)
    {
      rkk2[i] = h1*plus[i]; // k2
      tempnodes[i] = innodes[i] + (3./40)*rkk1[i] + (9./40)*rkk2[i];
    }
  getplus(plus, tempnodes, Ks, ks, andtrack);
  for(int i = 0; i < NUMNODES; i++)
    {
      rkk3[i] = h1*plus[i]; // k3
      tempnodes[i] = innodes[i] + (3./10)*rkk1[i] + (-9./10)*rkk2[i] + (6./5)*rkk3[i];
    }
  getplus(plus, tempnodes, Ks, ks, andtrack);
  for(int i = 0; i < NUMNODES; i++)
    {
      rkk4[i] = h1*plus[i]; // k4
      tempnodes[i] = innodes[i] + (-11./54)*rkk1[i] + (5./2)*rkk2[i] + (-70./27)*rkk3[i] + (35./27)*rkk4[i];
    }
  getplus(plus, tempnodes, Ks, ks, andtrack);
  for(int i = 0; i < NUMNODES; i++)
    {
      rkk5[i] = h1*plus[i]; // k5
      tempnodes[i] = innodes[i] + (1631./55296)*rkk1[i] + (175./512)*rkk2[i] + (575./13824)*rkk3[i] + (44275./110592)*rkk4[i] + (253./4096)*rkk5[i];
    }
  getplus(plus, tempnodes, Ks, ks, andtrack);
  for(int i = 0; i < NUMNODES; i++)
    {
      rkk6[i] = h1*plus[i]; // k6
      outnodes[i] = innodes[i] + (37./378)*rkk1[i] + (250./621)*rkk3[i] + (125./594)*rkk4[i] + (512./1771)*rkk6[i];
      yerr[i] = (37./378 - 2825./27648)*rkk1[i] + (250./621 - 18575./48384)*rkk3[i] + (125./594 - 13575./55296)*rkk4[i] 
	+ (-277./14336)*rkk5[i] + (512./1771 - 1./4)*rkk6[i];
    }
}

/* differential equations defined */
void getplus(vector<double> &plus, vector<double> &nodes, vector<double> &Ks,
	     vector<double> &ks, vector<int> &andtrack)
{
  for(int i = 0; i < NUMNODES; i++)
    plus[i] = 0;
  vector<double> conc(NUMNODES*(NUMNODES+1));

  for(int i = 0; i < NUMNODES; i++)
    {
      for(int j = 0; j < (NUMNODES+1); j++)
	if(ks[i*(NUMNODES+1)+j] > 0)
	  conc[i*(NUMNODES+1)+j] = TOTALNODE - nodes[i];
	else 	
	  conc[i*(NUMNODES+1)+j] = nodes[i];

      if(andtrack[i] == 0) // node is under OR logic
	{
	  for(int j = 0; j < (NUMNODES); j++)
	    {
	      if(fabs(ks[i*(NUMNODES+1)+j]) > 0)
		plus[i] += nodes[j]*ks[i*(NUMNODES+1)+j]*conc[i*(NUMNODES+1)+j]/(conc[i*(NUMNODES+1)+j] + Ks[i*(NUMNODES+1)+j] + nodes[j]);
	    }
	}
      else // node is under AND logic
	{
	  if(i != 0) // node is not input node
	    {
	      double in1 = 0, in2 = 0, in3 = 0;
	      int sign1 = 1, sign2 = 1, sign3 = 1;
	      if(Ks[i*(NUMNODES+1)+0] > 0) // if there's regulation by A on node i
		{
		  in1 = ks[i*(NUMNODES+1)+0]*nodes[0]/(conc[i*(NUMNODES+1)+0] + Ks[i*(NUMNODES+1)+0] + nodes[0]);
		  if(ks[i*(NUMNODES+1)+0] < 0)
                    sign1 = -1;
                }
	      if(ks[i*(NUMNODES+1)+0] * ks[i*(NUMNODES+1)+1] > 0) // if reg by A and by B have the same sign
		in1 *= ks[i*(NUMNODES+1)+1]*nodes[1]/(conc[i*(NUMNODES+1)+1] + Ks[i*(NUMNODES+1)+1] + nodes[1]);
	      else if(Ks[i*(NUMNODES+1)+1] > 0) // if reg by A and by B have different sign
		{
		  in2 = ks[i*(NUMNODES+1)+1]*nodes[1]/(conc[i*(NUMNODES+1)+1] + Ks[i*(NUMNODES+1)+1] + nodes[1]);
		  if(ks[i*(NUMNODES+1)+1] < 0)
                    sign2 = -1;
                }
	      if(ks[i*(NUMNODES+1)+0] * ks[i*(NUMNODES+1)+2] > 0) // if reg by A and by C have the same sign
		in1 *= ks[i*(NUMNODES+1)+2]*nodes[2]/(conc[i*(NUMNODES+1)+2] + Ks[i*(NUMNODES+1)+2] + nodes[2]);
	      else if(ks[i*(NUMNODES+1)+1] * ks[i*(NUMNODES+1)+2] > 0) // if reg by B and by C have the same sign
		in2 *= ks[i*(NUMNODES+1)+2]*nodes[2]/(conc[i*(NUMNODES+1)+2] + Ks[i*(NUMNODES+1)+2] + nodes[2]);
	      else if(Ks[i*(NUMNODES+1)+2] > 0) // if reg by C exists
		{
		  in3 = ks[i*(NUMNODES+1)+2]*nodes[2]/(conc[i*(NUMNODES+1)+2] + Ks[i*(NUMNODES+1)+2] + nodes[2]);
		  if(ks[i*(NUMNODES+1)+2] < 0)
                    sign3 = -1;
                }

              if(in1 > 0 && sign1 < 0)
		in1 *= sign1;
              if(in2 > 0 && sign2 < 0)
                in2 *= sign2;
              if(in3 > 0 && sign3 < 0)
                in3 *= sign3;

	      plus[i] += in1*conc[i*(NUMNODES+1)+0] + in2*conc[i*(NUMNODES+1)+1] + in3*conc[i*(NUMNODES+1)+2];
	    }
	  else // node is input node
	    {
	      double in0 = base_input*kIA/((TOTALNODE-nodes[0]) + KIA + base_input);
	      double in1 = 0, in2 = 0, in3 = 0;
	      int sign1 = 1, sign2 = 1, sign3 = 1;
	      if(ks[0*(NUMNODES+1)+0] > 0) // if reg by A has the same sign as input (positive)
		in0 *= ks[0*(NUMNODES+1)+0]*nodes[0]/(conc[0*(NUMNODES+1)+0] + Ks[0*(NUMNODES+1)+0] + nodes[0]);
	      else if(Ks[0*(NUMNODES+1)+0] > 0) // if there's regulation by A on node i
		{
		  in1 = ks[0*(NUMNODES+1)+0]*nodes[0]/(conc[0*(NUMNODES+1)+0] + Ks[0*(NUMNODES+1)+0] + nodes[0]);
		  if(ks[0*(NUMNODES+1)+0] < 0)
                    sign1 = -1;
                }
	      if(ks[0*(NUMNODES+1)+1] > 0) // if reg by B has the same sign as input
		in0 *= ks[0*(NUMNODES+1)+1]*nodes[1]/(conc[0*(NUMNODES+1)+1] + Ks[0*(NUMNODES+1)+1] + nodes[1]);
	      else if(ks[0*(NUMNODES+1)+0] * ks[0*(NUMNODES+1)+1] > 0) // if reg by A and by B have the same sign
		in1 *= ks[0*(NUMNODES+1)+1]*nodes[1]/(conc[0*(NUMNODES+1)+1] + Ks[0*(NUMNODES+1)+1] + nodes[1]);
	      else if(Ks[0*(NUMNODES+1)+1] > 0) // if reg by B exists
		{
		  in2 = ks[0*(NUMNODES+1)+1]*nodes[1]/(conc[0*(NUMNODES+1)+1] + Ks[0*(NUMNODES+1)+1] + nodes[1]);
		  if(ks[0*(NUMNODES+1)+1] < 0)
                    sign2 = -1;
                }
	      if(ks[0*(NUMNODES+1)+2] > 0) // if reg by C has the same sign as input
		in0 *= ks[0*(NUMNODES+1)+2]*nodes[2]/(conc[0*(NUMNODES+1)+2] + Ks[0*(NUMNODES+1)+2] + nodes[2]);
	      else if(ks[0*(NUMNODES+1)+0] * ks[0*(NUMNODES+1)+2] > 0) // if reg by A and by C have the same sign
		in1 *= ks[0*(NUMNODES+1)+2]*nodes[2]/(conc[0*(NUMNODES+1)+2] + Ks[0*(NUMNODES+1)+2] + nodes[2]);
	      else if(ks[0*(NUMNODES+1)+1] * ks[0*(NUMNODES+1)+2] > 0) // if reg by B and by C have the same sign
		in2 *= ks[0*(NUMNODES+1)+2]*nodes[2]/(conc[0*(NUMNODES+1)+2] + Ks[0*(NUMNODES+1)+2] + nodes[2]);
	      else if(Ks[0*(NUMNODES+1)+2] > 0) // if reg by C exists
		{
		  in3 = ks[0*(NUMNODES+1)+2]*nodes[2]/(conc[0*(NUMNODES+1)+2] + Ks[0*(NUMNODES+1)+2] + nodes[2]);
		  if(ks[0*(NUMNODES+1)+2] < 0)
                    sign3 = -1;
                }

              if(in1 > 0 && sign1 < 0)
		in1 *= sign1;
              if(in2 > 0 && sign2 < 0)
                in2 *= sign2;
              if(in3 > 0 && sign3 < 0)
                in3 *= sign3;

	      plus[0] += in0*(TOTALNODE - nodes[0]) + in1*conc[0*(NUMNODES+1)+0] + in2*conc[0*(NUMNODES+1)+1] + in3*conc[0*(NUMNODES+1)+2];
	    }
	}
      // effect of constitutive activator/deactivator
      if(ks[i*(NUMNODES+1)+3] != 0)
	plus[i] += FCONC*ks[i*(NUMNODES+1)+3]*conc[i*(NUMNODES+1)+3]/(conc[i*(NUMNODES+1)+3] + Ks[i*(NUMNODES+1)+3] + FCONC);
    }

  if(andtrack[0] == 0) // add effect of input if input node is using OR logic
    plus[0] += base_input*kIA*(TOTALNODE-nodes[0])/((TOTALNODE-nodes[0]) + KIA + base_input);
}

/* calculate timecourse metrics */
int get_metrics(vector<double> &timecourse, vector<double> &timevec, double &init, 
		double &final, double &peakheight, double &peaktime, double &halfup, 
		double &halfdown, double &sstime, double &integral, double &ten, 
		double &ninety, double &inttopeak, double &inttohalf)
{
  // initialize quantities to track
  double halfheight;
  bool pastpeak = false, pasthalf = false;
  init = timecourse[0]; // output initial value
  final = timecourse.back(); // output final value
  peakheight = 0; // maximum output amplitude
  peaktime = 0; // time maximum output amplitude was reached
  halfup = 0; // time 50% max output amplitude was reached, starting from low amplitude
  halfdown = sstime; // time 50% max output amplitude was reached, starting from high amplitude
  sstime = timecourse.size() - 1; // time output reached steady state
  integral = 0; // integral of output
  ten = 0; // time 10% max output amplitude was reached, starting from low amplitude
  ninety = 0; // time 90% max output amplitude was reached, starting from high amplitude
  inttopeak = 0; // integral to 50% max amplitude
  inttohalf = 0; // integral to maximum amplitude

  int sign = 1; // = (1, -1) if output turns (on, off) in response to input
  int ipeaktime = 1;
  vector<double> temp(timecourse.size());

  // make a copy of the timecourse, shifting initial value to zero and taking absolute value.
  for(unsigned int i = 0; i < timecourse.size(); i++)
    {
      temp[i] = fabs(timecourse[i] - init);
      if(timevec[i] < INTTIME)
	integral += temp[i];
    }

  for(unsigned int i = 0; i < temp.size(); i++)
    if(temp[i] > peakheight)
      {
	peakheight = temp[i];
	ipeaktime = i;
	peaktime = timevec[i];
      }
  halfheight = peakheight/2;
  for(unsigned int i = 0; i < temp.size(); i++)
    {
      if(ninety == 0 and temp[i] >= 0.9*peakheight)
	ninety = timevec[i];
      if(ten == 0 and temp[i] >= 0.1*peakheight)
	ten = timevec[i];
      if(halfup == 0 and temp[i] >= halfheight)
	{
	  halfup = timevec[i];
	  pasthalf = true;
	}
      if(pastpeak == false && temp[i] >= peakheight)
	{
	  pastpeak = true;
	  pasthalf = false;
	}
      if(pastpeak == true && pasthalf == false && temp[i] <= halfheight)
	{
	  halfdown = timevec[i];
	  break;
	}
    }
  for(unsigned int i = 0; i < temp.size(); i++)
    {
      if(temp[i] < halfheight)
	inttohalf += temp[i];
      if(temp[i] < peakheight)
	inttopeak += temp[i];
      else
	break;
    }

  if(timecourse[ipeaktime] < init)
    sign = -1;
  return sign;
}

/* initialize parameter vectors */
void start(vector<double> &ks, vector<double> &Ks, vector<int> &andtrack)
{
  // initialize kcat's
  ks[0*(NUMNODES+1)+0] = kAA;
  ks[0*(NUMNODES+1)+1] = kBA;
  ks[0*(NUMNODES+1)+2] = kCA;
  ks[0*(NUMNODES+1)+3] = kFAA;
  ks[1*(NUMNODES+1)+0] = kAB;
  ks[1*(NUMNODES+1)+1] = kBB;
  ks[1*(NUMNODES+1)+2] = kCB;
  ks[1*(NUMNODES+1)+3] = kFBB;
  ks[2*(NUMNODES+1)+0] = kAC;
  ks[2*(NUMNODES+1)+1] = kBC;
  ks[2*(NUMNODES+1)+2] = kCC;
  ks[2*(NUMNODES+1)+3] = kFCC;

  // initialize Km's
  Ks[0*(NUMNODES+1)+0] = KAA;
  Ks[0*(NUMNODES+1)+1] = KBA;
  Ks[0*(NUMNODES+1)+2] = KCA;
  Ks[0*(NUMNODES+1)+3] = KFAA;
  Ks[1*(NUMNODES+1)+0] = KAB;
  Ks[1*(NUMNODES+1)+1] = KBB;
  Ks[1*(NUMNODES+1)+2] = KCB;
  Ks[1*(NUMNODES+1)+3] = KFBB;
  Ks[2*(NUMNODES+1)+0] = KAC;
  Ks[2*(NUMNODES+1)+1] = KBC;
  Ks[2*(NUMNODES+1)+2] = KCC;
  Ks[2*(NUMNODES+1)+3] = KFCC;

  // initialize logic coder
  andtrack[0] = Aand;
  andtrack[1] = Band;
  andtrack[2] = Cand;
}

/* end */
