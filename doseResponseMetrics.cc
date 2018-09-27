/*****************************************************
 * doseResponseMetrics.cc
 * Jaline Gerardin 2010
 *
 * given a dose response, calculate steepness and EC50
 *
 * command line arguments:
 * [number of data points] [responses] [doses]
 *****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
using namespace std;

int main(int argc, char**argv)
{
  // number of data points in dose response
  int numpoints = atoi(argv[1]);

  // read in response values
  vector<double> ydata(numpoints);
  for(int i = 0; i < numpoints; i++)
    ydata[i] = atof(argv[i+2]);

  // find maximum and minimum response values
  double max = ydata[numpoints-1];
  double min = ydata[0];
  for(int i = 0; i < numpoints; i++)
    {
      if(ydata[i] > max)
	max = ydata[i];
      if(ydata[i] < min)
	min = ydata[i];
    }

  // calculate values for 90%, 50%, and 10% response
  double ninety = (max - min)*0.9 + min;
  double fifty = (max - min)*0.5 + min;
  double ten = (max - min)*0.1 + min;

  // read in dose values
  vector<double> xdata(numpoints);
  for(int i = 0; i < numpoints; i++)
    xdata[i] = atof(argv[i+2+numpoints]);

  // identify indexes bracketing 10%, 50%, and 90% response
  int tenm = 0, tenp = 1;
  int ninem = numpoints - 2, ninep = numpoints - 1;
  int fifm = 0, fifp = numpoints-1;
  for(int i = 0; i < numpoints; i++)
    {
      if(ydata[i] > ten)
	{
	  tenm = i-1;
	  tenp = i;
	  break;
	}
    }
  for(int i = 0; i < numpoints; i++)
    {
      if(ydata[i] > fifty)
	{
	  fifm = i-1;
	  fifp = i;
	  break;
	}
    }
  for(int i = numpoints-1; i >= 0; i--)
    {
      if(ydata[i] < ninety)
	{
	  ninem = i;
	  ninep = i+1;
	  break;
	}
    }

  // use linear interpolation between bracketing indexes to estimate
  // doses for 10%, 50%, and 90%
  ten -= ydata[tenm];
  ten /= (ydata[tenp] - ydata[tenm]);
  fifty -= ydata[fifm];
  fifty /= (ydata[fifp] - ydata[fifm]);
  ninety -= ydata[ninem];
  ninety /= (ydata[ninep] - ydata[ninem]);

  double tenint = (xdata[tenp] - xdata[tenm])*ten + xdata[tenm];
  double nineint = (xdata[ninep] - xdata[ninem])*ninety + xdata[ninem];
  double fifint = (xdata[fifp] - xdata[fifm])*fifty + xdata[fifm];

  // write to stdout
  cout << max << "\t" << tenint/nineint << "\t"
       << fifint << endl;

  return 0;
}
