/**
 * Document: MaxCompiler tutorial (maxcompiler-tutorial.pdf)
 * Chapter: 4      Example: 2      Name: Simple
 * MaxFile name: Simple
 * Summary:
 * 	 Takes a stream and for each value x calculates x^2 + x.
 */
#include <stdint.h>
#include <MaxSLiCInterface.h>
#include "Maxfiles.h"

int check(float *dataOut, float *expected, int size)
{
	int status = 0;
	for(int i=0; i < size; i++)
	{
		if(dataOut[i] != expected[i]) {
			fprintf(stderr, "Output data @ %d = %1.8g (expected %1.8g)\n",
				i, dataOut[i], expected[i]);
			status = 1;
		}else {
		   fprintf(stderr, "Test OK for Output data @ %d = %1.8g (expected %1.8g)\n",
				i, dataOut[i], expected[i]); 
		}
	}
	return status;
}

//dy/dx = 1 + y^2
static void derivs(float x, float *y, float *dydx, int n){
  int i = 0 ;
  for(i = 0 ; i < n ; i++){
    dydx[i] = 1.0 + y[i]*y[i];   
    //dydx[i] = 1.0 ;    

  }  
}

static void rk4_calc(const float y[], const float dydx[], int n, float x, float h, float yout[]){
  int i = 0 ;
  float xh, hh, h6, *dym, *dyt, *yt;
  
  dym = (float*)malloc(sizeof(float)*n);
  dyt = (float*)malloc(sizeof(float)*n);
  yt = (float*)malloc(sizeof(float)*n);
  
  hh = h*0.5;
  h6 = h/6.0;
  xh = x + hh;
  
  for(i = 0 ; i < n ; i++){
    yt[i] = y[i] + hh*dydx[i];    
  }
  derivs(xh, yt, dyt, n);

  for(i = 0 ; i < n ; i++){
    yt[i] = y[i] + hh*dyt[i];    
  }
  derivs(xh, yt, dym, n);
    
  for(i = 0 ; i < n ; i++){
    yt[i]  = y[i] + h*dym[i];
    dym[i] = dym[i] + dyt[i];
  }
  derivs(x+h, yt, dyt, n);

  for(i = 0 ; i < n ; i++){
    yout[i] = y[i] + h6*(dydx[i] + dyt[i] + 2.0*dym[i]);    
  }
}

const float y[8]          = { 1, 0, 2, 0, 4, 1, 8, 3 };
const float dydx[8]       = { 2, 1, 3, 1, 5, 2, 9, 4 };
float x = 1.0;
float h = 2.0;
float yout[8]       = { 0, 0, 0, 0, 0, 0, 0, 0 };

float expected[8]   = { 5, 3.5, 6.5, 3.5, 9.5, 5, 15.5, 8 };
const int n = 8;


int main()
{

	Simple(n, h, x, dydx, y, yout);

	printf("Running DFE RK4.\n");
	
	//check expected values
	rk4_calc(y, dydx, n, x, h, expected);
	int status = check(yout, expected, n);
	if (status)
		printf("Test failed.\n");
	else
		printf("Test passed OK!\n");
	return 0;
}
