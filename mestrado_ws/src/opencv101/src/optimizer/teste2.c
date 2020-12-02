#include <stdio.h>
#include "stdlib.h"
#include "dos.h"

//i: integer;
double Jmin=0;
double  result[6];
double  u_opt[5];
double  x_00[4]={0.01,0.01,0.01,0};
double  u_00[4]={0.5,0,0.5,0};
double  sc_0[12]={0,0,0,0,0,0,0,0,0,0,0,0};
double  sc_1[8]={0,0,0,0,0,0,0,0};
double  ufc_0[12]={0,0,0,0,0,0,0,0,0,0,0,0};
double  ufc_1[8]={0,0,0,0,0,0,0,0};
double  Parameters_0[18]={3,3,1,1,0.05,1,0.5,0,0,0,0,0,0,0,0,0,0,0}; 
double  Parameters_1[16]={3,3,1,1,0.05,1,0.5,0,0,0,0,0,0,0,0,0}; 

int main()
{
	while(1) {
	   
	   
	   optimizer(Parameters_0,x_00,u_00,sc_0,ufc_0,result);
	   
	   Jmin=result[0]; 
	   u_opt[0]=result[1];
	   u_opt[1]=result[2];
	   u_opt[2]=result[3];
	   u_opt[3]=result[4];
	   u_opt[4]=result[5];
	   
	   optimizer(Parameters_1,x_00,u_00,sc_1,ufc_1,result);
	   
	   Jmin=result[0]; 
	   u_opt[0]=result[1];
	   u_opt[1]=result[2];
	   u_opt[2]=result[3];
	   u_opt[3]=result[4];
	   u_opt[4]=result[5];
	   
	   optimizer(Parameters_1,x_00,u_00,sc_1,ufc_1,result);
	   
	   Jmin=result[0]; 
	   u_opt[0]=result[1];
	   u_opt[1]=result[2];
	   u_opt[2]=result[3];
	   u_opt[3]=result[4];
	   u_opt[4]=result[5];


	   //delay(20);
	   //end;
	   //writeln('o valor de Jmin aqui é ', Jmin);
	   //writeln(' ');
	   //writeln('o valor de u1 aqui é ', u_opt[0]);
	   //writeln(' ');
	   //writeln('o valor de u2 aqui é ', u_opt[1]);	
	   //writeln(' ');
	   //writeln('o valor de u3 aqui é ', u_opt[2]);
	   //writeln(' ');
	   //writeln('o valor de u4 aqui é ', u_opt[3]);
   }
	return 0;
}
