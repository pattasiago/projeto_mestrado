program teste;

{$mode objfpc}{$H+}


uses   unitoptimizer, SysUtils;


var
  //i: integer;
  Jmin: Double=0;
  result:array[0..5] of double;
  u_opt:array[0..4] of double;
  x_0: array[0..3] of double=(0.01,0.01,0.01,0);
  u_0: array[0..3] of double=(0.5,0,0.5,0);
  sc_0: array[0..11] of double=(0,0,0,0,0,0,0,0,0,0,0,0);
  sc_1: array[0..7] of double=(0,0,0,0,0,0,0,0);
  ufc_0: array[0..11] of double=(0,0,0,0,0,0,0,0,0,0,0,0);
  ufc_1: array[0..7] of double=(0,0,0,0,0,0,0,0);
  Parameters_0: array[0..17] of double=(3,3,1,1,0.05,1,0.5,0,0,0,0,0,0,0,0,0,0,0); 
  Parameters_1: array[0..15] of double=(3,3,1,1,0.05,1,0.5,0,0,0,0,0,0,0,0,0); 

BEGIN


		

 	while(true) do begin
	   
	   
	   optimizer(Parameters_0,x_0,u_0,sc_0,ufc_0,result);
	   
	   Jmin:=result[0]; 
	   u_opt[0]:=result[1];
	   u_opt[1]:=result[2];
	   u_opt[2]:=result[3];
	   u_opt[3]:=result[4];
	   u_opt[4]:=result[5];
	   
	   optimizer(Parameters_1,x_0,u_0,sc_1,ufc_1,result);
	   
	   Jmin:=result[0]; 
	   u_opt[0]:=result[1];
	   u_opt[1]:=result[2];
	   u_opt[2]:=result[3];
	   u_opt[3]:=result[4];
	   u_opt[4]:=result[5];
	   
	   optimizer(Parameters_1,x_0,u_0,sc_1,ufc_1,result);
	   
	   Jmin:=result[0]; 
	   u_opt[0]:=result[1];
	   u_opt[1]:=result[2];
	   u_opt[2]:=result[3];
	   u_opt[3]:=result[4];
	   u_opt[4]:=result[5];


	   //sleep(20);
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
	   end;
	   	
END.
