/* **************************************************************************** */
/*                                 user functions                               */
/* **************************************************************************** */
#include "o8para.h"
#define  X extern
#include "o8comm.h"
#include "o8fint.h"
#include "o8cons.h"
#include "o8fuco.h"
#undef   X
#include "matrix.h"
#include "stdio.h"
#include "stdlib.h"




int Hp, Hc; 
double x_0[100];
double ufc[100];
double sc[100];
double Fe=0;
int f1=0;
double Ts;
double Curv;
double CBp=0;
double f=0;
double f0=0;
//double xe0=0;
//double ye0=0;
//double alphae=0;
double qsip=0;
//double qsi=0;
double sp=0;
double s=0;
//double thetae0=0;
//double alphap=0;
double ur=0;
double uf=0;
double uf0=0;
double w_ref=0;
double H=0;
double demax=0;
double thetaemax=0;
double etamax=0;
int VISUAL = 0;
double alpha=0;
double p11=0;
double VPFMode = 0;

//double u1,u2,u3,u4,u5=0;
//double xe,ye,alphae,thetae,eta=0;
//double xe0,ye0,alphae0,thetae0,eta0=0;
//double xe1,ye1,alphae1,thetae1,eta1=0;

int id=0;
double k1;
double k2;
double k3;

double p12;
double p21;
double p13;
double p31;
double q;

double TC1=0;
double TC2=0;
double wu=0;
double ws=0;


int i = 0;
int ii = 0;

MAT   *Q;
MAT   *R;
MAT   *P;

MAT   *xb;
MAT   *xbcon;
MAT   *ub;
MAT   *ubcon;


void optimizer(double Params[], double x00[], double u00[], double sc0[], double ufc0[], double result[]){  
	
if (Params[20]==1){
VISUAL = 1;
}
else if (Params[20]==0){
VISUAL = 0;	
}	
	
	
void donlp2(void);
void solchk(void);
void error_model();
void calccoterms(int k); 


//Definição das constantes 
////Horizonte de predição
Hp=Params[0];
//Horizonte de controle
Hc=Params[1]; 
//Período de amostragem
Ts=Params[4]; 
//Curvatura no trecho de caminho
Curv=Params[5];
//	printf("\n");
//	printf("Curv %e",Curv);
//	printf("\n");
//Taxa de progressão do veículo virtual 
qsip=Params[6];
//Taxa de progressão do veículo virtual 
id=Params[7]; 
//Velocidade de navegação de referência 
ur=Params[8]; 
uf=Params[9];
w_ref=Params[10];

H=Params[19];



alpha=Params[22];

VPFMode=Params[23];

uf0=uf;
//	printf("\n");
//	printf("id %d",id);
//	printf("\n");

if (id==0){
	k1=Params[11];
	k2=Params[12];
	k3=Params[13];
	q=Params[14];
	p12=Params[15];
	p13=Params[16];
	ws=Params[17];
	wu=Params[18];
	//printf("\n");
	//printf("k1 = %e",k1);
	//printf("\n");
	//printf("k2 = %e",k2);
	//printf("\n");
	//printf("k3 = %e",k3);
	//printf("\n");
	//printf("q = %e",q);
	//printf("\n");
	//printf("p12 = %e",p12);
	//printf("\n");
	//printf("p13 = %e",p13);
	//printf("\n");
	//printf("ws = %e",ws);
	//printf("\n");
	//printf("wu = %e",wu);
	//printf("\n");

	for(i=0;i<3*(Hp+1);i++){
	ufc[i]=ufc0[i];
	sc[i]=sc0[i];
	//printf("\n");
	//printf("ufc %d = %e",id,ufc[i]);
	//printf("\n");		
	}
}
else if (id==1){
	
	k1=Params[11];
	k2=Params[12];
	q=Params[13];
	p21=Params[14];
	ws=Params[15];
	wu=Params[16];
	//printf("\n");
	//printf("k1 = %e",k1);
	//printf("\n");
	//printf("k2 = %e",k2);
	//printf("\n");
	//printf("q = %e",q);
	//printf("\n");
	//printf("p21 = %e",p21);
	//printf("\n");
	//printf("ws = %e",ws);
	//printf("\n");
	//printf("wu = %e",wu);
	//printf("\n");

	for(i=0;i<2*(Hp+1);i++){
	sc[i]=sc0[i];
	ufc[i]=ufc0[i];
	//printf("\n");
	//printf("ufc %d = %e",id,ufc[i]);
	//printf("\n");			
	}
}
else if (id==2){
	k1=Params[11];
	k3=Params[12];
	q=Params[13];
	p31=Params[14];
	ws=Params[15];
	wu=Params[16];
	//printf("\n");
	//printf("k1 = %e",k1);
	//printf("\n");
	//printf("k3 = %e",k3);
	//printf("\n");
	//printf("q = %e",q);
	//printf("\n");
	//printf("p31 = %e",p31);
	//printf("\n");
	//printf("ws = %e",ws);
	//printf("\n");
	//printf("wu = %e",wu);
	//printf("\n");
	for(i=0;i<2*(Hp+1);i++){
	ufc[i]=ufc0[i];
	sc[i]=sc0[i];
	//printf("\n");
	//printf("ufc %d = %e",id,ufc[i]);
	//printf("\n");		
	}
}
 

if (VISUAL==1){
x_0[0] = 0;
for(i=1;i<=2*(Hc+1);i++){
	x_0[i] = u00[i-1];
}
i=0;


//Inicialização das matrizes do modelo

//Matriz de ponderação do custo terminal
P=m_get(3,3);
m_ident(P);
sm_mlt(Params[21],P,P);
p11=P->me[0][0];
//m_output(P);

//Matriz de ponderação dos estados
Q=m_get(3,3);
m_ident(Q);
//sm_mlt(Params[2],Q,Q);
//Q->me[0][0]=(Params[2]*10);
//Q->me[1][1]=(Params[2]*100);
//Q->me[2][2]=(Params[2]*0.001);
//Q->me[2][2]=(Params[2]*0.001);
//m_output(Q);


//Matriz de ponderação das entradas
R=m_get(2,2);
m_ident(R);
sm_mlt(Params[3],R,R);
//R->me[1][1]=(Params[3]*100);
R->me[1][1]=(Params[3]*0.1);
//R->me[3][3]=(Params[3]/5);
//m_output(R);

//Matrizes de estados atuais e preditos
xb=m_get(3,Hp+1);
for(i=0;i<=Hp;i++){
xb->me[0][i]=x00[3*i];
xb->me[1][i]=x00[3*i+1];
xb->me[2][i]=x00[3*i+2];
}
//m_output(xb);

if (VPFMode==2){
if(i>=1){
	if(abs(xb->me[0][i])>=abs(xb->me[0][i-1])){
		demax=abs(xb->me[0][i]);
		}
	else{
	    demax=abs(xb->me[0][i-1]);
		}
	if(abs(xb->me[1][i])>=abs(xb->me[1][i-1])){
     	thetaemax=abs(xb->me[1][i]);
		}
	else{
	    thetaemax=abs(xb->me[1][i-1]);
		}
	if(abs(xb->me[2][i])>=abs(xb->me[2][i-1])){
     	etamax=abs(xb->me[2][i]);
		}
	else{
	    etamax=abs(xb->me[2][i-1]);
		}
}
}


//Matrizes de estados atuais e preditos
xbcon=m_get(3,Hp+1);
for(i=0;i<=Hp;i++){
xbcon->me[0][i]=x00[3*i];
xbcon->me[1][i]=x00[3*i+1];
xbcon->me[2][i]=x00[3*i+2];
}
//m_output(xbcon);

//Matrizes de entradas atuais
ub=m_get(2,Hc+1);
for(i=0;i<=Hc;i++){
ub->me[0][i]=u00[2*i];
ub->me[1][i]=u00[2*i+1];
}
//m_output(ub);

//Matrizes de entradas atuais
ubcon=m_get(2,Hc+1);
for(i=0;i<=Hc;i++){
ubcon->me[0][i]=u00[2*i];
ubcon->me[1][i]=u00[2*i+1];
}
//m_output(ubcon);

i=0;

//chamada função otimizador
donlp2();

TC1=0;
TC2=0;
result[0] = fx;
result[1] = x[1];
result[2] = x[2];
//result[3] = x[3];
//result[4] = x[4];
//result[5] = x[5];
f=0;

//getchar();
//res[0] = 111;0

//M_FREE(A);
//M_FREE(B);
//M_FREE(C);
//m_free(Q);
//m_free(R);
//m_free(ub);
//m_free(xb);
//m_free(ubcon);
//m_free(xbcon);
}

else if (VISUAL==0){
x_0[0] = 0;
for(i=1;i<=2*(Hc+1);i++){
	x_0[i] = u00[i-1];
}
i=0;


//Inicialização das matrizes do modelo

//Matriz de ponderação dos estados
Q=m_get(2,2);
m_ident(Q);
sm_mlt(Params[2],Q,Q);
Q->me[0][0]=(Params[2]*0.01);
//Q->me[0][0]=(Params[2]*0.01);
//m_output(Q);

//Matriz de ponderação das entradas
R=m_get(2,2);
m_ident(R);
sm_mlt(Params[3],R,R);
R->me[0][0]=(Params[3]*0.01);
//m_output(R);

//Matrizes de estados atuais e preditos
xb=m_get(2,Hp+1);
for(i=0;i<=Hp;i++){
xb->me[0][i]=x00[2*i];
xb->me[1][i]=x00[2*i+1];
}
//m_output(xb);


//Matrizes de estados atuais e preditos
xbcon=m_get(2,Hp+1);
for(i=0;i<=Hp;i++){
xbcon->me[0][i]=x00[2*i];
xbcon->me[1][i]=x00[2*i+1];

}
//m_output(xbcon);

//Matrizes de entradas atuais
ub=m_get(2,Hc+1);
for(i=0;i<=Hc;i++){
ub->me[0][i]=u00[2*i];
ub->me[1][i]=u00[2*i+1];
}
//m_output(ub);

//Matrizes de entradas atuais
ubcon=m_get(2,Hc+1);
for(i=0;i<=Hc;i++){
ubcon->me[0][i]=u00[2*i];
ubcon->me[1][i]=u00[2*i+1];
}
//m_output(ubcon);

i=0;

//chamada função otimizador
donlp2();

TC1=0;
TC2=0;
result[0] = fx;
result[1] = x[1];
result[2] = x[2];
f=0;

//getchar();
//res[0] = 111;0

//M_FREE(A);
//M_FREE(B);
//M_FREE(C);
//m_free(Q);
//m_free(R);
//m_free(ub);
//m_free(xb);
//m_free(ubcon);
//m_free(xbcon);

}

return;
}

//Resolução das equações diferenciais do modelo - Método de Euller
void error_model(){


int kk=0;
double u1,u2,u3=0;
double xe0,ye0,thetae0,eta0,beta0, de0=0;
double xe1,ye1,thetae1,eta1,beta1, de1=0;

if (VISUAL==0){
for (kk=0;kk<=Hp;kk++){
	
u1 = ub->me[0][kk];
u2 = ub->me[1][kk];

eta0 = xb->me[0][kk];
beta0 = xb->me[1][kk];


eta1=eta0+Ts*(u1);
beta1=beta0+Ts*(u2);

xbcon->me[0][kk] = eta1;
xbcon->me[1][kk] = beta1;

//printf("\nxe=%e\n",xe1);
//printf("\nye=%e\n",ye1);
//printf("\nalphae=%e\n",alphae1);
//printf("\nthetae=%e\n",thetae1);
}
}

else if (VISUAL==1){
	
for (kk=0;kk<=Hp;kk++){
	
u1 = ub->me[0][kk];
u2 = ub->me[1][kk];

de0 = xb->me[0][kk];
thetae0 = xb->me[1][kk];
eta0 = xb->me[2][kk];


//ATENÇÃO!!!


//printf("\nCurv=%e\n",Curv);
//printf("\nur=%e\n",ur);
eta1=eta0+Ts*(u2);  
thetae1=thetae0+Ts*(u1);

w_ref=(u1*cos(thetae1)+Curv*(eta1+ur))/(cos(thetae1)-Curv*de0);
de1=de0+Ts*(H*w_ref+((eta1+ur)+de0*w_ref)*tan(thetae1));




xbcon->me[0][kk] = de1;
xbcon->me[1][kk] = thetae1;
xbcon->me[2][kk] = eta1;


//printf("\nxe=%e\n",xe1);
//printf("\nye=%e\n",ye1);
//printf("\nalphae=%e\n",alphae1);
//printf("\nthetae=%e\n",thetae1);

}	
	
	
}
	

	
return;	
	
}

void calccoterms(int k){

double sc1,sc2,sc3=0;
double ufc1,ufc2,ufc3=0;
double xe,ye,thetae,eta,de=0;
double u1,u2,u3=0;



if (VISUAL==0){

de=5;
de=de+de;

CBp=Curv/(1+q*Curv);

u1 = ub->me[0][k];
u2 = ub->me[1][k];


eta = xb->me[0][k];
beta = xb->me[1][k];


if (id==0){
sc1=sc[k];	
sc2=sc[k+(Hp+1)];
sc3=sc[k+(2*(Hp+1))];	
ufc1=ufc[k];	
ufc2=ufc[k+(Hp+1)];
ufc3=ufc[k+(2*(Hp+1))];

//	printf("\n");
//	printf("sc1 %e",sc1);
//	printf("sc2 %e",sc2);
//  printf("sc3 %e",sc3);
//	printf("\n"); 

//ufc1=ufc1+Ts*(u3);	
eta=eta+Ts*(u1);

ufc1=eta+ur;




TC1=wu*pow(k1*ufc1-k2*ufc2,2)+wu*pow(k1*ufc1-k3*ufc3,2);

sp=qsip/(1-CBp*q); 
sc1=sc1+Ts*sp;


TC2=ws*pow(sc1-sc2+p12,2)+ws*pow(sc1-sc3+p13,2);
Fe=((sc1-sc2+p12)+(sc1-sc3+p13))/2;
	
	//printf("\n");
	//printf("TC2 %e",TC2);
	//printf("\n");

}
else if (id==1){
	

	
sc1=sc[k];	
sc2=sc[k+(Hp+1)];
ufc1=ufc[k];	
ufc2=ufc[k+(Hp+1)];

//ufc2=ufc2+Ts*(u3);

eta=eta+Ts*(u1);

ufc2=eta+ur;


TC1=wu*pow(k2*ufc2-k1*ufc1,2);


sp=qsip/(1-CBp*q); 
sc2=sc2+Ts*sp;
TC2=ws*pow(sc2-sc1+p21,2);
Fe=sc2-sc1+p21;
	
}
else if (id==2){
sc1=sc[k];	
sc3=sc[k+(Hp+1)];
ufc1=ufc[k];	
ufc3=ufc[k+(Hp+1)];


//ufc3=ufc3+Ts*(u3);	

eta=eta+Ts*(u1);

ufc3=eta+ur;


TC1=wu*pow(k3*ufc3-k1*ufc1,2);

//qsip=(ufc3)*cos(thetae)-u1;
sp=qsip/(1+CBp*q); //TRAZER q1 e Cpath
sc3=sc3+Ts*sp;

TC2=ws*pow(sc3-sc1+p31,2);	
Fe=sc3-sc1+p31;
}
	
	
}

else if (VISUAL==1){
	
xe=5;
ye=3;
ye=xe+ye;
	
CBp=Curv/(1+q*Curv);

u1 = ub->me[0][k];
u2 = ub->me[1][k];
//u3 = ub->me[2][k];
//u4 = ub->me[3][k];
//u5 = ub->me[4][k];


//qsip = (alphap-u2)/Curv;

de = xb->me[0][k];
thetae = xb->me[1][k];
eta = xb->me[2][k];


if (id==0){
sc1=sc[k];	
sc2=sc[k+(Hp+1)];
sc3=sc[k+(2*(Hp+1))];	
ufc1=ufc[k];	
ufc2=ufc[k+(Hp+1)];
ufc3=ufc[k+(2*(Hp+1))];

//	printf("\n");
//	printf("sc1 %e",sc1);
//	printf("sc2 %e",sc2);
//  printf("sc3 %e",sc3);
//	printf("\n"); 

//ufc1=ufc1+Ts*(u3);	
eta=eta+Ts*(u2);

ufc1=eta+ur;



thetae=thetae+Ts*(u1);


TC1=wu*pow(k1*ufc1-k2*ufc2,2)+wu*pow(k1*ufc1-k3*ufc3,2);


w_ref=(u1*cos(thetae)+Curv*ufc1)/(cos(thetae)-Curv*de);

qsip=(ufc1+w_ref*de)/cos(thetae);
sp=qsip/(1-CBp*q); 
sc1=sc1+Ts*sp;


TC2=ws*pow(sc1-sc2+p12,2)+ws*pow(sc1-sc3+p13,2);
Fe=((sc1-sc2+p12)+(sc1-sc3+p13))/2;
	
	//printf("\n");
	//printf("TC2 %e",TC2);
	//printf("\n");

}
else if (id==1){
sc1=sc[k];	
sc2=sc[k+(Hp+1)];
ufc1=ufc[k];	
ufc2=ufc[k+(Hp+1)];

//ufc2=ufc2+Ts*(u3);

eta=eta+Ts*(u2);

ufc2=eta+ur;

thetae=thetae+Ts*(u1);

TC1=wu*pow(k2*ufc2-k1*ufc1,2);

w_ref=(u1*cos(thetae)+Curv*ufc1)/(cos(thetae)-Curv*de);

qsip=(ufc1+w_ref*de)/cos(thetae);
sp=qsip/(1-CBp*q); //TRAZER q1 e Cpath
sc2=sc2+Ts*sp;
TC2=ws*pow(sc2-sc1+p21,2);
Fe=sc2-sc1+p21;
	
}
else if (id==2){
sc1=sc[k];	
sc3=sc[k+(Hp+1)];
ufc1=ufc[k];	
ufc3=ufc[k+(Hp+1)];


//ufc3=ufc3+Ts*(u3);	

eta=eta+Ts*(u2);

ufc3=eta+ur;

thetae=thetae+Ts*(u1);

TC1=wu*pow(k3*ufc3-k1*ufc1,2);

w_ref=(u1*cos(thetae)+Curv*ufc1)/(cos(thetae)-Curv*de);

qsip=(ufc1+w_ref*de)/cos(thetae);
sp=qsip/(1+CBp*q); //TRAZER q1 e Cpath
sc3=sc3+Ts*sp;

TC2=ws*pow(sc3-sc1+p31,2);	
Fe=sc3-sc1+p31;
}
	
	
}	
	
	
	
	
return;	
	
	
	
}




/* **************************************************************************** */
/*                              donlp2-intv size initialization                 */
/* **************************************************************************** */
void user_init_size(void){
	if (VISUAL==0){
    n      = 2*(Hc+1);
    //nlin   =  1;
    //nonlin =  5*(Hc+1)+1;
    nonlin =  2*(Hc+1);
    //
    //nonlin =  0;
    iterma = 4000;
    nstep = 20;
}
    else if (VISUAL==1){
	n      = 2*(Hc+1);
    nlin   =  0;
    if (VPFMode==1){
    nonlin =  3*(Hc+1)+4;
    }
    else{ 
    nonlin =  3*(Hc+1);
     }
    iterma = 4000;
    nstep = 20;			
	}
}

/* **************************************************************************** */
/*                              donlp2-intv standard setup                           */
/* **************************************************************************** */



void user_init(void) {
    //static INTEGER  i,j;
    //static double   xst0[4];
                                  
    /* name is ident of the example/user and can be set at users will       */
    /* the first static character must be alphabetic. 40 characters maximum */

    if (id==0){
    strcpy(name,"robot1");
	}
	else if (id==1){
    strcpy(name,"robot2");
	}
	else if (id==2){
    strcpy(name,"robot3");
	}
	
//	printf("\nVPFMode=%e\n",VPFMode);
//printf("\nVISUAL=%d\n",VISUAL);
   
    
    /* x is initial guess and also holds the current solution */
    /* problem dimension n = dim(x), nlin=number of linear constraints
       nonlin = number of nonlinear constraints  */
    
    analyt = FALSE;
    epsdif = 1.e-7;   /* gradients exact to machine precision */
    /* if you want numerical differentiation being done by donlp2 then:*/
    epsfcn   = 1.e-8; /* function values exact to machine precision */
    //epsdif = 1.e-14;   /* gradients exact to machine precision */
    /* if you want numerical differentiation being done by donlp2 then:*/
    //epsfcn   = 1.e-16; /* function values exact to machine precision */

    /*  bloc    = TRUE; */
    /* if one wants to evaluate all functions  in an independent process */
    /* difftype = 3; *//* the most accurate and most expensive choice */
    difftype = 1;    
    bloc = FALSE;
    nreset = n;
    
    
    del0 = 0.2e0;
    tau0 = 1.e0;
    tau  = 0.1e0;
    
    
    
    for (i = 0 ; i <= n ; i++) {
       x[i] = x_0[i]; 
    }
    i=0;
    
      // printf("\n\n o valor de x[3] é: %e \n\n", x[3]);  
    

    /*  set lower and upper bounds */
    big = 1.e20;
    
i=0;

if (VISUAL==0){
//Restrições nas entradas
for (i=0;i<=Hp;i++){
low[2*i+1]=-5;
low[2*i+2]=-5;
up[2*i+1]=5;
up[2*i+2]=5;
}
i=0;

//Restrições na saída ?? ANALISAR ??
for (i=0;i<=Hp;i++){
low[2*(Hp+1)+2*i+1]=-0.01;
low[2*(Hp+1)+2*i+2]=-0.2;
up[2*(Hp+1)+2*i+1]=0.01;
up[2*(Hp+1)+2*i+2]=0.2;
}
//low[5*(Hp+1)+5*Hp+6]=-0.05;
//up[5*(Hp+1)+5*Hp+6]=0.05;
i=0;
return;
}
 else if (VISUAL==1) 
 {
//Restrições nas entradas
for (i=0;i<=Hp;i++){
	if(VPFMode==2){
	low[2*i+1]=-1000;
	low[2*i+2]=-1000;
	up[2*i+1]=1000;	
    up[2*i+2]=1000;	
		
	}
	else{
	low[2*i+1]=-5;
	low[2*i+2]=-5;
	up[2*i+1]=5;	
    up[2*i+2]=5;		
	}
}
i=0;


//Restrições na saída ?? ANALISAR ??

for (i=0;i<=Hp;i++){
	if (VPFMode==2){
		low[2*(Hp+1)+3*i+1]=-0.05-abs(demax);
		low[2*(Hp+1)+3*i+2]=-0.3-abs(thetaemax);
		low[2*(Hp+1)+3*i+3]=-0.02-abs(etamax);
		up[2*(Hp+1)+3*i+1]=0.05+abs(demax);
		up[2*(Hp+1)+3*i+2]=0.3+abs(thetaemax);
		up[2*(Hp+1)+3*i+3]=0.02+abs(etamax);
		}
		else{
		low[2*(Hp+1)+3*i+1]=-0.05;
		low[2*(Hp+1)+3*i+2]=-0.3;
		low[2*(Hp+1)+3*i+3]=-0.01;
		up[2*(Hp+1)+3*i+1]=0.05;
		up[2*(Hp+1)+3*i+2]=0.3;
		up[2*(Hp+1)+3*i+3]=0.01;
		}
}

if (VPFMode==1){
low[2*(Hp+1)+3*Hp+4]=0;
low[2*(Hp+1)+3*Hp+5]=-50;
low[2*(Hp+1)+3*Hp+6]=-5;
low[2*(Hp+1)+3*Hp+7]=-5;
up[2*(Hp+1)+3*Hp+4]=0.025;
up[2*(Hp+1)+3*Hp+5]=0;
up[2*(Hp+1)+3*Hp+6]=5;
up[2*(Hp+1)+3*Hp+7]=5;
i=0;
}     
}

}

/* **************************************************************************** */
/*                                 special setup                                */
/* **************************************************************************** */
void setup(void) {
	//silent=TRUE; //RESOLVER PROBLEMA DE TRAVAMENTO QUANTO TRUE !!!!
    //intakt = TRUE;
    te0=FALSE;
    te1=FALSE;
    te2=FALSE;
    te3=FALSE;
    cold=FALSE;

    return;
}

/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk(void) {
	
//void donlp2(void);	
	//if (abs(Fe)>=0.05){
	//printf("\n");
	//printf("ERRO GRANDE !!!!");
	//printf("\n");
	//ws=ws+0.001;	
	//wu=wu+0.001;
	//donlp2();
	
	//}
    return;
}


/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */
void ef(double x[],double *fx) {



if (VISUAL==0){
//variáveis auxiliares
int k;
VEC *xb0;
VEC *ub0;
VEC *xb0TQ;
VEC *ub0TR;


xb0=v_get(2);
ub0=v_get(2);
xb0TQ=v_get(2);
ub0TR=v_get(2);	
	
for (k=0;k<=Hp;k++){
ub->me[0][k]=x[2*k+1];
ub->me[1][k]=x[2*k+2];
}
//m_output(ub);
k=0;

error_model();

f=0;
for (k=0;k<=Hp;k++){
//printf("\nk=%d\n",k);

	get_col(xbcon,k,xb0);
	get_col(ub,k,ub0);
	
	calccoterms(k);
	
	
	xb0TQ = vm_mlt(Q,xb0,xb0TQ);
	ub0TR = vm_mlt(R,ub0,ub0TR);
	
	f0 = in_prod(xb0TQ,xb0) + in_prod(ub0TR,ub0) + TC1 + TC2;
	//f0 = in_prod(xb0TQ,xb0) + in_prod(ub0TR,ub0);
	
	//printf("\nf0=%e\n",f0);
	f = f + Ts*f0;
	//printf("\nf=%e\n",f);
}
k=0;
*fx=f;
//getchar();

//printf("\nfx=%e\n",*fx);


v_free(xb0);
v_free(ub0);
v_free(xb0TQ);
v_free(ub0TR);
//getchar();
}

else if (VISUAL==1){
	//variáveis auxiliares
int k;
VEC *xb0;
VEC *ub0;
VEC *xb0TQ;
VEC *ub0TR;

double f0;
double TC=0;

VEC *xb0F;
VEC *xb0TP;


xb0=v_get(3);
ub0=v_get(2);
xb0TQ=v_get(3);
ub0TR=v_get(2);	

xb0F=v_get(3);
xb0TP=v_get(3);

	
for (k=0;k<=Hp;k++){
ub->me[0][k]=x[2*k+1];
ub->me[1][k]=x[2*k+2];
}
//m_output(ub);
k=0;

error_model();

f=0;
for (k=0;k<=Hp;k++){
//printf("\nk=%d\n",k);

	get_col(xbcon,k,xb0);
	get_col(ub,k,ub0);
	
	calccoterms(k);

	xb0TQ = vm_mlt(Q,xb0,xb0TQ);
	ub0TR = vm_mlt(R,ub0,ub0TR);
	
	f0 = in_prod(xb0TQ,xb0) + in_prod(ub0TR,ub0) + TC1 + TC2;
	//f0 = in_prod(xb0TQ,xb0) + in_prod(ub0TR,ub0);
	
	//printf("\nf0=%e\n",f0);
	f = f + Ts*f0;
	//printf("\nf=%e\n",f);
}
get_col(xbcon,Hp,xb0F);

xb0TP = vm_mlt(P,xb0F,xb0TP);

TC=in_prod(xb0TP,xb0F)/2;


if(VPFMode==1){
*fx=f+TC;
}
else{
*fx=f;
}
//getchar();

//printf("\nfx=%e\n",*fx);

//printf("\nxe=%e\n",xbcon->me[0][0]);
//printf("\nye=%e\n",xbcon->me[1][0]);
//printf("\nalphae=%e\n",xbcon->me[2][0]);
//printf("\nthetae=%e\n",xbcon->me[3][0]);


v_free(xb0);
v_free(xb0F);
v_free(ub0);
v_free(xb0TQ);
v_free(xb0TP);
v_free(ub0TR);

//getchar();	
	
	
	
	
	
}
return;
}

/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf(DOUBLE x[],DOUBLE gradf[]) {
    return;
}

/* **************************************************************************** */
/*  compute nonlinear constraints */
/* **************************************************************************** */
void econ(INTEGER type ,INTEGER liste[], DOUBLE x[],DOUBLE con[],
             LOGICAL err[]) {

int k;
	
	
if (VISUAL==0){
for (k=0;k<=Hp;k++){
ub->me[0][k]=x[2*k+1];
ub->me[1][k]=x[2*k+2];
}	
		 
error_model();
k=0;
for (k=0;k<=Hp;k++){
//printf("\nk=%d\n",k);
//calccoterms(k);
con[2*k+1]=xbcon->me[0][k];
con[2*k+2]=xbcon->me[1][k];
//con[5*k+5]=xbcon->me[4][k];
//err[1]=TRUE;
//calccoterms(k);

}

//printf("\nFe=%e\n",Fe);
//con[5*Hp+6]=Fe;

//if (id==0){
//con[5*Hp+6]=TC2;
//}
//else if (id==1){
//con[5*Hp+6]=Fe[2];
//} 
//else if (id==1){
//con[5*Hp+6]=Fe[3];
//} 

}
else if (VISUAL==1){
	
for (k=0;k<=Hp;k++){
ub->me[0][k]=x[2*k+1];
ub->me[1][k]=x[2*k+2];

}	

error_model();
			 
//constraints();

for (k=0;k<=Hp;k++){
//printf("\nk=%d\n",k);
con[3*k+1]=xbcon->me[0][k];
con[3*k+2]=xbcon->me[1][k];
con[3*k+3]=xbcon->me[2][k];
//err[2]=TRUE;
}
if(VPFMode==1){
con[3*Hp+4]=abs(xbcon->me[1][Hp])-abs(xbcon->me[0][Hp]);
con[3*Hp+5]=p11*xbcon->me[0][Hp]*(-w_ref*H-(eta+ur)*tan(xbcon->me[1][Hp]));
con[3*Hp+6]=(Curv*(eta+ur)-alpha*xbcon->me[1][Hp]*cos(xbcon->me[1][Hp]))/(cos(xbcon->me[1][Hp])-Curv*xbcon->me[0][Hp]);
con[3*Hp+7]=(Curv*(eta+ur)-alpha*xbcon->me[1][Hp]*cos(xbcon->me[1][Hp]))/(cos(xbcon->me[1][Hp])-Curv*xbcon->me[0][Hp]);

}


//printf("\nxe=%e\n",con[1]);
//printf("\nye=%e\n",con[2]);
//printf("\nalphae=%e\n",con[3]);
//printf("\nthetae=%e\n",con[4]);
}

return;

}

/* **************************************************************************** */
/*          compute the gradient of the  nonlinear constraints               */
/* **************************************************************************** */
void econgrad(INTEGER liste[] ,INTEGER shift ,  DOUBLE x[], DOUBLE **grad) {
	return;
}

/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern(INTEGER mode) {

return;

}
