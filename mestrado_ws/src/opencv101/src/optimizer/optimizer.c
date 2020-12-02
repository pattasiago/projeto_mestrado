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
double resultado[100];
double Fe=0;
int f1=0;
double Ts;
double Curv;
double CBp=0;
double f=0;
double f0=0;
double qsip=0;
double sp=0;
double s=0;
double ur=0;
double uf=0;
double uf0=0;
double w_ref=0;
double w_0=0;
double w_i=0;
double H=0;
double Hr1=0;
double Hpp=0;
double Hpp0=0;
double xep=0;
double yep=0;
double demax=0;
double thetaemax=0;
double etamax=0;
double lijemax=0;
double psiijemax=0;
double betaijemax=0;
int VISUAL = 0;
double alpha=0;
double beta=0;
double p11=0;
double p22=0;
double q11=0;
double q22=0;
double r11=0;
double VPFMode = 0;
double thetael = 0;
double thetae0_r1 = 0;
double uopt0_r1 = 0;
double thetael0 = 0;
double dep0 = 0;
double de00 = 0;
double L = 0;
double v_0 = 0;
double vj,wj=0;


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

double lijd=0;
double psiijd=0;


int i = 0;
int ii = 0;

MAT   *Q;
MAT   *R;
MAT   *P;

MAT   *xb;
MAT   *xbcon;
MAT   *ub;
MAT   *ubcon;

VEC   *lij000;
VEC   *psiij000;
VEC   *betaij000;
VEC   *gammaij000;




void optimizer(double Params[], double x00[], double u00[], double lij00[], double psiij00[], double betaij00[], double gammaij00[], double result[]){  
	
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
Hr1=Params[18];
H=Params[19];
alpha=Params[22];
VPFMode=Params[23];
L=Params[24];
v_0=Params[25];
w_0=Params[26];
beta=Params[27];
psiijd=Params[28];
vj=Params[29];
wj=Params[30];

//	printf("\n");
//	printf("id %d",id);
//	printf("\n");

//Leader control
if (VISUAL==1){
x_0[0] = 0;
for(i=1;i<=(Hc+1);i++){
	x_0[i] = u00[i-1];
}
i=0;



//Inicialização das matrizes do modelo

//Matriz de ponderação do custo terminal
P=m_get(2,2);
m_ident(P);
sm_mlt(Params[21],P,P);
P->me[0][0]=(Params[21]*20);
P->me[1][1]=(Params[21]*2);
p11=P->me[0][0];
p22=P->me[0][0];
//m_output(P);

//Matriz de ponderação dos estados
Q=m_get(2,2);
m_ident(Q);
sm_mlt(Params[2],Q,Q);
Q->me[0][0]=(Params[2]*100);
Q->me[1][1]=(Params[2]*1);
q11=Q->me[0][0];
q22=Q->me[0][0];
//m_output(Q);



//Matriz de ponderação das entradas
R=m_get(1,1);
m_ident(R);
sm_mlt(Params[3],R,R);
R->me[0][0]=(Params[3]);
//m_output(R);

//Matrizes de estados atuais e preditos
xb=m_get(2,Hp+1);
for(i=0;i<=Hp;i++){
xb->me[0][i]=x00[2*i];
xb->me[1][i]=x00[2*i+1];
//m_output(xb);
//printf("\nVPFMode=%e\n",VPFMode);
if (VPFMode==2){
if(i>=1){
//	 printf("\nAQUIn");
	if(fabs(xb->me[0][i])>=fabs(xb->me[0][i-1])){
		demax=fabs(xb->me[0][i]);
		}
	else{
	    demax=fabs(xb->me[0][i-1]);
		}
	if(fabs(xb->me[1][i])>=fabs(xb->me[1][i-1])){
     	thetaemax=fabs(xb->me[1][i]);
		}
	else{
	    thetaemax=fabs(xb->me[1][i-1]);
		}
}
}
}
//m_output(xb);
//printf("\ndemaxe=%e\n",demax);
//printf("\ndthetaemaxe=%e\n",thetaemax);
//printf("\ndetamaxe=%e\n",etamax);


//////////////////////
//Matrizes de estados atuais e preditos
xbcon=m_get(2,Hp+1);
for(i=0;i<=Hp;i++){
xbcon->me[0][i]=x00[2*i];
xbcon->me[1][i]=x00[2*i+1];
}
//m_output(xbcon);

//Matrizes de entradas atuais
ub=m_get(1,Hc+1);
for(i=0;i<=Hc;i++){
ub->me[0][i]=u00[i];
}
//m_output(ub);

//Matrizes de entradas atuais
ubcon=m_get(1,Hc+1);
for(i=0;i<=Hc;i++){
ubcon->me[0][i]=u00[i];
}
//m_output(ubcon);

i=0;

//chamada função otimizador
donlp2();

TC1=0;
TC2=0;
//for(i=0;(i<=1);i++){
result[0] = resultado[0];
result[1] = resultado[1];

//}
f=0;

//M_FREE(A);
//M_FREE(B);
//M_FREE(C);
//m_free(Q);
//m_free(R);
//m_free(ub);
//m_free(xb);
//m_free(ubcon);
//m_free(xbcon);
//return;
}

//Followers control
else if (VISUAL==0){
x_0[0] = 0;
for(i=1;i<=2*(Hc+1);i++){
	x_0[i] = u00[i-1];
}
i=0;


//Matriz de ponderação do custo terminal
P=m_get(3,3);
m_ident(P);
sm_mlt(Params[21],P,P);
P->me[0][0]=(Params[21]*0.02);
P->me[1][1]=(Params[21]*0.001);
P->me[2][2]=(Params[21]*0.0001);
p11=P->me[0][0];
p22=P->me[1][1];
//p33=P->me[2][2];
//m_output(P);

//Inicialização das matrizes do modelo

//Matriz de ponderação dos estados
Q=m_get(3,3);
m_ident(Q);
sm_mlt(Params[2],Q,Q);
Q->me[0][0]=(Params[2]*0.02);
Q->me[1][1]=(Params[2]*0.001);
Q->me[2][2]=(Params[2]*0.0001);
//m_output(Q);

//Matriz de ponderação das entradas
R=m_get(2,2);
m_ident(R);
sm_mlt(Params[3],R,R);
R->me[0][0]=(Params[3]*0.000001);
R->me[1][1]=(Params[3]*0.000001);
//m_output(R);

//Matrizes de estados atuais e preditos
xb=m_get(3,Hp+1);
for(i=0;i<=Hp;i++){
xb->me[0][i]=x00[3*i];
xb->me[1][i]=x00[3*i+1];
xb->me[2][i]=x00[3*i+2];
if (VPFMode==2){
if(i>=1){
//	 printf("\nAQUIn");
	if(fabs(xb->me[0][i])>=fabs(xb->me[0][i-1])){
		lijemax=fabs(xb->me[0][i]);
		}
	else{
	    lijemax=fabs(xb->me[0][i-1]);
		}
	if(fabs(xb->me[1][i])>=fabs(xb->me[1][i-1])){
     	psiijemax=fabs(xb->me[1][i]);
		}
	else{
	    psiijemax=fabs(xb->me[1][i-1]);
		}
	if(fabs(xb->me[2][i])>=fabs(xb->me[2][i-1])){
     	betaijemax=fabs(xb->me[2][i]);
		}
	else{
	    betaijemax=fabs(xb->me[2][i-1]);
		}
}
}

}
//m_output(xb);


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



lij000=v_get(Hp+1);
for(i=0;i<=Hp;i++){
lij000->ve[i]=lij00[i];
}

psiij000=v_get(Hp+1);
for(i=0;i<=Hp;i++){
psiij000->ve[i]=psiij00[i];
}

betaij000=v_get(Hp+1);
for(i=0;i<=Hp;i++){
betaij000->ve[i]=betaij00[i];
}

gammaij000=v_get(Hp+1);
for(i=0;i<=Hp;i++){
gammaij000->ve[i]=gammaij00[i];
}

i=0;

//chamada função otimizador
donlp2();

TC1=0;
TC2=0;
for(i=0;(i<=2);i++){
result[i] = resultado[i];
}
f=0;



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
double u1,u2=0;
double de0, thetae0;
double de1, thetae1=0;
double lij0,psiij0,gammaij0=0;
//double lij0,psiij0,gammaij0=0;
//double lij1,psiij1,betaij1=0;
//double lijp,psiijp,betaijp,gammaijp=0;
//double lijp,psiijp,betaijp=0;

double lije0,psiije0,betaije0=0;
double lije1,psiije1,betaije1=0;

if (VISUAL==0){
for (kk=0;kk<=Hp;kk++){
	
u1 = ub->me[0][kk];
u2 = ub->me[1][kk];

lije0=xb->me[0][kk];
psiije0=xb->me[1][kk]; 
betaije0=xb->me[2][kk];

lij0 = lij000->ve[kk]; 
psiij0 = psiij000->ve[kk]; 
//betaij0 = betaij000->ve[kk];
gammaij0 = gammaij000->ve[kk];  


//lijp=-u1+L*wj*sin(gammaij0);
//betaijp=u2;
////////wj=w_0-u2;
////////vj=-(u1+L*wj*sin(gammaij0))/cos(gammaij0);
//psiijp=(-vj*sin(gammaij0)/lij0)+(L*wj*cos(gammaij0)/lij0)+(v_0*sin(psiij0)/lij0)-w_0;

//gammaijp=betaijp+psiijp;

//lij1=lij0+Ts*(lijp);
//psiij1=psiij0+Ts*(psiijp);
//betaij1=betaij0+Ts*(betaijp);
//gammaij1=gammaij0+Ts*(gammaijp);

//////lijp=u1*cos(gammaij1)+L*u2*sin(gammaij1)-v_0*cos(psiij1);
//////psiijp=(-u1*sin(gammaij1)/lij1)+(L*u2*cos(gammaij1)/lij1)+(v_0*sin(psiij1)/lij1)-w_0;

//lije1=lije0+Ts*(-lijp);
//psiije1=psiije0+Ts*(-psiijp);
////betaije1=betaije0+Ts*(betaijp);

wj=w_0-u2;
vj=(v_0*cos(psiij0)-u1)/cos(gammaij0);
lije1=lije0+Ts*(u1-L*wj*sin(gammaij0));
psiije1=psiije0+Ts*((vj*sin(gammaij0)/lij0)-(L*wj*cos(gammaij0)/lij0)-(v_0*sin(psiij0)/lij0)+w_0);
betaije1=betaije0+Ts*(u2);
  

xbcon->me[0][kk] = lije1;
xbcon->me[1][kk] = psiije1;
xbcon->me[2][kk] = betaije1;




//printf("\nu1=%e\n",u1);
//printf("\nu1=%e\n",u2);
//printf("\npsiije1=%e\n",psiije1);
//printf("\nbetaije1=%e\n",betaije1);
//printf("\nbetaije1=%e\n",w_0);
//printf("\nthetae=%e\n",thetae1);
}
//gammaij000->ve[Hp]=gammaij1;
}

else if (VISUAL==1){
	
for (kk=0;kk<=Hp;kk++){
	
u1 = ub->me[0][kk];

de0 = xb->me[0][kk];
thetae0 = xb->me[1][kk];


thetae1=thetae0+Ts*(u1);

w_ref=(u1*cos(thetae1)+Curv*(ur))/(cos(thetae1)-Curv*de0);
de1=de0+Ts*(H*w_ref+((ur)+de0*w_ref)*tan(thetae1));

xbcon->me[0][kk] = de1;
xbcon->me[1][kk] = thetae1;

}





//printf("\nxe=%e\n",xe1);
//printf("\nye=%e\n",ye1);
//printf("\nalphae=%e\n",alphae1);
//printf("\nthetae=%e\n",thetae1);

}	
	
return;	
}
	

	

	




/* **************************************************************************** */
/*                              donlp2-intv size initialization                 */
/* **************************************************************************** */
void user_init_size(void){
	if (VISUAL==0){
    n      = 2*(Hc+1);
    if (VPFMode==1){
    nonlin =  3*(Hc+1)+5;
    }
    else{ 
    nonlin =  3*(Hc+1);
     }
    iterma = 4000;
    nstep = 20;
}
    else if (VISUAL==1){
	n      = 1*(Hc+1);
    nlin   =  0;
    if (VPFMode==1){
    nonlin =  2*(Hc+1)+4;
    }
    else{ 
    nonlin =  2*(Hc+1);
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
if(VPFMode==2){
low[2*i+1]=-1000;
low[2*i+2]=-1000;
up[2*i+1]=1000;
up[2*i+2]=1000;
}
else{
low[2*i+1]=-0.3;
low[2*i+2]=-5;
up[2*i+1]=0.3;
up[2*i+2]=5;		
}
}
i=0;

//Restrições na saída ?? ANALISAR ??
for (i=0;i<=Hp;i++){
		if (VPFMode==2){
		low[2*(Hp+1)+3*i+1]=-0.05-fabs(lijemax);
		low[2*(Hp+1)+3*i+2]=-0.5-fabs(psiijemax);
		low[2*(Hp+1)+3*i+3]=-0.5-fabs(betaijemax);
		up[2*(Hp+1)+3*i+1]=0.05+fabs(lijemax);
		up[2*(Hp+1)+3*i+2]=0.5+fabs(psiijemax);
		up[2*(Hp+1)+3*i+3]=0.5+fabs(betaijemax);
		}
		else{
		low[2*(Hp+1)+3*i+1]=-0.05;
		low[2*(Hp+1)+3*i+2]=-0.5;
		low[2*(Hp+1)+3*i+3]=-1.5;
		up[2*(Hp+1)+3*i+1]=0.05;
		up[2*(Hp+1)+3*i+2]=0.5;
		up[2*(Hp+1)+3*i+3]=1.5;
	}
}
if (VPFMode==1){
low[2*(Hp+1)+3*Hp+4]=0;
low[2*(Hp+1)+3*Hp+5]=0;
low[2*(Hp+1)+3*Hp+6]=0;
low[2*(Hp+1)+3*Hp+7]=-0.6/alpha;
low[2*(Hp+1)+3*Hp+8]=-(5+w_0)/beta;
up[2*(Hp+1)+3*Hp+4]=2;
up[2*(Hp+1)+3*Hp+5]=0.5;
up[2*(Hp+1)+3*Hp+6]=2;
up[2*(Hp+1)+3*Hp+7]=0.6/alpha;
up[2*(Hp+1)+3*Hp+8]=(5-w_0)/beta;
i=0;
}   

i=0;
return;
}

else if (VISUAL==1) 
 //printf("\nVPFMode=%e\n",VPFMode);
 
 {
//Restrições nas entradas
for (i=0;i<=Hp;i++){
	if(VPFMode==2){
	low[i+1]=-1000;
	up[i+1]=1000;	
	}
	else{
	low[i+1]=-5;
	up[i+1]=5;		
	}
}
i=0;



for (i=0;i<=Hp;i++){
	if (VPFMode==2){
		low[1*(Hp+1)+2*i+1]=-0.05-fabs(demax);
		low[1*(Hp+1)+2*i+2]=-0.3-fabs(thetaemax);
		up[1*(Hp+1)+2*i+1]=0.05+fabs(demax);
		up[1*(Hp+1)+2*i+2]=0.3+fabs(thetaemax);
		}
		else{
		low[1*(Hp+1)+2*i+1]=-0.10;
		low[1*(Hp+1)+2*i+2]=-0.10;
		up[1*(Hp+1)+2*i+1]=0.10;
		up[1*(Hp+1)+2*i+2]=0.10;
		}
}
if (VPFMode==1){
low[(Hp+1)+2*Hp+3]=0;
low[(Hp+1)+2*Hp+4]=0;
low[(Hp+1)+2*Hp+5]=0;
low[(Hp+1)+2*Hp+6]=-50;
up[(Hp+1)+2*Hp+3]=10;
up[(Hp+1)+2*Hp+4]=1;
up[(Hp+1)+2*Hp+5]=1;
up[(Hp+1)+2*Hp+6]=50;

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
	
if(VISUAL==1){	
resultado[0] = fx;
//for(i=1;(i<=1);i++){
resultado[1] = x[1];
//printf("\n");
//printf("result %d =  %e",i,resultado[i]);
//printf("\n");
//}
}
else if (VISUAL==0){
resultado[0] = fx;
for(i=1;(i<=2);i++){
resultado[i] = x[i];
//printf("\n");
//printf("result %d =  %e",i,resultado[i]);
//printf("\n");
}	
	
}
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
	
	
	xb0TQ = vm_mlt(Q,xb0,xb0TQ);
	ub0TR = vm_mlt(R,ub0,ub0TR);
	
    f0 = in_prod(xb0TQ,xb0) + in_prod(ub0TR,ub0);
	
	f = f + Ts*f0;

}
k=0;
get_col(xbcon,Hp,xb0F);

xb0TP = vm_mlt(P,xb0F,xb0TP);

TC=in_prod(xb0TP,xb0F)/2;


if(VPFMode==1){
*fx=f+TC;
}
else{
*fx=f;
}
;


v_free(xb0);
v_free(ub0);
v_free(xb0TQ);
v_free(ub0TR);

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


xb0=v_get(2);
ub0=v_get(1);
xb0TQ=v_get(2);
ub0TR=v_get(1);	

xb0F=v_get(2);
xb0TP=v_get(2);

	
for (k=0;k<=Hp;k++){
ub->me[0][k]=x[k+1];
}
//m_output(ub);
k=0;

error_model();

f=0;
for (k=0;k<=Hp;k++){
//printf("\nk=%d\n",k);

	get_col(xbcon,k,xb0);
	get_col(ub,k,ub0);
	

	xb0TQ = vm_mlt(Q,xb0,xb0TQ);
	ub0TR = vm_mlt(R,ub0,ub0TR);
	
	f0 = in_prod(xb0TQ,xb0) + in_prod(ub0TR,ub0);
	
	//printf("\nf0=%e\n",f0);
	f = f + Ts*f0;
	//printf("\nf=%e\n",f);
}
k=0;
get_col(xbcon,Hp,xb0F);

xb0TP = vm_mlt(P,xb0F,xb0TP);

TC=in_prod(xb0TP,xb0F)/2;


if(VPFMode==1){
*fx=f+TC;
}
else{
*fx=f;
}


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
con[3*k+1]=xbcon->me[0][k];
con[3*k+2]=xbcon->me[1][k];
con[3*k+3]=xbcon->me[2][k];
}

//printf("AQUI");


if(VPFMode==1){
con[3*Hp+4]=fabs(xbcon->me[2][Hp])-fabs(xbcon->me[1][Hp]);
con[3*Hp+5]=xbcon->me[0][Hp]*wj*sin(gammaij000->ve[Hp]);
con[3*Hp+6]=(xbcon->me[1][Hp])*(-vj*sin(gammaij000->ve[Hp])+wj*cos(gammaij000->ve[Hp])+v_0*sin(psiij000->ve[Hp])-w_0);
con[3*Hp+7]=xbcon->me[0][Hp];
con[3*Hp+8]=xbcon->me[2][Hp];
}

}
else if (VISUAL==1){
	
for (k=0;k<=Hp;k++){
ub->me[0][k]=x[k+1];
}	

error_model();
			 
//constraints();
for (k=0;k<=Hp;k++){
//printf("\nk=%d\n",k);
con[2*k+1]=xbcon->me[0][k];
con[2*k+2]=xbcon->me[1][k];
}

//w_ref=(ub->me[0][0]*cos(xbcon->me[1][Hp])+Curv*(ur))/(cos(xbcon->me[1][Hp])-Curv*xbcon->me[0][Hp]);


if(VPFMode==1){
con[2*Hp+3]=p22*alpha-q22-r11*(alpha*alpha)-p11*w_ref*tan(xbcon->me[1][Hp])-q11;
con[2*Hp+4]=fabs(xbcon->me[1][Hp])-fabs(xbcon->me[0][Hp]);
con[2*Hp+5]=-p11*xbcon->me[0][Hp]*(w_ref*H+v_0*tan(xbcon->me[1][Hp]));
con[2*Hp+6]=cos(xbcon->me[1][Hp])*(alpha*xbcon->me[1][Hp]-w_ref);
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
