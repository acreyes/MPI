#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

#define send_data_tag 2001
#define return_data_tag 2002

double gam = 5./3.;//adiabatic index
double thet = 1.5;//slope limiter

typedef struct Vars{//conserved variables
  double mass;
  double xvelocity;
  double yvelocity;
  double energy;
  double press;
} Vars;

typedef struct Flux{
  double rhov;
  double momx;
  double momy;
  double energy;
} Flux;


void Print(int Nx, int Ny, Vars * U, double dx, double dy);

double Et(double rho, double xvel, double yvel, double P){//calculates the energy from pressure and other conserved variables
  return(.5*rho*(pow(xvel/rho,2)+pow(yvel/rho,2)) + P/(gam-1.));
}

double sgn(double x){//returns sign of the number
  if(x < 0) return(-1.);
  else if(x >0) return(1.);
  else return(0.);
}

double min(double x, double y, double z){//returns minium of set
  double num = x;
  if(num > y) num = y;
  if(num > z) num = z;
  return(num);
}
double minmod(double x, double y, double z){//slope limiter
  double minm = .25*fabs(sgn(x) + sgn(y))*(sgn(x) + sgn(z))*min(fabs(x), fabs(y), fabs(z));
  return(minm);
}

Vars U_L(Vars Ui, Vars Uimo, Vars Uipo){//interpolates the conserved variables at interface for the left state
  Vars UL;
  UL.mass = Ui.mass + 0.5*minmod(thet*(Ui.mass - Uimo.mass), 0.5*(Uipo.mass - Uimo.mass),thet*(Uipo.mass-Ui.mass));
  UL.xvelocity = Ui.xvelocity + 0.5*minmod(thet*(Ui.xvelocity - Uimo.xvelocity), 0.5*(Uipo.xvelocity - Uimo.xvelocity),thet*(Uipo.xvelocity-Ui.xvelocity));
  UL.yvelocity = Ui.yvelocity + 0.5*minmod(thet*(Ui.yvelocity - Uimo.yvelocity), 0.5*(Uipo.yvelocity - Uimo.yvelocity),thet*(Uipo.yvelocity-Ui.yvelocity));
  UL.press = Ui.press + 0.5*minmod(thet*(Ui.press - Uimo.press), 0.5*(Uipo.press - Uimo.press),thet*(Uipo.press-Ui.press));
  UL.energy = Et(UL.mass, UL.xvelocity, UL.yvelocity, UL.press);
  return(UL);
}

Vars U_R(Vars Ui, Vars Uipo, Vars Uipt){//does U_L for the right state
  Vars UR;
  UR.mass = Uipo.mass - 0.5*minmod(thet*(Uipo.mass-Ui.mass),0.5*(Uipt.mass-Ui.mass),thet*(Uipt.mass-Uipo.mass));
  UR.xvelocity = Uipo.xvelocity - 0.5*minmod(thet*(Uipo.xvelocity-Ui.xvelocity),.5*(Uipt.xvelocity-Ui.xvelocity),thet*(Uipt.xvelocity-Uipo.xvelocity));
  UR.yvelocity = Uipo.yvelocity - 0.5*minmod(thet*(Uipo.yvelocity-Ui.yvelocity),.5*(Uipt.yvelocity-Ui.yvelocity),thet*(Uipt.yvelocity-Uipo.yvelocity));
  UR.press = Uipo.press - 0.5*minmod(thet*(Uipo.press-Ui.press), 0.5*(Uipt.press-Ui.press), thet*(Uipt.press-Uipo.press));
  UR.energy = Et(UR.mass, UR.xvelocity, UR.yvelocity, UR.press);
  return(UR);
}

void init_sys(int Nx, int Ny, Vars * U, double dx, double dy, int type){//various initial conditions
  int i, j, randpx, randmx, randpy, randmy;
  if(type == 0){//Kelvin-Helmholtz instabilities
    srand(time(NULL));
    for(i=0;i<Nx;++i){
      for(j=0;j<Ny;++j){
	randpx = rand() % 50;//positive random number between 0 and 50
	randmx = -1*(rand() % 50);//negative random number between 0 and 50
	randpy = rand() % 50;//positive random number between 0 and 50
	randmy = -1*(rand() % 50);//negative random number between 0 and 50
	if(fabs(-.5+j*dy) < .2){
	  U[i*Ny+j].mass = 2.;
	  U[i*Ny+j].xvelocity = 0.5*U[i*Ny+j].mass + .01*(randpx + randmx)/50.;
	  U[i*Ny+j].yvelocity = .01*(randpy + randmy)/50.;
	  U[i*Ny+j].press = 2.5;
	}
	else{
	  U[i*Ny+j].mass = 1.;
	  U[i*Ny+j].xvelocity = -0.5*U[i*Ny+j].mass + .01*(randpx + randmx)/50.;
	  U[i*Ny+j].yvelocity = .01*(randpy + randmy)/50.;
	  U[i*Ny+j].press = 2.5;
	}
      }
    }
  }
  else if(type == 1){//blast
    for(i=0;i<Nx;++i){
      for(j=0;j<Ny;++j){
	if(pow(-.5+i*dx,2) < .005 && pow(-.5+j*dx,2) < .005){
	  U[i*Ny+j].mass = 1.;
	  U[i*Ny+j].xvelocity = 0.;
	  U[i*Ny+j].yvelocity = 0.;
	  U[i*Ny+j].press = 10.;
	}
	else{
	  U[i*Ny+j].mass = 1.;
	  U[i*Ny+j].xvelocity = 0.;
	  U[i*Ny+j].yvelocity = 0.;
	  U[i*Ny+j].press = .1;
	}
      }
    }
  }
  else if(type == 2){//shock tube
    for(i=0;i<Nx;++i){
      for(j=0;j<Ny;++j){
	if(j < Ny/2){
	  U[i*Ny+j].mass = 1.;
	  U[i*Ny+j].xvelocity = 0.;
	  U[i*Ny+j].yvelocity = 0.;
	  U[i*Ny+j].press = 1.;
	}
	else{
	  U[i*Ny+j].mass = .1;
	  U[i*Ny+j].xvelocity = 0.;
	  U[i*Ny+j].yvelocity = 0.;
	  U[i*Ny+j].press = .125;
	}
      }
    }
  }
  else if(type == 3){//shock tube
    for(i=0;i<Nx;++i){
      for(j=0;j<Ny;++j){
	if(i < Ny/2){
	  U[i*Ny+j].mass = 1.;
	  U[i*Ny+j].xvelocity = 0.;
	  U[i*Ny+j].yvelocity = 0.;
	  U[i*Ny+j].press = 1.;
	}
	else{
	  U[i*Ny+j].mass = .1;
	  U[i*Ny+j].xvelocity = 0.;
	  U[i*Ny+j].yvelocity = 0.;
	  U[i*Ny+j].press = .125;
	}
      }
    }
  }
  for(i=0;i<Nx;++i){
    for(j=0;j<Ny;++j){
      U[i*Ny+j].energy = Et(U[i*Ny+j].mass, U[i*Ny+j].xvelocity, U[i*Ny+j].yvelocity, U[i*Ny+j].press);
    }
  }

}

double lambda(double vel, Vars U, double pm){
  double c = sqrt(gam*U.press/U.mass);
  return(vel + pm*c);
}
double MAX(double a, double b, double c){
  double max = a;
  if(max < b) max = b;
  if(max < c) max = c;
  return(max);
}

double alpha(Vars UL, Vars UR, double pm, int xy){
  double alph;
  if(xy == 1){// x integration
    alph = MAX(0., pm*lambda(UR.xvelocity, UR, pm), pm*lambda(UL.xvelocity, UL,pm));
  }
  else if(xy == -1){// y integration
    alph = MAX(0., pm*lambda(UR.yvelocity, UR, pm), pm*lambda(UL.yvelocity, UL,pm));
  }
  return(alph);
}

Flux get_F(Vars U){//get the flux from the conserved variables
  Flux F;  
  F.rhov = U.xvelocity;
  F.momx = pow(U.xvelocity,2)/U.mass + U.press;
  F.momy = U.xvelocity*U.yvelocity/U.mass;
  F.energy = (U.energy + U.press)*U.xvelocity/U.mass;
  return(F);
}

Flux get_G(Vars U){
  Flux G;
  G.rhov = U.yvelocity;
  G.momx = U.yvelocity*U.xvelocity/U.mass;
  G.momy = pow(U.yvelocity,2)/U.mass + U.press;
  G.energy = (U.energy + U.press)*U.yvelocity/U.mass;
  return(G);
}

double Maxalpha(Vars Uimo, Vars Ui, Vars Uipo, Vars Uipt, int type){//calculates maxalpha for stability condition
  double maxalph = 0.;
  Vars UL = U_L(Ui, Uimo, Uipo); Vars UR = U_R(Ui, Uipo, Uipt); //Calculates the conserved variables at the interfaces
  double alphap = alpha(UL, UR, 1., type);
  double alpham = alpha(UL, UR, -1., type);
  maxalph = MAX(maxalph, alphap, alpham);
  return(maxalph);
}

Flux Fhll( Vars Uimo, Vars Ui, Vars Uipo, Vars Uipt){//calculates the HLL flux for a given interface in x
  Flux F_HLL;
  Vars UL = U_L(Ui, Uimo, Uipo); Vars UR = U_R(Ui, Uipo, Uipt); //Calculates the conserved variables at the interfaces
  double alphap = alpha(UL, UR, 1., 1);
  double alpham = alpha(UL, UR, -1., 1);
  Flux FL = get_F(UL); Flux FR = get_F(UR); //Calculates the Fluxes from the left and right at the interface
  F_HLL.rhov = (alphap*FL.rhov + alpham*FR.rhov - alphap*alpham*(UR.mass - UL.mass))/(alphap + alpham);
  F_HLL.momx = (alphap*FL.momx + alpham*FR.momx - alphap*alpham*(UR.xvelocity - UL.xvelocity))/(alphap + alpham);
  F_HLL.momy = (alphap*FL.momy + alpham*FR.momy - alphap*alpham*(UR.yvelocity - UL.yvelocity))/(alphap + alpham);
  F_HLL.energy = (alphap*FL.energy + alpham*FR.energy - alphap*alpham*(UR.energy - UL.energy))/(alphap + alpham);
  return(F_HLL);
}

Flux Ghll( Vars Uimo, Vars Ui, Vars Uipo, Vars Uipt){//calculates the HLL flux for a given interface in y
  Flux G_HLL;
  Vars UL = U_L(Ui, Uimo, Uipo); Vars UR = U_R(Ui, Uipo, Uipt); //Calculates the conserved variables at the interfaces
  double alphap = alpha(UL, UR, 1., -1);
  double alpham = alpha(UL, UR, -1., -1);
  Flux FL = get_G(UL); Flux FR = get_G(UR); //Calculates the Fluxes from the left and right at the interface
  G_HLL.rhov = (alphap*FL.rhov + alpham*FR.rhov - alphap*alpham*(UR.mass - UL.mass))/(alphap + alpham);
  G_HLL.momx = (alphap*FL.momx + alpham*FR.momx - alphap*alpham*(UR.xvelocity - UL.xvelocity))/(alphap + alpham);
  G_HLL.momy = (alphap*FL.momy + alpham*FR.momy - alphap*alpham*(UR.yvelocity - UL.yvelocity))/(alphap + alpham);
  G_HLL.energy = (alphap*FL.energy + alpham*FR.energy - alphap*alpham*(UR.energy - UL.energy))/(alphap + alpham);
  return(G_HLL);
}

Flux L(Flux FL, Flux FR,  Flux GL, Flux GR, double dy, double dx){
  Flux LU;
  LU.rhov = -1.*(FR.rhov - FL.rhov)/dx - 1.*(GR.rhov - GL.rhov)/dy;
  LU.momx = -1.*(FR.momx - FL.momx)/dx - 1.*(GR.momx - GL.momx)/dy;
  LU.momy = -1.*(FR.momy - FL.momy)/dx - 1.*(GR.momy - GL.momy)/dy;
  LU.energy = -1.*(FR.energy - FL.energy)/dx - 1.*(GR.energy - GL.energy)/dy;
  return(LU);
}

Vars make_U1(double dt, Vars U, Flux LU){
  Vars U1;
  U1.mass = U.mass + dt*LU.rhov;
  U1.xvelocity = U.xvelocity + dt*LU.momx;
  U1.yvelocity = U.yvelocity + dt*LU.momy;
  U1.energy = U.energy + dt*LU.energy;
  U1.press = (gam-1.)*(U1.energy - .5*(pow(U1.yvelocity,2) + pow(U1.xvelocity,2))/U1.mass);
  return(U1);
}

Vars make_U2(double dt, Vars U, Vars U1, Flux LU1){
  Vars U2;
  U2.mass = 0.75*U.mass + 0.25*U1.mass + 0.25*dt*LU1.rhov;
  U2.xvelocity = 0.75*U.xvelocity + 0.25*U1.xvelocity + 0.25*dt*LU1.momx;
  U2.yvelocity = 0.75*U.yvelocity + 0.25*U1.yvelocity + 0.25*dt*LU1.momy;
  U2.energy = 0.75*U.energy + 0.25*U1.energy + 0.25*dt*LU1.energy;
  U2.press = (gam-1.)*(U2.energy - .5*(pow(U2.yvelocity,2) + pow(U2.xvelocity,2))/U2.mass);
  return(U2);
}

Vars make_UN(double dt, Vars U, Vars U2, Flux LU2){
  Vars UN;
  UN.mass = U.mass/3. + 2.*U2.mass/3. + 2.*dt*LU2.rhov/3.;
  UN.xvelocity = U.xvelocity/3. + 2.*U2.xvelocity/3. + 2.*dt*LU2.momx/3.;
  UN.yvelocity = U.yvelocity/3. + 2.*U2.yvelocity/3. + 2.*dt*LU2.momy/3.;
  UN.energy = U.energy/3. + 2.*U2.energy/3. + 2.*dt*LU2.energy/3.;
  UN.press = (gam-1.)*(UN.energy - .5*(pow(UN.yvelocity,2) + pow(UN.xvelocity,2))/UN.mass);
  return(UN);
}


void Set_B2A(int Nx, int Ny, Vars * A, Vars * B){
  int i, j, count; count = 0;
  for(i=0;i<Nx;++i){
    for(j=0;j<Ny;++j){
      B[i*Ny+j] = A[i*Ny+j];
      //count += 1;
    }
  }
  //printf("%d\n", count);
}

void set_boundary(Vars * U, int Nx, int Ny, int type){
  int i,j; Vars Ul1, Ul2;

  if(type == 0){
    for(j=0;j<Ny;++j){//periodic boundary conditions for symmetry about x=1/2
      U[1*Ny+j] = U[0*Ny+j];
      U[0*Ny+j] = U[(Nx-1)*Ny+j];
      U[(Nx-1)*Ny+j] = U[(Nx-2)*Ny+j];
      U[(Nx-2)*Ny+j] = U[(Nx-3)*Ny+j];
    }
    
    for(i=0;i<Nx;++i){
      U[i*Ny+1] = U[i*Ny+0];
      U[i*Ny+0] = U[i*Ny+Ny-1];
      U[i*Ny+Ny-1] = U[i*Ny+Ny-2];
      U[i*Ny+Ny-2] = U[i*Ny+Ny-3];     
    }
    
  }
  else if(type == 1){
    for(j=0;j<Ny;++j){//periodic boundary conditions for symmetry about y=1/2
      U[1*Ny+j] = U[0*Ny+j];
      U[0*Ny+j] = U[(Nx-1)*Ny+j];
      U[(Nx-1)*Ny+j] = U[(Nx-2)*Ny+j];
      U[(Nx-2)*Ny+j] = U[(Nx-3)*Ny+j];
    }
  }
  else if(type == 2){
    for(i=0;i<Ny;++i){//periodic boundary conditions for symmetry about x=1/2
      U[i*Ny+1] = U[i*Ny+0];
      U[i*Ny+0] = U[i*Ny+Ny-1];
      U[i*Ny+Ny-1] = U[i*Ny+Ny-2];
      U[i*Ny+Ny-2] = U[i*Ny+Ny-3];
    }
  }
}

double advance_system(int Nx, int Ny, Vars * U, double dx, double dy, int my_id, int zones_to_do, int num_procs, MPI_Datatype mpi_Vars){
  int  N = Nx;  
  Vars * U_n, * Un,  * U_temp; int i, j, k, ierr, inst[2], bound;
  U_n = malloc(N*N*sizeof(Vars));
  Un = malloc(N*N*sizeof(Vars));
  U_temp = malloc(zones_to_do*(N-4)*sizeof(Vars));
  bound = 0;
  
  MPI_Status status;
    
  Flux FL, FR, GL, GR;
  double maxalphax, maxalphay, dt;
 
  if(my_id == 0){
    int start, end;
    start = 2;
    for(k=num_procs-1;k>0;--k){
      /*loop over all processors*/
      end = start + zones_to_do - 1;
      inst[0] = start; inst[1] = end;
      /*Send Start, End indices*/
      ierr = MPI_Send( &inst[0], 2, MPI_INT, k, send_data_tag, MPI_COMM_WORLD);
      start = end + 1;
    }
    end = start + zones_to_do -1;
    /*Calculate dt*/
    /*start calculating maxalpha in my zones*/
    double my_maxalpha = 0.; double maxalpha;
    for(i=2;i<Nx-2;++i){
      for(j=2;j<Nx-2;++j){
	maxalphax = Maxalpha(U[(i-2)*Ny+j],U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j], 1);
	maxalphay = Maxalpha(U[(i-2)*Ny+j],U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j], -1);
	my_maxalpha = MAX(my_maxalpha, maxalphax, maxalphay);
      }
    }
    /* /\*recieve maxalphas from other processors and compare to mine*\/ */
    /* for(k=1;k<num_procs;++k){ */
    /*   /\*recieve start, end for what the process did*\/ */
    /*   ierr = MPI_Recv( inst, 2, MPI_INT, k, send_data_tag, MPI_COMM_WORLD, &status); */
    /*   /\*recieve alpha*\/ */
    /*   ierr = MPI_Recv( &maxalpha, 1, MPI_DOUBLE, k, send_data_tag, MPI_COMM_WORLD, &status); */
    /*   /\*Add values to U_n*\/ */
    /*   //printf("dt = %f\n", maxalpha); */
    /*   if(maxalpha > my_maxalpha) my_maxalpha = maxalpha; */
    /* } */

    dt = .5*dx/my_maxalpha;
    //printf("dt = %f\n", my_maxalpha);
    /*Broadcast dt*/
    ierr = MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Set_B2A(Nx, Ny, U, U_n);

    /////////////////////////////////////////////////
    /*Calculate U1*/
    for(i=2;i<Nx-2;++i){
      for(j=start;j<end+1;++j){
	FL = Fhll(U[(i-2)*Ny+j],U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j]);
	FR = Fhll(U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j],U[(i+2)*Ny+j]);
	GL = Ghll(U[i*Ny+j-2],U[i*Ny+j-1],U[i*Ny+j],U[i*Ny+j+1]);
	GR = Ghll(U[i*Ny+j-1],U[i*Ny+j],U[i*Ny+j+1],U[i*Ny+j+2]);
	U_n[i*Ny+j] = make_U1(dt, U[i*Ny+j], L(FL, FR, GL, GR, dy, dx));
      }
    }

    /*recieve U_temps from other procs to complete U1*/
    for(k=1;k<num_procs;++k){
      /*recieve start, end for what the process did*/
      ierr = MPI_Recv( inst, 2, MPI_INT, k, send_data_tag, MPI_COMM_WORLD, &status);
      /*recieve U_temp*/
      ierr = MPI_Recv( U_temp, zones_to_do*(N-4), mpi_Vars, k, send_data_tag, MPI_COMM_WORLD, &status);
      /*Add values to U_n*/
      for(i=2;i<Nx-2;++i){
	for(j=inst[0];j<inst[1]+1;++j){
	  U_n[i*Ny+j] = U_temp[(i-2)*zones_to_do+j-inst[0]]; 
	}
      }
    }

    //set_boundary(U_n, Nx, Ny, bound);

    Set_B2A(Nx, Ny, U_n, Un);
    
    /*Broadcast U_n*/
    ierr = MPI_Bcast(U_n, N*N, mpi_Vars, 0, MPI_COMM_WORLD);
    
    //something wrong between here...//
    /////////////////////////////////////////////////
    /*Start Calculating U2*/
    for(i=2;i<Nx-2;++i){
      for(j=start;j<end+1;++j){
    	FL = Fhll(U_n[(i-2)*Ny+j],U_n[(i-1)*Ny+j],U_n[i*Ny+j],U_n[(i+1)*Ny+j]);
    	FR = Fhll(U_n[(i-1)*Ny+j],U_n[i*Ny+j],U_n[(i+1)*Ny+j],U_n[(i+2)*Ny+j]);
    	GL = Ghll(U_n[i*Ny+j-2],U_n[i*Ny+j-1],U_n[i*Ny+j],U_n[i*Ny+j+1]);
    	GR = Ghll(U_n[i*Ny+j-1],U_n[i*Ny+j],U_n[i*Ny+j+1],U_n[i*Ny+j+2]);
    	Un[i*Ny+j] = make_U2(dt, U[i*Ny+j], U_n[i*Ny+j], L(FL, FR, GL, GR, dy, dx));
      }
    }

    /*recieve U_temps from other procs to complete U2*/
    for(k=1;k<num_procs;++k){
      /*recieve start, end for what the process did*/
      ierr = MPI_Recv( inst, 2, MPI_INT, k, send_data_tag, MPI_COMM_WORLD, &status);
      /*recieve U_temp*/
      ierr = MPI_Recv( U_temp, zones_to_do*(N-4), mpi_Vars, k, send_data_tag, MPI_COMM_WORLD, &status);
      /*Add values to U_n*/
      for(i=2;i<Nx-2;++i){
    	for(j=inst[0];j<inst[1]+1;++j){
    	  Un[i*Ny+j] = U_temp[(i-2)*zones_to_do + j - inst[0]];
    	}
      }
    }
    //set_boundary(Un, Nx, Ny, bound);  
    Set_B2A(Nx, Ny, Un, U_n);
    /*Broadcast Un*/
    ierr = MPI_Bcast(Un, N*N, mpi_Vars, 0, MPI_COMM_WORLD);

    /////////////////////////////////////////////////
    /*Start calculating U at next timestep into U_n*/
    for(i=2;i<Nx-2;++i){
      for(j=start;j<end+1;++j){
    	FL = Fhll(Un[(i-2)*Ny+j],Un[(i-1)*Ny+j],Un[i*Ny+j],Un[(i+1)*Ny+j]);
    	FR = Fhll(Un[(i-1)*Ny+j],Un[i*Ny+j],Un[(i+1)*Ny+j],Un[(i+2)*Ny+j]);
    	GL = Ghll(Un[i*Ny+j-2],Un[i*Ny+j-1],Un[i*Ny+j],Un[i*Ny+j+1]);
    	GR = Ghll(Un[i*Ny+j-1],Un[i*Ny+j],Un[i*Ny+j+1],Un[i*Ny+j+2]);
    	U_n[i*Ny+j] = make_UN(dt, U[i*Ny+j], Un[i*Ny+j], L(FL, FR, GL, GR, dy, dx));
      }
    }

    /*recieve U_temps from other procs to complete U1*/
    for(k=1;k<num_procs;++k){
      /*recieve start, end for what the process did*/
      ierr = MPI_Recv( inst, 2, MPI_INT, k, send_data_tag, MPI_COMM_WORLD, &status);
      /*recieve U_temp*/
      ierr = MPI_Recv( U_temp, zones_to_do*(N-4), mpi_Vars, k, send_data_tag, MPI_COMM_WORLD, &status);
      /*Add values to U_n*/
      for(i=2;i<Nx-2;++i){
    	for(j=inst[0];j<inst[1]+1;++j){
    	  U_n[i*Ny+j] = U_temp[(i-2)*zones_to_do+j-inst[0]];
    	}
      }
    }
    ///and here///

    
    Set_B2A(Nx, Ny, U_n, U);
    set_boundary(U, Nx, Ny, bound);
  }
  else{
    /*I am a Slave*/
    /*Recieve my start and end indices*/
    
    ierr = MPI_Recv( inst, 2, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD, &status);
    /* /\*Calculate maxalpha in my zones*\/ */
    /* double maxalpha = 0.; */
    /* for(i=2;i<Nx-2;++i){ */
    /*   for(j=inst[0];j<inst[1]+1;++j){ */
    /* 	maxalphax = Maxalpha(U[(i-2)*Ny+j],U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j], 1); */
    /* 	maxalphay = Maxalpha(U[(i-2)*Ny+j],U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j], -1); */
    /* 	maxalpha = MAX(maxalpha, maxalphax, maxalphay); */
    /* 	//printf("slave alpha = %f %d %d\n", maxalpha, i, j); */
    /*   } */
    /* } */
    /* /\*send my result*\/ */


    /* ierr = MPI_Send( &inst[0], 2, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD); */
    /* ierr = MPI_Send( &maxalpha, 1, MPI_DOUBLE, 0, send_data_tag, MPI_COMM_WORLD); */
    /* /\*recieve broadcast of dt from root*\/ */
    ierr = MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /////////////////////////////////////////////////
    /*Calculate my part of U1, put it in U_temp*/
    for(i=2;i<Nx-2;++i){
      for(j=inst[0];j<inst[1]+1;++j){
	FL = Fhll(U[(i-2)*Ny+j],U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j]);
	FR = Fhll(U[(i-1)*Ny+j],U[i*Ny+j],U[(i+1)*Ny+j],U[(i+2)*Ny+j]);
	GL = Ghll(U[i*Ny+j-2],U[i*Ny+j-1],U[i*Ny+j],U[i*Ny+j+1]);
	GR = Ghll(U[i*Ny+j-1],U[i*Ny+j],U[i*Ny+j+1],U[i*Ny+j+2]);
	U_temp[(i-2)*zones_to_do+j-inst[0]] = make_U1(dt, U[i*Ny+j], L(FL, FR, GL, GR, dy, dx));
      }
    }
    /*Send my part of U1 back to Root*/
    ierr = MPI_Send( &inst[0], 2, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD);
    ierr = MPI_Send( &U_temp[0],(N-4)*zones_to_do , mpi_Vars, 0, send_data_tag, MPI_COMM_WORLD);

    /*Recieve Broadcast of U_n from Root*/
    ierr = MPI_Bcast(U_n, N*N, mpi_Vars, 0, MPI_COMM_WORLD);
    
    /////////////////////////////////////////////////
    /*Calculate my part of U2, put it in U_temp*/
    for(i=2;i<Nx-2;++i){
      for(j=inst[0];j<inst[1]+1;++j){
    	FL = Fhll(U_n[(i-2)*Ny+j],U_n[(i-1)*Ny+j],U_n[i*Ny+j],U_n[(i+1)*Ny+j]);
	FR = Fhll(U_n[(i-1)*Ny+j],U_n[i*Ny+j],U_n[(i+1)*Ny+j],U_n[(i+2)*Ny+j]);
	GL = Ghll(U_n[i*Ny+j-2],U_n[i*Ny+j-1],U_n[i*Ny+j],U_n[i*Ny+j+1]);
	GR = Ghll(U_n[i*Ny+j-1],U_n[i*Ny+j],U_n[i*Ny+j+1],U_n[i*Ny+j+2]);
    	U_temp[(i-2)*zones_to_do+j-inst[0]] = make_U2(dt, U[i*Ny+j], U_n[i*Ny+j], L(FL, FR, GL, GR, dy, dx));
      }
    }
    /*Send my part of U2 back to Root*/
    ierr = MPI_Send( &inst[0], 2, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD);
    ierr = MPI_Send( &U_temp[0], zones_to_do*(N-4), mpi_Vars, 0, send_data_tag, MPI_COMM_WORLD);

    /*Recieve Broadcast of Un from Root*/
    ierr = MPI_Bcast(Un, N*N, mpi_Vars, 0, MPI_COMM_WORLD);
    /////////////////////////////////////////////////
    /*Calculate my part of U at nex timestep, put it in U_temp*/
    for(i=2;i<Nx-2;++i){
      for(j=inst[0];j<inst[1]+1;++j){
    	FL = Fhll(Un[(i-2)*Ny+j],Un[(i-1)*Ny+j],Un[i*Ny+j],Un[(i+1)*Ny+j]);
    	FR = Fhll(Un[(i-1)*Ny+j],Un[i*Ny+j],Un[(i+1)*Ny+j],Un[(i+2)*Ny+j]);
    	GL = Ghll(Un[i*Ny+j-2],Un[i*Ny+j-1],Un[i*Ny+j],Un[i*Ny+j+1]);
    	GR = Ghll(Un[i*Ny+j-1],Un[i*Ny+j],Un[i*Ny+j+1],Un[i*Ny+j+2]);
    	U_temp[(i-2)*zones_to_do+j-inst[0]] = make_UN(dt, U[i*Ny+j], Un[i*Ny+j], L(FL, FR, GL, GR, dy, dx));
      }
    }
    /*Send my part of UN back to Root*/
    ierr = MPI_Send( &inst[0], 2, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD);
    ierr = MPI_Send( &U_temp[0], zones_to_do*(N-4), mpi_Vars, 0, send_data_tag, MPI_COMM_WORLD);
  }
  free(U_n); free(Un); free(U_temp);
  return(dt);
}

void Write_Cons(int Nx, int Ny, Vars * U, double dx, double dy, FILE * fid){
  int i, j; for(i=0;i<Nx;++i){
    for(j=2;j<Ny-2;++j){
      // printf("%d %d\n", i, j);
      fprintf(fid, "%e %e  %e  %e  %e  %e  %e\n", i*dx, j*dy, U[i*Ny+j].mass, U[i*Ny+j].xvelocity/U[i*Ny+j].mass, U[i*Ny+j].yvelocity/U[i*Ny+j].mass, U[i*Ny+j].energy, U[i*Ny+j].press);
    }
    fprintf(fid,"\n");
  }
}
void Print(int Nx, int Ny, Vars * U, double dx, double dy){
  int i, j; for(i=0;i<Nx;++i){
    for(j=0;j<Ny;++j){
      //if(i*dx == .5) printf("%f\n", U[i*Ny+j].mass);
      printf( "%e %e  %e  %e  %e  %e  %e\n", i*dx, j*dy, U[i*Ny+j].mass, U[i*Ny+j].xvelocity/U[i*Ny+j].mass, U[i*Ny+j].yvelocity/U[i*Ny+j].mass, U[i*Ny+j].energy, U[i*Ny+j].press);
    }
    printf("\n");
  }
}
int main(int argc, char ** argv){
  int my_id, root, ierr, num_procs;
  MPI_Status status;

  ierr = MPI_Init(&argc, &argv);//Creat processes
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  /*Make MPI data type for Vars*/
  const int nitems=5;
  int blocklengths[5] = {1, 1, 1, 1, 1};
  MPI_Datatype types[5] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Datatype mpi_Vars;
  MPI_Aint offsets[5];

  offsets[0] = offsetof(Vars, mass);
  offsets[1] = offsetof(Vars, xvelocity);
  offsets[2] = offsetof(Vars, yvelocity);
  offsets[3] = offsetof(Vars, energy);
  offsets[4] = offsetof(Vars, press);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_Vars);
  MPI_Type_commit(&mpi_Vars);
  /*start the program*/

  
  int N, type; N = num_procs*100;
  type = 1;
  int zones_to_do = N/num_procs;
  double dt; int count = 0;char str[80];

  FILE *fid, *finit;

  double dx = 1./(double)N;
  double t, T; t = 0.; T = .2;
  int num = 30;
  Vars * U = malloc((N+4)*(N+4)*sizeof(Vars)); init_sys(N+4, N+4, U, dx, dx, 1);
  if(my_id == 0){
    /*I am root*/
    
    finit = fopen("2Dinit.dat","w");
    Write_Cons(N+4, N+4, U, dx, dx, finit);
    fclose(finit);
    int count = 0;
    
  }
  while(t<T){
    //printf("before\n");
    dt = advance_system(N+4, N+4, U, dx, dx, my_id, zones_to_do, num_procs, mpi_Vars);
    t+=dt;    
    //break; 
    //printf("what time is it = %f\n", dt);
    /*Broadcast U*/
    ierr = MPI_Bcast(U, (N+4)*(N+4), mpi_Vars, 0, MPI_COMM_WORLD);
    /*
    if(my_id == 0){ 
      if( count % 1 == 0){
	sprintf(str, "T_%d.dat", count);
	fid = fopen(str, "w");
	Write_Cons(N+4, N+4, U, dx, dx, fid);
	fclose(fid);
	//printf("T=%f\n", t);
      }
      count += 1;
      }*/
  }
  if(my_id == 0){
    /*I am Root*/
    printf("%d\n", count);
    fid = fopen("22data.dat","w");
    Write_Cons(N+4, N+4, U, dx, dx, fid);
    fclose(fid);
  }
  free(U);
  MPI_Finalize();
}

