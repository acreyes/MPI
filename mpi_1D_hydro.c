#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>

#define send_data_tag 2001
#define return_data_tag 2002

double gam = 1.4;//adiabatic index
double thet = 1.5;//slope limiter

typedef struct Vars{//conserved variables
  double mass;
  double velocity;
  double energy;
  double press;
} Vars;

typedef struct Flux{
  double rhov;
  double mom;
  double energy;
} Flux;

double Et(double rho, double vel, double P){//calculates the energy from pressure and other conserved variables
  return(.5*rho*pow(vel/rho,2) + P/(gam-1.));
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
  UL.velocity = Ui.velocity + 0.5*minmod(thet*(Ui.velocity - Uimo.velocity), 0.5*(Uipo.velocity - Uimo.velocity),thet*(Uipo.velocity-Ui.velocity));
  UL.press = Ui.press + 0.5*minmod(thet*(Ui.press - Uimo.press), 0.5*(Uipo.press - Uimo.press),thet*(Uipo.press-Ui.press));
  UL.energy = Et(UL.mass, UL.velocity, UL.press);
  return(UL);
}

Vars U_R(Vars Ui, Vars Uipo, Vars Uipt){//does U_L for the right state
  Vars UR;
  UR.mass = Uipo.mass - 0.5*minmod(thet*(Uipo.mass-Ui.mass),0.5*(Uipt.mass-Ui.mass),thet*(Uipt.mass-Uipo.mass));
  UR.velocity = Uipo.velocity - 0.5*minmod(thet*(Uipo.velocity-Ui.velocity),.5*(Uipt.velocity-Ui.velocity),thet*(Uipt.velocity-Uipo.velocity));
  UR.press = Uipo.press - 0.5*minmod(thet*(Uipo.press-Ui.press), 0.5*(Uipt.press-Ui.press), thet*(Uipt.press-Uipo.press));
  UR.energy = Et(UR.mass, UR.velocity, UR.press);
  return(UR);
}

void init_sys(int N, Vars * U, double dx, int type){//various initial conditions
  if(type == 1){//regular shock tube
    int i;
    for(i=0;i<N;++i){
      if(i < N/2){
	U[i].mass = 1.0;
	U[i].velocity = 0.0;
	U[i].press = 1.0;
      }
      else{
	U[i].mass = 0.1;
	U[i].velocity = 0.0;
	U[i].press = 0.125;
	U[i].energy = Et(U[i].mass, U[i].velocity, U[i].press);
      }
    }
    for(i=0;i<N;++i) U[i].energy = Et(U[i].mass, U[i].velocity, U[i].press);
  }
  else if(type == 2){//shock tube with sin perubation on the right side
    int i;
    for(i=0;i<N;++i){
      double x = (double)i*dx;
      if(x < 1./4.){
	U[i].mass = 3.856;
	U[i].velocity = 5.629;
	U[i].press = 10.33;
      }
      else{
	U[i].mass = 1. - 0.2*sin(32*x);
	U[i].velocity = 0.0;
	U[i].press = 1.;
	U[i].energy = Et(U[i].mass, U[i].velocity, U[i].press);
      }
    }
    for(i=0;i<N;++i) U[i].energy = Et(U[i].mass, U[i].velocity, U[i].press);
  }
}

double lambda(Vars U, double pm){
  double c = sqrt(gam*U.press/U.mass);
  return(U.velocity/U.mass + pm*c);
}
double MAX(double a, double b, double c){
  double max = a;
  if(max < b) max = b;
  if(max < c) max = c;
  return(max);
}
double alpha(Vars UL, Vars UR, double pm){
  double alph = MAX(0., pm*lambda(UR, pm), pm*lambda(UL,pm));
  return(alph);
}

Flux get_F(Vars U){//get the flux from the conserved variables
  Flux F;  
  F.rhov = U.velocity;
  F.mom = pow(U.velocity,2)/U.mass + U.press;
  F.energy = (U.energy + U.press)*U.velocity/U.mass;
  return(F);
}

double Maxalpha(int N, Vars * U){//calculates maxalpha for stability condition
  int i; double maxalph = 0.;
  for(i=2;i<N+3;++i){
    Vars UL = U_L(U[i-1], U[i-2], U[i]); Vars UR = U_R(U[i-1], U[i], U[i+1]); //Calculates the conserved variables at the interfaces
    double alphap = alpha(UL, UR, 1.);
    double alpham = alpha(UL, UR, -1.);
    maxalph = MAX(maxalph, alphap, alpham);
  }
  return(maxalph);
}

Flux hll( Vars * U, int i){//calculates the HLL flux for a given interface
  Flux F_HLL;
  Vars UL = U_L(U[i-1], U[i-2], U[i]); Vars UR = U_R(U[i-1], U[i], U[i+1]); //Calculates the conserved variables at the interfaces
  double alphap = alpha(UL, UR, 1.);
  double alpham = alpha(UL, UR, -1.);
  Flux FL = get_F(UL); Flux FR = get_F(UR); //Calculates the Fluxes from the left and right at the interface
  //Flux FL = get_F(U[i-1]); Flux FR = get_F(U[i]); //Calculates the Fluxes from the left and right at the interface
  F_HLL.rhov = (alphap*FL.rhov + alpham*FR.rhov - alphap*alpham*(UR.mass - UL.mass))/(alphap + alpham);
  F_HLL.mom = (alphap*FL.mom + alpham*FR.mom - alphap*alpham*(UR.velocity - UL.velocity))/(alphap + alpham);
  F_HLL.energy = (alphap*FL.energy + alpham*FR.energy - alphap*alpham*(UR.energy - UL.energy))/(alphap + alpham);
  return(F_HLL);
}

Flux L(Flux FL, Flux FR, double dx){
  Flux LU;
  LU.rhov = -1.*(FR.rhov - FL.rhov)/dx;
  LU.mom = -1.*(FR.mom - FL.mom)/dx;
  LU.energy = -1.*(FR.energy - FL.energy)/dx;
  return(LU);
}

Vars make_U1(double dt, Vars U, Flux LU){
  Vars U1;
  U1.mass = U.mass + dt*LU.rhov;
  U1.velocity = U.velocity + dt*LU.mom;
  U1.energy = U.energy + dt*LU.energy;
  U1.press = (gam-1.)*(U1.energy - .5*pow(U1.velocity,2)/U1.mass);
  return(U1);
}

Vars make_U2(double dt, Vars U, Vars U1, Flux LU1){
  Vars U2;
  U2.mass = 0.75*U.mass + 0.25*U1.mass + 0.25*dt*LU1.rhov;
  U2.velocity = 0.75*U.velocity + 0.25*U1.velocity + 0.25*dt*LU1.mom;
  U2.energy = 0.75*U.energy + 0.25*U1.energy + 0.25*dt*LU1.energy;
  U2.press = (gam-1.)*(U2.energy - .5*pow(U2.velocity,2)/U2.mass);
  return(U2);
}

Vars make_UN(double dt, Vars U, Vars U2, Flux LU2){
  Vars UN;
  UN.mass = U.mass/3. + 2.*U2.mass/3. + 2.*dt*LU2.rhov/3.;
  UN.velocity = U.velocity/3. + 2.*U2.velocity/3. + 2.*dt*LU2.mom/3.;
  UN.energy = U.energy/3. + 2.*U2.energy/3. + 2.*dt*LU2.energy/3.;
  UN.press = (gam-1.)*(UN.energy - .5*pow(UN.velocity,2)/UN.mass);
  return(UN);
}
void Write_Cons(int N, Vars * U, double dx, FILE * fid){
  int i; for(i=0;i<N;++i){
    fprintf(fid, "%e %e  %e  %e  %e\n", -0.5 + i*dx, U[i].mass, U[i].velocity/U[i].mass, U[i].energy, U[i].press);
  }
}

void advance_system(int N, Vars * U, double dx, double dt, int type, int my_id, int zones_to_do, int num_procs, MPI_Datatype mpi_Vars){
  int ierr, inst[2], i, j;
  Flux FL, FR; Vars * U_n = malloc(N*sizeof(Vars));
  Vars * Un = malloc(N*sizeof(Vars));
  Vars * U_temp = malloc(zones_to_do*sizeof(Vars));
  MPI_Status status;

  if(my_id == 0){
    int start, end;

    for(i=0;i<N;++i){
      U_n[i] = U[i];
    }
    start = 2;
    for(j=num_procs-1;j>0;--j){
      /*loop over all processors*/
      end = start + zones_to_do - 1;
      inst[0] = start; inst[1] = end;
      ierr = MPI_Send( &inst[0], 2, MPI_INT, j, send_data_tag, MPI_COMM_WORLD);
      start = end + 1;
    }
    end = start + zones_to_do -1;

    /*begin calculating U at next timestep, U1*/
    for(i=start;i<end+1;++i){
      FL = hll(U, i); FR = hll(U, i+1);
      U_n[i] = make_U1(dt, U[i], L(FL, FR, dx));
    }
    
    /*recieve U_temps from other procs*/
    for(j=1;j<num_procs;++j){
      /*recieve start, end for what the process did*/
      ierr = MPI_Recv( inst, 2, MPI_INT, j, send_data_tag, MPI_COMM_WORLD, &status);
      /*recieve U_temp*/
      ierr = MPI_Recv( U_temp, zones_to_do, mpi_Vars, j, send_data_tag, MPI_COMM_WORLD, &status);
      /*Add values to U_n*/
      for(i=inst[0];i<inst[1]+1;++i){
	U_n[i] = U_temp[i-inst[0]];
      }
    }
    /*copy U_n to Un*/
    for(i=0;i<N;++i){
      Un[i] = U_n[i];
    }

    /*Broadcast U_n*/
    ierr = MPI_Bcast(U_n, N, mpi_Vars, 0, MPI_COMM_WORLD);

    /*Now do U2*/
    for(i=start;i<end+1;++i){
      FL = hll(U_n, i); FR = hll(U_n, i+1);
      Un[i] = make_U2(dt, U[i], U_n[i], L(FL, FR, dx));
    }
    /*recieve U_temps from other procs*/
    for(j=1;j<num_procs;++j){
      /*recieve start, end for what the process did*/
      ierr = MPI_Recv( inst, 2, MPI_INT, j, send_data_tag, MPI_COMM_WORLD, &status);
      /*recieve U_temp*/
      ierr = MPI_Recv( U_temp, zones_to_do, mpi_Vars, j, send_data_tag, MPI_COMM_WORLD, &status);
      /*Add values to U_n*/
      for(i=inst[0];i<inst[1]+1;++i){
	Un[i] = U_temp[i-inst[0]];
      }
    }
    /*copy Un to U_n*/
    for(i=0;i<N;++i){
      U_n[i] = Un[i];
    }

    /*Broadcast Un*/
    ierr = MPI_Bcast(Un, N, mpi_Vars, 0, MPI_COMM_WORLD);

    /*Now do UN*/
    for(i=start;i<end+1;++i){
      FL = hll(Un, i); FR = hll(Un, i+1);
      U_n[i] = make_UN(dt, U[i], Un[i], L(FL, FR, dx));
    }
    /*recieve U_temps from other procs*/
    for(j=1;j<num_procs;++j){
      /*recieve start, end for what the process did*/
      ierr = MPI_Recv( inst, 2, MPI_INT, j, send_data_tag, MPI_COMM_WORLD, &status);
      /*recieve U_temp*/
      ierr = MPI_Recv( U_temp, zones_to_do, mpi_Vars, j, send_data_tag, MPI_COMM_WORLD, &status);
      /*Add values to U_n*/
      for(i=inst[0];i<inst[1]+1;++i){
	U_n[i] = U_temp[i-inst[0]];
      }
    }
    /*copy Un to U_n*/
    for(i=2;i<N-2;++i){
      U[i] = U_n[i];
    }


    
  }
  else{
    /*I am a Slave*/
    ierr = MPI_Recv( inst, 2, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD, &status);

    /*calculate U at next timestep*/
    for(i=inst[0];i<inst[1]+1;++i){
      FL = hll(U, i); FR = hll(U, i+1);
      U_temp[i-inst[0]] = make_U1(dt, U[i], L(FL, FR, dx));
    }
    ierr = MPI_Send( &inst[0], 2, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD);
    ierr = MPI_Send( &U_temp[0], zones_to_do, mpi_Vars, 0, send_data_tag, MPI_COMM_WORLD);

    /*Recieve Broadcast of U_n from Root*/
    ierr = MPI_Bcast(U_n, N, mpi_Vars, 0, MPI_COMM_WORLD);

    /*Now do U2*/
    for(i=inst[0];i<inst[1]+1;++i){
      FL = hll(U_n, i); FR = hll(U_n, i+1);
      U_temp[i-inst[0]] = make_U2(dt, U[i], U_n[i], L(FL, FR, dx));
    }
    ierr = MPI_Send( &inst[0], 2, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD);
    ierr = MPI_Send( &U_temp[0], zones_to_do, mpi_Vars, 0, send_data_tag, MPI_COMM_WORLD);

    /*Recieve Broadcast of Un from Root*/
    ierr = MPI_Bcast(Un, N, mpi_Vars, 0, MPI_COMM_WORLD);

    /*Now do UN*/
    for(i=inst[0];i<inst[1]+1;++i){
      FL = hll(Un, i); FR = hll(Un, i+1);
      U_temp[i-inst[0]] = make_UN(dt, U[i], Un[i], L(FL, FR, dx));
    }
    ierr = MPI_Send( &inst[0], 2, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD);
    ierr = MPI_Send( &U_temp[0], zones_to_do, mpi_Vars, 0, send_data_tag, MPI_COMM_WORLD);

  }
  
}



int main(int argc, char ** argv){
  int my_id, root, ierr, num_procs;
  MPI_Status status;

  ierr = MPI_Init(&argc, &argv);//Creat processes
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  /*Make MPI data type for Vars*/
  const int nitems=4;
  int blocklengths[4] = {1,1,1,1};
  MPI_Datatype types[4] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Datatype mpi_Vars;
  MPI_Aint offsets[4];

  offsets[0] = offsetof(Vars, mass);
  offsets[1] = offsetof(Vars, velocity);
  offsets[2] = offsetof(Vars, energy);
  offsets[3] = offsetof(Vars, press);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_Vars);
  MPI_Type_commit(&mpi_Vars);
  /*start the program*/

  
<<<<<<< HEAD
  int N, type; N = num_procs*100;
=======
  int N, type; N = num_procs*50;
>>>>>>> 679565316e4d55272e8fe48d200ca7b01a54c1cc
  type = 2;
  int zones_to_do = N/num_procs;
  double dt = 0.001; int count = 0;char str[80];

  FILE *fid, *finit;
  int num,counts,dat_count;
  num = 3; counts = 0; dat_count=0;


  double dx = 1./(double)N;
<<<<<<< HEAD
  double t, T; t = 0.; T = .5;
=======
  double t, T; t = 0.; T = .6;
>>>>>>> 679565316e4d55272e8fe48d200ca7b01a54c1cc
  Vars * U = malloc((N+4)*sizeof(Vars)); init_sys(N+4, U, dx, type);
  if(my_id == 0){
    /*I am root*/
    
    finit = fopen("init.dat","w");
    Write_Cons(N+4, U, dx, finit);
  }
  while(t<T){
    advance_system(N+4, U, dx, dt, type, my_id, zones_to_do, num_procs, mpi_Vars);
    t+=dt;
    /*Broadcast U*/
    ierr = MPI_Bcast(U, N+4, mpi_Vars, 0, MPI_COMM_WORLD);
    //break;
    /*    if(my_id == 0){//i am root
      if(counts % num == 0){//write data
	char name[800];
	sprintf(name,"T_%d.dat",dat_count);
	fid = fopen(name, "w");
	Write_Cons(N+4,U,dx,fid);
	dat_count += 1;
      }
      counts+=1;
      }*/
    if(my_id == 0){//i am root
      printf("T=%f\n", t);
    }
  }
  if(my_id == 0){
    /*I am Root*/
    fid = fopen("data.dat","w");
    Write_Cons(N+4, U, dx, fid);
  }
  MPI_Finalize();
}
