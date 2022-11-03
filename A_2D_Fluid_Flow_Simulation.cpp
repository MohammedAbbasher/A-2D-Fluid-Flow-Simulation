#include <iostream>
#include<cmath>
#include <fstream>
using namespace std;



int N= 256;                                      //If you make this variable a global variable, there will be no need to pass it inside functions
int Iterations = 4;


class FluidFlow{
    public:

        //This function is map to indexing, it helps us using a big 1D array instead of a 2D array
        int Index(int x, int y){      
            return x+y*N;
        }

        //This function is to guarantee that the fluid will not flow out of the box and will get back 
        void SetBounds(int b, float x []){                                       
        //First part for any point at the boundaries but not in the corners
            //This is for the two vertical lines in the 2D box
            for(int j = 1; j < N - 1; j++) {
                x[Index(0  , j)] = b == 1 ? -x[Index(1  , j)] : x[Index(1  , j)];
                x[Index(N-1, j)] = b == 1 ? -x[Index(N-2, j)] : x[Index(N-2, j)];
            }  
            //This part for the two horizontal lines in the 2D box
            for(int i = 1; i < N - 1; i++) {
                x[Index(i, 0 )] = b == 2 ? -x[Index(i, 1)] : x[Index(i, 1)];
                x[Index(i, N-1)] = b == 2 ? -x[Index(i, N-2)] : x[Index(i, N-2)];
            }
        //Second part is for the corners
            x[Index(0, 0)]       = 0.5f * (x[Index(1, 0)]  + x[Index(0, 1)]);     //the constant 0.5 is a double, when we put f we declare it as a float                                              
            x[Index(0, N-1)]     = 0.5f * (x[Index(1, N-1)]+ x[Index(0, N-2)]);
            x[Index(N-1, 0)]     = 0.5f * (x[Index(N-2, 0)] + x[Index(N-1, 1)]);                     
            x[Index(N-1, N-1)]   = 0.5f * (x[Index(N-2, N-1)]+ x[Index(N-1, N-2)]);
        }

        //This is a linear function describes the diffusion phenomena, we take the diffusion to be a linear combination of the values of the nearby cells
        void LinearSolver(int b, float x[], float x0[], float a, float c){    
            float cRecip = 1.0 / c;
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[Index(i, j)] = (x0[Index(i, j)] + a*(x[Index(i+1, j)] +x[Index(i-1, j)] +x[Index(i  , j+1)]+x[Index(i  , j-1)])) * cRecip;     //At the end we divide everything by c
                }
            }    
            SetBounds(b, x);
        }

        void Diffuse (int b, float x[], float x0[], float diff, float Dt){
            float a = Dt * diff * (N - 2) * (N - 2);          
            LinearSolver(b, x, x0, a, 1 + 6 * a);
        }

        void Project(float Velocity_x[], float Velocity_y[], float p[], float div[]){      
            for (int k = 1; k < N - 1; k++) {
                for (int j = 1; j < N - 1; j++) {
                    for (int i = 1; i < N - 1; i++) {
                        div[Index(i, j)] = -0.5f*(Velocity_x[Index(i+1, j)] -Velocity_x[Index(i-1, j)] +Velocity_y[Index(i  , j+1)]-Velocity_y[Index(i  , j-1)])/N;
                        p[Index(i, j)] = 0;
                    }
                }
            }
            SetBounds(0, div); 
            SetBounds(0, p);
            LinearSolver(0, p, div, 1, 6);
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    Velocity_x[Index(i, j)] -= 0.5f * N*(p[Index(i+1, j)]-p[Index(i-1, j)]);
                    Velocity_y[Index(i, j)] -= 0.5f * N*(p[Index(i, j+1)]-p[Index(i, j-1)]);      
                }
            }   
            SetBounds(1, Velocity_x);
            SetBounds(2, Velocity_y);
        }

        void Advect(int b, float d[], float d0[],  float Velocity_x[], float Velocity_y[], float Dt){
            float i0;
            float i1;
            float j0;
            float j1;
            float s0;
            float s1;
            float t0;
            float t1;
            float tmp1; 
            float tmp2;
            float x;
            float y;
            float Dtx = Dt * (N - 2);
            float Dty = Dt * (N - 2);
            float Nfloat = N; 
            float ifloat, jfloat;
    
            for(int j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
                for(int i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                    tmp1 = Dtx * Velocity_x[Index(i, j)];
                    tmp2 = Dty * Velocity_y[Index(i, j)];
                    x    = ifloat - tmp1; 
                    y    = jfloat - tmp2;
               
                    if(x < 0.5f) {
                    x = 0.5f;} 
                    if(x > Nfloat + 0.5f){
                    x = Nfloat + 0.5f; } 
                    i0 = floorf(x); 
                    i1 = i0 + 1.0f;
                    if(y < 0.5f){
                    y = 0.5f;}  
                    if(y > Nfloat + 0.5f){
                        y = Nfloat + 0.5f;}  
                    j0 = floorf(y);
                    j1 = j0 + 1.0f;    
                    s1 = x - i0; 
                    s0 = 1.0f - s1; 
                    t1 = y - j0; 
                    t0 = 1.0f - t1;                            
                    int i0i = i0;
                    int i1i = i1;
                    int j0i = j0;
                    int j1i = j1;
               
                    d[Index(i, j)] =  s0 * ( t0 * d0[Index(i0i, j0i)] + t1* d0[Index(i0i, j1i)]) +s1 * ( t0 * d0[Index(i1i, j0i)] + t1* d0[Index(i1i, j1i)]);
                }
            }
    
            SetBounds(b, d);
        }
};


class FluidBox: public FluidFlow {        
    private:   
                                         
        int Size;
        float Dt;                                           //time step
        float Diffusivity;                                  //diffusion constant
        float Viscousity;                                   //Viscousity of the fluid

        //Arrays need to be of size N*N, to do so we need to use dynamic memory allocation
        float S[256*256];                                   //previous density
        float Density[256*256];                             //current density
        float U[256*256];                                   //Current velocity in the x and y directions 
        float V[256*256];                                   
        float U0[256*256];                                  //previous velocities 
        float V0[256*256];

    public:

        //This Constructor
        FluidBox (float Dt, float Diffusivity, float Viscosity){                                    
            Size = N;
            Dt = Dt;
            Diffusivity = Diffusivity;
            Viscosity = Viscosity;
        }

        //We are talking about an incompressible fluid so change in density will actually change the properties of the fluid 
        void AddDensity(int x, int y, float amount){
            Density[Index(x, y)] += amount;                   //We are going to add this density at one point like an input but we still have not decided where exactly
        }

        void AddVelocity(int x, int y, float amountX, float amountY){
            U[Index(x, y)] += amountX;       
            V[Index(x, y)] += amountY;
        }

        void Procees() {                 

            //Firts Diffuse the velocities in x and y directions and then project
            Diffuse(1, U0, U, Viscousity, Dt);
            Diffuse(2, V0, V, Viscousity, Dt);
            Project(U0, V0, U, V);

            //Second Advect the two velocities and then project
            Advect(1, U, U0, V0, V0, Dt);
            Advect(2, V, V0, U0, V0, Dt);
            Project(U, V, U0, V0);
    
            //Third Diffuse and Advect the density
            Diffuse(0, S, Density, Diffusivity, Dt);
            Advect(0, Density, S, U, V, Dt);

       
        }   
};



int main(){
//The idea is to build the box, add some velocity and density and then watch that velocity and density advects and diffuses
    FluidBox Fluid(1,2,2) ;                               //We have chosen the density,viscousity and time step to be the choice of the user and the rest will stay secret   
    Fluid.AddDensity(0.5f*N,0.5f*N,10);    
    Fluid.AddVelocity(0.5f*N,0.5f*N,2,10);        

    return 0;
}

