#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;
int main()
{
double dx,dy;
double e,en,ed,emax= pow(10,-6); 
int x=101,y=101,it=0; 
double Re,U=1; 
int i,j;
cout<<"Enter Reynolds number :"<<endl;
cin>>Re;
dx=pow((x-1),-1); 
dy=pow((y-1),-1);
double beta=dx/dy; 
ofstream stream("stream.dat");   //for stream function
ofstream fp("u_centerline.dat");    //for centerline u plot
ofstream v_line("v_centerline.dat");    //for centerline u plot
ofstream fp_uvec("u_vectors.dat");         //for u vector plot
ofstream fp_vvec("v_vectors.dat");          // for v vector plot
 stream<<"TITLE = STREAM FUNCTION"<<endl<<"VARIABLES = \"x\", \"y\", \"phi\""<<endl;
    stream<<"Zone T = \"BLOCK1\", i=101, j=101, F=POINT\n\n"<<endl;
    fp<<"TITLE = u_centerline"<<endl<<"VARIABLES = \"y\", \"u\""<<endl;
    fp<<"Zone T = \"BLOCK1\", F=POINT\n\n"<<endl;
    v_line<<"TITLE = v_centerline"<<endl<<"VARIABLES = \"x\", \"v\""<<endl;
    v_line<<"Zone T = \"BLOCK1\", F=POINT\n\n"<<endl;
    fp_uvec<<"TITLE = u_vectors"<<endl<<"VARIABLES = \"x\", \"y\", \"U\""<<endl;
    fp_uvec<<"Zone T = \"BLOCK1\", i=101, j=101, F=POINT\n\n"<<endl;
    fp_vvec<<"TITLE = v_vectors"<<endl<<"VARIABLES = \"x\", \"y\", \"V\""<<endl;
    fp_vvec<<"Zone T = \"BLOCK1\", i=101, j=101, F=POINT\n\n"<<endl;

double w[y][x],w1[y][x],S[y][x],Sp[y][x]; 
double u[y][x],v[y][x]; 
for(i=0;i<x;i++) 
{ 
for(j=0;j<y;j++)
{
S[j][i]=0;
w[j][i]=0;
w1[j][i]=0;
u[j][i]=0;
v[j][i]=0;
}
}
for(i=1;i<y-1;i++) 
{
w[y-1][i]=-(2*U/dy);
w1[y-1][i]=-(2*U/dy);
}
do
{
en=0; 
ed=0;  
for(j=0;j<y;j++)
{
for(i=0;i<x;i++)
{
w1[j][i]=w[j][i]; 
}
}
for(j=1;j<y-1;j++)
{
for(i=1;i<x-1;i++)
{
S[j][i]=(S[j][i+1]+S[j][i-1]+(beta*beta)*(S[j+1][i]+S[j-1][i])+(dx*dx)*w[j][i])/(2*(1+pow(beta,2)));
}
}
for(j=1;j<x-1;j++) 
{
w[j][0]=-(2*S[j][1])/pow(dx,2); 
w[j][x-1]=-(2*S[j][x-2])/pow(dx,2);
}
for(i=1;i<x-1;i++)
{
w[0][i]=-(2*S[1][i])/pow(dy,2);
w[y-1][i]=-(2*S[y-2][i]+2*dy*U)/pow(dy,2);
}
for(j=1;j<y-1;j++)
{
for(i=1;i<x-1;i++)
{
w[j][i]=(w[j][i+1]+w[j][i-1]+(beta*beta)*(w[j+1][i]+w[j-1][i])-(0.25*Re*beta*(w[j][i+1]-
w[j][i-1]) *(S[j+1][i]-S[j-1][i]))
 +(0.25*Re*beta*(w[j+1][i]-w[j-1][i])*(S[j][i+1]-S[j][i-1])))/(2*(1+pow(beta,2)));
}
}
for(j=0;j<y;j++)
{
for(i=0;i<x;i++)
{
en=en+fabs(w[j][i]-w1[j][i]); 
ed=ed+fabs(w1[j][i]); 
}
}
e=en/ed ;
it++;
}while(e>emax);
cout<<"number of iterations =  "<<it<<"\n";
for(j=1;j<y-1;j++)
{
for(i=1;i<x-1;i++)
{
u[j][i]=(S[j+1][i]-S[j-1][i])/(2*dx); 
v[j][i]=-1*(S[j][i+1]-S[j][i-1])/(2*dy);
}
}
for(i=0;i<y;i++)
{
u[y-1][i]=U;
}
for(j=0;j<y;j++)
{
for(i=0;i<x;i++)
{
stream<<j*dx<<"\t"<<i*dy<<"\t"<<S[i][j]<<"\n";
fp_uvec<<j*dx<<"\t"<<i*dy<<"\t"<<u[i][j]<<"\n";
fp_vvec<<j*dx<<"\t"<<i*dy<<"\t"<<v[i][j]<<"\n";
}
}
for(i=0;i<x;i++)
{
fp<<u[i][50]<<"\t"<<i*dy<<"\n";
v_line<<i*dx<<"\t"<<v[50][i]<<"\n";
}
return 0;
}
