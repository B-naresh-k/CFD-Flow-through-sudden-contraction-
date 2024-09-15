#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
	int i,j;  
	int m,n; 

	m=75;n=30;
	double U=0.2; 
	printf("M=%d\tN=%d\n",m,n);
	
	double h,l;
	
	l=15.0;h=2.0;
	printf("L=%lf\tH=%lfU=%lf\n",l,h,U); 
	
	double Re;
	Re=100.0;
	
	double z=0.9;
	double w[m][n], w_old[m][n],psi[m][n], psi_old[m][n],u[m][n],v[m][n];  
	double dx=l/m, dy=h/n, beta=dx/dy ;
	printf("dx=%lf\tdy=%lf\n\tbeta=%lf\n",dx,dy,beta);
	
	int iteration=0;
	double error_psi, error_w;
	double p=(m)*(n);

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			psi[i][j]=0.0;
			w[i][j]=0.0;
	
		}
	}
	
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==(m-1))
			{
				u[i][j]=1.0;
				v[i][j]=0.0;
				//psi[i][j+1]=dy*psi[i][j];
				if(j<n-1)
				{
				            psi[i][j+1]=psi[i][j]+1.0*dy;
				 }
				else{
						psi[m-1][n-1]=1.0*h;
				}
			
			
		
				
			}
			else if(j==0)
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			
			
			}
			else if(i==0)
			{
				if(j<=14)
				{
				
					u[i][j]=u[i+1][j];
					v[i][j]=v[i+1][j];
				    psi[i][j]=psi[i+1][j];
		       	}
				else if(j>14)
			    {
					u[i][j]=0.0;
					v[i][j]=0.0;
					psi[i][j]=1.0*h;
			   	}
					
			}		
			
			else if(j==(n-1))
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=psi[0][j];
				/*if(j<n-1)
				{
				            psi[i][j+1]=psi[i][j]+1.0*dy;
				 }
				else{
						psi[m-1][j]=1.0*h;
				}
			
			*/
				
			}
			
			
			else
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			}
			
		}
		
	
	}
	

		for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
				if(i==(m-1))
				{
					w[i][j]=-2.0*(psi[m-2][j]-psi[m-1][j])/(pow(dx,2));
					
				}
				else if(j==0)
				{
					w[i][j]=-2.0*(psi[i][1]-psi[i][0])/(pow(dy,2)); 
			
				}
				else if(i==0)
				{
					if(j<=14)
					{
							w[i+1][j]=w[i][j];
					}
					else if(j>14)
					{
							w[i][j]=-2.0*(psi[1][j]-psi[0][j])/(pow(dx,2));
									}				
				 
				
				
				}
				else if(j==n-1)
				{
					w[i][j]=-2.0*(psi[i][n-2]-psi[i][n-1])/(pow(dy,2)); 
				
				}
				else
				{
					w[i][j]=0.0;
				}
			
			}		
		}
	
	do
	{
		for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
				psi_old[i][j]=psi[i][j];
				w_old[i][j]=w[i][j];
			}
		}
		
		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				psi[i][j]= (0.5/(1.0+pow(beta,2)))*(dx*dx*w[i][j]+(pow(beta,2))*(psi[i][j+1]+psi[i][j-1])+psi[i+1][j]+psi[i-1][j]);
			
			}
		}
		
		
		
		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2.0*dy);
				v[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2.0*dx);
			
			}
		}
		
		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				//w[i][j]=(1.0-u[i][j]*dx*Re/2.0)*w[i+1][j]+(1.0+u[i][j]*dx*Re/2.0)*w[i-1][j]+(1.0-v[i][j]*dy*Re/2.0)*w[i][j+1]+(1.0+v[i][j]*dy*Re/2.0)*w[i][j-1];
			/*	w[i][j]=(0.5/(1.0+pow(beta,2)))*( (1.0-((psi[i][j+1]-psi[i][j-1]))*((beta*Re/4.0)))*w[i+1][j]
				+ (1.0+((psi[i][j+1]-psi[i][j-1])*((beta*Re)/4.0)))*w[i-1][j]
				+ ((1.0+((psi[i+1][j])*(Re/(4.0*beta))))*(pow(beta,2)*w[i][j+1]))
				+ ((1.0-((psi[i+1][j]-psi[i-1][j])*(Re/(4.0*beta))))*(pow(beta,2)*w[i][j-1])) );
				*/
				
				/*w[i][j]=(1-z)*((0.5/(1+pow(beta,2)))*((1.0-u[i][j]*dx*Re/2.0)*w[i+1][j]+(1.0+u[i][j]*dx*Re/2.0)*w[i-1][j]
				+ (1.0-v[i][j]*dy*Re/2.0)*w[i][j+1]*pow(beta,2)+(1.0+v[i][j]*dy*Re/2.0)*w[i][j-1]*pow(beta,2)))+z*w_old[i][j];
			
			
		*/
		
		
			w[i][j]=(1.0-z)*((0.5/(1.0+pow(beta,2)))*((1.0-((psi[i][j+1]-psi[i][j-1])*((beta*Re)/4.0)))*w[i+1][j]
				+(1.0+((psi[i][j+1]-psi[i][j-1])*((beta*Re)/4.0)))*w[i-1][j]
				+((1.0+((psi[i+1][j]-psi[i-1][j])*(Re/(4.0*beta))))*(pow(beta,2)*w[i][j+1]))
				+((1.0-((psi[i+1][j]-psi[i-1][j])*(Re/(4.0*beta))))*(pow(beta,2)*w[i][j-1]))))+z*w_old[i][j]; 
		
		
			}
		}
		
		
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
				if(i==(m-1))
			{
				u[i][j]=1.0;
				v[i][j]=0.0;
			
				if(j<n-1)
				{
		            psi[i][j+1]=psi[i][j]+1.0*dy;
				 }
				else{
						psi[m-1][n-1]=1.0*h;
				}
			
			
		
				
			}
			else if(j==0)
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			
			
			}
			else if(i==0)
			{
				if(j<=14)
				{
				
					u[i][j]=u[i+1][j];
					v[i][j]=v[i+1][j];
				    psi[i][j]=psi[i+1][j];
		       	}
				else if(j>14)
			    {
					u[i][j]=0.0;
					v[i][j]=0.0;
					psi[i][j]=1.0*h;
			   	}
					
			}		
			
			else if(j==n-1)
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=psi[0][j];
			
			}
				
			}
		}
			
			for(i=0;i<m;i++)
			{
			for(j=0;j<n;j++)
			{
					if(i==(m-1))
				{
					w[i][j]=-2.0*(psi[m-2][j]-psi[m-1][j])/(pow(dx,2));
					
				}
				else if(j==0)
				{
					w[i][j]=-2.0*(psi[i][1]-psi[i][0])/(pow(dy,2)); 
			
				}
				else if(i==0)
				{
					if(j<=14)
					{
							w[i+1][j]=w[i][j];
					}
					else if(j>14)
					{
							w[i][j]=-2.0*(psi[1][j]-psi[0][j])/(pow(dx,2));
									}				
				 
				
				
				}
				else if(j==n-1)
				{
					w[i][j]=-2.0*(psi[i][n-2]-psi[i][n-1])/(pow(dy,2)); 
				
				}
			
			}		
		}
		
		
	
		error_psi=0.0;
		error_w=0.0;
		printf("error_psi=%lf\terror_w=%lf\t",error_psi,error_w);
		for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
				error_psi=error_psi+pow(psi[i][j]-psi_old[i][j],2);
				error_w=error_w+pow(w[i][j]-w_old[i][j],2);
			}
		}
		error_psi=sqrt(error_psi/p);
		error_w=sqrt(error_w/p);
		
		printf("iteration=%d\terror_psi=%lf\terror_w=%lf\n",iteration,error_psi,error_w);
		iteration++;
	

	}while(error_psi>1.0e-6 || error_w>1.0e-6);
	

	double x=0.0,y=0.0;
	FILE *fp;
	fp=fopen("output.dat","w");
	fprintf(fp,"ZONE I=%d, J=%d\n",m,n);
	for(j=0;j<n;j++)
	{
		y=j*dy;
		for(i=0;i<m;i++)
		{
			x=i*dx;
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,u[i][j],v[i][j],psi[i][j],w[i][j]);
		}
		
	}
	fclose(fp);
	
	
	printf("......................final......................\n");
	for(i=0;i<m-1;i++)
	{
		for(j=0;j<n-1;j++)
		{
			psi[i][j]=0.0;
			w[i][j]=0.0;
	
		}
	
	}
	
	return 0;
	
}

