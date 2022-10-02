#include<stdio.h>
#include<math.h>

//MP1 By Litana & Yap
// Delta z = 0.2 m
// Delta t=0.5 years
double thomas(double a[][20], double b[][200], int var , int m, int n, double C2[][200]);
int main()
{
	double dz=.2,dt=.5; //pre-defined
	double V=0.0063,R,D,a,b,c,sum=0,t;
	int sz, st, i,j,k,dummy=1, check;

	sz=3/dz+1;
	st=60/dt+1;

	if (dz<0.25)
        sz=3.0000000/dz+2;

	double A[sz-1][sz-1],C[sz][st];
   	double A1[20][20], A2[20][20], C2[20][200], dum[20][20], dum2[20][200], dum3[20];

   	printf("\n==================================== MP 1: ONE-DIMENSIONAL POLLUTANT MIGRATION ====================================\n\n");
        printf("Given:\n> Delta z = 0.2 m\n> Delta t = 0.5 years\n");
        //printf("Please check the Executive report for the full solution of the mathematical model.");
        printf("[NOTE]: Please open the file 'output.csv' to check the summary of calculated data for ALL methods. ");
    printf("\n===================================================================================================================\n\n");
        char ch;
            printf("Please input anything to continue: ");
            scanf("%c", &ch);
            getchar();
            printf("\n\n\n");

	FILE *ptr;
        ptr = fopen("output.csv", "w");
        fclose(ptr);


	for(dummy=1;dummy<7;dummy++)
    {
        //initializing alpha beta gamma
        switch (dummy)
        {
            case 1:
            R=1.770;
            D=0.018;
            break;
            case 2:
            R=1.24;
            D=0.008;
            break;
            case 3:
            R=6;
            D=0.012;
            break;
            case 4:
            R=1;
            D=0.018;
            break;
            case 5:
            R=1.45;
            D=0.014;
            break;
            case 6:
            R=13.5;
            D=0.012;
            break;
        }

        a=1-(2.0*D/dz/dz/R*dt);
        b=(V/2/dz/R*dt)+(D/dz/dz/R*dt);
        c=0.0-(V/2/dz/R*dt)+(D/dz/dz/R*dt);

//boundary condition C[0][t]
        for(j=1;j<st;j++)
        {
            t=dt*j;
            if(t>0&&t<=10)
                C[0][j]=10000.0-100*pow((t-10),2);
            else if(t>10&&t<=18)
                C[0][j]=10000.0;
            else if(t>18&&t<=20)
                C[0][j]=0.0-5000*t+100000;
            else
                C[0][j]=0;

            C[0][j]*=c;
        }

//initial condition C[z][0]
        for(i=0;i<sz;i++)
        {
            C[i][0]=0;
        }


//Constants matrix
        for(i=0;i<sz-1;i++)
        {
            for(j=0;j<sz-1;j++)
            {
                if (i==j)
                    A[i][j]=a;
                else if (i==(j-1)&&j!=(sz))
                    A[i][j]=b;
                else if (i==(j+1))
                    A[i][j]=c;
                else
                    A[i][j]=0;
                if (i==(sz-2)&&j==(sz-3))
                    A[i][j]+=b;
            }
        }

        for(j=0;j<st;j++)
        {
//Matrix Multiplication
			for (i = 1; i < sz; i++)
            {
                for (k = 1; k < sz; k++)
                {
                  sum = sum + A[i-1][k-1]*C[k][j];
                }
//BC matrix
            if(i==1)
				C[i][j+1] = sum + C[0][j];
			else
				C[i][j+1] = sum;

            sum = 0;
            }
        }
     for(i=0;i<sz;i++)
        {
            C[i][0]=0;
        }
//Blackscreen & CSV Print

        FILE *ptr;
        ptr = fopen("output.csv", "a");

        if (dummy==1){
        printf("\nFTCS METHOD\n\n");
        fprintf(ptr,"\nFTCS METHOD\n\n");
        }

        switch (dummy)
        {
            case 1:
            printf("ACETONE:\n\t");
            fprintf(ptr, "ACETONE:\n\t ,");
            break;
            case 2:
            printf("DCM:\n\t");
            fprintf(ptr, "DCM:\n\t ,");
            break;
            case 3:
            printf("CALCIUM:\n\t");
            fprintf(ptr, "CALCIUM:\n\t ,");
            break;
            case 4:
            printf("CHOLRIDE:\n\t");
            fprintf(ptr, "CHOLRIDE:\n\t ,");
            break;
            case 5:
            printf("SODIUM:\n\t");
            fprintf(ptr, "SODIUM:\n\t ,");
            break;
            case 6:
            printf("MAGNESIUM:\n\t");
            fprintf(ptr, "MAGNESIUM:\n\t ,");
            break;
        }
        for(j=0;j<st;j++)
        {
            printf("%.1lf \t",dt*j );
            fprintf(ptr,"%.2lf,",dt*j );
        }
        printf("\n      ",dt*j );
        for(j=0;j<st;j++)
        {
            printf("________",dt*j );
        }

        for(i=0;i<sz;i++)
        {

            printf("\n%.1lf m|\t",dz*i);
            fprintf(ptr,"\n%.2lf,",dz*i);
            for(j=0;j<st;j++)
            {
                printf("%.2lf\t",C[i][j]);
                fprintf(ptr,"%.15lf,",C[i][j]);
            }
        }
        printf("\n\n");

        fprintf(ptr,"\n\n");
        fclose(ptr);
    }

    //CRANK-NICHOLSON METHOD

    for(dummy=1;dummy<7;dummy++)
    {
        //initializing alpha beta gamma
        switch (dummy)
        {
            case 1:
            R=1.770;
            D=0.018;
            break;
            case 2:
            R=1.24;
            D=0.008;
            break;
            case 3:
            R=6;
            D=0.012;
            break;
            case 4:
            R=1;
            D=0.018;
            break;
            case 5:
            R=1.45;
            D=0.014;
            break;
            case 6:
            R=13.5;
            D=0.012;
            break;
        }

        a=R/dt;
        b=V/4/dt;
        c=(D/2/dz/dz);

//boundary condition C[0][t]
        for(j=1;j<st;j++)
        {
            t=dt*j;
            if(t>0&&t<=10)
                C2[0][j]=10000.0-100*pow((t-10),2);
            else if(t>10&&t<=18)
                C2[0][j]=10000.0;
            else if(t>18&&t<=20)
                C2[0][j]=0.0-5000*t+100000;
            else
                C2[0][j]=0;

            C2[0][j]*=(c-b);
        }

//initial condition C[z][0]
        for(i=0;i<sz;i++)
        {
            C2[i][0]=0;
        }


//1st Constants matrix
        for(i=0;i<sz-1;i++)
        {
            for(j=0;j<sz-1;j++)
            {
                if (i==(sz-2)&&j==(sz-3))
                    A1[i][j]=-2.0*c;
				else if (i==j)
                    A1[i][j]=a+(2.0*c);
                else if (i==(j-1)&&j!=(sz))
                    A1[i][j]=0.0-c-b;
                else if (i==(j+1))
                    A1[i][j]=b-c;
                else
                    A1[i][j]=0;
            }
        }

//2nd Constants matrix
        for(i=0;i<sz-1;i++)
        {
            for(j=0;j<sz-1;j++)
            {
                if (i==(sz-2)&&j==(sz-3))
                    A2[i][j]=2.0*c;
				else if (i==j)
                    A2[i][j]=a-(2.0*c);
                else if (i==(j-1)&&j!=(sz))
                    A2[i][j]=c+b;
                else if (i==(j+1))
                    A2[i][j]=c-b;
                else
                    A2[i][j]=0;

            }
        }
        for(j=0;j<st;j++)
        {
//Multiplication
			for (i = 1; i < sz; i++)
            {
                for (k = 1; k < sz; k++)
                {
                  sum = sum + A2[i-1][k-1]*C2[k][j];
                }
//BC matrix
            if(i==1)
				{C2[i][j+1] = sum + C2[0][j+1]+C2[0][j];} //added current and prev bound matrix only
			else{
                C2[i][j+1] = sum ;
			}
            sum = 0;
            }
            for(i=0; i<sz;i++){
                dum2[i][j]=C2[i+1][j+1];
            }

    //THOMAS ALGORITHM

            thomas(A1, dum2, sz-1, i, j, C2);

        }
     for(i=0;i<sz;i++)
        {
            C2[i][0]=0;
        }

//Blackscreen & CSV Print

        FILE *ptr;
        ptr = fopen("output.csv", "a");

        if (dummy==1){
        printf("\nCRANK NICHOLSON METHOD\n\n");
        fprintf(ptr,"\nCRANK NICHOLSON METHOD\n\n");
        }

        switch (dummy)
        {
            case 1:
            printf("ACETONE:\n\t");
            fprintf(ptr, "ACETONE:\n\t ,");
            break;
            case 2:
            printf("DCM:\n\t");
            fprintf(ptr, "DCM:\n\t ,");
            break;
            case 3:
            printf("CALCIUM:\n\t");
            fprintf(ptr, "CALCIUM:\n\t ,");
            break;
            case 4:
            printf("CHOLRIDE:\n\t");
            fprintf(ptr, "CHOLRIDE:\n\t ,");
            break;
            case 5:
            printf("SODIUM:\n\t");
            fprintf(ptr, "SODIUM:\n\t ,");
            break;
            case 6:
            printf("MAGNESIUM:\n\t");
            fprintf(ptr, "MAGNESIUM:\n\t ,");
            break;
        }
        for(j=0;j<st;j++)
        {
            printf("%.1lf \t",dt*j );
            fprintf(ptr,"%.2lf,",dt*j );
        }
        printf("\n      ",dt*j );
        for(j=0;j<st;j++)
        {
            printf("________",dt*j );
        }
        for(i=0;i<sz;i++)
        {
            printf("\n%.1lf m|\t",dz*i);
            fprintf(ptr,"\n%.1lf,",dz*i);
            for(j=0;j<st;j++)
            {
                printf("%.2lf\t",C2[i][j]);
                fprintf(ptr,"%.15lf,",C2[i][j]);
            }
        }
        printf("\n\n");

        fprintf(ptr,"\n\n");
        fclose(ptr);
    }
    printf("\n[NOTE]: Please open the file 'output.csv' to check the summary of calculated data for ALL methods. \n\n\n");
}

double thomas(double a[20][20], double b[20][200], int var, int m, int n, double C2[20][200]){
    int i,j,k, h, maxcount;
    double x[20]={0};
    double constants[20]={0}, u[20][20];
    double dummy, max;
    double a1, b1, c1, p;
    double gamma[40]={0};


    //Tri-diagonal Operations
    for(i=0; i<var; i++){
        for(j=0; j<var; j++){
            if(i==j){
                u[i][j]=1;
                b1=a[i][j];}
            if(i!=0 && j==(i-1)){
                a1=a[i][j];}
            if(i!=var && j==(i+1)){
                c1=a[i][j];}
            if (i!=j){
                u[i][j]=0;}
            if(j==(i+1)){
                    if(i==0){
                        u[i][j]= c1/b1;
                        constants[i]=b[i][n]/a[i][j-1];
                        gamma[i]= c1/b1;
                        p=constants[i];
                    }else if (i!=0){
                        u[i][j]= c1/(b1-a1*(gamma[i-1]));
                        constants[i]=(b[i][n]-(a1*p))/(b1-(a1*gamma[i-1]));
                        p=constants[i];
                        gamma[i]=u[i][j];
                    }}}
                if (i==var-1){
                constants[i]=(b[i][n]-(a1*p))/(b1-(a1*gamma[i-1]));
            }}

    //BACK SUBSTITUTION
    for(i=var-1; i>=0; i--){
            if(i==(var-1)){
                x[i]=constants[i];
                dummy=x[i];
            }else{
                x[i]=constants[i]-gamma[i]*dummy;
                dummy=x[i];
            }}

    for(i=1; i<var+1; i++){
        C2[i][n+1]= x[i-1];

    }
    return(C2[20][200]);
}



