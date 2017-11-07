

// pfCPU.cpp : 定义控制台应用程序的入口点。
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <windows.h>

#include <string.h>
#include "nicslu.h"

//#ifdef _DLL //动态连接
//
//#pragma comment(lib,"mkl_intel_lp64_dll.lib")
//#pragma comment(lib,"mkl_intel_thread_dll.lib")
//#pragma comment(lib,"mkl_core_dll.lib")
//
//#else //静态连接
//
//#pragma comment(lib,"mkl_intel_lp64.lib")
//#pragma comment(lib,"mkl_intel_thread.lib")
//#pragma comment(lib,"mkl_core.lib")
//
//#endif // _DLL



double tsolve=0,tpf=0,tY=0,tformJ0=0,tLU=0,tcnt=0;
double tpf0,tsolve0;
#define M_PI 3.141592653589793
const int  N=13659;
int formJ(int *Ji,int *Jj,double *J,int *from,int *to,double *G,double *B,double *V,double *ang,double *P,double *Q,int n,int r,int lineN,int *NodetoFuncP,int *NodetoFuncQ,int *FunctoNode,int *type){
	int nnzJ=-1;
	double value;
	for(int i=0;i<lineN;i++){
		if((type[from[i]]!=1)&&(type[to[i]]!=1)){
			//H中两个非零元素
			value=V[from[i]]*V[to[i]]*B[i];//*cos(ang[from[i]]-ang[to[i]])-G[i]*sin(ang[from[i]]-ang[to[i]]);
			//if(abs(value)>0.000000001){
				nnzJ++;
				Ji[nnzJ]=NodetoFuncP[from[i]];
				Jj[nnzJ]=NodetoFuncP[to[i]];
				J[nnzJ]=value;
				//if(nnzJ==985)
				//	printf("//");
			//}
			value=V[from[i]]*V[to[i]]*B[i];//*cos(ang[to[i]]-ang[from[i]])-G[i]*sin(ang[to[i]]-ang[from[i]]);
			//if(abs(value)>0.000000001){
				nnzJ++;
				Ji[nnzJ]=NodetoFuncP[to[i]];
				Jj[nnzJ]=NodetoFuncP[from[i]];
				J[nnzJ]=value;
				//if(nnzJ==985)
				//	printf("//");
			//}
			//L中两个非零元素
			if((type[from[i]]==3)&&(type[to[i]]==3)){
				value=V[from[i]]*V[to[i]]*B[i];//*cos(ang[from[i]]-ang[to[i]])-G[i]*sin(ang[from[i]]-ang[to[i]]);
				//if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncQ[from[i]];
					Jj[nnzJ]=NodetoFuncQ[to[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				//}
				value=V[from[i]]*V[to[i]]*B[i];//*cos(ang[to[i]]-ang[from[i]])-G[i]*sin(ang[to[i]]-ang[from[i]]);
				//if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncQ[to[i]];
					Jj[nnzJ]=NodetoFuncQ[from[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				//}
			}
			//N中两个非零元素
			if(type[to[i]]==3){
				value=-V[from[i]]*V[to[i]]*G[i];//*cos(ang[from[i]]-ang[to[i]])-B[i]*sin(ang[from[i]]-ang[to[i]]);
	/*			if(abs(value)>0.000000001){*/
					nnzJ++;
					Ji[nnzJ]=NodetoFuncP[from[i]];
					Jj[nnzJ]=NodetoFuncQ[to[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				//}
			}
			if(type[from[i]]==3){
				value=-V[from[i]]*V[to[i]]*G[i];//*cos(ang[to[i]]-ang[from[i]])-B[i]*sin(ang[to[i]]-ang[from[i]]);
				//if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncP[to[i]];
					Jj[nnzJ]=NodetoFuncQ[from[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				//}
			}
			//M中两个非零元素
			if(type[from[i]]==3){
				value=V[from[i]]*V[to[i]]*G[i];//*cos(ang[from[i]]-ang[to[i]])+B[i]*sin(ang[from[i]]-ang[to[i]]);
				//if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncQ[from[i]];
					Jj[nnzJ]=NodetoFuncP[to[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				//}
			}
			if(type[to[i]]==3){
				value=V[from[i]]*V[to[i]]*G[i];//*cos(ang[to[i]]-ang[from[i]])+B[i]*sin(ang[to[i]]-ang[from[i]]);
				//if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncQ[to[i]];
					Jj[nnzJ]=NodetoFuncP[from[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				//}
			}
		}
	}
	for(int i=0;i<n;i++){//H对角线元素
		nnzJ++;
		Ji[nnzJ]=i;
		Jj[nnzJ]=i;
		J[nnzJ]=V[FunctoNode[i]]*V[FunctoNode[i]]*(B[FunctoNode[i]+lineN]+Q[i]/(V[FunctoNode[i]]*V[FunctoNode[i]]));
		//if(nnzJ==985)
		//			printf("//");
	}
	for(int i=0;i<n-r;i++){//L对角线元素
		nnzJ++;
		Ji[nnzJ]=i+n;
		Jj[nnzJ]=i+n;
		J[nnzJ]=V[FunctoNode[i+n]]*V[FunctoNode[i+n]]*(B[FunctoNode[i+n]+lineN]-Q[NodetoFuncP[FunctoNode[i+n]]]/(V[FunctoNode[i+n]]*V[FunctoNode[i+n]]));
		//if(nnzJ==985)
		//			printf("//");
	}
	for(int i=0;i<n-r;i++){//N和M对角线元素
		//if(type[FunctoNode[i]]==3){
		nnzJ++;
		Ji[nnzJ]=NodetoFuncP[FunctoNode[i+n]];
		Jj[nnzJ]=i+n;
		J[nnzJ]=-V[FunctoNode[i+n]]*V[FunctoNode[i+n]]*(G[FunctoNode[i+n]+lineN]-P[NodetoFuncP[FunctoNode[i+n]]]/(V[FunctoNode[i+n]]*V[FunctoNode[i+n]]));
		//if(nnzJ==985)
		//			printf("//");
		nnzJ++;

		Ji[nnzJ]=i+n;
		Jj[nnzJ]=NodetoFuncP[FunctoNode[i+n]];
		J[nnzJ]=V[FunctoNode[i+n]]*V[FunctoNode[i+n]]*(G[FunctoNode[i+n]+lineN]-P[NodetoFuncP[FunctoNode[i+n]]]/(V[FunctoNode[i+n]]*V[FunctoNode[i+n]]));

		//if(nnzJ==985)
		//			printf("//");
	}
	//for(int i=0;i<n+1+lineN;i++)
	//	printf("%d %f %f\n",i,G[i],B[i]);
	//for(int i=0;i<nnzJ;i++)
	//	printf("%d %d %f\n",Ji[i],Jj[i],J[i]);
	return nnzJ+1;
}

double *calc_fx(_handle_t solver,double *V,double *angle,int *type,int *FunctoNode,int *NodetoFuncP,int *NodetoFuncQ,double *G,double *B,double *Pg,double *Pl,double *Qg,double *Ql,int lineN,int nfunc,int *from,int *to){
	int n=N-1;
	double *cntI=(double *)malloc(nfunc*sizeof(double));
	double *pf=(double *)malloc(nfunc*sizeof(double));
	double *fx=(double *)malloc(nfunc*sizeof(double));
	double deltaP,deltaQ;
	double count;
	for(int i=0;i<nfunc;i++){
		cntI[i]=0;
	}
	LARGE_INTEGER tstart,tend,tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&tstart);
	for(int i=0;i<lineN+N;i++){
		if(i<lineN){
			if(type[from[i]]==2){
				deltaP=V[to[i]]*(G[i]*cos(angle[from[i]]-angle[to[i]])+B[i]*sin(angle[from[i]]-angle[to[i]]));
				cntI[NodetoFuncP[from[i]]]+=deltaP;
			}
			if(type[from[i]]==3){
				deltaP=V[to[i]]*(G[i]*cos(angle[from[i]]-angle[to[i]])+B[i]*sin(angle[from[i]]-angle[to[i]]));
				cntI[NodetoFuncP[from[i]]]+=deltaP;
				deltaQ=V[to[i]]*(G[i]*sin(angle[from[i]]-angle[to[i]])-B[i]*cos(angle[from[i]]-angle[to[i]]));
				cntI[NodetoFuncQ[from[i]]]+=deltaQ;
			}
			if(type[to[i]]==2){
				deltaP=V[from[i]]*(G[i]*cos(angle[to[i]]-angle[from[i]])+B[i]*sin(angle[to[i]]-angle[from[i]]));
				cntI[NodetoFuncP[to[i]]]+=deltaP;
			}
			if(type[to[i]]==3){
				deltaP=V[from[i]]*(G[i]*cos(angle[to[i]]-angle[from[i]])+B[i]*sin(angle[to[i]]-angle[from[i]]));
				cntI[NodetoFuncP[to[i]]]+=deltaP;
				deltaQ=V[from[i]]*(G[i]*sin(angle[to[i]]-angle[from[i]])-B[i]*cos(angle[to[i]]-angle[from[i]]));
				cntI[NodetoFuncQ[to[i]]]+=deltaQ;
			}
		}
		else{
			int j=i-lineN;
			if(type[j]==2){
				cntI[NodetoFuncP[j]]+=V[j]*G[i];
			}
			if(type[j]==3){
				cntI[NodetoFuncP[j]]+=V[j]*G[i];
				cntI[NodetoFuncQ[j]]+=-V[j]*B[i];
			}
		}

	}
	for(int i=0;i<nfunc;i++){
		int node=FunctoNode[i];//node为实际节点号
		if(i<n)
			pf[i]=-(Pg[node]-Pl[node]-V[node]*cntI[i]);
		else 
			pf[i]=-(Qg[node]-Ql[node]-V[node]*cntI[i]);
		 	//printf("pf[%d]=%f\n",i,pf[i]);
	}
	QueryPerformanceCounter(&tend);
	tpf0+=(tend.QuadPart - tstart.QuadPart)*1000.0/tc.QuadPart;

	//QueryPerformanceCounter(&tstart);
	NicsLU_Solve(solver, pf, fx);
	//QueryPerformanceCounter(&tend);
	//tsolve0+=(tend.QuadPart - tstart.QuadPart)*1000.0/tc.QuadPart;
	//for(int i=0;i<nfunc;i++){
	//	count=0;
	//	for(int j=0;j<nfunc;j++){
	//		count+=-inv_J[i*nfunc+j]*pf[j];
	//	}
	//	fx[i]=count;
	//}
	//for(int i=0;i<nfunc;i++)
	//		printf("fx[%d]=%f\n",i,fx[i]);
	free(cntI);
	free(pf);
	return fx;
}
void sort(int *col_idx, double *a, int start, int end)
{
	int i, j, it;
	double dt;

	for (i=end-1; i>start; i--)
		for(j=start; j<i; j++)
			if (col_idx[j] > col_idx[j+1]){

				if (a){
					dt=a[j]; 
					a[j]=a[j+1]; 
					a[j+1]=dt;
				}
				it=col_idx[j]; 
				col_idx[j]=col_idx[j+1]; 
				col_idx[j+1]=it;

			}
}
void coo2csr(int n, int nz, double *a, int *i_idx, int *j_idx,
			 double *csr_a, int *col_idx, int *row_start)
{
	int i, l;

	for (i=0; i<=n; i++) row_start[i] = 0;

	/* determine row lengths */
	for (i=0; i<nz; i++) row_start[i_idx[i]+1]++;


	for (i=0; i<n; i++) row_start[i+1] += row_start[i];


	/* go through the structure  once more. Fill in output matrix. */
	for (l=0; l<nz; l++){
		i = row_start[i_idx[l]];
		csr_a[i] = a[l];
		col_idx[i] = j_idx[l];
		row_start[i_idx[l]]++;
	}
	/* shift back row_start */
	for (i=n; i>0; i--) row_start[i] = row_start[i-1];

	row_start[0] = 0;

	for (i=0; i<n; i++){
		sort (col_idx, csr_a, row_start[i], row_start[i+1]);
	}
}
int mergeY(int *from,int *to,double *G,double *B,int lineN,int N){
	int i=0;
	while (i<lineN){
		for(int j=0;j<i;j++){
			if(((from[i]==from[j])&&(to[i]==to[j]))||((from[i]==to[j])&&(to[i]==from[j]))){
				G[j]+=G[i];
				B[j]+=B[i];
	//			printf("%d %d ",i,j);
				for(int k=i;k<lineN-1;k++){
					from[k]=from[k+1];
					to[k]=to[k+1];
					G[k]=G[k+1];
					B[k]=B[k+1];
				}
				for(int k=lineN-1;k<lineN+N-1;k++){
					G[k]=G[k+1];
					B[k]=B[k+1];
				}
				lineN--;
				i--;
			}
		}
		i++;
	}
	return lineN;
}
int main()
{
	//struct busstation
	//{
	//	double V,ang,Pg,Qg,Pl,Ql;
	//	int type;
	//}bus[N];
	int n=N-1;
	int bushu,total=0;

	
	//double deltat=0;
	//int minit=10000;
	//double mindeltat=0;
	int iteration=1;
	for(int ite=0;ite<iteration;ite++){
		//printf("第%d次:\n",ite);
		double t=0,tmax=50;
		double deltat=0.2;
		total++;
		//int k;
		//	double R[N*N],X[N*N],C[N*N]={0},tr[N*N],shift[N*N];
		double *R = (double*)malloc(5*N*sizeof(double));
		double *X = (double*)malloc(5*N*sizeof(double));
		double *C = (double*)malloc(5*N*sizeof(double));
		double *tr = (double*)malloc(5*N*sizeof(double));
		double *shift = (double*)malloc(5*N*sizeof(double));
		//	int from[N*N],to[N*N];
		int *from = (int*)malloc(5*N*sizeof(int));
		int *to = (int*)malloc(5*N*sizeof(int));
		//double inix[2*N];
		double *V = (double*)malloc(N*sizeof(double));
		double *ang = (double*)malloc(N*sizeof(double));
		double *Pg = (double*)malloc(N*sizeof(double));
		double *Qg = (double*)malloc(N*sizeof(double));
		double *Pl = (double*)malloc(N*sizeof(double));
		double *Ql = (double*)malloc(N*sizeof(double));
		double *GG = (double*)malloc(N*sizeof(double));
		double *BB = (double*)malloc(N*sizeof(double));
		int *type = (int*)malloc(N*sizeof(int));
		//double V[N],ang[N],Pg[N],Qg[N],Pl[N],Ql[N],GG[N],BB[N];
		//int type[N];
		long int *node = (long int*)malloc(N*sizeof(long int));
		//	int from[N*N],to[N*N];

		//double inix[2*N];
		int *FunctoNode = (int*)malloc(2*N*sizeof(int));
		int *NodetoFuncP = (int*)malloc(N*sizeof(int));
		int *NodetoFuncQ = (int*)malloc(N*sizeof(int));
		double cstV,cstth;
		//   for(long int i=0;i<N*N;i++){ 
		//       R[i]=1.0e308;
		//	X[i]=1.0e308;
		//	tr[i]=1;
		//}

		FILE *fp;
		if((fp=fopen("net_13659ill2.txt","rt+"))==NULL){
			printf("Cannot open file strike any key exit!");
			getch();
			exit(1);
		}
		int nPV=0,nPQ=0;

		for (int i=0;i<N;i++){
			//fscanf(fp,"%d ",&node);
			fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %d\n",&node[i],&V[i],&ang[i],&Pg[i],&Qg[i],&Pl[i],&Ql[i],&GG[i],&BB[i],&type[i]);
			ang[i]=ang[i]*M_PI/180;	 
			if(type[i]==2){ //PV节点
				//inix[nPQ+nPV]=ang[node-1];
				ang[i]=0;
				FunctoNode[nPQ+nPV]=i;
				NodetoFuncP[i]=nPQ+nPV;
				nPV++;
			}
			if(type[i]==3){ //PQ节点
				//inix[nPQ+N-1]=V[node-1];
				ang[i]=0;
				V[i]=1;
				FunctoNode[nPQ+N-1]=i;
				//inix[nPQ+nPV]=ang[node-1];
				FunctoNode[nPQ+nPV]=i;
				NodetoFuncP[i]=nPQ+nPV;
				NodetoFuncQ[i]=nPQ+N-1;
				nPQ++;
			}
			//if(type[node-1]==1){ //参考节点
			// cstV=V[node-1];
			// cstth=ang[node-1];
			//}
		}
		//printf("%d ",nPV);
		int lineN=0;
		long int fromNode,toNode;
		while (!feof(fp)){
			//fscanf(fp,"%d %d ",&from,&to);
			fscanf(fp,"%d %d %lf %lf %lf %lf %lf\n",&fromNode,&toNode,&R[lineN],&X[lineN],&C[lineN],&tr[lineN],&shift[lineN]);
			for(int i=0;i<N;i++){
				if (node[i]==fromNode) from[lineN]=i;
				if (node[i]==toNode) to[lineN]=i;
			}
			lineN++;
		}	
		fclose(fp);
		double deltaP,deltaQ;
		int nfunc=2*n-nPV;


		double *G=(double *)malloc((lineN+N)*sizeof(double));
		double *B=(double *)malloc((lineN+N)*sizeof(double));
		
		double *fx1=(double *)malloc(nfunc*sizeof(double));
		double *fx2=(double *)malloc(nfunc*sizeof(double));
		//double *fx4=(double *)malloc(nfunc*sizeof(double));
		
		//for(int j=0;j<100;j++){
		LARGE_INTEGER t1,t2,tc;
		QueryPerformanceFrequency(&tc);
		QueryPerformanceCounter(&t1);
		for(int i=0;i<lineN+N;i++){
			double cntg,cntb;
			if(i<lineN){
				G[i]=-(R[i]/(R[i]*R[i]+X[i]*X[i]))/tr[i];
				B[i]=-(-X[i]/(R[i]*R[i]+X[i]*X[i]))/tr[i];
			}
			else
				if(i<lineN+N){
					int j=i-lineN;
					cntg=0;cntb=0;
					for(int k=0;k<lineN;k++){
						if(from[k]==j){
							cntg=cntg+(R[k]/(R[k]*R[k]+X[k]*X[k]))/(tr[k]*tr[k]);
							cntb=cntb+(-X[k]/(R[k]*R[k]+X[k]*X[k])+0.5*C[k])/(tr[k]*tr[k]);
						}
						if(to[k]==j){
							cntg=cntg+R[k]/(R[k]*R[k]+X[k]*X[k]);
							cntb=cntb-X[k]/(R[k]*R[k]+X[k]*X[k])+0.5*C[k];
						}
					}
					G[i]=cntg+GG[j];
					B[i]=cntb+BB[j];
				}
		}   
		QueryPerformanceCounter(&t2);
		double tY0=(t2.QuadPart - t1.QuadPart)*1000.0/tc.QuadPart;
		lineN=mergeY(from,to,G,B,lineN,N);
		//printf("%d\n",lineN);
		//for(int i=0;i<lineN+N;i++)
		//	if(i<lineN)
		//		printf("%d %d %f %f\n",from[i],to[i],G[i],B[i]);
		//	else
		//		printf("%d %d %f %f\n",i-lineN,i-lineN,G[i],B[i]);
        double *J=(double*)malloc((N*N)*sizeof(double));
		int *Ji=(int*)malloc((N*N)*sizeof(int));
		int *Jj=(int*)malloc((N*N)*sizeof(int));

		double *Ptot=(double*)malloc(N*sizeof(double));
		double *Qtot=(double*)malloc(N*sizeof(double));

		
		double *cntI=(double *)malloc(2*n*sizeof(double));
		LARGE_INTEGER t1cnt,t2cnt;
		QueryPerformanceCounter(&t1);
		for(int i=0;i<lineN+N;i++){
			if(i<lineN){
				if(type[from[i]]!=1){
					deltaP=V[to[i]]*(G[i]*cos(ang[from[i]]-ang[to[i]])+B[i]*sin(ang[from[i]]-ang[to[i]]));
					cntI[NodetoFuncP[from[i]]]+=deltaP;
					deltaQ=V[to[i]]*(G[i]*sin(ang[from[i]]-ang[to[i]])-B[i]*cos(ang[from[i]]-ang[to[i]]));
					cntI[NodetoFuncP[from[i]]+n]+=deltaQ;
				}
				if(type[to[i]]!=1){
					deltaP=V[from[i]]*(G[i]*cos(ang[to[i]]-ang[from[i]])+B[i]*sin(ang[to[i]]-ang[from[i]]));
					cntI[NodetoFuncP[to[i]]]+=deltaP;
					deltaQ=V[from[i]]*(G[i]*sin(ang[to[i]]-ang[from[i]])-B[i]*cos(ang[to[i]]-ang[from[i]]));
					cntI[NodetoFuncP[to[i]]+n]+=deltaQ;
				}
			}
		else{
			int j=i-lineN;
			if(type[j]!=1){
				cntI[NodetoFuncP[j]]+=V[j]*G[i];
				cntI[NodetoFuncP[j]+n]+=-V[j]*B[i];
			}
		}
	}

		for(int i=0;i<2*(N-1);i++){
			if(i<N-1)
				Ptot[i]=V[FunctoNode[i]]*cntI[i];
			else 
				Qtot[i-N+1]=V[FunctoNode[i-N+1]]*cntI[i];
		}
		//for(int i=0;i<N-1;i++)
		//	printf("%d %f %f\n",i,Ptot[i],Qtot[i]);
		int nnzJ=formJ(Ji,Jj,J,from,to,G,B,V,ang,Ptot,Qtot,N-1,nPV,lineN,NodetoFuncP,NodetoFuncQ,FunctoNode,type);
		//for(int i=0;i<nnzJ;i++)
		//	if((Ji[i]==70) && (Jj[i]==156))
		//		printf("%d %d %f\n",Ji[i],Jj[i],J[i]);
		//for(int j=i+1;j<nnzJ;j++)
		//      if((Ji[i]==Ji[j])&&(Jj[i]==Jj[j]))
		//       printf("%d %d %f\n",Ji[i],Jj[i],J[i]);

		double* val = (double *)malloc(nnzJ*sizeof(double));

		int* colindJ = (int *)malloc(nnzJ*sizeof(int));
		int* rowptrJ = (int *)malloc((nfunc+1)*sizeof(int));
		coo2csr(nfunc,nnzJ,J,Ji,Jj,val,colindJ,rowptrJ);
		//QueryPerformanceCounter(&t2);
		//double tformJ00=(t2.QuadPart - t1.QuadPart)*1000.0/tc.QuadPart;

		//int iteration=2;
		//for(int ite=0;ite<iteration;ite++){
		_handle_t solver = NULL;
		_uint_t i;
		_double_t *cfg;
		const _double_t *stat;


		if (__FAIL(NicsLU_Initialize(&solver, &cfg, &stat)))
		{
			printf("Failed to initialize\n");
			system("pause");
			return -1;
		}
		//    printf("Version %.0lf\nLicense to %.0lf\n", stat[31], stat[29]);
		//for(int i=0;i<52;i++)
		//   printf("%f ",cfg[i]);
		unsigned int* rowptr = (unsigned int *)malloc((nfunc+1)*sizeof(unsigned int));
		unsigned int* colind = (unsigned int *)malloc(nnzJ*sizeof(unsigned int));

		for(int i=0;i<nnzJ;i++)
			colind[i]=(unsigned int)colindJ[i];
		for(int i=0;i<nfunc+1;i++)
			rowptr[i]=(unsigned int)rowptrJ[i];
		//for(int i=0;i<nnzJ;i++)
		//	printf("%f ",val[i]);
		//printf("\n");
		//for(int i=0;i<nnzJ;i++)
		//	printf("%d ",colind[i]);
		//printf("\n");
		//for(int i=0;i<nfunc+1;i++)
		//	printf("%d ",rowptr[i]);
		//QueryPerformanceCounter(&t1);
		int cao=NicsLU_Analyze(solver, nfunc, val, colind, rowptr, MATRIX_ROW_REAL, NULL, NULL);
		int cao2=NicsLU_Factorize(solver, val, 1);
		//QueryPerformanceCounter(&t2);
		//double tLU0=(t2.QuadPart - t1.QuadPart)*1000.0/tc.QuadPart;
		//printf("%d %d",cao,cao2);
		//int lnz=(int) stat[9],unz=(int) stat[10];
		//double * lx = (double *)malloc((lnz+10)*sizeof(double));
		//double * ux = (double *)malloc((unz+10)*sizeof(double));
		////	_double_t sci[30]={0},scj[30]={0};
		//unsigned int * li = (unsigned int *)malloc((lnz+10)*sizeof(unsigned int));
		//unsigned int * ui = (unsigned int *)malloc((unz+10)*sizeof(unsigned int));
		//unsigned int * pi = (unsigned int *)malloc(nfunc*sizeof(unsigned int));
		//unsigned int * pj = (unsigned int *)malloc(nfunc*sizeof(unsigned int));
		//size_t * lp = (size_t *)malloc((nfunc+10)*sizeof(size_t));
		//size_t * up = (size_t *)malloc((nfunc+10)*sizeof(size_t));
		//NicsLU_GetFactors(solver,lx,li,lp,ux,ui,up,1,pi,pj,NULL,NULL);


		bushu=0;

		//QueryPerformanceCounter(&t2);
		//double fx1[500],fx2[500],fx3[500],fx4[500];
		double anglelast=ang[FunctoNode[0]];
		tpf0=0;tsolve0=0;
		while(t<tmax){
			bushu++;
			fx1=calc_fx(solver,V,ang,type,FunctoNode,NodetoFuncP,NodetoFuncQ,G,B,Pg,Pl,Qg,Ql,lineN,nfunc,from,to);
			//for(int i=0;i<nfunc;i++)
			//	printf("fx[%d]=%f\n",i,fx1[i]);
			for(int i=0;i<nfunc;i++){
				if(i<N-1) 
					ang[FunctoNode[i]]+=deltat*fx1[i];
				else 
					V[FunctoNode[i]]+=deltat*fx1[i]*V[FunctoNode[i]];
				//if(i==0)
				//	printf("第%d步:angle[%d]=%f\n",bushu,i,ang[FunctoNode[i]]);
			}
			fx2=calc_fx(solver,V,ang,type,FunctoNode,NodetoFuncP,NodetoFuncQ,G,B,Pg,Pl,Qg,Ql,lineN,nfunc,from,to);
			for(int i=0;i<nfunc;i++){
				if(i<N-1) 
					ang[FunctoNode[i]]=ang[FunctoNode[i]]-deltat*fx1[i]+0.5*deltat*(fx2[i]+fx1[i]);
				else if(i<nfunc)
					V[FunctoNode[i]]=V[FunctoNode[i]]-deltat*fx1[i]+(V[FunctoNode[i]]-deltat*fx1[i])*0.5*deltat*(fx2[i]+fx1[i]);
				if(i==0)
					printf("第%d步:angle[%d]=%f\n",bushu,i,ang[FunctoNode[i]]);
			}
			t=t+deltat;
			if((abs(ang[FunctoNode[0]]-anglelast)<0.0001)) break;
			anglelast=ang[FunctoNode[0]];
		}
		QueryPerformanceCounter(&t2);
		double tcnt0=(t2.QuadPart - t1.QuadPart)*1000.0/tc.QuadPart;
		//printf("Use Time:%f\n",tcnt0);
		/*double *V_s = (double*)malloc(N*sizeof(double));
		double *ang_s = (double*)malloc(N*sizeof(double));*/
		double V_s,ang_s;
		FILE *sta;
		sta=fopen("re13659_2.txt","rt+");
		double maxdV=0,maxdang=0;
		for(int i=0;i<N;i++){
			fscanf(sta,"%lf",&V_s);
			if(abs(V_s-V[i])>maxdV)
				maxdV=abs(V_s-V[i]);
		}
		int inde=0;
		for(int i=0;i<N;i++){
			fscanf(sta,"%lf",&ang_s);
			double pian=0;
			if(abs(ang_s+pian-ang[i])>maxdang){//标幺值为pi/6
				maxdang=abs(ang_s+pian-ang[i]);
				if (i==38) printf("%f %f\n",ang_s,ang[i]);
				inde=i;
			}
		}
		printf("%d %f %f\n",inde,maxdV,maxdang);
		NicsLU_Free(solver);
		free(G);
		free(B);
		free(fx1);
		free(fx2);
		//}
		free(J);
		free(Ji);
		free(Jj);
		free(Ptot);
		free(Qtot);
		free(val);
		free(colind);
		free(rowptr);
		free(colindJ);
		free(rowptrJ);
		free(R);
		free(X);
		free(C);
		free(tr);
		free(shift);
		free(from);
		free(to);
		free(V);
		free(ang);
		free(Pg);
		free(Qg);
		free(Pl);
		free(Ql);
		free(GG);
		free(BB);
		free(type);
		free(node);
		free(FunctoNode);
		free(NodetoFuncP);
		free(NodetoFuncQ);
		//tpf+=tpf0;
		//double inix[2*N];
		//if(bushu==17){
		////	//tY+=tY0;
		////	//tformJ0+=tformJ00;
		////	//tLU+=tLU0;
		//tpf+=tpf0;
		////	//tsolve+=tsolve0;
		////	tcnt+=tcnt0;
		//}
		//else
		//	ite--;
	}
	//int iteration=1;
	printf("iteration times:%d\n",bushu);
	printf("formY time: %f\n",tY/iteration);
	printf("formJ0 time: %f\n",tformJ0/iteration);
	printf("LU time: %f\n",tLU/iteration);
	printf("calculate pf time:%f ms\n",tpf/iteration);
	printf("solve linear eqa time:%f ms\n",tsolve/iteration);
	printf("total time: %f ms\n",tcnt/iteration);
	//getchar();
	return 0;
}

