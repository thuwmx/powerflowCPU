// pfCPU.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <windows.h>

#include <string.h>
#include <mkl.h>

#include <mkl_pardiso.h>
#include "mkl_types.h"

//#ifdef _DLL //��̬����
//
//#pragma comment(lib,"mkl_intel_lp64_dll.lib")
//#pragma comment(lib,"mkl_intel_ilp64_dll.lib")
//#pragma comment(lib,"mkl_intel_thread_dll.lib")
//#pragma comment(lib,"mkl_core_dll.lib")
//
//#else //��̬����
//
//#pragma comment(lib,"mkl_intel_lp64.lib")
//#pragma comment(lib,"mkl_intel_ilp64.lib")
//#pragma comment(lib,"mkl_intel_thread.lib")
//#pragma comment(lib,"mkl_core.lib")
//#pragma comment(lib,"mkl_lapack95_lp64.lib")
//#pragma comment(lib,"mkl_lapack95_ilp64.lib")
//#endif // _DLL

//#ifdef _DLL //��̬����
//
//#pragma comment(lib,"mkl_intel_lp64_dll.lib")
//#pragma comment(lib,"mkl_intel_thread_dll.lib")
//#pragma comment(lib,"mkl_core_dll.lib")
//
//#else //��̬����
//
#pragma comment(lib,"mkl_intel_ilp64.lib")
#pragma comment(lib,"mkl_intel_thread.lib")
#pragma comment(lib,"mkl_core.lib")
#pragma comment(lib,"mkl_lapack95_ilp64.lib")
#pragma comment(lib,"mkl_sequential.lib")
//
//#endif // _DLL



double tsolve=0,tpf=0,tY=0,tformJ=0;

#define M_PI 3.141592653589793
const int  N=39;
int mergeY(int *from,int *to,double *G,double *B,int lineN,int N){
	int i=0;
	while (i<lineN){
		for(int j=0;j<i;j++){
			if(((from[i]==from[j])&&(to[i]==to[j]))||((from[i]==to[j])&&(to[i]==from[j]))){
				G[j]+=G[i];
				B[j]+=B[i];
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
int formJ(int *Ji,int *Jj,double *J,int *from,int *to,double *G,double *B,double *V,double *ang,double *P,double *Q,int n,int r,int lineN,int *NodetoFuncP,int *NodetoFuncQ,int *FunctoNode,int *type){
	int nnzJ=-1;
	double value;
	for(int i=0;i<lineN;i++){
		if((type[from[i]]!=1)&&(type[to[i]]!=1)){
			//H����������Ԫ��
			value=V[from[i]]*(B[i]*cos(ang[from[i]]-ang[to[i]])-G[i]*sin(ang[from[i]]-ang[to[i]]))*V[to[i]];
			if(abs(value)>0.000000001){
				nnzJ++;
				Ji[nnzJ]=NodetoFuncP[from[i]];
				Jj[nnzJ]=NodetoFuncP[to[i]];
				J[nnzJ]=value;
				//if(nnzJ==985)
				//	printf("//");
			}
			value=V[to[i]]*(B[i]*cos(ang[to[i]]-ang[from[i]])-G[i]*sin(ang[to[i]]-ang[from[i]]))*V[from[i]];
			if(abs(value)>0.000000001){
				nnzJ++;
				Ji[nnzJ]=NodetoFuncP[to[i]];
				Jj[nnzJ]=NodetoFuncP[from[i]];
				J[nnzJ]=value;
				//if(nnzJ==985)
				//	printf("//");
			}
			//L����������Ԫ��
			if((type[from[i]]==3)&&(type[to[i]]==3)){
				value=V[from[i]]*(B[i]*cos(ang[from[i]]-ang[to[i]])-G[i]*sin(ang[from[i]]-ang[to[i]]))*V[to[i]];
				if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncQ[from[i]];
					Jj[nnzJ]=NodetoFuncQ[to[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				}
				value=V[to[i]]*(B[i]*cos(ang[to[i]]-ang[from[i]])-G[i]*sin(ang[to[i]]-ang[from[i]]))*V[from[i]];
				if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncQ[to[i]];
					Jj[nnzJ]=NodetoFuncQ[from[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				}
			}
			//N����������Ԫ��
			if(type[to[i]]==3){
				value=V[from[i]]*(-G[i]*cos(ang[from[i]]-ang[to[i]])-B[i]*sin(ang[from[i]]-ang[to[i]]))*V[to[i]];
				if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncP[from[i]];
					Jj[nnzJ]=NodetoFuncQ[to[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				}
			}
			if(type[from[i]]==3){
				value=V[to[i]]*(-G[i]*cos(ang[to[i]]-ang[from[i]])-B[i]*sin(ang[to[i]]-ang[from[i]]))*V[from[i]];
				if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncP[to[i]];
					Jj[nnzJ]=NodetoFuncQ[from[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				}
			}
			//M����������Ԫ��
			if(type[from[i]]==3){
				value=V[from[i]]*(G[i]*cos(ang[from[i]]-ang[to[i]])+B[i]*sin(ang[from[i]]-ang[to[i]]))*V[to[i]];
				if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncQ[from[i]];
					Jj[nnzJ]=NodetoFuncP[to[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				}
			}
			if(type[to[i]]==3){
				value=V[to[i]]*(G[i]*cos(ang[to[i]]-ang[from[i]])+B[i]*sin(ang[to[i]]-ang[from[i]]))*V[from[i]];
				if(abs(value)>0.000000001){
					nnzJ++;
					Ji[nnzJ]=NodetoFuncQ[to[i]];
					Jj[nnzJ]=NodetoFuncP[from[i]];
					J[nnzJ]=value;
					//if(nnzJ==985)
					//printf("//");
				}
			}
		}
	}
	for(int i=0;i<n;i++){//H�Խ���Ԫ��
		nnzJ++;
		Ji[nnzJ]=i;
		Jj[nnzJ]=i;
		J[nnzJ]=V[FunctoNode[i]]*V[FunctoNode[i]]*B[FunctoNode[i]+lineN]+Q[i];
		//if(nnzJ==985)
		//			printf("//");
	}
	for(int i=0;i<n-r;i++){//L�Խ���Ԫ��
		nnzJ++;
		Ji[nnzJ]=i+n;
		Jj[nnzJ]=i+n;
		J[nnzJ]=V[FunctoNode[i+n]]*V[FunctoNode[i+n]]*B[FunctoNode[i+n]+lineN]-Q[NodetoFuncP[FunctoNode[i+n]]];
		//if(nnzJ==985)
		//			printf("//");
	}
	for(int i=0;i<n-r;i++){//N��M�Խ���Ԫ��
		//if(type[FunctoNode[i]]==3){
		nnzJ++;
		Ji[nnzJ]=NodetoFuncP[FunctoNode[i+n]];
		Jj[nnzJ]=i+n;
		J[nnzJ]=-V[FunctoNode[i+n]]*V[FunctoNode[i+n]]*G[FunctoNode[i+n]+lineN]-P[NodetoFuncP[FunctoNode[i+n]]];
		//if(nnzJ==985)
		//			printf("//");
		nnzJ++;

		Ji[nnzJ]=i+n;
		Jj[nnzJ]=NodetoFuncP[FunctoNode[i+n]];
		J[nnzJ]=V[FunctoNode[i+n]]*V[FunctoNode[i+n]]*G[FunctoNode[i+n]+lineN]-P[NodetoFuncP[FunctoNode[i+n]]];

		//if(nnzJ==985)
		//			printf("//");
	}
	//for(int i=0;i<n+1+lineN;i++)
	//	printf("%d %f %f\n",i,G[i],B[i]);
	//for(int i=0;i<nnzJ;i++)
	//	printf("%d %d %f\n",Ji[i],Jj[i],J[i]);
	return nnzJ+1;
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
double *calc_fx(void *pt,MKL_INT *iparm,double *V,double *angle,int *type,int *FunctoNode,int *NodetoFuncP,int *NodetoFuncQ,double *G,double *B,double *Pg,double *Pl,double *Qg,double *Ql,int lineN,int nfunc,int *from,int *to,double *ang,int nPV){
	int n=N-1;
	double *cntI=(double *)malloc(2*n*sizeof(double));
	double *pf=(double *)malloc(nfunc*sizeof(double));
	double *fx=(double *)malloc(nfunc*sizeof(double));
	double deltaP,deltaQ;
	double count;
	for(int i=0;i<2*n;i++){
		cntI[i]=0;
	}
	LARGE_INTEGER tstart,tend,tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&tstart);
	for(int i=0;i<lineN+N;i++){
		if(i<lineN){
			if(type[from[i]]!=1){
				deltaP=V[to[i]]*(G[i]*cos(angle[from[i]]-angle[to[i]])+B[i]*sin(angle[from[i]]-angle[to[i]]));
				cntI[NodetoFuncP[from[i]]]+=deltaP;
				deltaQ=V[to[i]]*(G[i]*sin(angle[from[i]]-angle[to[i]])-B[i]*cos(angle[from[i]]-angle[to[i]]));
				cntI[NodetoFuncP[from[i]]+n]+=deltaQ;
			}
			if(type[to[i]]!=1){
				deltaP=V[from[i]]*(G[i]*cos(angle[to[i]]-angle[from[i]])+B[i]*sin(angle[to[i]]-angle[from[i]]));
				cntI[NodetoFuncP[to[i]]]+=deltaP;
				deltaQ=V[from[i]]*(G[i]*sin(angle[to[i]]-angle[from[i]])-B[i]*cos(angle[to[i]]-angle[from[i]]));
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
	double *Ptot=(double*)malloc(N*sizeof(double));
	double *Qtot=(double*)malloc(N*sizeof(double));
	for(int i=0;i<2*(N-1);i++){
		if(i<N-1)
			Ptot[i]=V[FunctoNode[i]]*cntI[i];
		else
			Qtot[i-N+1]=V[FunctoNode[i-N+1]]*cntI[i];
	}
	for(int i=0;i<nfunc;i++){
		int node=FunctoNode[i];//nodeΪʵ�ʽڵ��
		if(i<n)
			pf[i]=-(Pg[node]-Pl[node]-Ptot[i]);
		else 
			pf[i]=-(Qg[node]-Ql[node]-Qtot[NodetoFuncP[node]]);
		// 	printf("pf[%d]=%f\n",i,pf[i]);
	}
	QueryPerformanceCounter(&tend);
	tpf+=(tend.QuadPart - tstart.QuadPart)*1000.0/tc.QuadPart;



	double *J=(double*)malloc((N*N)*sizeof(double));
	int *Ji=(int*)malloc((N*N)*sizeof(int));
	int *Jj=(int*)malloc((N*N)*sizeof(int));

	QueryPerformanceCounter(&tstart);

	int* rowptrJ = (int *)malloc((nfunc+2)*sizeof(int));



	int nnzJ=formJ(Ji,Jj,J,from,to,G,B,V,ang,Ptot,Qtot,N-1,nPV,lineN,NodetoFuncP,NodetoFuncQ,FunctoNode,type);
	double* val = (double *)malloc((nnzJ+1)*sizeof(double));
	int* colindJ = (int *)malloc((nnzJ+1)*sizeof(int));
	coo2csr(nfunc,nnzJ,J,Ji,Jj,val,colindJ,rowptrJ);
	QueryPerformanceCounter(&tend);
	tformJ+=(tend.QuadPart - tstart.QuadPart)*1000.0/tc.QuadPart;

	MKL_INT mtype = 11; /*real and nonsymmetric*/
	/* RHS and solution vectors. */
	MKL_INT nrhs = 1; /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	
	/* Pardiso control parameters. */
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum; /* Double dummy */
	MKL_INT idum; /* Integer dummy. */
	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 1; /* Print statistical information in file */
	error = 0; /* Initialize error flag */
	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	///* -------------------------------------------------------------------- */
	for(int i=nnzJ; i>0;i--)
		val[i]=val[i-1];
	for(int i=nnzJ; i>0;i--){
		colindJ[i]=colindJ[i-1];
		colindJ[i]++;
	}
	for(int i=nfunc+1; i>0;i--){
     	rowptrJ[i]=rowptrJ[i-1];
		rowptrJ[i]++;
	}
	for(int i=1;i<nnzJ+1;i++)
		printf("%f ",val[i]);
	printf("\n");
	for(int i=1;i<nnzJ+1;i++)
		printf("%d ",colindJ[i]);
	printf("\n");
	for(int i=1;i<=nfunc+1;i++)
		printf("%d ",rowptrJ[i]);
	//rowptrJ[nfunc]--;
	phase = 11;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&nfunc, val, rowptrJ, colindJ, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&nfunc, val, rowptrJ, colindJ, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
	/* Set right hand side to one. */

	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&nfunc, val, rowptrJ, colindJ, &idum, &nrhs,
		iparm, &msglvl, pf, fx, &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1; /* Release internal memory. */
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		&nfunc, &ddum, rowptrJ, colindJ, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	//for (i = 0; i < n; i++) {
	//	if(i==0)
	//	  printf("\n fx [%d] = % f", i, fx[i] );
	//}
	printf ("\n");

	//for(int i=0;i<20;i++)
	//   printf("%f ",stat[i]);


	QueryPerformanceCounter(&tstart);
	QueryPerformanceCounter(&tend);
	tsolve+=(tend.QuadPart - tstart.QuadPart)*1000.0/tc.QuadPart;
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
	free(Ptot);
	free(Qtot);
	free(J);
	free(Ji);
	free(Jj);
	free(rowptrJ);
	free(val);
	free(colindJ);

	return fx;
}

int main()
{
	//struct busstation
	//{
	//	double V,ang,Pg,Qg,Pl,Ql;
	//	int type;
	//}bus[N];

	int iteration=1;
	for(int ite=0;ite<iteration;ite++){
		int n=N-1;
		int k;
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
		if((fp=fopen("net_39.txt","rt+"))==NULL){
			printf("Cannot open file strike any key exit!");
			getch();
			exit(1);
		}
		int nPV=0,nPQ=0;

		for (int i=0;i<N;i++){
			//fscanf(fp,"%d ",&node);
			fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %d\n",&node[i],&V[i],&ang[i],&Pg[i],&Qg[i],&Pl[i],&Ql[i],&GG[i],&BB[i],&type[i]);
			ang[i]=ang[i]*M_PI/180;	 
			if(type[i]==2){ //PV�ڵ�
				//inix[nPQ+nPV]=ang[node-1];
				ang[i]=0;
				FunctoNode[nPQ+nPV]=i;
				NodetoFuncP[i]=nPQ+nPV;
				nPV++;
			}
			if(type[i]==3){ //PQ�ڵ�
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
			//if(type[node-1]==1){ //�ο��ڵ�
			// cstV=V[node-1];
			// cstth=ang[node-1];
			//}
		}
		int nfunc=2*n-nPV;
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

		double *G=(double *)malloc((lineN+N)*sizeof(double));
		double *B=(double *)malloc((lineN+N)*sizeof(double));
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
		tY+=(t2.QuadPart - t1.QuadPart)*1000.0/tc.QuadPart;
		lineN=mergeY(from,to,G,B,lineN,N);

		//if((fp=fopen("csrJ39.txt","rt+"))==NULL){
		//     printf("Cannot open file strike any key exit!");
		//     getch();
		//     exit(1);
		//   }
		int nnz;
		//fscanf(fp,"%d %d",&nfunc,&nnz);

		//for(int i=0;i<nnz;i++)
		//	fscanf(fp,"%lf",&val[i]);
		//for(int i=0;i<nnz;i++)
		//	fscanf(fp,"%d",&colind[i]);
		//for(int i=0;i<nfunc+1;i++)
		//	fscanf(fp,"%d",&rowptr[i]);
		//fclose(fp);

		//int iteration=2;
		//for(int ite=0;ite<iteration;ite++){

		//    printf("Version %.0lf\nLicense to %.0lf\n", stat[31], stat[29]);
		//for(int i=0;i<52;i++)
		//   printf("%f ",cfg[i]);
		void *pt[65];
		for (int i = 0; i <=64; i++) {
			pt[i] = 0;
		}
		MKL_INT iparm[65];
		for (int i = 0; i <=64; i++) {
			iparm[i] = 0;
		}
		iparm[1] = 0; /* No solver default */
		//iparm[1] = 2; /* Fill-in reordering from METIS */
		///* Numbers of processors, value of OMP_NUM_THREADS */
		//iparm[2] = 1;
		//iparm[3] = 0; /* No iterative-direct algorithm */
		//iparm[4] = 0; /* No user fill-in reducing permutation */
		//iparm[5] = 0; /* Write solution into x */
		//iparm[6] = 0; /* Not in use */
		//iparm[7] = 2; /* Max numbers of iterative refinement steps */
		//iparm[8] = 0; /* Not in use */
		//iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
		//iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
		//iparm[11] = 0; /* Not in use */
		//iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
		//iparm[13] = 0; /* Output: Number of perturbed pivots */
		//iparm[14] = 0; /* Not in use */
		//iparm[15] = 0; /* Not in use */
		//iparm[16] = 0; /* Not in use */
		//iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
		//iparm[18] = -1; /* Output: Mflops for LU factorization */
		//iparm[19] = 0; /* Output: Numbers of CG Iterations */
		//iparm[34] = 1;//zero-based CSR
		double t=0,tmax=50;
		double deltat=1;
		double *fx1=(double *)malloc(nfunc*sizeof(double));
		//double *fx2=(double *)malloc(nfunc*sizeof(double));
		//double *fx3=(double *)malloc(nfunc*sizeof(double));
		//double *fx4=(double *)malloc(nfunc*sizeof(double));
		//double tcnt=0;
		//for(int j=0;j<100;j++){

		//FILE *fpB;
		//fpB=fopen("39_B.txt","w+");
		//fprintf(fpB,"%d\n",lineN+N);
		//for(int i=0;i<lineN;i++){
		//	fprintf(fpB,"%d ",from[i]);
		//	fprintf(fpB,"%d ",to[i]);
		//	fprintf(fpB,"%f\n",B[i]);
		//}
		//for(int i=0;i<N;i++){
		//	fprintf(fpB,"%d ",i);
		//	fprintf(fpB,"%d ",i);
		//	fprintf(fpB,"%f\n",B[i+lineN]);
		//}
		//   fclose(fpB);
		int bushu=0;


		//QueryPerformanceCounter(&t2);
		//double fx1[500],fx2[500],fx3[500],fx4[500];
		double anglelast=ang[FunctoNode[0]];
		while(t<tmax){
			bushu++;
			fx1=calc_fx(pt,iparm,V,ang,type,FunctoNode,NodetoFuncP,NodetoFuncQ,G,B,Pg,Pl,Qg,Ql,lineN,nfunc,from,to,ang,nPV);
			//for(int i=0;i<nfunc;i++)
			//	printf("fx[%d]=%f\n",i,fx1[i]);
			for(int i=0;i<nfunc;i++){
				if(i<N-1) 
					ang[FunctoNode[i]]+=deltat*fx1[i];
				else 
					V[FunctoNode[i]]+=deltat*fx1[i];
				//if(i==38) printf("V=%f",V[nodeV[i]]);
				if(i==0)
					printf("��%d��:angle[%d]=%f\n",bushu,i,ang[FunctoNode[i]]);
			}

			//fx3=calc_fx(inv_J,V,ang,type,FunctoNode,NodetoFuncP,NodetoFuncQ,G,B,Pg,Pl,Qg,Ql,lineN,nfunc,from,to);
			//for(int i=0;i<nfunc;i++){
			//  if(i<N-1) 
			//       ang[FunctoNode[i]]+=deltat*fx3[i]-0.5*deltat*fx2[i];
			//     else 
			//    V[FunctoNode[i]]+=deltat*fx3[i]-0.5*deltat*fx2[i];
			//}
			//fx4=calc_fx(inv_J,V,ang,type,FunctoNode,NodetoFuncP,NodetoFuncQ,G,B,Pg,Pl,Qg,Ql,lineN,nfunc,from,to);
			//for(int i=0;i<nfunc;i++){
			//  if(i<N-1) 
			//        ang[FunctoNode[i]]=ang[FunctoNode[i]]-deltat*fx3[i]+deltat/6*(fx1[i]+2*fx2[i]+2*fx3[i]+fx4[i]);
			//     else 
			//     V[FunctoNode[i]]=V[FunctoNode[i]]-deltat*fx3[i]+deltat/6*(fx1[i]+2*fx2[i]+2*fx3[i]+fx4[i]);
			//if(i==0)
			//          printf("��%d��:angle[%d]=%f\n",bushu,i,ang[FunctoNode[i]]);
			//	}
			t=t+deltat;
			if(abs(ang[FunctoNode[0]]-anglelast)<0.0001) break;
			anglelast=ang[FunctoNode[0]];
		}
		//QueryPerformanceCounter(&t2);
		//tcnt+=(t2.QuadPart - t1.QuadPart)*1.0/tc.QuadPart*1000;
		//printf("Use Time:%f\n",(t2.QuadPart - t1.QuadPart)*1.0/tc.QuadPart*1000);

		free(G);
		free(B);
		free(fx1);
		//free(fx2);
		//}

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
	}
	//int iteration=1;
	printf("formY time: %f ms\n",tY/iteration);
	printf("calculate pf time:%f ms\n",tpf/iteration);
	printf("solve linear eqa time:%f ms\n",tsolve/iteration);

	getchar();
	return 0;
}

