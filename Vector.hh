struct vector_double
{
	int dim;
	double* coef;
};

typedef struct vector_double* vvector;

vvector vector_create(int n) 
{
	vvector r=(vvector)malloc(sizeof(struct vector_double));
	r->dim=n;
	r->coef=(double*)malloc(n*sizeof(double));
	return r;
}

void vector_free(vvector v) 
{
	free(v->coef);
	free(v);
}

vvector vector_null(int n) 
{
	int i;
	vvector r=vector_create(n);
	for(i=0;i<r->dim;i++) 
	{
	  r->coef[i]=0.;
	}
	return r;
}

void vector_print(vvector v) 
{
	int i;
	printf("[");
	for(i=0;i<v->dim;i++) 
	{
		printf("%f ",v->coef[i]);
	}
	printf("]");
}

double vector_norm(vvector v) 
{
	int i;
	double x=0.;
	for(i=0;i<v->dim;i++) 
	{
		x+=v->coef[i]*v->coef[i];
	}
	return sqrt(x);
}