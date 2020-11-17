#pragma once

template <class A>
class Fast_Sort
{
public:
	int *SORT_STACK;    // stack for Fast_Sort // default 1000 //__070321__//
	A *SORT_TEMP;
	A *SORT_TEMP2;

public:
	Fast_Sort(int len=1000);
	~Fast_Sort(void);

	void fast_sort_2(A **Array,int totnum,int main_patial,int head_patial=0,int real_patial=0);   // 2d (descending order)
	void fast_sort_2(A *Array, int totnum,int main_patial,int head_patial=0,int real_patial=0);   // 2d (descending order)
	void fast_sort_1(A *Array, int totnum,int head_patial=0);                                     // 1d (descending order)
	void fast_sort_2(A **Array,int *ix,int totnum,int main_patial,int head_patial=0,int real_patial=0);  // 2d (descending order)[index]
	void fast_sort_2(A *Array, int *ix,int totnum,int main_patial,int head_patial=0,int real_patial=0);  // 2d (descending order)[index]
	void fast_sort_1(A *Array, int *ix,int totnum,int head_patial=0);                                    // 1d (descending order)[index]

	void fast_sort_2up(A **Array,int totnum,int main_patial,int head_patial=0,int real_patial=0);   // 2d (ascending order)
	void fast_sort_2up(A *Array, int totnum,int main_patial,int head_patial=0,int real_patial=0);   // 2d (ascending order)
	void fast_sort_1up(A *Array, int totnum,int head_patial=0);                                     // 1d (ascending order)
	void fast_sort_2up(A **Array,int *ix,int totnum,int main_patial,int head_patial=0,int real_patial=0);  // 2d (ascending order)[index]
	void fast_sort_2up(A *Array, int *ix,int totnum,int main_patial,int head_patial=0,int real_patial=0);  // 2d (ascending order)[index]
	void fast_sort_1up(A *Array, int *ix,int totnum,int head_patial=0);                                    // 1d (ascending order)[index]
};




template <class A>
Fast_Sort<A>::Fast_Sort(int len)
{
	SORT_STACK=new int[len];
	SORT_TEMP=new A[len];
	SORT_TEMP2=new A[len];
}
template <class A>
Fast_Sort<A>::~Fast_Sort(void)
{
	delete [] SORT_STACK;
	delete [] SORT_TEMP;
	delete [] SORT_TEMP2;
}


template <class A>
void Fast_Sort<A>::fast_sort_2(A **Array,int totnum,int main_patial,int head_patial,int real_patial) // (descending order)
{
	int top;
	int i;
	int k;
	int j;
	int p,q;

	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			for(k=0;k<main_patial;k++) SORT_TEMP2[k]=Array[p+head_patial][k];

			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;
					if(Array[i+head_patial][real_patial]<=SORT_TEMP2[real_patial])break;
				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					if(Array[j+head_patial][real_patial]>=SORT_TEMP2[real_patial])break;
				}

				if(i<j)  // SWAP_I
				{
					for(k=0;k<main_patial;k++) SORT_TEMP[k]= Array[i+head_patial][k];
					for(k=0;k<main_patial;k++) Array[i+head_patial][k]=Array[j+head_patial][k];
					for(k=0;k<main_patial;k++) Array[j+head_patial][k]=SORT_TEMP[k];
				}
				else break;
			}

			// SWAP_II
			{
				for(k=0;k<main_patial;k++) Array[p+head_patial][k]=Array[j+head_patial][k];
				for(k=0;k<main_patial;k++) Array[j+head_patial][k]=SORT_TEMP2[k];
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}

template <class A>
void Fast_Sort<A>::fast_sort_2(A *Array,int totnum,int main_patial,int head_patial,int real_patial)  // (descending order)
{
	int top;
	int i;
	int k;
	int j;
	int p,q;

	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			for(k=0;k<main_patial;k++) SORT_TEMP2[k]=Array[main_patial*(p+head_patial)+k];

			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;
					if(Array[main_patial*(i+head_patial)+real_patial]<=SORT_TEMP2[real_patial])break;
				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					if(Array[main_patial*(j+head_patial)+real_patial]>=SORT_TEMP2[real_patial])break;
				}

				if(i<j)  // SWAP_I
				{
					for(k=0;k<main_patial;k++) SORT_TEMP[k]= Array[main_patial*(i+head_patial)+k];
					for(k=0;k<main_patial;k++) Array[main_patial*(i+head_patial)+k]=Array[main_patial*(j+head_patial)+k];
					for(k=0;k<main_patial;k++) Array[main_patial*(j+head_patial)+k]=SORT_TEMP[k];
				}
				else break;
			}

			// SWAP_II
			{
				for(k=0;k<main_patial;k++) Array[main_patial*(p+head_patial)+k]=Array[main_patial*(j+head_patial)+k];
				for(k=0;k<main_patial;k++) Array[main_patial*(j+head_patial)+k]=SORT_TEMP2[k];
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}


template <class A>
void Fast_Sort<A>::fast_sort_1(A *Array,int totnum,int head_patial)  // (descending order)
{
	int top;
	A v;
	int temp;
	int i;
	int j;
	int p,q;

	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			v=Array[p+head_patial];
			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;
					if(Array[i+head_patial]<=v)break;
				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					if(Array[j+head_patial]>=v)break;
				}

				if(i<j)  // SWAP_I
				{
					temp=Array[i+head_patial];
					Array[i+head_patial]=Array[j+head_patial];
					Array[j+head_patial]=temp;
				}
				else break;
			}

			// SWAP_II
			{
				Array[p+head_patial]=Array[j+head_patial];
				Array[j+head_patial]=v;
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}

template <class A>
void Fast_Sort<A>::fast_sort_2(A **Array,int *ix,int totnum,int main_patial,int head_patial,int real_patial)  // (descending order)[index]
{
	int top;
	A v;
	int temp;
	int i;
	int j;
	int p,q;
	int index;
	int index_v;


	//init//
	for(i=0;i<totnum;i++) ix[i]=i+head_patial; 


	//start//
	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			index_v=ix[p];
			v=Array[index_v][real_patial];


			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;			
					index=ix[i];
					if(Array[index][real_patial]<=v)break;

				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					index=ix[j];
					if(Array[index][real_patial]>=v)break;

				}

				if(i<j)  // SWAP_I
				{
					temp=ix[i];
					ix[i]=ix[j];
					ix[j]=temp;
				}
				else break;
			}

			// SWAP_II
			{
				ix[p]=ix[j];
				ix[j]=index_v;
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}

template <class A>
void Fast_Sort<A>::fast_sort_2(A *Array,int *ix,int totnum,int main_patial,int head_patial,int real_patial)  // (descending order)[index]
{
	int top;
	int i;
	A v;
	int temp;
	int j;
	int p,q;
	int index;
	int index_v;


	//init//
	for(i=0;i<totnum;i++) ix[i]=i+head_patial; 


	//start//
	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			index_v=ix[p];
			v=Array[main_patial*index_v+real_patial];

			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;				
					index=ix[i];
					if(Array[main_patial*index+real_patial]<=v)break;
				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					index=ix[j];
					if(Array[main_patial*index+real_patial]>=v)break;
				}

				if(i<j)  // SWAP_I
				{
					temp=ix[i];
					ix[i]=ix[j];
					ix[j]=temp;
				}
				else break;
			}

			// SWAP_II
			{
				ix[p]=ix[j];
				ix[j]=index_v;
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}

template <class A>
void Fast_Sort<A>::fast_sort_1(A *Array,int *ix,int totnum,int head_patial)  // (descending order)[index]
{
	int top;
	A v;
	int temp;
	int i;
	int j;
	int p,q;
	int index;
	int index_v;


	//init//
	for(i=0;i<totnum;i++) ix[i]=i+head_patial; 


	//start//
	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			index_v=ix[p];
			v=Array[index_v];

			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;				
					index=ix[i];
					if(Array[index]<=v)break;

				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					index=ix[j];
					if(Array[index]>=v)break;
				}

				if(i<j)  // SWAP_I
				{
					temp=ix[i];
					ix[i]=ix[j];
					ix[j]=temp;
				}
				else break;
			}

			// SWAP_II
			{
				ix[p]=ix[j];
				ix[j]=index_v;
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}




//__071123__//
template <class A>
void Fast_Sort<A>::fast_sort_2up(A **Array,int totnum,int main_patial,int head_patial,int real_patial)  // (ascending order)
{
	int top;
	int i;
	int k;
	int j;
	int p,q;

	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			for(k=0;k<main_patial;k++) SORT_TEMP2[k]=Array[p+head_patial][k];

			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;
					if(Array[i+head_patial][real_patial]>=SORT_TEMP2[real_patial])break;
				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					if(Array[j+head_patial][real_patial]<=SORT_TEMP2[real_patial])break;
				}

				if(i<j)  // SWAP_I
				{
					for(k=0;k<main_patial;k++) SORT_TEMP[k]= Array[i+head_patial][k];
					for(k=0;k<main_patial;k++) Array[i+head_patial][k]=Array[j+head_patial][k];
					for(k=0;k<main_patial;k++) Array[j+head_patial][k]=SORT_TEMP[k];
				}
				else break;
			}

			// SWAP_II
			{
				for(k=0;k<main_patial;k++) Array[p+head_patial][k]=Array[j+head_patial][k];
				for(k=0;k<main_patial;k++) Array[j+head_patial][k]=SORT_TEMP2[k];
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}

template <class A>
void Fast_Sort<A>::fast_sort_2up(A *Array,int totnum,int main_patial,int head_patial,int real_patial)  // (ascending order)
{
	int top;
	int i;
	int k;
	int j;
	int p,q;

	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			for(k=0;k<main_patial;k++) SORT_TEMP2[k]=Array[main_patial*(p+head_patial)+k];

			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;
					if(Array[main_patial*(i+head_patial)+real_patial]>=SORT_TEMP2[real_patial])break;
				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					if(Array[main_patial*(j+head_patial)+real_patial]<=SORT_TEMP2[real_patial])break;
				}

				if(i<j)  // SWAP_I
				{
					for(k=0;k<main_patial;k++) SORT_TEMP[k]= Array[main_patial*(i+head_patial)+k];
					for(k=0;k<main_patial;k++) Array[main_patial*(i+head_patial)+k]=Array[main_patial*(j+head_patial)+k];
					for(k=0;k<main_patial;k++) Array[main_patial*(j+head_patial)+k]=SORT_TEMP[k];
				}
				else break;
			}

			// SWAP_II
			{
				for(k=0;k<main_patial;k++) Array[main_patial*(p+head_patial)+k]=Array[main_patial*(j+head_patial)+k];
				for(k=0;k<main_patial;k++) Array[main_patial*(j+head_patial)+k]=SORT_TEMP2[k];
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}

template <class A>
void Fast_Sort<A>::fast_sort_1up(A *Array,int totnum,int head_patial)  // (ascending order)
{
	int top;
	A v;
	int temp;
	int i;
	int j;
	int p,q;

	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			v=Array[p+head_patial];
			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;
					if(Array[i+head_patial]>=v)break;
				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					if(Array[j+head_patial]<=v)break;
				}

				if(i<j)  // SWAP_I
				{
					temp=Array[i+head_patial];
					Array[i+head_patial]=Array[j+head_patial];
					Array[j+head_patial]=temp;
				}
				else break;
			}

			// SWAP_II
			{
				Array[p+head_patial]=Array[j+head_patial];
				Array[j+head_patial]=v;
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}

template <class A>
void Fast_Sort<A>::fast_sort_2up(A **Array,int *ix,int totnum,int main_patial,int head_patial,int real_patial)  // (ascending order)[index]
{
	int top;
	A v;
	int temp;
	int i;
	int j;
	int p,q;
	int index;
	int index_v;


	//init//
	for(i=0;i<totnum;i++) ix[i]=i+head_patial; 


	//start//
	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			index_v=ix[p];
			v=Array[index_v][real_patial];


			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;				
					index=ix[i];
					if(Array[index][real_patial]>=v)break;

				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					index=ix[j];
					if(Array[index][real_patial]<=v)break;

				}

				if(i<j)  // SWAP_I
				{
					temp=ix[i];
					ix[i]=ix[j];
					ix[j]=temp;
				}
				else break;
			}

			// SWAP_II
			{
				ix[p]=ix[j];
				ix[j]=index_v;
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}

template <class A>
void Fast_Sort<A>::fast_sort_2up(A *Array,int *ix,int totnum,int main_patial,int head_patial,int real_patial)  // (ascending order)[index]
{
	int top;
	int i;
	A v;
	int temp;
	int j;
	int p,q;
	int index;
	int index_v;


	//init//
	for(i=0;i<totnum;i++) ix[i]=i+head_patial; 


	//start//
	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			index_v=ix[p];
			v=Array[main_patial*index_v+real_patial];

			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;				
					index=ix[i];
					if(Array[main_patial*index+real_patial]>=v)break;
				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					index=ix[j];
					if(Array[main_patial*index+real_patial]<=v)break;
				}

				if(i<j)  // SWAP_I
				{
					temp=ix[i];
					ix[i]=ix[j];
					ix[j]=temp;
				}
				else break;
			}

			// SWAP_II
			{
				ix[p]=ix[j];
				ix[j]=index_v;
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}

template <class A>
void Fast_Sort<A>::fast_sort_1up(A *Array,int *ix,int totnum,int head_patial)  // (ascending order)[index]
{
	int top;
	A v;
	int temp;
	int i;
	int j;
	int p,q;
	int index;
	int index_v;


	//init//
	for(i=0;i<totnum;i++) ix[i]=i+head_patial; 


	//start//
	p=0;
	q=totnum-1;
	top=0;	
	for(;;)
	{
		while(p<q)
		{
			j=q+1;
			index_v=ix[p];
			v=Array[index_v];

			i=p;
			for(;;)
			{
				for(;;)
				{
					if(i>=j-1) break;			
					i++;					
					index=ix[i];
					if(Array[index]>=v)break;

				}

				for(;;)
				{
					if(j<=p) break;			
					j--;
					index=ix[j];
					if(Array[index]<=v)break;
				}

				if(i<j)  // SWAP_I
				{
					temp=ix[i];
					ix[i]=ix[j];
					ix[j]=temp;
				}
				else break;
			}

			// SWAP_II
			{
				ix[p]=ix[j];
				ix[j]=index_v;
			}
			//==PATITION==OVER==//


			if(j-p<q-j)
			{
				SORT_STACK[top+1]=j+1;
				SORT_STACK[top+2]=q;
				q=j-1;
			}
			else
			{
				SORT_STACK[top+1]=p;
				SORT_STACK[top+2]=j-1;
				p=j+1;
			}
			top=top+2;
		}
		if(top==0) return;
		q=SORT_STACK[top];
		p=SORT_STACK[top-1];
		top=top-2;
	}
}


