#include "iftSort.h"

#include "ift/core/dtypes/FHeap.h"
#include "ift/core/dtypes/GQueue.h"
#include "ift/core/io/Stream.h"

/**
@file
@brief A description is missing here
*/
void iftBucketSort(int *value, int *index, int nelems, uchar order)
{
  int i,j,maxval= IFT_INFINITY_INT_NEG, *sval=NULL, *sind=NULL;
  iftGQueue *Q=NULL;

  for(i=0; i < nelems; i++) 
    if (value[i] > maxval)
      maxval = value[i]; 

 
  Q    = iftCreateGQueue(maxval+1,nelems,value);
  sval = iftAllocIntArray(nelems);
  sind = iftAllocIntArray(nelems);

  if (order == IFT_DECREASING){
    iftSetRemovalPolicy(Q,IFT_MAXVALUE);
  }

  for(i=0; i < nelems; i++)
    iftInsertGQueue(&Q,i);

  j = 0;
  while(!iftEmptyGQueue(Q)) {
    i = iftRemoveGQueue(Q);
    sval[j]  = value[i];
    sind[j]  = index[i];
    j++;
  }

  iftDestroyGQueue(&Q);

  for(i=0; i < nelems; i++){
    value[i]=sval[i];
    index[i]=sind[i];
  }
  iftFree(sval);
  iftFree(sind);
}

void iftFHeapSort(float *value, int *index, int nelems, uchar order)
{
  int i,j, *sind=NULL;
  float *sval=NULL;
  iftFHeap *H=NULL;

  H    = iftCreateFHeap(nelems,value);
  sval = iftAllocFloatArray(nelems);
  sind = iftAllocIntArray(nelems);

  if (order == IFT_DECREASING){
    iftSetRemovalPolicyFHeap(H, IFT_MAXVALUE);
  }

  for(i=0; i < nelems; i++)
    iftInsertFHeap(H, i);

  j = 0;
  while(!iftEmptyFHeap(H)) {
    i = iftRemoveFHeap(H);
    sval[j]  = value[i];
    sind[j]  = index[i];
    j++;
  }

  iftDestroyFHeap(&H);

  for(i=0; i < nelems; i++){
    value[i]=sval[i];
    index[i]=sind[i];
  }
  iftFree(sval);
  iftFree(sind);
}

void iftQuickSort( int *value, int *index, int i0, int i1, uchar order ) 
{
  int m, d;
  

  if( i0 < i1 ) {
    /* random index to avoid bad pivots.*/
    // d = iftRandomInteger( i0, i1 );
    // to guarantee the same behavior on ties 
    d = (i0 + i1) / 2;
    iftSwap( value[ d ], value[ i0 ] );
    iftSwap( index[ d ], index[ i0 ] );
    m = i0;

    if(order == IFT_INCREASING ) {
      for( d = i0 + 1; d <= i1; d++ ) {
	if( value[ d ] < value[ i0 ] ) {
	  m++;
	  iftSwap( value[ d ], value[ m ] );
	  iftSwap( index[ d ], index[ m ] );
	}
      }
    }
    else {
      for( d = i0 + 1; d <= i1; d++ ) {
	if( value[ d ] > value[ i0 ] ) {
	  m++;
	  iftSwap( value[ d ], value[ m ] );
	  iftSwap( index[ d ], index[ m ] );
	}
      }
    }
    iftSwap( value[ m ], value[ i0 ] );
    iftSwap( index[ m ], index[ i0 ] );

    iftQuickSort( value, index, i0, m - 1, order );
    iftQuickSort( value, index, m + 1, i1, order );
  }
}

void iftSQuickSort( char **value, int *index, int i0, int i1, uchar order, int size ) 
{
  int m, d;
  

  if( i0 < i1 ) {
    /* random index to avoid bad pivots.*/
    // d = iftRandomInteger( i0, i1 );
    // to guarantee the same behavior on ties 
    d = (i0 + i1) / 2;
      iftSwapString(value[d], value[i0]);
    iftSwap( index[ d ], index[ i0 ] );
    m = i0;

    if(order == IFT_INCREASING ) {
      for( d = i0 + 1; d <= i1; d++ ) {
	if( strcmp( value[ d ], value[ i0 ] ) < 0 ){
	  m++;
        iftSwapString(value[d], value[m]);
	  iftSwap( index[ d ], index[ m ] );
	}
      }
    }
    else {
      for( d = i0 + 1; d <= i1; d++ ) {
	if( strcmp( value[ d ] , value[ i0 ] ) > 0 ) {
	  m++;
        iftSwapString(value[d], value[m]);
	  iftSwap( index[ d ], index[ m ] );
	}
      }
    }
      iftSwapString(value[m], value[i0]);
    iftSwap( index[ m ], index[ i0 ] );
    iftSQuickSort( value, index, i0, m - 1, order , size);
    iftSQuickSort( value, index, m + 1, i1, order , size);
  }
}

void iftFQuickSort( float *value, int *index, int i0, int i1, uchar order ) 
{
  int m, d;


	if( i0 < i1 ) {
		/* random index to avoid bad pivots.*/
		// d = iftRandomInteger( i0, i1 );
    // to guarantee the same behavior on ties 
    d = (i0 + i1) / 2;
		iftSwap( value[ d ], value[ i0 ] );
		iftSwap( index[ d ], index[ i0 ] );
		m = i0;

		if(order == IFT_INCREASING ) {
			for( d = i0 + 1; d <= i1; d++ ) {
				if( value[ d ] < value[ i0 ] ) {
					m++;
					iftSwap( value[ d ], value[ m ] );
					iftSwap( index[ d ], index[ m ] );
				}
			}
		}
		else {
			for( d = i0 + 1; d <= i1; d++ ) {
				if( value[ d ] > value[ i0 ] ) {
					m++;
					iftSwap( value[ d ], value[ m ] );
					iftSwap( index[ d ], index[ m ] );
				}
			}
		}
		iftSwap( value[ m ], value[ i0 ] );
		iftSwap( index[ m ], index[ i0 ] );
		iftFQuickSort( value, index, i0, m - 1, order );
		iftFQuickSort( value, index, m + 1, i1, order );
	}
}

void iftDQuickSort( double *value, int *index, int i0, int i1, uchar order ) 
{
  int m, d;
  

  if( i0 < i1 ) {
    /* random index to avoid bad pivots.*/
    // d = iftRandomInteger( i0, i1 );
    // to guarantee the same behavior on ties 
    d = (i0 + i1) / 2;
    iftSwap( value[ d ], value[ i0 ] );
    iftSwap( index[ d ], index[ i0 ] );
    m = i0;

    if(order == IFT_INCREASING ) {
      for( d = i0 + 1; d <= i1; d++ ) {
	if( value[ d ] < value[ i0 ] ) {
	  m++;
	  iftSwap( value[ d ], value[ m ] );
	  iftSwap( index[ d ], index[ m ] );
	}
      }
    }
    else {
      for( d = i0 + 1; d <= i1; d++ ) {
	if( value[ d ] > value[ i0 ] ) {
	  m++;
	  iftSwap( value[ d ], value[ m ] );
	  iftSwap( index[ d ], index[ m ] );
	}
      }
    }
    iftSwap( value[ m ], value[ i0 ] );
    iftSwap( index[ m ], index[ i0 ] );
    iftDQuickSort( value, index, i0, m - 1, order );
    iftDQuickSort( value, index, m + 1, i1, order );
  }
}

int iftSelectTheKthElem(int *value, int nelems, int k)
{  
  if (nelems <= k)
      iftError("k cannot be larger than the size of the array", "iftSelectTheKthElem");
 
  int first = 0, last = nelems-1;
 
  // if first == last then we reached the kth element
  while (first < last) {
    int i1 = first, i2 = last;
    int pivot  = value[(i1 + i2) / 2];
 
    while (i1 < i2) {
      if (value[i1] >= pivot) { // put the large values at the end
	int tmp = value[i2];
	value[i2] = value[i1];
	value[i1] = tmp;
	i2--;
      } else { // the value is smaller than the pivot, skip
	i1++;
      }
    }
 
    // if we stepped up (i1++) we need to step one down
    if (value[i1] > pivot)
      i1--;
 
    // the i1 pointer is at the end of the first k elements
    if (k <= i1) {
      last = i1;
    } else {
      first = i1 + 1;
    }
  }
 
  return value[k];
}


int iftSelectIndexOfTheKthElem(int *value, int *index, int nelems, int k)
{  
  if (nelems <= k)
      iftError("k cannot be larger than the size of the array", "iftSelectIndexOfTheKthElem");
 
  int first = 0, last = nelems-1;
 
  // if first == last then we reached the kth element
  while (first < last) {
    int i1 = first, i2 = last;
    int pivot  = value[(i1 + i2) / 2];
 
    while (i1 < i2) {
      if (value[i1] >= pivot) { // put the large values at the end
	int tmp = value[i2];
	value[i2] = value[i1];
	value[i1] = tmp;
	tmp = index[i2];
	index[i2] = index[i1];
	index[i1] = tmp;
	i2--;
      } else { // the value is smaller than the pivot, skip
	i1++;
      }
    }
 
    // if we stepped up (i1++) we need to step one down
    if (value[i1] > pivot)
      i1--;
 
    // the i1 pointer is at the end of the first k elements
    if (k <= i1) {
      last = i1;
    } else {
      first = i1 + 1;
    }
  }
 
  return index[k];
}
