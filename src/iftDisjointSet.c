#include "iftDisjointSet.h"

#include "ift/core/io/Stream.h"


iftDisjointSet* iftCreateDisjointSet(int n){
	iftDisjointSet *dset = (iftDisjointSet*)iftAlloc(1,sizeof(iftDisjointSet));

	dset->n = n;

	dset->parent = iftAllocIntArray(dset->n);
	dset->rank = iftAllocIntArray(dset->n);

	int i;
	for(i = 0; i < n; i++)
		dset->parent[i] = i;

	return dset;
}

void iftDestroyDisjointSet(iftDisjointSet **dset)
{
	if(dset != NULL && *dset != NULL)
	{
		iftFree((*dset)->parent);
		iftFree((*dset)->rank);

		iftFree(*dset);

		*dset = NULL;
	}
}

int iftDisjointSetFind(iftDisjointSet* dset, int elem){
	if(elem < 0 || elem >= dset->n)
        iftError("Invalid element", "iftDisjointSetFind");

	if(dset->parent[elem] != elem)
		dset->parent[elem] = iftDisjointSetFind(dset, dset->parent[elem]);

	return dset->parent[elem];
}

int iftDisjointSetUnion(iftDisjointSet* dset, int elem1, int elem2){
	int root1 = iftDisjointSetFind(dset, elem1);
	int root2 = iftDisjointSetFind(dset, elem2);

	if (root1 == root2)
		return root1;

	if( dset->rank[root1] < dset->rank[root2]){
		dset->parent[root1] = root2;
		return root2;
	}

	if (dset->rank[root1] > dset->rank[root2]){
		dset->parent[root2] = root1;
		return root1;
	}

	dset->parent[root2] = root1;
	dset->rank[root1]++;

	return root1;
}
