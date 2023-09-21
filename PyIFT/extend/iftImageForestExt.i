void SetPathval(iftImage *in)
{
    iftImageForest *fst = ($self);
    iftDestroyImage(&fst->pathval);
    fst->pathval = in;
}

void SetLabel(iftImage *in)
{
    iftImageForest *fst = ($self);
    iftDestroyImage(&fst->label);
    fst->label = in;
}

void SetRoot(iftImage *in)
{
    iftImageForest *fst = ($self);
    iftDestroyImage(&fst->root);
    fst->root = in;
}

void SetMarker(iftImage *in)
{
    iftImageForest *fst = ($self);
    iftDestroyImage(&fst->marker);
    fst->marker = in;
}

void SetPred(iftImage *in)
{
    iftImageForest *fst = ($self);
    iftDestroyImage(&fst->pred);
    fst->pred = in;
}

void SetImage(iftImage *in)
{
    iftImageForest *fst = ($self);
    iftDestroyImage(&fst->img);
    fst->img = in;
}

iftImage *GetPathval(void)
{
    iftImageForest *fst = ($self);
    return fst->pathval;
}

iftImage *GetLabel(void)
{
    iftImageForest *fst = ($self);
    return fst->label;
}

iftImage *GetRoot(void)
{
    iftImageForest *fst = ($self);
    return fst->root;
}

iftImage *GetMarker(void)
{
    iftImageForest *fst = ($self);
    return fst->marker;
}

iftImage *GetPred(void)
{
    iftImageForest *fst = ($self);
    return fst->pred;
}

iftImage *GetImage(void)
{
    iftImageForest *fst = ($self);
    return fst->img;
}