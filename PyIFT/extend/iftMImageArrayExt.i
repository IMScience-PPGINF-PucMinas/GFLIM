iftMImage * __getitem__(int index)
{
    iftMImageArray *arr = ($self);

    if (index < 0 || index >= arr->n) {
        SWIG_Error(0, "Index out of bonds.");
        return NULL;
    }

    return iftCopyMImage(arr->val[index]);
}