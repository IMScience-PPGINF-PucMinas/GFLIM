iftLabeledSet* __add__(iftLabeledSet* lset){
    iftLabeledSet *copy = iftCopyOrderedLabeledSet(($self));
    iftConcatLabeledSet(&copy,&lset);
    return copy;
}

iftLabeledSet* __sub__(iftLabeledSet* lset){
    iftLabeledSet *copy = iftCopyOrderedLabeledSet(($self));
    iftRemoveSubsetLabeledSet(&copy,&lset);
    return copy;
}

PyObject* AsDict(){
    iftLabeledSet* s = ($self);

    PyObject *dict = PyDict_New();
    while(s != NULL){
        PyDict_SetItem(dict, PyInt_FromLong(s->elem), PyInt_FromLong(s->label));
        s = s->next;
    }

    return dict;
}

void Print(){
    iftLabeledSet *ptr = ($self);

    while(ptr != NULL){
        printf("%d, %d\n",ptr->elem, ptr->label);
        ptr = ptr->next;
    }
}