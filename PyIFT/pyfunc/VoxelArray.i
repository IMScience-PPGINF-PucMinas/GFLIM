//! swig(newobject)
iftVoxelArray *CreateVoxelArrayFromNumPy(PyObject *input)
{
    int is_new = 0;
    PyArrayObject *ary = obj_to_array_allow_conversion(input, NPY_LONG, &is_new);

    if (ary == NULL)
        return NULL;

    if (!require_contiguous(ary))
        SWIG_Error(0, "Input numpy array is not contiguous");

    int n_dims = array_numdims(ary);

    if (n_dims != 2) {
        char error[256];
        sprintf(error, "Voxel Numpy Array must have 2 dimensions, %d found.", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    int n_voxels = array_size(ary, 0);
    int v_dim = array_size(ary, 1);

    if (v_dim < 1 || v_dim > 3) {
        char error[256];
        sprintf(error, "Voxel Numpy Array must have length between (1, 3) on its second dimension, %d found.", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    long *ptr = array_data(ary);
    iftVoxelArray *out = iftCreateVoxelArray(n_voxels);

    for (int i = 0; i < n_voxels; i++)
    {
        int idx = i * v_dim;
        out->val[i].x = ptr[idx];
        if (v_dim > 1)
            out->val[i].y = ptr[idx + 1];
        if (v_dim > 2)
            out->val[i].z = ptr[idx + 2];
    }

    return out;
}


static int _get_direction(const iftImage *obj, const iftAdjRel *A, int src)
{
    iftVoxel u = iftGetVoxelCoord(obj, src);

    // find object pixel
    int j = 1;
    for (; j < A->n; j++) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, j);
        if (iftValidVoxel(obj, v)) {
            int q = iftGetVoxelIndex(obj, v);
            if (obj->val[q])
                break;
        }
    }

    for (int i = 0; i < A->n; i++, j = (j + 1) % A->n) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, j);
        if (iftValidVoxel(obj, v)) {
            int q = iftGetVoxelIndex(obj, v);
            if (!obj->val[q])
                break;
        } else break;
    }

    return j;
}


//! swig(newobject)
PyObject *GetPaths(const iftVoxelArray *voxel_array, const iftImage *contour, const iftImage *obj)
{
    int init_index = IFT_NIL;
    for (int i = 0; i < voxel_array->n && init_index == IFT_NIL; i++) {
        int p = iftGetVoxelIndex(contour, voxel_array->val[i]);
        if (contour->val[p])
            init_index = i;
    }

    if (init_index < 0) {
        SWIG_Error(0, "No anchor in contour found.");
        return NULL;
    }

    iftImage *dist = iftCreateImage(contour->xsize, contour->ysize, contour->zsize);
    iftLIFO *Q = iftCreateLIFO(contour->n);
    iftAdjRel *A = iftClockCircular(sqrtf(2.0f));

    iftImage *pred = iftCreateImage(contour->xsize, contour->ysize, contour->zsize);

    for (int i = 0; i < contour->n; i++)
        pred->val[i] = IFT_NIL;

    int src = iftGetVoxelIndex(contour, voxel_array->val[init_index]);
    iftInsertLIFO(Q, src);

    int ori_ref = _get_direction(obj, A, src);

    while (!iftEmptyLIFO(Q))
    {
        int p = iftRemoveLIFO(Q);
        iftVoxel u = iftGetVoxelCoord(contour, p);
        int j = ori_ref;
        if (pred->val[p] != IFT_NIL) {
            iftVoxel t = iftGetVoxelCoord(pred, pred->val[p]);
            j = ift6NeighborsClockFromLeft(&t, &u);
        }

        iftSet *adj = NULL;
        for (int i = 0; i < A->n; i++, j = (j + 1) % A->n)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, j);
            if (!iftValidVoxel(contour, v))
                continue;
            int q = iftGetVoxelIndex(contour, v);
            if (contour->val[q] != 0 && Q->color[q] != IFT_BLACK &&
                dist->val[q] <= dist->val[p])
            {
                pred->val[q] = p;
                dist->val[q] = dist->val[p] + 1;
                if (Q->color[q] == IFT_WHITE)
                    iftInsertSet(&adj, q);
            }
        }

        for (iftSet *s = adj; s; )
        {
            iftInsertLIFO(Q, s->elem);
            iftSet *aux = s;
            s = s->next;
            iftFree(aux);
        }
    }

    int p = src;
    for (int i = 1; i < A->n; i++) {
        iftVoxel v = iftGetAdjacentVoxel(A, voxel_array->val[init_index], i);
        if (iftValidVoxel(contour, v)) {
            int q = iftGetVoxelIndex(contour, v);
            if (dist->val[p] < dist->val[q])
                p = q;
        }
    }

    iftDestroyLIFO(&Q);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&dist);

    PyObject *list = PyList_New(voxel_array->n);
    for (int i = 0; i < voxel_array->n; i++)
        PyList_SetItem(list, i, PyList_New(0));

    int cur_vox = init_index, cur_p;
    do {
        cur_vox--;
        if (cur_vox < 0)
            cur_vox = voxel_array->n - 1;
        cur_p = iftGetVoxelIndex(contour, voxel_array->val[cur_vox]);
    } while (contour->val[cur_p] == 0);

    PyObject *cur_l = PyList_GetItem(list, cur_vox);

    for (; p != IFT_NIL; p = pred->val[p])
    {
        PyList_Append(cur_l, PyLong_FromLong((long) p));
        if (cur_p == p) {
            do {
                cur_vox--;
                if (cur_vox == -1) // cycle
                    cur_vox = voxel_array->n - 1;
                cur_p = iftGetVoxelIndex(contour, voxel_array->val[cur_vox]);
                cur_l = PyList_GetItem(list, cur_vox);
            } while (contour->val[cur_p] == 0); // to avoid error when voxel is not in contour
        }
    }

    iftDestroyImage(&pred);

    // sorting to the correct order
    for (int i = 0; i < PyList_Size(list); i++) {
        PyObject *l = PyList_GetItem(list, i);
        int info = PyList_Reverse(l);
        if (info) {
            SWIG_Error(0, "It was not possible to sort the path pixels on the correct order");
        }
    }

    return list;
}


int AnchorIndex(PyObject *paths, int pixel_index)
{
    int path_idx = -1;
    PyObject *path = NULL;
    for (int i = 0; i < PyList_Size(paths) && path_idx == -1; i++) {
        path = PyList_GetItem(paths, i);
        for (int j = 0; j < PyList_Size(path); j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            if (pos == pixel_index)
                path_idx = i;
        }
    }
    return path_idx;
}


int InsertAnchorAndSplitPaths(iftVoxelArray *anchors, PyObject *paths, const iftImage *ref, int x, int y)
{
    if (!PyList_CheckExact(paths))
        SWIG_Error(0, "Input \"paths\" must be a list.");

    iftVoxel k; k.x = x; k.y = y; k.z = 0;
    int pt = iftGetVoxelIndex(ref, k);

    int path_idx = -1, pos_idx = -1;
    PyObject *path = NULL;
    for (int i = 0; i < PyList_Size(paths) && path_idx == -1; i++) {
        path = PyList_GetItem(paths, i);
        for (int j = 0; j < PyList_Size(path) && pos_idx == -1; j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            if (pos == pt) {
                path_idx = i;
                pos_idx = j;
            }
        }
    }
    if (path_idx == -1)
        return NULL;

    anchors->n++;
    anchors->val = iftRealloc(anchors->val, anchors->n * sizeof *anchors->val);
    for (int i = (anchors->n-2); i > path_idx; i--)
        iftCopyVoxel(&anchors->val[i], &anchors->val[i + 1]);
    iftCopyVoxel(&k, &anchors->val[path_idx + 1]);

    PyObject *back_path = PyList_GetSlice(path, 0, pos_idx);
    PyObject *front_path = PyList_GetSlice(path, pos_idx, PyList_Size(path));

    PyList_SetItem(paths, path_idx, back_path);
    if (path_idx == (PyList_Size(paths) - 1)) {
        PyList_Append(paths, front_path);
    } else {
        PyList_Insert(paths, path_idx + 1, front_path);
    }
//    Py_DECREF(path);

    return path_idx + 1;
}


iftVoxelArray *JoinPaths(PyObject *paths, const iftImage *ref)
{
    if (!PyList_CheckExact(paths))
        SWIG_Error(0, "Input \"paths\" must be a list.");

    iftVoxelArray *arr = iftCreateVoxelArray(ref->n);

    int idx = 0;
    for (int i = 0; i < PyList_Size(paths); i++) {
        PyObject *path = PyList_GetItem(paths, i);
        for (int j = 0; j < PyList_Size(path); j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            iftVoxel v = iftGetVoxelCoord(ref, pos);
            arr->val[idx].x = v.x; arr->val[idx].y = v.y; arr->val[idx].z = v.z;
            idx++;
        }
    }
    arr->n = idx;
    arr->val = iftRealloc(arr->val, idx * sizeof *arr->val);

    return arr;
}


int iftFindAnchorInPaths(PyObject *paths, int index, bool split_path)
{
    if (!PyList_CheckExact(paths))
        SWIG_Error(0, "Input \"paths\" must be a list.");

    int path_idx = IFT_NIL, pos_idx = IFT_NIL;
    PyObject *path = NULL;
    for (int i = 0; i < PyList_Size(paths) && path_idx == IFT_NIL; i++) {
        path = PyList_GetItem(paths, i);
        for (int j = 0; j < PyList_Size(path) && pos_idx == IFT_NIL; j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            if (index == pos) {
                path_idx = i;
                pos_idx = j;
            }
        }
    }

    if (pos_idx == IFT_NIL)
        return IFT_NIL;

    if (split_path) {
        PyObject *back_path = PyList_GetSlice(path, 0, pos_idx);
        PyObject *front_path = PyList_GetSlice(path, pos_idx, PyList_Size(path));

        PyList_SetItem(paths, path_idx, back_path);
        if (path_idx == (PyList_Size(paths) - 1)) {
            PyList_Append(paths, front_path);
        } else {
            PyList_Insert(paths, path_idx + 1, front_path);
        }
    }

    return path_idx;
}

int FindAnchorSuggestionByArea(PyObject *paths,
                               const iftImage *gradient,
                               const iftImage *source,
                               const iftImage *target,
                               const iftImage *forbidden,
                               float window)
{
    iftImage *err = iftXor(source, target);

    window = iftMin(1.0, window);
    window = iftMax(0.0, window);

    iftAdjRel *A = iftIs3DImage(gradient) ? iftSpheric(1.75) : iftCircular(M_SQRT2);
    iftImage *largest = iftSelectLargestComp(err, A);
    iftDestroyImage(&err);
    err = iftDilate(largest, A, NULL);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&largest);

    int path_idx = IFT_NIL, max_sum = 0;
    PyObject *path = NULL;
    for (int i = 0; i < PyList_Size(paths); i++) {
        path = PyList_GetItem(paths, i);
        int path_sum = 0;
        for (int j = 0; j < PyList_Size(path); j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            path_sum += (err->val[pos] && !forbidden->val[pos]);
        }
        if (path_sum > max_sum) {
            max_sum = path_sum;
            path_idx = i;
        }
    }

    iftDestroyImage(&err);

    if (path_idx == IFT_NIL)
        return IFT_NIL;

    path = PyList_GetItem(paths, path_idx);
    int size = PyList_Size(path);
    int mid = size / 2;
    int lower = mid - (size * window / 2), upper = mid + (size * window / 2);
    upper = iftMin(size, upper);
    lower = iftMax(0, lower);
    while (lower > 0 && forbidden->val[lower])
        lower--;
    int weakest = PyLong_AsLong(PyList_GetItem(path, lower));

    for (int i = lower; i < upper; i++) {
        long p = PyLong_AsLong(PyList_GetItem(path, i));
        if (!forbidden->val[p] && gradient->val[p] < gradient->val[weakest])
            weakest = p;
    }

    return weakest;
}


int SelectAnchorByLargestError(PyObject *paths,
                               const iftImage *source,
                               const iftImage *target,
                               const iftVoxelArray *src_anchors,
                               const iftVoxelArray *tgt_anchors)
{
    iftImage *err = iftXor(source, target);

    iftAdjRel *A = iftIs3DImage(source) ? iftSpheric(1.75) : iftCircular(M_SQRT2);
    iftImage *aux = iftSelectLargestComp(err, A);
    iftDestroyImage(&err);
    err = iftDilate(aux, A, NULL);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&aux);

    int size = PyList_Size(paths);
    int *count = iftAlloc(size, sizeof *count);
    PyObject *path = NULL;
    for (int i = 0; i < size; i++) {
        int prev_id = i - 1;
        if (prev_id < 0)
            prev_id = size - 1;
        path = PyList_GetItem(paths, i);
        bool found = 0;
        for (int j = 0; j < PyList_Size(path); j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            found |= err->val[pos] != 0;
        }
        count[prev_id] |= found;
        count[i] |= found;
    }

    iftDestroyImage(&err);

    int *dist = iftAlloc(size, sizeof *dist);
    int most = 0;
    for (int i = 0; i < size; i++)
    {
        if (count[i]) {
            dist[i] = iftSquaredVoxelDistance(src_anchors->val[i], tgt_anchors->val[i]);
            if (dist[i] > dist[most])
                most = i;
        }
    }

    int most_dist = dist[most];
    iftFree(dist);
    iftFree(count);
    if (most_dist == 0)
        return IFT_NIL;

    return most;
}


//! swig(newobject)
iftImage *iftPathsImage(PyObject *paths, const iftImage *ref)
{
    if (!PyList_CheckExact(paths))
        SWIG_Error(0, "Input \"paths\" must be a list.");

    iftImage *image = iftCreateImageFromImage(ref);

    PyObject *path = NULL;
    for (int i = 0; i < PyList_Size(paths); i++) {
        path = PyList_GetItem(paths, i);
        int count = 0;
        for (int j = 0; j < PyList_Size(path); j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
           image->val[pos] = ++count;
        }
    }

    return image;
}


void iftImageFillPath(iftImage *image,
                      PyObject *paths,
                      int value)
{
    PyObject *path = NULL;
    for (int i = 0; i < PyList_Size(paths); i++) {
        path = PyList_GetItem(paths, i);
        for (int j = 0; j < PyList_Size(path); j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            image->val[pos] = value;
        }
    }
}

void ImageFillVoxelArray(iftImage *image,
                         iftVoxelArray *arr,
                         int value)
{
    for (int i = 0; i < arr->n; i++)
    {
        int p = iftGetVoxelIndex(image, arr->val[i]);
        image->val[p] = value;
    }
}

void FImageFillVoxelArray(iftFImage *image,
                          iftVoxelArray *arr,
                          float value)
{
    for (int i = 0; i < arr->n; i++)
    {
        int p = iftGetVoxelIndex(image, arr->val[i]);
        image->val[p] = value;
    }
}

/*
void iftAddSecurityAnchors(const iftImage *source,
                           const iftImage *target,
                           const iftImage *src_contour,
                           const iftImage *tgt_contour,
                           PyObject *src_paths,
                           PyObject *tgt_paths,
                           iftVoxelArray *src_anchors,
                           iftVoxelArray *tgt_anchors)
{
    iftAdjRel *A = iftIs3DImage(source) ? iftSpheric(1.75) : iftCircular(M_SQRT2);
    iftImage *err = iftXor(source, target);

    iftImage *aux = iftSelectCompAboveArea(err, A, 100);
    err = iftDilate(aux, A, NULL);
    iftDestroyImage(&aux);

    iftImage *both = iftAnd(src_contour, tgt_contour);
    aux = iftAnd(err, both);
    iftDestroyImage(&err);
    iftDestroyImage(&both);

    const int min_dist = 5;
    for (int i = 0; i < PyList_Size(src_paths); i++) {
        PyObject *path = PyList_GetItem(src_paths, i);
        int size = PyList_Size(path);
        for (int j = 0; j < size && j < min_dist; j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            aux->val[pos] = 0;
        }

        for (int j = iftMax(0, size - min_dist); j < size; j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            aux->val[pos] = 0;
        }
    }

    iftWriteImageByExt(iftNormalize(aux, 0, 255), "aux.png");
}
*/

bool iftAddSecurityAnchors(const iftImage *source,
                           const iftImage *target,
                           const iftImage *src_contour,
                           const iftImage *tgt_contour,
                           int index,
                           PyObject *src_paths,
                           PyObject *tgt_paths,
                           iftVoxelArray *src_anchors,
                           iftVoxelArray *tgt_anchors)
{
    iftAdjRel *A = iftIs3DImage(source) ? iftSpheric(1.75) : iftCircular(M_SQRT2);
    iftImage *err = iftXor(source, target);

    iftImage *aux = iftSelectCompAboveArea(err, A, 100);
    err = iftDilate(aux, A, NULL);
    iftDestroyImage(&aux);

    iftImage *both = iftAnd(src_contour, tgt_contour);
    aux = iftAnd(err, both);
    iftDestroyImage(&err);
    iftDestroyImage(&both);

    PyObject *path = NULL;
    int size = PyList_Size(tgt_paths);
    for (int i = 0; i < size; i++) {
        path  = PyList_GetItem(tgt_paths, i);
        long pos = PyLong_AsLong(PyList_GetItem(path, 0));
        aux->val[pos] = 0;
        pos = PyLong_AsLong(PyList_GetItem(path, PyList_Size(path) - 1));
        aux->val[pos] = 0;
    }

    int next_index = index + 1;
    const int min_dist = 5;
    int next_cross = IFT_NIL;
    int cross_index = IFT_NIL;

    path = PyList_GetItem(src_paths, index);
    if (!path) return NULL;
    size = PyList_Size(path);
    for (int j = min_dist; j < (size - min_dist); j++) {
        long pos = PyLong_AsLong(PyList_GetItem(path, j));
        if (aux->val[pos]) {
            next_cross = pos;
            cross_index = j;
            break;
        }
    }

    if (next_cross != IFT_NIL)
    {
        path = PyList_GetItem(tgt_paths, index);
        if (!path) return NULL;
        size = PyList_Size(path);
        for (int j = 1; j < size - 1; j++) {
            long pos = PyLong_AsLong(PyList_GetItem(path, j));
            if (pos == next_cross)
            {
                iftVoxel u = iftGetVoxelCoord(source, pos);
                iftInsertVoxel(src_anchors, next_index, &u);
                PyObject *src_p = PyList_GetItem(src_paths, index);
                PyObject *prev_slice = PyList_GetSlice(src_p, 0, cross_index);
                PyObject *next_slice = PyList_GetSlice(src_p, cross_index, PyList_Size(src_p));

                PyList_SetItem(src_paths, index, prev_slice);
                if (next_index == size)
                    PyList_Append(src_paths, next_slice);
                else
                    PyList_Insert(src_paths, next_index, next_slice);

                iftInsertVoxel(tgt_anchors, next_index, &u);
                prev_slice = PyList_GetSlice(path, 0, j);
                next_slice = PyList_GetSlice(path, j, size);

                PyList_SetItem(tgt_paths, index, prev_slice);
                if (next_index == size)
                    PyList_Append(tgt_paths, next_slice);
                else
                    PyList_Insert(tgt_paths, next_index, next_slice);
                break;
            }
        }
    }

    if (PyList_Size(tgt_paths) != PyList_Size(src_paths))
    {
        PyErr_SetString(PyExc_RuntimeError, "Target and source paths size does not match.");
        return NULL;
    }

    int prev_index = index - 1;
    if (prev_index < 0)
        prev_index = src_anchors->n - 1;

    int prev_cross = IFT_NIL;
    cross_index = IFT_NIL;
    path = PyList_GetItem(src_paths, prev_index);
    if (!path) return NULL;
    size = PyList_Size(path);
    for (int j = iftMax(0, size - min_dist); j >= min_dist; j--) {
        long pos = PyLong_AsLong(PyList_GetItem(path, j));
        if (aux->val[pos]) {
            prev_cross = pos;
            cross_index = j;
            break;
        }
    }

    if (prev_cross == IFT_NIL)
        return false;

    path = PyList_GetItem(tgt_paths, prev_index);
    if (!path) return NULL;
    size = PyList_Size(path);
    bool added = false;
    for (int j = 1; j < size - 1; j++) {
        long pos = PyLong_AsLong(PyList_GetItem(path, j));
        if (pos == prev_cross)
        {
            iftVoxel u = iftGetVoxelCoord(source, pos);
            iftInsertVoxel(src_anchors, index, &u);
            PyObject *src_p = PyList_GetItem(src_paths, prev_index);
            PyObject *prev_slice = PyList_GetSlice(src_p, 0, cross_index);
            PyObject *next_slice = PyList_GetSlice(src_p, cross_index, PyList_Size(src_p));

            PyList_SetItem(src_paths, prev_index, prev_slice);
            PyList_Insert(src_paths, index, next_slice);

            iftInsertVoxel(tgt_anchors, index, &u);
            prev_slice = PyList_GetSlice(path, 0, j);
            next_slice = PyList_GetSlice(path, j, size);

            PyList_SetItem(tgt_paths, prev_index, prev_slice);
            PyList_Insert(tgt_paths, index, next_slice);
            added = true;
            break;
        }
    }

    if (PyList_Size(tgt_paths) != PyList_Size(src_paths))
    {
        PyErr_SetString(PyExc_RuntimeError, "Target and source paths size does not match.");
        return NULL;
    }
    return added;
}



//! swig(newobject)
PyObject *iftInsideOutSideNeighbour(const iftImage *object,
                                    int index)
{
    iftAdjRel *A = iftClockCircular(sqrtf(2.0f));
    iftAdjRel *L = iftLeftSide(A, 3.0), *R = iftRightSide(A, 3.0);
    int j = _get_direction(object, A, index);


    iftVoxel u = iftGetVoxelCoord(object, index);
    for (int i = 0; i < A->n; i++, j = (j + 1) % A->n) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, j);
        if (iftValidVoxel(object, v))
        {
            int p = iftGetVoxelIndex(object, v);
            if (object->val[p])
                break;
        }
    }

    iftVoxel left = iftGetAdjacentVoxel(L, u, j);
    iftVoxel right = iftGetAdjacentVoxel(R, u, j);

    int l = index, r = index;
    if (iftValidVoxel(object, left))
        l = iftGetVoxelIndex(object, left);
    if (iftValidVoxel(object, right))
        r = iftGetVoxelIndex(object, right);

    iftDestroyAdjRel(&A);
    iftDestroyAdjRel(&L);
    iftDestroyAdjRel(&R);
    return Py_BuildValue("(ii)", l, r);
}