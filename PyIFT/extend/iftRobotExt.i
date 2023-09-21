void setOrig(const iftImage *image)
{
    iftRobot *bot = ($self);
    iftDestroyImage(&bot->orig);
    bot->orig = iftCopyImage(image);
}

void setGt(const iftImage *image)
{
    iftRobot *bot = ($self);
    iftDestroyImage(&bot->gt_img);
    bot->gt_img = iftCopyImage(image);
}

void setSegm(const iftImage *image)
{
    iftRobot *bot = ($self);
    iftDestroyImage(&bot->segm);
    bot->segm = iftCopyImage(image);
}

void setError(const iftImage *image)
{
    iftRobot *bot = ($self);
    iftDestroyImage(&bot->error);
    bot->error = iftCopyImage(image);
}

void *setSeeds(iftLabeledSet *seeds)
{
    iftRobot *bot = ($self);
    iftDestroyLabeledSet(&bot->seeds);
    bot->seeds = iftCopyLabeledSet(seeds);
}

void setMarkerRadius(const iftAdjRel *A)
{
    iftRobot *bot = ($self);
    iftDestroyAdjRel(&bot->mrk_radius);
    bot->mrk_radius = iftCopyAdjacency(A);
}

void setMarkerMaxRadius(const iftAdjRel *A)
{
    iftRobot *bot = ($self);
    iftDestroyAdjRel(&bot->limit_mrk_rad);
    bot->limit_mrk_rad = iftCopyAdjacency(A);
}

