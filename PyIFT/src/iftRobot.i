%include "iftImage.i"
%include "iftRepresentation.i"
%include "iftSimilarity.i"



typedef struct {

    int stop_iter; /* default */
    int cur_iter; /* default */

    float stop_acur; /* default */
    float cur_acur; /* default */

    int max_seed_per_iter; /* default */
    bool converged; /* default */

    iftImage *orig; /* default */
    iftImage *gt_img; /* default */
    iftImage *segm; /* default */
    iftImage *error; /* default */

    iftLabeledSet* seeds; /* default */

    iftAdjRel *mrk_radius; /* default */
    iftAdjRel *limit_mrk_rad; /* default */

} iftRobot;

%extend iftRobot {

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
	
	

	~iftRobot() {
		iftRobot* ptr = ($self);
		iftDestroyRobot(&ptr);
	}
};

%newobject iftCreateRobot;
%feature("autodoc", "2");
iftRobot *iftCreateRobot(const char *orig_path, const char *gt_path, int stop_iter, float stop_acur);

%feature("autodoc", "2");
void iftResetRobotSegmentation(iftRobot *bot);

%feature("autodoc", "2");
void iftRobotUpdateError(iftRobot *bot);

%feature("autodoc", "2");
void iftRobotFindSeedsMSSkel(iftRobot *bot);

%feature("autodoc", "2");
void iftRobotFindSeedsCenterOfMass(iftRobot *bot);

%feature("autodoc", "2");
void iftRobotPrintInfo(iftRobot *bot);

%feature("autodoc", "2");
bool iftRobotContinue(iftRobot *bot);

