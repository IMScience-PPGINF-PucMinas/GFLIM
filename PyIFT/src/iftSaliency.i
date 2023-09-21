%include "iftDataSet.i"
%include "iftImage.i"
%include "iftIGraph.i"



typedef struct {
    int itself_iterations;
    int number_superpixels;
    int obj_seeds_number;
    float oisf_gamma;
    int oisf_iterations;
    float oisf_beta;
    float oisf_alpha;
    float oisf_dist_penalty;
    float superpixel_increase_ratio;
    float query_importance;
    float normalization_value;
    float integration_lambda;
    int integration_iteration;
    int enhancement_type;
    iftITSELFPriors prior_params;
}iftITSELFParameters;

%extend iftITSELFParameters {

	~iftITSELFParameters() {
		iftITSELFParameters* ptr = ($self);
		iftDestroyITSELFParameters(&ptr);
	}
};

%newobject iftInitializeITSELFParametersByDefaultU2Net;
%feature("autodoc", "2");
iftITSELFParameters *iftInitializeITSELFParametersByDefaultU2Net();

%newobject iftInitializeITSELFParametersByDefaultScribbles;
%feature("autodoc", "2");
iftITSELFParameters *iftInitializeITSELFParametersByDefaultScribbles();

%newobject iftITSELFParametersSetParam;
%feature("autodoc", "2");
void iftITSELFParametersSetParam(iftITSELFParameters *params, char *parameter_name, float value);

%newobject iftSESS;
%feature("autodoc", "2");
iftImage *iftSESS(iftImage *original_image, iftImage *initial_saliency, iftITSELFParameters *params, iftMImage *features);

