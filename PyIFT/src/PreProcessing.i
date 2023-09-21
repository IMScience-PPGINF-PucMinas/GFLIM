%include "iftImage.i"
%include "iftPlane.i"

%newobject iftBrainMRIPreProcessing;
%feature("autodoc", "2");
iftImage *iftBrainMRIPreProcessing(const iftImage *mri, int nbits, const iftPlane *msp_in, const iftImage *mask,
                                   const iftImage *ref_mri, const iftImage *ref_mask, bool skip_n4,
                                   bool skip_median_filter, bool skip_msp_alignment, bool skip_hist_matching,
                                   iftPlane **msp_out);

