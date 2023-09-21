#include "iftReconstruction.h"


#include "ift/core/io/Stream.h"


/*********************** PRIVATE ********************************************************************/

int iftDFour1(double data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			iftSwap(data[j], data[i]);
			iftSwap(data[j+1], data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
	return(1);
}

int iftGammln(double xx, double *yy)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	*yy = -tmp+log(2.5066282746310005*ser/x);
	return(1);
}

//////////////////////////////////////////////////////////////////////////////////////////
// Generate backprojection filter (see Kak-Slaney Chapter 3 Eq.123 (page 91)
// and if desired window the filter with XXXX filter
/////////////////////////////////////////////////////////////////////////////////////////

int iftGenerateTomoFilter(iftComplex* H, int nFilterPts, double dSignalXInc, int bWindow)
{
	int i;

	for (i = 0; i < nFilterPts; i++) 
	{
		H[i].r = 0.0; //real
		H[i].i = 0.0; //imag
	}

	iftComplex* Hi = H+1;
	iftComplex* Hf = &H[nFilterPts - 1];

	int mid = nFilterPts / 2;

	for(i = 1; i < nFilterPts/2 - 1; i += 2)
	{
		(*Hi).r =  -0.5 / ((double)(i * i) * IFT_PI * IFT_PI * dSignalXInc * dSignalXInc);
		(*Hf).r = (*Hi).r;
		Hi += 2; Hf -= 2;
	}

	//zero position
	(*H).r = 1.0f / (8.0 * dSignalXInc * dSignalXInc);

	if(mid%2 != 0) //check to see if midpoint is even
		(*(H+mid)).r = -0.5 / ((double)(mid * mid) * IFT_PI * IFT_PI * dSignalXInc * dSignalXInc);

	//convert to frequency domain
	double *F = (double *)H;
	iftDFour1(F-1, nFilterPts, 1);

	////////////////////////////////////////////////////////////////////////////////
	//generate window filter in frequency space and multiple it by the backprojection
	//  filter (imaginary component is zero)
	////////////////////////////////////////////////////////////////////////////////

	if(bWindow)
	{
		double con = 2.0 * IFT_PI / nFilterPts;
		Hi = H+1;
		Hf = &H[nFilterPts - 1];
		mid = nFilterPts / 2;

		for(i = 1; i < nFilterPts/2 - 1; i++)
		{
			(*Hi).r = (*Hf).r *= 0.54 + 0.46 * cos(con * (double)i);
			Hi++; Hf--;
		}

		//zero point
		(*H).r *= 1.0f;

		(*(H+mid)).r *= 0.54 + 0.46 * cos(con * (double)mid);
	}
	return(1);
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Extend the projection on either end by replicating the first/last signal value
///////////////////////////////////////////////////////////////////////////////////////////////

int iftExtendProjections_Constant(double *pProj, iftComplex* R, int nSignalXPts, int extFactor)
{
	int j;
	iftComplex* pR = R;

	int ind = extFactor * nSignalXPts;

 	for(j = 0; j < ind; j++)
 		(*pR++).r = pProj[0];

 	for(j = 0; j < nSignalXPts; j++)
 		(*pR++).r = pProj[j];

 	for(j = 0; j < ind; j++)
 		(*pR++).r = pProj[nSignalXPts-1];
	return(1);
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Generates the filter for lambda tomography using the filter proposed in
//		1)A.Faridani, E.Ritman, and K.Smith, "Local Tomography," SIAM(Soc. Ind. Appl. Math.)
//		  J. Appl. Math. 52, 459-484 (1992)
//      2)P. Huabsomboom, "3D filtered backprojection algorithm for local tomography",
//        Master's Thesis, Oregon State University, Dept. of Math, 2000.
//////////////////////////////////////////////////////////////////////////////////////////////

int iftGenerateLambdaTomoFilter(iftComplex* H, int nFilterPts, double dSignalInc, int scaleFactor)
{
	int i;
	double gam1, gam2;

	for (i = 0; i < nFilterPts; i++) 
	{
		H[i].r = 0.0; //real
		H[i].i = 0.0; //imag
	}
		
	double alpha = 11.4174;

	//PSF radius = omega
	double omega = dSignalInc * (double)scaleFactor * sqrt((2.0*alpha + 1.0) / 3.0);
	double con = 50.848111446;

	//calculates ln Gamma(x)
	iftGammln(12.4174, &gam2);
	iftGammln(13.9174, &gam1);

	con = 2.0 * exp(gam1 - gam2) /sqrt(IFT_PI);
	con /= (omega*omega*omega);

	iftComplex* Hi = H+1;
	iftComplex* Hf = &H[nFilterPts - 1];

	for(i = 1; i < nFilterPts/2 - 1; i++)
	{
		double rad = (double)i * dSignalInc / omega;
		if(rad < 1.0)
		{
			(*Hi).r = con * pow( (1.0-rad*rad),alpha-1.0 ) * (1.0 - (2*alpha+1)*rad*rad);
			(*Hf).r = (*Hi).r;
		}
		else
		{
			(*Hi).r = 0.0;
			(*Hf).r = (*Hi).r;
		}
		Hi++; Hf--;
	}

	//zero position
	(*H).r = con;

	//convert to frequency domain
	double *F = (double *)H;
	iftDFour1(F-1, nFilterPts, 1);
	return(1);
} 







/**************    PUBLIC  **************************************************************************/	

////////////////////////////////////////////////////////////////////////////////////////////////
// Filtered Back Projection
// output 3D voxel grid corresponds to the detector being projected back to the origin
//                 (imaginary detector through origin of phantom)
// uses projective geometry
// voxel volume centered around origin
// short scan
// subsampling 
///////////////////////////////////////////////////////////////////////////////////////////////

int iftFilteredBackProjection(
			      double *proj,                    //projection data
			      int numProj,                     //number of projections
			      int nSignalXPts,                 //pixels detector
			      int nSignalYPts,                 //lines detector
			      double dSignalXInc,              //pixel spacing
			      double dSignalYInc,              //line spacing
			      double AOR,			    //column on detector that corresponds to the axis of rotation 
			      double midPlane,                 //row on detector that corresponds to the midplane
			      double dFocalLength,             //source to actual AOR distance
			      double dSourceDetectorLength,    //source to detector distance
			      double *theta,                   //angle data
			      double *vox,                     //slice data
			      int numSlices,                   //number of slices,
			      int pixels,		       //cols per slice
			      int lines,                       //rows per slice
			      double subSample,                //sub sampling factor
			      int extProj,                     //extend projection by extProj*nSignalXPts 
			      int bWindow                      //use smoothing window?
			      )
{
	
  int i, iProj, iSlice, iX, iY, l, ind;
  int nFilterPts;
  iftComplex    *H, * R;
  double *deltaProj;			
  int subSampleX = (int)subSample;
  int subSampleY = (int)subSample;
  int subSampleZ = (int)subSample;
  
  // scale signalInc  to adjust for imaginary detector through origin of phantom 
  dSignalXInc *= dFocalLength / dSourceDetectorLength;
  dSignalYInc *= dFocalLength / dSourceDetectorLength;
	
  //center in pixel coordinates
  /* double centerX = (double)AOR;          */
  /* double centerY = (double)midPlane; */
  /* double centerZ = centerY; */
	
  //increase nSignalXPts and change centerX (0.0 pt) to account for
  //extended projections

  int nSignalXPts_ext = (2 * extProj + 1) * nSignalXPts;
  double centerX_ext  = extProj * nSignalXPts + AOR;
  
  //////////////////////////////////////////////////////////////////////////
  //determine filter length with appropriate zero padding
  // increase array size to 2^n larger than 2*m_nSignalPts
  //////////////////////////////////////////////////////////////////////////
  
  double dLogBase2 = log(2.0 * (double)nSignalXPts_ext) / log(2.0);
  int    iLogBase2 = (int)(floor (dLogBase2));
  if(dLogBase2 == (double)iLogBase2)
    nFilterPts = 1 << iLogBase2;
  else
    nFilterPts = 1 << (iLogBase2 + 1);
  
  /////////////////////////////////////////////////////////////////////////////
  //generate backprojection filter in frequency domain
  //(see Kak-Slaney Chapter 3 Eq.123 (page 91)
  /////////////////////////////////////////////////////////////////////////////
  
  H = iftAllocComplexArray(nFilterPts);
  
  double *F;
  
  iftGenerateTomoFilter(H, nFilterPts, dSignalXInc, bWindow);
  
  R = iftAllocComplexArray(nFilterPts);
  double disFromCenterX, disFromCenterY, disFromCenterY2;
  
  double P[3][4]; // projection matrix
  double Vx[3], Vy[3], Vz[3], PM0[3], M0[3];
  double p, eta, p1, p2, con;
  int lo_p, lo_eta;
  double frac_p, frac_eta;
  double *pProj, *pVox;
  double ccos, ssin;
  double fac;
  double dFocalLength_pix  = dFocalLength / dSignalXInc; //focal length in pixel units
  double dFocalLength2_pix = dFocalLength_pix * dFocalLength_pix;
	
  int nProj = numProj;
  
  double dBeta = 2.0 * IFT_PI / (double)(numProj);
  
  //array to hold error projection used for linear interpolation
  
  deltaProj = iftAllocDoubleArray(nSignalXPts*nSignalYPts);
  
  for(iProj = 0; iProj < nProj; iProj++)
    {
      printf("\rprojection # %d \n",iProj);
      fflush(stdout);
      
      ////////////////////////////////////////////////////////////////////////////////////
      //generate filtered backprojection for one slanted fan at a time
      // y direction of projection image : slanted fan position
      // x direction of projection image is the projection data for the slanted fan
      ////////////////////////////////////////////////////////////////////////////////////
      for (l = 0; l < nSignalYPts; l++) //one slanted fan at a time (y direction of projection image)
	{
	  ///////////////////////////////////
	  //generate modified input signal 
	  ///////////////////////////////////
	  disFromCenterY = ((double)l - midPlane ) * dSignalYInc;
	  disFromCenterY2 = disFromCenterY * disFromCenterY;
	  pProj = proj + iProj * nSignalYPts * nSignalXPts + l*nSignalXPts;
	  
	  //extend projections
	  
	  iftExtendProjections_Constant(pProj, R, nSignalXPts, extProj);
	  
	  for (i = 0; i < nSignalXPts_ext; i++) 
	    {
	      R[i].i = 0.0f;     //imag
	      disFromCenterX = ((double)i - centerX_ext ) * dSignalXInc;
	      fac = dFocalLength * dFocalLength + disFromCenterX * disFromCenterX + disFromCenterY2;
	      R[i].r *= dFocalLength / sqrt(fac); // real
	    } 
	  //pad with zeros
	  for (i = nSignalXPts_ext; i < nFilterPts; i++)
	    {
	      R[i].r = 0.0f;  
	      R[i].i = 0.0f;
	    }
	  
	  //convert to modified projection data to frequency domain
	  F = (double *)R;
	  iftDFour1(F-1, nFilterPts, 1);
	  
	  //mulitply modified projection data by filter data
	  double real, imag;
	  for(i = 0; i < nFilterPts; i++)
	    {
	      real = R[i].r * H[i].r - R[i].i * H[i].i;
	      imag = R[i].i * H[i].r + R[i].r * H[i].i;
	      R[i].r = real;
	      R[i].i = imag;	      
	    }
	  
	  //convert filtered projection to spatial domain
	  F = (double *)R;
	  iftDFour1(F-1, nFilterPts, -1);
	  
		
	  //normalize by nFilterPts and put filtered projection back into original projection image
	  int ind = extProj * nSignalXPts;
	  for(i = 0; i < nSignalXPts; i++)
	    {
	      pProj[i] = R[i + ind].r /(double)nFilterPts;
	    }
	
	  /* precalculate difference for linear interpolation
	     and normalize by nFilterPts and put filtered
	     projection back into original and put filtered
	     projection back into original projection image */
	  double *ptr = deltaProj + l * nSignalXPts;
	  for(i = 0; i < nSignalXPts - 1; i++)
	    ptr[i] = (pProj[i+1] - pProj[i]);
	  
	  ptr[nSignalXPts - 1] = 0.0;
	}
    
      ///////////////////////////////////////////////////////////////////////////////
      //perform weighted backprojection of filtered projection along corresponding
      //    ray of the titled fan
      /////////////////////////////////////////////////////////////////////////////////
      ccos = cos(theta[iProj]);
      ssin = sin(theta[iProj]);
      
      //calculate projection matrix for given view angle
      
      P[0][0] = AOR * ccos - dFocalLength_pix * ssin * subSampleX;
      P[0][1] = AOR * ssin + dFocalLength_pix * ccos * subSampleY;
      P[0][2] = 0.0;
      P[0][3] = dFocalLength_pix * AOR;
      P[1][0] = midPlane * ccos;
      P[1][1] = midPlane * ssin;
      P[1][2] = dFocalLength_pix * subSampleZ;
      P[1][3] = dFocalLength_pix * midPlane;
      P[2][0] = ccos;
      P[2][1] = ssin;
      P[2][2] = 0.0;
      P[2][3] = dFocalLength_pix;
      
      P[0][0] = AOR * ccos + dFocalLength_pix * ssin * subSampleX;
      P[0][1] = -AOR * ssin + dFocalLength_pix * ccos * subSampleY;
      P[0][2] = 0.0;
      P[0][3] = dFocalLength_pix * AOR;
      P[1][0] = midPlane * ccos;
      P[1][1] = -midPlane * ssin;
      P[1][2] = dFocalLength_pix * subSampleZ;
      P[1][3] = dFocalLength_pix * midPlane;
      P[2][0] = ccos;
      P[2][1] = -ssin;
      P[2][2] = 0.0;
      P[2][3] = dFocalLength_pix;
      
      pProj = proj + iProj * nSignalYPts * nSignalXPts;
      pVox = vox;
	    
      M0[0] = -0.5 * (double)(pixels-1);
      M0[1] = -0.5 * (double)(lines-1);
      M0[2] = -0.5 * (double)(numSlices-1);
      
      PM0[0] = P[0][0] * M0[0] + P[0][1] * M0[1] + P[0][2] * M0[2] + P[0][3];
      PM0[1] = P[1][0] * M0[0] + P[1][1] * M0[1] + P[1][2] * M0[2] + P[1][3];
      PM0[2] = P[2][0] * M0[0] + P[2][1] * M0[1] + P[2][2] * M0[2] + P[2][3];
      
      Vz[0] = PM0[0];
      Vz[1] = PM0[1];
      Vz[2] = PM0[2];
      
      for(iSlice = 0; iSlice < numSlices; iSlice++) //loop over slices of reconstructed image
	{
	  Vy[0] = Vz[0];
	  Vy[1] = Vz[1];
	  Vy[2] = Vz[2];
	  for(iY = 0; iY < lines; iY++) //loop over lines of reconstructed image
	    {
	      Vx[0] = Vy[0];
	      Vx[1] = Vy[1];
	      Vx[2] = Vy[2];
	      for(iX = 0; iX < pixels; iX++)  //loop over pixes of reconstructed image
		{
		  
		  con = 1.0 / Vx[2];
		  //determine p (in pixels) location in distance on virtual detector for the ray going thru pt(x,y,z)
		  p = Vx[0] * con;
		  //determine eta location in distance on virtual detector for the ray going thru pt(x,y,z)
		  eta = Vx[1] * con;
		  
		  //use bilinear interpolation to get projection value
		  lo_p = (int)(p + 1.0) - 1;
		  lo_eta = (int)(eta + 1.0) - 1;
		  if(lo_p >= 0 && lo_p < nSignalXPts-1 && lo_eta >= 0 && lo_eta < nSignalYPts-1)
		    {
		      frac_p = p - (double)lo_p;
		      frac_eta = eta - (double)lo_eta;
		      
		      ind = lo_eta*nSignalXPts + lo_p;
		      
		      p1 = pProj[ind] + frac_p * deltaProj[ind];
		      
		      ind += nSignalXPts;
		      p2 =  pProj[ind] + frac_p * deltaProj[ind];
		      
		      *pVox += (p1 + frac_eta * (p2 - p1)) * (con * con);
		      
		    }

		  pVox++;
			
		  Vx[0] += P[0][0];
		  Vx[1] += P[1][0];
		  Vx[2] += P[2][0];
		}
	      Vy[0] += P[0][1];
	      Vy[1] += P[1][1];
	      Vy[2] += P[2][1];
	    }
	  Vz[0] += P[0][2];
	  Vz[1] += P[1][2];
	  Vz[2] += P[2][2];
	}
      
    }
	
  //normalize image
  double norm = dBeta;
  norm *= dSignalXInc; //ds
  norm *= dFocalLength2_pix;  //
  pVox = vox;
  for(i = 0; i < pixels * lines * numSlices; i++)
    {
      *pVox *= norm;
      pVox++;
    }
  
  iftFree(H);
  iftFree(R);
  iftFree(deltaProj);
  
  return(1);
}

