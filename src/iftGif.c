#include "iftGif.h"

#include "ift/core/io/Stream.h"


const int kGifTransIndex = 0;

void iftGifGetClosestPaletteColor(iftGifPalette* pPal, int r, int g, int b, int* bestInd, int* bestDiff,int treeRoot)
{

    // base case, reached the bottom of the tree
    if(treeRoot > (1<<pPal->bitDepth)-1)
    {
        int ind = treeRoot-(1<<pPal->bitDepth);
        if(ind == kGifTransIndex) return;

        // check whether this color is better than the current winner
        int r_err = r - ((int32_t)pPal->r[ind]);
        int g_err = g - ((int32_t)pPal->g[ind]);
        int b_err = b - ((int32_t)pPal->b[ind]);
        int diff = iftGifIAbs(r_err)+iftGifIAbs(g_err)+iftGifIAbs(b_err);

        if(diff < (*bestDiff))
        {
            (*bestInd) = ind;
            (*bestDiff) = diff;
        }

        return;
    }

    // take the appropriate color (r, g, or b) for this node of the k-d tree
    int comps[3]; comps[0] = r; comps[1] = g; comps[2] = b;
    int splitComp = comps[pPal->treeSplitElt[treeRoot]];

    int splitPos = pPal->treeSplit[treeRoot];
    if(splitPos > splitComp)
    {
        // check the left subtree
        iftGifGetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2);
        if( (*bestDiff) > splitPos - splitComp )
        {
            // cannot prove there's not a better value in the right subtree, check that too
            iftGifGetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2+1);
        }
    }
    else
    {
        iftGifGetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2+1);
        if( (*bestDiff) > splitComp - splitPos )
        {
            iftGifGetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2);
        }
    }
}

void iftGifSwapPixels(uint8_t* image, int pixA, int pixB)
{
    uint8_t rA = image[pixA*4];
    uint8_t gA = image[pixA*4+1];
    uint8_t bA = image[pixA*4+2];
    uint8_t aA = image[pixA*4+3];

    uint8_t rB = image[pixB*4];
    uint8_t gB = image[pixB*4+1];
    uint8_t bB = image[pixB*4+2];
    uint8_t aB = image[pixA*4+3];

    image[pixA*4] = rB;
    image[pixA*4+1] = gB;
    image[pixA*4+2] = bB;
    image[pixA*4+3] = aB;

    image[pixB*4] = rA;
    image[pixB*4+1] = gA;
    image[pixB*4+2] = bA;
    image[pixB*4+3] = aA;
}

int iftGifPartition(uint8_t* image, const int left, const int right, const int elt, int pivotIndex)
{
    const int pivotValue = image[(pivotIndex)*4+elt];
    iftGifSwapPixels(image, pivotIndex, right-1);
    int storeIndex = left;
    bool split = 0;
    for(int ii=left; ii<right-1; ++ii)
    {
        int arrayVal = image[ii*4+elt];
        if( arrayVal < pivotValue )
        {
            iftGifSwapPixels(image, ii, storeIndex);
            ++storeIndex;
        }
        else if( arrayVal == pivotValue )
        {
            if(split)
            {
                iftGifSwapPixels(image, ii, storeIndex);
                ++storeIndex;
            }
            split = !split;
        }
    }
    iftGifSwapPixels(image, storeIndex, right-1);
    return storeIndex;
}

// Perform an incomplete sort, finding all elements above and below the desired median
void iftGifPartitionByMedian(uint8_t* image, int left, int right, int com, int neededCenter)
{
    if(left < right-1)
    {
        int pivotIndex = left + (right-left)/2;

        pivotIndex = iftGifPartition(image, left, right, com, pivotIndex);

        // Only "sort" the section of the array that contains the median
        if(pivotIndex > neededCenter)
            iftGifPartitionByMedian(image, left, pivotIndex, com, neededCenter);

        if(pivotIndex < neededCenter)
            iftGifPartitionByMedian(image, pivotIndex+1, right, com, neededCenter);
    }
}

// Builds a palette by creating a balanced k-d tree of all pixels in the image
void iftGifSplitPalette(uint8_t* image, int numPixels, int firstElt, int lastElt, int splitElt, int splitDist, int treeNode, bool buildForDither, iftGifPalette* pal)
{
    if(lastElt <= firstElt || numPixels == 0)
        return;

    // base case, bottom of the tree
    if(lastElt == firstElt+1)
    {
        if(buildForDither)
        {
            // Dithering needs at least one color as dark as anything
            // in the image and at least one brightest color -
            // otherwise it builds up error and produces strange artifacts
            if( firstElt == 1 )
            {
                // special case: the darkest color in the image
                uint32_t r=255, g=255, b=255;
                for(int ii=0; ii<numPixels; ++ii)
                {
                    r = iftGifIMin(r, image[ii*4+0]);
                    g = iftGifIMin(g, image[ii*4+1]);
                    b = iftGifIMin(b, image[ii*4+2]);
                }

                pal->r[firstElt] = r;
                pal->g[firstElt] = g;
                pal->b[firstElt] = b;

                return;
            }

            if( firstElt == (1 << pal->bitDepth)-1 )
            {
                // special case: the lightest color in the image
                uint32_t r=0, g=0, b=0;
                for(int ii=0; ii<numPixels; ++ii)
                {
                    r = iftGifIMax(r, image[ii*4+0]);
                    g = iftGifIMax(g, image[ii*4+1]);
                    b = iftGifIMax(b, image[ii*4+2]);
                }

                pal->r[firstElt] = r;
                pal->g[firstElt] = g;
                pal->b[firstElt] = b;

                return;
            }
        }

        // otherwise, take the average of all colors in this subcube
        uint64_t r=0, g=0, b=0;
        for(int ii=0; ii<numPixels; ++ii)
        {
            r += image[ii*4+0];
            g += image[ii*4+1];
            b += image[ii*4+2];
        }

        r += numPixels / 2;  // round to nearest
        g += numPixels / 2;
        b += numPixels / 2;

        r /= numPixels;
        g /= numPixels;
        b /= numPixels;

        pal->r[firstElt] = (uint8_t)r;
        pal->g[firstElt] = (uint8_t)g;
        pal->b[firstElt] = (uint8_t)b;

        return;
    }

    // Find the axis with the largest range
    int minR = 255, maxR = 0;
    int minG = 255, maxG = 0;
    int minB = 255, maxB = 0;
    for(int ii=0; ii<numPixels; ++ii)
    {
        int r = image[ii*4+0];
        int g = image[ii*4+1];
        int b = image[ii*4+2];

        if(r > maxR) maxR = r;
        if(r < minR) minR = r;

        if(g > maxG) maxG = g;
        if(g < minG) minG = g;

        if(b > maxB) maxB = b;
        if(b < minB) minB = b;
    }

    int rRange = maxR - minR;
    int gRange = maxG - minG;
    int bRange = maxB - minB;

    // and split along that axis. (incidentally, this means this isn't a "proper" k-d tree but I don't know what else to call it)
    int splitCom = 1;
    if(bRange > gRange) splitCom = 2;
    if(rRange > bRange && rRange > gRange) splitCom = 0;

    int subPixelsA = numPixels * (splitElt - firstElt) / (lastElt - firstElt);
    int subPixelsB = numPixels-subPixelsA;

    iftGifPartitionByMedian(image, 0, numPixels, splitCom, subPixelsA);

    pal->treeSplitElt[treeNode] = splitCom;
    pal->treeSplit[treeNode] = image[subPixelsA*4+splitCom];

    iftGifSplitPalette(image,              subPixelsA, firstElt, splitElt, splitElt-splitDist, splitDist/2, treeNode*2,   buildForDither, pal);
    iftGifSplitPalette(image+subPixelsA*4, subPixelsB, splitElt, lastElt,  splitElt+splitDist, splitDist/2, treeNode*2+1, buildForDither, pal);
}

// Finds all pixels that have changed from the previous image and
// moves them to the fromt of th buffer.
// This allows us to build a palette optimized for the colors of the
// changed pixels only.
int iftGifPickChangedPixels( const uint8_t* lastFrame, uint8_t* frame, int numPixels )
{
    int numChanged = 0;
    uint8_t* writeIter = frame;

    for (int ii=0; ii<numPixels; ++ii)
    {
        if(lastFrame[0] != frame[0] ||
           lastFrame[1] != frame[1] ||
           lastFrame[2] != frame[2])
        {
            writeIter[0] = frame[0];
            writeIter[1] = frame[1];
            writeIter[2] = frame[2];
            ++numChanged;
            writeIter += 4;
        }
        lastFrame += 4;
        frame += 4;
    }

    return numChanged;
}


// Creates a palette by placing all the image pixels in a k-d tree and then averaging the blocks at the bottom.
// This is known as the "modified median split" technique
void iftGifMakePalette( const uint8_t* lastFrame, const uint8_t* nextFrame, uint32_t width, uint32_t height, int bitDepth, bool buildForDither, iftGifPalette* pPal )
{
    pPal->bitDepth = bitDepth;

    // SplitPalette is destructive (it sorts the pixels by color) so
    // we must create a copy of the image for it to destroy
    int imageSize = width*height*4*sizeof(uint8_t);
    uint8_t* destroyableImage = (uint8_t*)malloc(imageSize);
    memcpy(destroyableImage, nextFrame, imageSize);

    int numPixels = width*height;
    if(lastFrame)
        numPixels = iftGifPickChangedPixels(lastFrame, destroyableImage, numPixels);

    const int lastElt = 1 << bitDepth;
    const int splitElt = lastElt/2;
    const int splitDist = splitElt/2;

    iftGifSplitPalette(destroyableImage, numPixels, 1, lastElt, splitElt, splitDist, 1, buildForDither, pPal);

    free(destroyableImage);

    // add the bottom node for the transparency index
    pPal->treeSplit[1 << (bitDepth-1)] = 0;
    pPal->treeSplitElt[1 << (bitDepth-1)] = 0;

    pPal->r[0] = pPal->g[0] = pPal->b[0] = 0;
}

// Implements Floyd-Steinberg dithering, writes palette value to alpha
void iftGifDitherImage( const uint8_t* lastFrame, const uint8_t* nextFrame, uint8_t* outFrame, uint32_t width, uint32_t height, iftGifPalette* pPal )
{
    int numPixels = width*height;

    // quantPixels initially holds color*256 for all pixels
    // The extra 8 bits of precision allow for sub-single-color error values
    // to be propagated
    int32_t* quantPixels = (int32_t*)malloc(sizeof(int32_t)*numPixels*4);

    for( int ii=0; ii<numPixels*4; ++ii )
    {
        uint8_t pix = nextFrame[ii];
        int32_t pix16 = ((int32_t)pix)  * 256;
        quantPixels[ii] = pix16;
    }

    for( uint32_t yy=0; yy<height; ++yy )
    {
        for( uint32_t xx=0; xx<width; ++xx )
        {
            int32_t* nextPix = quantPixels + 4*(yy*width+xx);
            const uint8_t* lastPix = lastFrame? lastFrame + 4*(yy*width+xx) : NULL;

            // Compute the colors we want (rounding to nearest)
            int32_t rr = (nextPix[0] + 127) / 256;
            int32_t gg = (nextPix[1] + 127) / 256;
            int32_t bb = (nextPix[2] + 127) / 256;

            // if it happens that we want the color from last frame, then just write out
            // a transparent pixel
            if( lastFrame &&
                lastPix[0] == rr &&
                lastPix[1] == gg &&
                lastPix[2] == bb )
            {
                nextPix[0] = rr;
                nextPix[1] = gg;
                nextPix[2] = bb;
                nextPix[3] = kGifTransIndex;
                continue;
            }

            int32_t bestDiff = 1000000;
            int32_t bestInd = kGifTransIndex;

            // Search the palete
            int treeRoot = 1;
            iftGifGetClosestPaletteColor(pPal, rr, gg, bb, &bestInd, &bestDiff,treeRoot);

            // Write the result to the temp buffer
            int32_t r_err = nextPix[0] - ((int32_t)pPal->r[bestInd]) * 256;
            int32_t g_err = nextPix[1] - ((int32_t)pPal->g[bestInd]) * 256;
            int32_t b_err = nextPix[2] - ((int32_t)pPal->b[bestInd]) * 256;

            nextPix[0] = pPal->r[bestInd];
            nextPix[1] = pPal->g[bestInd];
            nextPix[2] = pPal->b[bestInd];
            nextPix[3] = bestInd;

            // Propagate the error to the four adjacent locations
            // that we haven't touched yet
            int quantloc_7 = (yy*width+xx+1);
            int quantloc_3 = (yy*width+width+xx-1);
            int quantloc_5 = (yy*width+width+xx);
            int quantloc_1 = (yy*width+width+xx+1);

            if(quantloc_7 < numPixels)
            {
                int32_t* pix7 = quantPixels+4*quantloc_7;
                pix7[0] += iftGifIMax( -pix7[0], r_err * 7 / 16 );
                pix7[1] += iftGifIMax( -pix7[1], g_err * 7 / 16 );
                pix7[2] += iftGifIMax( -pix7[2], b_err * 7 / 16 );
            }

            if(quantloc_3 < numPixels)
            {
                int32_t* pix3 = quantPixels+4*quantloc_3;
                pix3[0] += iftGifIMax( -pix3[0], r_err * 3 / 16 );
                pix3[1] += iftGifIMax( -pix3[1], g_err * 3 / 16 );
                pix3[2] += iftGifIMax( -pix3[2], b_err * 3 / 16 );
            }

            if(quantloc_5 < numPixels)
            {
                int32_t* pix5 = quantPixels+4*quantloc_5;
                pix5[0] += iftGifIMax( -pix5[0], r_err * 5 / 16 );
                pix5[1] += iftGifIMax( -pix5[1], g_err * 5 / 16 );
                pix5[2] += iftGifIMax( -pix5[2], b_err * 5 / 16 );
            }

            if(quantloc_1 < numPixels)
            {
                int32_t* pix1 = quantPixels+4*quantloc_1;
                pix1[0] += iftGifIMax( -pix1[0], r_err / 16 );
                pix1[1] += iftGifIMax( -pix1[1], g_err / 16 );
                pix1[2] += iftGifIMax( -pix1[2], b_err / 16 );
            }
        }
    }

    // Copy the palettized result to the output buffer
    for( int ii=0; ii<numPixels*4; ++ii )
    {
        outFrame[ii] = quantPixels[ii];
    }

    free(quantPixels);
}

// Picks palette colors for the image using simple thresholding, no dithering
void iftGifThresholdImage( const uint8_t* lastFrame, const uint8_t* nextFrame, uint8_t* outFrame, uint32_t width, uint32_t height, iftGifPalette* pPal )
{
    uint32_t numPixels = width*height;
    for( uint32_t ii=0; ii<numPixels; ++ii )
    {
        // if a previous color is available, and it matches the current color,
        // set the pixel to transparent
        if(lastFrame &&
           lastFrame[0] == nextFrame[0] &&
           lastFrame[1] == nextFrame[1] &&
           lastFrame[2] == nextFrame[2])
        {
            outFrame[0] = lastFrame[0];
            outFrame[1] = lastFrame[1];
            outFrame[2] = lastFrame[2];
            outFrame[3] = kGifTransIndex;
        }
        else
        {
            // palettize the pixel
            int32_t bestDiff = 1000000;
            int32_t bestInd = 1;
            int treeRoot =1;
            iftGifGetClosestPaletteColor(pPal, nextFrame[0], nextFrame[1], nextFrame[2], &bestInd, &bestDiff,treeRoot);

            // Write the resulting color to the output buffer
            outFrame[0] = pPal->r[bestInd];
            outFrame[1] = pPal->g[bestInd];
            outFrame[2] = pPal->b[bestInd];
            outFrame[3] = bestInd;
        }

        if(lastFrame) lastFrame += 4;
        outFrame += 4;
        nextFrame += 4;
    }
}

// insert a single bit
void iftGifWriteBit( iftGifBitStatus* stat, uint32_t bit )
{
    bit = bit & 1;
    bit = bit << stat->bitIndex;
    stat->byte |= bit;

    ++stat->bitIndex;
    if( stat->bitIndex > 7 )
    {
        // move the newly-finished byte to the chunk buffer
        stat->chunk[stat->chunkIndex++] = stat->byte;
        // and start a new byte
        stat->bitIndex = 0;
        stat->byte = 0;
    }
}

// write all bytes so far to the file
void iftGifWriteChunk( FILE* f, iftGifBitStatus* stat )
{
    fputc(stat->chunkIndex, f);
    fwrite(stat->chunk, 1, stat->chunkIndex, f);

    stat->bitIndex = 0;
    stat->byte = 0;
    stat->chunkIndex = 0;
}

void iftGifWriteCode( FILE* f, iftGifBitStatus* stat, uint32_t code, uint32_t length )
{
    for( uint32_t ii=0; ii<length; ++ii )
    {
        iftGifWriteBit(stat, code);
        code = code >> 1;

        if( stat->chunkIndex == 255 )
        {
            iftGifWriteChunk(f, stat);
        }
    }
}

// write a 256-color (8-bit) image palette to the file
void iftGifWritePalette( const iftGifPalette* pPal, FILE* f )
{
    fputc(0, f);  // first color: transparency
    fputc(0, f);
    fputc(0, f);

    for(int ii=1; ii<(1 << pPal->bitDepth); ++ii)
    {
        uint32_t r = pPal->r[ii];
        uint32_t g = pPal->g[ii];
        uint32_t b = pPal->b[ii];

        fputc(r, f);
        fputc(g, f);
        fputc(b, f);
    }
}

// write the image header, LZW-compress and write out the image
void iftGifWriteLzwImage(FILE* f, uint8_t* image, uint32_t left, uint32_t top,  uint32_t width, uint32_t height, uint32_t delay, iftGifPalette* pPal)
{
    // graphics control extension
    fputc(0x21, f);
    fputc(0xf9, f);
    fputc(0x04, f);
    fputc(0x05, f); // leave prev frame in place, this frame has transparency
    fputc(delay & 0xff, f);
    fputc((delay >> 8) & 0xff, f);
    fputc(kGifTransIndex, f); // transparent color index
    fputc(0, f);

    fputc(0x2c, f); // image descriptor block

    fputc(left & 0xff, f);           // corner of image in canvas space
    fputc((left >> 8) & 0xff, f);
    fputc(top & 0xff, f);
    fputc((top >> 8) & 0xff, f);

    fputc(width & 0xff, f);          // width and height of image
    fputc((width >> 8) & 0xff, f);
    fputc(height & 0xff, f);
    fputc((height >> 8) & 0xff, f);

    //fputc(0, f); // no local color table, no transparency
    //fputc(0x80, f); // no local color table, but transparency

    fputc(0x80 + pPal->bitDepth-1, f); // local color table present, 2 ^ bitDepth entries
    iftGifWritePalette(pPal, f);

    const int minCodeSize = pPal->bitDepth;
    const uint32_t clearCode = 1 << pPal->bitDepth;

    fputc(minCodeSize, f); // min code size 8 bits

    iftGifLzwNode* codetree = (iftGifLzwNode*)malloc(sizeof(iftGifLzwNode)*4096);

    memset(codetree, 0, sizeof(iftGifLzwNode)*4096);
    int32_t curCode = -1;
    uint32_t codeSize = minCodeSize+1;
    uint32_t maxCode = clearCode+1;

    iftGifBitStatus* stat = (iftGifBitStatus*)calloc(1,sizeof(iftGifBitStatus));
    stat->byte = 0;
    stat->bitIndex = 0;
    stat->chunkIndex = 0;

    iftGifWriteCode(f, stat, clearCode, codeSize);  // start with a fresh LZW dictionary

    for(uint32_t yy=0; yy<height; ++yy)
    {
        for(uint32_t xx=0; xx<width; ++xx)
        {
            uint8_t nextValue = image[(yy*width+xx)*4+3];

            // "loser mode" - no compression, every single code is followed immediately by a clear
            //WriteCode( f, stat, nextValue, codeSize );
            //WriteCode( f, stat, 256, codeSize );

            if( curCode < 0 )
            {
                // first value in a new run
                curCode = nextValue;
            }
            else if( codetree[curCode].m_next[nextValue] )
            {
                // current run already in the dictionary
                curCode = codetree[curCode].m_next[nextValue];
            }
            else
            {
                // finish the current run, write a code
                iftGifWriteCode( f, stat, curCode, codeSize );

                // insert the new run into the dictionary
                codetree[curCode].m_next[nextValue] = ++maxCode;

                if( maxCode >= (1ul << codeSize) )
                {
                    // dictionary entry count has broken a size barrier,
                    // we need more bits for codes
                    codeSize++;
                }
                if( maxCode == 4095 )
                {
                    // the dictionary is full, clear it out and begin anew
                    iftGifWriteCode(f, stat, clearCode, codeSize); // clear tree

                    memset(codetree, 0, sizeof(iftGifLzwNode)*4096);
                    curCode = -1;
                    codeSize = minCodeSize+1;
                    maxCode = clearCode+1;
                }

                curCode = nextValue;
            }
        }
    }

    // compression footer
    iftGifWriteCode( f, stat, curCode, codeSize );
    iftGifWriteCode( f, stat, clearCode, codeSize );
    iftGifWriteCode( f, stat, clearCode+1, minCodeSize+1 );

    // write out the last partial chunk
    while( stat->bitIndex ) iftGifWriteBit(stat, 0);
    if( stat->chunkIndex ) iftGifWriteChunk(f, stat);

    fputc(0, f); // image block terminator
    free(codetree);
    free(stat);
}

// Creates a gif file.
// The input GIFWriter is assumed to be uninitialized.
// The delay value is the time between frames in hundredths of a second - note that not all viewers pay much attention to this value.
bool iftGifBegin( iftGifWriter* writer, const char* filename, uint32_t width, uint32_t height, uint32_t delay, int32_t bitDepth, bool dither)
{
#if _MSC_VER >= 1400
    writer->f = 0;
    fopen_s(&writer->f, filename, "wb");
#else
    writer->f = fopen(filename, "wb");
#endif
    if(!writer->f) return false;

    writer->firstFrame = true;

    // allocate
    writer->oldImage = (uint8_t*)malloc(width*height*4);

    fputs("GIF89a", writer->f);

    // screen descriptor
    fputc(width & 0xff, writer->f);
    fputc((width >> 8) & 0xff, writer->f);
    fputc(height & 0xff, writer->f);
    fputc((height >> 8) & 0xff, writer->f);

    fputc(0xf0, writer->f);  // there is an unsorted global color table of 2 entries
    fputc(0, writer->f);     // background color
    fputc(0, writer->f);     // pixels are square (we need to specify this because it's 1989)

    // now the "global" palette (really just a dummy palette)
    // color 0: black
    fputc(0, writer->f);
    fputc(0, writer->f);
    fputc(0, writer->f);
    // color 1: also black
    fputc(0, writer->f);
    fputc(0, writer->f);
    fputc(0, writer->f);

    if( delay != 0 )
    {
        // animation header
        fputc(0x21, writer->f); // extension
        fputc(0xff, writer->f); // application specific
        fputc(11, writer->f); // length 11
        fputs("NETSCAPE2.0", writer->f); // yes, really
        fputc(3, writer->f); // 3 bytes of NETSCAPE2.0 data

        fputc(1, writer->f); // JUST BECAUSE
        fputc(0, writer->f); // loop infinitely (byte 0)
        fputc(0, writer->f); // loop infinitely (byte 1)

        fputc(0, writer->f); // block terminator
    }

    return true;
}

// Writes out a new frame to a GIF in progress.
// The GIFWriter should have been created by GIFBegin.
// AFAIK, it is legal to use different bit depths for different frames of an image -
// this may be handy to save bits in animations that don't change much.
bool iftGifWriteFrame( iftGifWriter* writer, const uint8_t* image, uint32_t width, uint32_t height, uint32_t delay, int bitDepth, bool dither)
{
    if(!writer->f) return false;

    const uint8_t* oldImage = writer->firstFrame? NULL : writer->oldImage;
    writer->firstFrame = false;

    iftGifPalette pal;
    iftGifMakePalette((dither? NULL : oldImage), image, width, height, bitDepth, dither, &pal);

    if(dither)
        iftGifDitherImage(oldImage, image, writer->oldImage, width, height, &pal);
    else
        iftGifThresholdImage(oldImage, image, writer->oldImage, width, height, &pal);

    iftGifWriteLzwImage(writer->f, writer->oldImage, 0, 0, width, height, delay, &pal);

    return true;
}

// Writes the EOF code, closes the file handle, and frees temp memory used by a GIF.
// Many if not most viewers will still display a GIF properly if the EOF code is missing,
// but it's still a good idea to write it out.
bool iftGifEnd( iftGifWriter* writer )
{
    if(!writer->f) return false;

    fputc(0x3b, writer->f); // end of file
    fclose(writer->f);
    free(writer->oldImage);

    writer->f = NULL;
    writer->oldImage = NULL;

    return true;
}

void iftWriteGif(iftImage *frames, int alpha, int delay, const char* filename){
    iftGifWriter* gifWriter = NULL;
    uint8_t* imageVec = NULL;

    float maximumIntensity = iftMaximumValue(frames);
    if(frames->Cb == NULL){
        //grayscale
        iftSetCbCr(frames,(iftMaxImageRange(iftImageDepth(frames))+1)/2);
    }

    if(maximumIntensity > 255){
        int value = 256;
        while(value-1 < maximumIntensity){
            value <<= 1;
        }
        for (int z = 0; z < frames->zsize; ++z) {
            for (int y = 0; y < frames->ysize; ++y) {
                for (int x = 0; x < frames->xsize; ++x) {
                    iftImgVal(frames,x,y,z) = (((float)iftImgVal(frames,x,y,z))/(value-1.0f))*255;
                }
            }
        }
    }
    int channelDepth = 8;
    bool useDither = false;
    gifWriter = (iftGifWriter*)calloc(1,sizeof(iftGifWriter));
    iftGifBegin(gifWriter,filename, frames->xsize,frames->ysize,delay,channelDepth,useDither);
    for (int z = 0; z < frames->zsize; ++z) {
        imageVec = convertImageYcbCr2IntergerArray8bits(frames,alpha,z);
        iftGifWriteFrame(gifWriter,imageVec,frames->xsize,frames->ysize,delay,channelDepth,useDither);
        free(imageVec);
    }
    iftGifEnd(gifWriter);
}

//uint8_t* convertImageYcbCr2IntergerArray8bits(iftImage* image, int alpha, int frameIndex){
//
//    uint8_t* output = NULL;
//    if (image == NULL){
//        return NULL;
//    }
//
//    //RGBA
//    output = (uint8_t*)calloc(4*image->xsize*image->ysize ,sizeof(uint8_t));
//    long index = 0;
//    iftColor RGBA;
//    iftColor yCbCrA;
//
//    index = 0;
//    for (int y = 0; y < image->ysize; ++y) {
//        for (int x = 0; x < image->xsize; ++x) {
//            yCbCrA.val[0] = iftImgVal(image,x,y,frameIndex);
//            yCbCrA.val[1] = iftImgCb(image,x,y,frameIndex);
//            yCbCrA.val[2] = iftImgCr(image,x,y,frameIndex);
//            yCbCrA.alpha = alpha;
//            RGBA = iftYCbCrtoRGB(yCbCrA,255);
//            output[index] = RGBA.val[0];
//            index++;
//            output[index] = RGBA.val[1];
//            index++;
//            output[index] = RGBA.val[2];
//            index++;
//            output[index] = RGBA.alpha;
//            index++;
//        }
//    }
//
//    return output;
//}

int iftGifIMax(int l, int r) {
    return l>r?l:r;
}

int iftGifIMin(int l, int r) {
    return l<r?l:r;
}

int iftGifIAbs(int i) {
    return i<0?-i:i;
}
