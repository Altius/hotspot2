#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <string>

using namespace std;

struct Position {
  long x;  // genomic coordinate
  float y; // score (e.g. normalized density)
  int idx; // index into a vector of (contiguous) Positions
};

// As the name implies, Range is a generic range of contiguous Positions,
// bounded on the left by "beg" and on the right by "end."
// (Each Range is inclusive of "beg" and "end," i.e.,
// a 0-based Range gets output as (beg-1, end].)
// In practice, a Range is used to store boundaries of local maxima
// and peaks (which begin life as local maxima).
struct Range {
  string chrom;
  string id;
  Position beg;
  Position end;
  float maxY;
  long x_of_maxY;
  long inputSummit; // wavelet summit
  long summitSeparation;
};

// Each contiguous "island" of bp in the input file gets called a "Region."
// Each Region is made up of Positions, some (or all) of which define local maxima.
// A local maximum usually gets identified as a peak, but sometimes it's a shoulder of a peak.
struct Region {
  string chrom;
  string id;
  long posSummit; // wavelet summit
  vector<Position> posns;
  list<Range> localMaxima;
};

// Function for sorting Ranges in descending order by height.
// (The height is maxY, but when this function is actually called,
// all sites within each Range have the same height.)
// If two Ranges have the same height, they get sorted in genomic order.
bool Height_GT(const Range& a, const Range& b);
bool Height_GT(const Range& a, const Range& b)
{
  if (fabs(a.maxY - b.maxY) < 0.0001) // == on floating-point numbers
    {
      if (a.beg.x == b.beg.x)
	return a.end.x < b.end.x;
      return a.beg.x < b.beg.x;
    }
  return a.maxY > b.maxY;
}

// Function for sorting Ranges in genomic order.
bool GenomicOrder_LT(const Range& a, const Range& b);
bool GenomicOrder_LT(const Range& a, const Range& b)
{
  if (a.chrom == b.chrom)
    {
      if (a.beg.x == b.beg.x)
	{
	  if (a.end.x == b.end.x)
	    return a.inputSummit < b.inputSummit;
	  return a.end.x < b.end.x;
	}
      return a.beg.x < b.beg.x;
    }
  return a.chrom < b.chrom;
}

// This function finds the boundaries of the peak
// whose local maximum is the contiguous region bounded by [reg.posns[idxL], reg.posns[idxR]].
// A boundary point is the first site before the first site where the height dips below half-maximum,
// or the corresponding local minimum if half-maximum is not attained,
// or the end of the Region, whichever comes first.
// The slots within the vector of Positions, idxL and idxR,
// will correspond to the locations of those bounds when this function returns.
// The return value will be true unless the peak found is < minWidth bp wide or doesn't enclose the current posSummit.
bool getPeakBoundaries(const Region& r, int& idxL, int& idxR, const long& minWidth);
bool getPeakBoundaries(const Region& r, int& idxL, int& idxR, const long& minWidth)
{
  float maxY(r.posns[idxL].y), halfmaxY(r.posns[idxL].y/2.0), minY(maxY);
  const float secondaryCutoff(999999. * maxY); // ensures we always return the local min when we fail to achieve half-max
  const int maxIdx = static_cast<int>(r.posns.size()) - 1; // in English: "maximum index," not "idx of max"
  int idxOfLocalMin(idxL);

  if (idxL > 0)
    idxL--; // begin our search at the first site to the left of the local maximum
  while (idxL > 0 && (r.posns[idxL].y - halfmaxY) > 0.0001) // the L boundary of the Region has index = 0
    {
      if (r.posns[idxL].y < minY || r.posns[idxL].y - minY < 0.0001) // i.e., "if y <= minY"
	{
	  minY = r.posns[idxL].y;
	  idxOfLocalMin = idxL;
	}
      if (r.posns[idxL].y - maxY > -0.0001) // i.e., "if y >= maxY"
	{
	  // While searching for the half-max point to the left of the local max,
	  // we climbed above the local max, so this is a subset of another peak...
	  // but we want to give it a chance to get called as a peak, so we backtrack to the local min.
	  idxL = idxOfLocalMin;
	  if (secondaryCutoff - r.posns[idxL].y > 0.0001)
	    break;
	  else
	    return false; // we'll actually never hit this line unless we change the algorithm
	}
      else
	idxL--;
    }
  if (halfmaxY - r.posns[idxL].y > -0.0001) // i.e., "if halfmaxY >= y"
    idxL++;
  else
    {
      if (0 == idxL && r.posns[idxL].y > minY)
	idxL = idxOfLocalMin;
    }

  idxOfLocalMin = idxR;
  minY = maxY;
  if (idxR < maxIdx)
    idxR++; // begin our search at the first site to the right of the local maximum
  while (idxR < maxIdx && (r.posns[idxR].y - halfmaxY) > 0.0001) // the R boundary of the Region has index = maxIdx
    {
      if (r.posns[idxR].y < minY || r.posns[idxR].y - minY < 0.0001) // i.e., "if y <= minY"
	{
	  minY = r.posns[idxR].y;
	  idxOfLocalMin = idxR;
	}
      if (r.posns[idxR].y - maxY > - 0.0001) // i.e., "if y >= maxY"
	{
	  // While searching for the half-max point to the left of the local max,
	  // we climbed above the local max, so this is a subset of another peak...
	  // but we want to give it a chance to get called as a peak, so we backtrack to the local min.
	  idxR = idxOfLocalMin;
	  if (secondaryCutoff - r.posns[idxR].y > 0.0001)
	    break;
	  else
	    return false; // we'll actually never hit this line unless we change the algorithm
	}
      else
	idxR++;
    }
  if (halfmaxY - r.posns[idxR].y > -0.0001) // i.e., "if halfmax >= y"
    idxR--;
  else
    {
      if (maxIdx == idxR && r.posns[idxR].y > minY)
	idxR = idxOfLocalMin;
    }

  // Require the peak to contain the posSummit received for this region.
  if (r.posSummit < r.posns[idxL].x || r.posSummit > r.posns[idxR].x)
    return false;
  if (idxR - idxL + 1 >= minWidth) // 0-based sites (100,120] are the 20 sites from 101-120; could have, e.g., idxL = 0 and idxR = 19 in this case
    return true;
  return false;
}

// This function finds full-width-at-half-maximum peaks for all (filtered) local maxima within the contiguous Region,
// If half-maximum is not attained, the local minimum is used in its place.
// After this function removes peaks overlapped by taller ones,
// it finds the peak whose summit is nearest to the input posSummit for the region,
// and appends it to the vector "peaks," which in practice will contain some overlaps that another function will resolve.
// A peak must be at least minWidth bp wide, and must contain the input posSummit,  in order to be reported.
// This function clears or "flushes" the contents of "reg" at the end.
void processRegionAndFlushIt(Region& reg, const long& minWidth, vector<Range>& peaks);
void processRegionAndFlushIt(Region& reg, const long& minWidth, vector<Range>& peaks)
{
  if (reg.localMaxima.empty())
    {
      if (!reg.posns.empty())
	reg.posns.clear();
      return;
    }

  reg.localMaxima.sort(Height_GT);
  // Local maxima are now sorted in descending order by y-value (height), if there are 2 or more of them.
  // Local maxima with equal y-values, if any, are sorted from L to R in genomic order.
  // We traverse the list of local maxima in this order, tracing the shape of each down to half of its maximum.
  list<Range>::iterator it1 = reg.localMaxima.begin(), it2;
  while (it1 != reg.localMaxima.end())
    {
      int idxL(it1->beg.idx), idxR(it1->end.idx);
 
      it1->maxY = it1->beg.y;
      it1->x_of_maxY = static_cast<int>(floor(0.5*static_cast<double>(it1->beg.x + it1->end.x))); // important to use double, since x can have many digits
      if (!getPeakBoundaries(reg, idxL, idxR, minWidth)) // idxL and idxR will be changed upon return
	{
	  // The criteria for being called a "peak" were not met by this local maximum:
	  // it's too narrow (narrower than midWidth), or it doesn't contain reg.posSummit.
	  it2 = it1;
	  it1++;
	  reg.localMaxima.erase(it2);
	  continue;
	}
      it1->beg = reg.posns[idxL];
      it1->end = reg.posns[idxR];
      it1->inputSummit = reg.posSummit;
      it1->summitSeparation = labs(it1->x_of_maxY - it1->inputSummit);
      it1++;
    }

  // We now have 0, 1, or more peaks that are at least minWidth bp wide
  // and contain reg.posSummit.
  if (reg.localMaxima.empty())
    {
      reg.posns.clear();
      return;
    }

  long minSeparation = reg.localMaxima.begin()->summitSeparation;
  it1 = it2 = reg.localMaxima.begin();
  it2++;
  while (it2 != reg.localMaxima.end())
    {
      if (it2->summitSeparation < minSeparation)
	{
	  minSeparation = it2->summitSeparation;
	  it1 = it2;
	}
      it2++;
    }

  it1->chrom = reg.chrom;
  it1->id = reg.id;
  peaks.push_back(*it1);
  reg.posns.clear();
  reg.localMaxima.clear();
  return;
}

// This function resolves overlapping peaks by assigning the overlapping region
// to either the left or right peak in an overlapping pair, thereby making them
// adjacent rather than overlapping.
// It writes the results to stdout.
void cleanUpAnyOverlapsAndWriteOutput(vector<Range>& peaks, const long& minWidth);
void cleanUpAnyOverlapsAndWriteOutput(vector<Range>& peaks, const long& minWidth)
{
  if (peaks.empty())
    return;

  int idxL(0), idxR(1);
  bool triedSwapping(false);

  sort(peaks.begin(), peaks.end(), GenomicOrder_LT);

  while (idxR < static_cast<int>(peaks.size()))
    {
    TopOfLoop:
      if (peaks[idxR].beg.x <= peaks[idxL].end.x) // then these two peaks overlap
	{
	  if (peaks[idxL].x_of_maxY < peaks[idxR].beg.x) // L's summit is left of R's beg, so shorten the L peak
	    peaks[idxL].end.x = peaks[idxR].beg.x - 1;
	  else
	    {
	      if (peaks[idxR].x_of_maxY > peaks[idxL].end.x) // R's summit is right of L's end, so shorten the R peak
		peaks[idxR].beg.x = peaks[idxL].end.x + 1;
	      else // these two overlapping peaks are confusing; swap (L,R) --> (R,L) and try again
		{
		  if (!triedSwapping)
		    {
		      Range temp = peaks[idxR];
		      peaks[idxR] = peaks[idxL];
		      peaks[idxL] = temp;
		      triedSwapping = true;
		      goto TopOfLoop;
		    }
		  else
		    {
		      if (peaks[idxL].beg.x == peaks[idxR].beg.x && peaks[idxL].end.x == peaks[idxR].end.x)
			{
			  // Two adjacent wavelet peaks yielded the same FWHM peak.
			  // Delete the one whose wavelet summit (coordinate, x-value) is farther from the FWHM summit.
			  // Do this by enforcing the idxL version to be the one we keep
			  // and deleting the idxR version.
			  long distL(labs(peaks[idxL].inputSummit - peaks[idxL].x_of_maxY)),
			    distR(labs(peaks[idxR].inputSummit - peaks[idxR].x_of_maxY));
			  if (distL > distR)
			    peaks[idxL] = peaks[idxR];
			  peaks.erase(peaks.begin() + idxR);
			}
		      else
			{
			  cerr << peaks[idxL].chrom << ':' << peaks[idxL].beg.x-1 << '-' << peaks[idxL].end.x
			       << ", waveletSummit = " << peaks[idxL].inputSummit
			       << ", FWHM summit = " << peaks[idxL].x_of_maxY
			       << ",\nunsure how to resolve overlap with "
			       << peaks[idxR].chrom << ':' << peaks[idxR].beg.x-1 << '-' << peaks[idxR].end.x
			       << ", waveletSummit = " << peaks[idxR].inputSummit
			       << ", FWHM summit = " << peaks[idxR].x_of_maxY << '.' << endl << endl;
			  exit(2);
			}
		    }
		}
	    }
	}
      // We may have shortened peaks[idxL], so ensure it's sufficiently wide before reporting it.
      if (peaks[idxL].end.x - (peaks[idxL].beg.x - 1) >= minWidth)
	cout << peaks[idxL].chrom << '\t' << peaks[idxL].beg.x - 1 << '\t'
	     << peaks[idxL].end.x << '\t' << peaks[idxL].id << '\t'
	     << peaks[idxL].maxY << '\t' << peaks[idxL].x_of_maxY
	     << '\t' << peaks[idxL].inputSummit << endl;
      idxL++;
      idxR++;
      triedSwapping = false;
    }
  // Write the last element in the vector.
  if (peaks[idxL].end.x - (peaks[idxL].beg.x - 1) >= minWidth)
    cout << peaks[idxL].chrom << '\t' << peaks[idxL].beg.x - 1 << '\t'
	 << peaks[idxL].end.x << '\t' << peaks[idxL].id << '\t'
	 << peaks[idxL].maxY << '\t' << peaks[idxL].x_of_maxY
	 << '\t' << peaks[idxL].inputSummit << endl;
}


// This function parses the input and calls the functions
// that process it and write the output to stdout.
// pExeName is simply the name of the program;
// we include it in an error message if something goes wrong.
bool parseInputFindPeaksWriteOutput(const char* pExeName, const long& minWidth);
bool parseInputFindPeaksWriteOutput(const char* pExeName, const long& minWidth)
{
  const int BUFSIZE(1000);
  char buf[BUFSIZE], *p;
  Region region;
  string curChrom;
  vector<Range> peaksPossiblyWithOverlaps;
  Position curPosn, prevPosn;
  Range curRange;
  short fieldnum;
  long begCoord, curPosSummit, prevPosSummit(0);
  bool ascending(true);
  long linenum(0);

  prevPosn.x = 0;
  region.chrom = string("xxxNONExxx");
  
  while (cin.getline(buf, BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      bool newChrom(false);
      if (!(p = strtok(buf,"\t")) || !*p)
	{
	MissingField:
	  cerr << "Error:  Failed to find required field "
	       << fieldnum << " on line " << linenum
	       << " of the input to program " << pExeName
	       << '.' << endl << endl;
	  return false;
	}
      if (0 != strcmp(p, region.chrom.c_str()))
	{
	  curChrom = string(p);
	  newChrom = true;
	}
      fieldnum++;
      if (!(p = strtok(NULL,"\t")))
	goto MissingField;
      begCoord = atol(p); // used only for sanity check below
      fieldnum++;
      if (!(p = strtok(NULL,"\t")))
	goto MissingField;
      curPosn.x = atol(p);
      if (curPosn.x != begCoord + 1)
	{
	  cerr << "Error:  Interval on line " << linenum
	       << " of the input to program " << pExeName << " is not 1bp wide as expected; "
	       << "interval is " << curChrom << ':' << begCoord
	       << '-' << curPosn.x << '.' << endl << endl;
	  return false;
	}
      fieldnum++;
      if (!(p = strtok(NULL,"\t")))
	goto MissingField;
      if (1 == linenum) // we currently assume/require the same id in all lines of input
	region.id = string(p);
      fieldnum++;
      if (!(p = strtok(NULL,"\t")))
	goto MissingField;
      curPosn.y = atof(p);
      fieldnum++;
      if (!(p = strtok(NULL,"\t")))
	goto MissingField;
      curPosSummit = atol(p);
      fieldnum++;
      if ((p = strtok(NULL,"\t")))
	{
	  cerr << "Error:  Expected exactly " << fieldnum - 1 << " columns of input for program "
	       << pExeName << ", but encountered at least " << fieldnum
	       << " on line " << linenum << '.' << endl << endl;
	  return false;
	}
      
      if (curPosSummit != prevPosSummit || (!newChrom && curPosn.x != prevPosn.x + 1))
	{
	  if (linenum > 1)
	    {
	      if (curPosSummit == prevPosSummit && curPosn.x != prevPosn.x + 1)
		{
		  cerr << "Error:  In the input to program " << pExeName << ", line " << linenum-1
		       << " contains " << region.chrom << ':'
		       << prevPosn.x-1 << '-' << prevPosn.x
		       << " for summit " << prevPosSummit
		       << ", but line " << linenum << " contains "
		       << curChrom << ':' << curPosn.x-1
		       << '-' << curPosn.x << " for the same summit.\n"
		       << "All entries for a given summit must be in ascending order by column 3."
		       << endl << endl;
		  return false;
		}
	      // Check whether curRange, at the right boundary of region,
	      // is a local maximum.
	      if (region.localMaxima.empty() || ascending)
		region.localMaxima.push_back(curRange);
	    }
	  processRegionAndFlushIt(region, minWidth, peaksPossiblyWithOverlaps); // this will simply return if region.localMaxima is empty
	  region.posSummit = curPosSummit;
	  if (newChrom)
	    region.chrom = curChrom;
	  curPosn.idx = 0;
	  curRange.beg = curPosn;
	  curRange.maxY = curPosn.y;
	  ascending = true; // if the next y value is lower, we want this one stored as a local maximum
	}
      else
	curPosn.idx++;
      region.posns.push_back(curPosn);
      if (fabs(curPosn.y - curRange.beg.y) < 0.0001)
	curRange.end = curPosn;
      else
	{
	  if (curRange.beg.y - curPosn.y > 0.0001)
	    {
	      if (ascending)
		region.localMaxima.push_back(curRange);
	      ascending = false;
	    }
	  else
	    ascending = true;
	  curRange.beg = curRange.end = curPosn;
	  curRange.maxY = curPosn.y;
	}
      prevPosn = curPosn;
      prevPosSummit = curPosSummit;
    }

  // Process the final region.
  // Check whether curRange, at the right boundary of region,
  // is a local maximum.
  if (region.localMaxima.empty() || curRange.beg.y - region.posns[curRange.beg.idx - 1].y > 0.0001)
    region.localMaxima.push_back(curRange);
  processRegionAndFlushIt(region, minWidth, peaksPossiblyWithOverlaps);

  cleanUpAnyOverlapsAndWriteOutput(peaksPossiblyWithOverlaps, minWidth);
  
  return true;
}

int main(int argc, const char* argv[])
{
  if (2 != argc)
    {
      cerr << "Usage:  [stdin] | " << argv[0] << " minWidth | [stdout],\n"
	   << "where minWidth is the minimum width (in bp) that a variable-width peak must have to be reported.\n"
	   << "6-column input is required, NOT in sort-bed order, but instead, sorted by column 6, then by column 3.\n"
	   << "The input columns are:  chr, pos-1, pos, ID, score, posSummit,\n"
	   << "where \"score\" is the score (e.g., normalized density) at site \"pos\" on chromosome \"chr\"."
	   << "Rows MUST be grouped by column 6 (i.e., sorted by posSummit);\n"
	   << "the sort order of posSummit can be anything (typicaly numeric, but alphabetic would also work).\n"
	   << "All rows with the same posSummit value (column 6) MUST be sorted in ascending order by pos (column 3).\n"
	   << "There will be 7 columns of output:\n"
	   << "chr, beg, end, ID, maxScore, posWithMaxScore, posSummit,\n"
	   << "where (beg, end] are the coordinates of each output variable-width peak (0-based),\n"
	   << "ID and posSummit echo their input values, posWithMaxScore is the central bp within (beg,end]\n"
	   << "where the maximum score was found, and maxScore is that score.\n"
	   << "**YOU MUST CALL sort-bed ON THE OUTPUT BEFORE USING IT IN bedmap/bedops OPERATIONS**."
	   << endl << endl;
      return -1;
    }

  ios_base::sync_with_stdio(false); // calling this static method in this way turns off checks, and can speed up I/O
  
  if (!parseInputFindPeaksWriteOutput(argv[0], atol(argv[1])))
    return -1;

  return 0;
}

