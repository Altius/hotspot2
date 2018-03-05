// To compile this code into an executable,
// simply enter the command
//
// $ g++ -O3 hotspot2_part1.cpp -o hotspot2_part1
//
// or substitute any desired name for the executable for the last argument.
// The argument -O3 (capital "oh") generates optimized code;
// it can be omitted if desired.
// Any C++ compiler can be used in place of g++.
//
#include "hotspot2_version.h" // for versioning
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime> // for seeding the random number generator
#include <deque>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits> // for epsilon()
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility> // for pair
#include <vector>

using namespace std;

map<string, string*> interned;
map<string*, int> chromAsInt;
const long double CHANGE_OF_SCALE(1000.);

string* intern(string s)
{
  map<string, string*>::iterator it = interned.find(s);
  if (it == interned.end())
    {
      interned[s] = new string(s);
      chromAsInt[interned[s]] = static_cast<int>(interned.size());
      return interned[s];
    }
  return it->second;
}

int idxFromChrom(string* ptr)
{
  map<string*, int>::const_iterator it = chromAsInt.find(ptr);
  if (chromAsInt.end() == it)
    {
      cerr << "Coding error:  Line " << __LINE__ << ", failed to find \""
	   << *ptr << "\" in the lookup table." << endl << endl;
      exit(2);
    }
  return it->second;
}

struct SiteRange {
  string* chrom;
  string* ID;
  long double pval;
  //  long double qval;
  long begPos;
  long endPos;
  int count;
  int negLog10P_scaled; // -log10(P-value), scaled and rounded to limit the number of sig. digits
  bool hasPval; // whether a P-value has been computed for the site range
#ifdef DEBUG
  bool sampled;
#endif
};

long double nextProbNegativeBinomial(const int& k, const long double& prevVal, const vector<long double>& params);
long double nextProbNegativeBinomial(const int& k, const long double& prevVal, const vector<long double>& params)
{
  if (params.size() < 2)
    {
      cerr << "Error:  nextProbNegativeBinomial() received an incorrect parameter vector; expected m, r."
           << endl
           << endl;
      exit(1);
    }
  if (0 == k)
    {
      cerr << "Error:  nextProbNegativeBinomial() received k = 0, which is invalid (require k > 0)."
           << endl
           << endl;
      exit(1);
    }
  const long double &m(params[0]), &r(params[1]), kk(static_cast<long double>(k));
  if (r + m == 0)
    {
      cerr << "Error:  nextProbNegativeBinomial() received m = " << m
           << " and r = " << r << ", which is invalid (require r+m > 0)."
           << endl
           << endl;
      exit(1);
    }
  return prevVal * m * (r + kk - 1.) / ((r + m) * kk);
}

long double nextProbBinomial(const int& k, const long double& prevVal, const vector<long double>& params);
long double nextProbBinomial(const int& k, const long double& prevVal, const vector<long double>& params)
{
  if (params.size() < 3)
    {
      cerr << "Error:  nextProbBinomial() received an incorrect parameter vector; expected m, v, n."
           << endl
           << endl;
      exit(1);
    }
  if (0 == k)
    {
      cerr << "Error:  nextProbBinomial() received k = 0, which is invalid (require k > 0)."
           << endl
           << endl;
      exit(1);
    }
  const long double &m(params[0]), &v(params[1]), &nn(params[2]), kk(static_cast<long double>(k));
  if (kk > nn + 0.5) // really "if kk > nn," but we could have "kk == nn" with nn=15.99999 and kk=16.000001, hence the 0.5
    return 0.;
  if (v < 1.0e-8)
    {
      cerr << "Error:  nextProbBinomial() received v = 0, which is invalid (require variance > 0)."
           << endl
           << endl;
      exit(1);
    }
  return prevVal * (m - v) * (nn + 1. - kk) / (v * kk);
}

long double nextProbPoisson(const int& k, const long double& prevVal, const vector<long double>& params);
long double nextProbPoisson(const int& k, const long double& prevVal, const vector<long double>& params)
{
  if (params.empty())
    {
      cerr << "Error:  nextProbPoisson() received an incorrect parameter vector; expected m."
           << endl
           << endl;
      exit(1);
    }
  if (0 == k)
    {
      cerr << "Error:  nextProbPoisson() received k = 0, which is invalid (require k > 0)."
           << endl
           << endl;
      exit(1);
    }
  const long double &m(params[0]), kk(static_cast<long double>(k));
  return prevVal * m / kk;
}

class SiteManager {
public:
  SiteManager(ofstream& ofsJustPvals) : m_ofsJustNegLog10PscaledAndNumOccs(ofsJustPvals) {};
  void addSite(const SiteRange& s);
  void processPvalue(const long double& pval
#ifdef DEBUG
		     , const bool& sampled
#endif
		     );
  void writeLastUnreportedSite();

private:
  SiteManager(); // require the above constructor to be used
  SiteManager(const SiteManager&); // deny use of the copy constructor
  //  void initialize(ofstream& ofsJustPvals);
  deque<SiteRange> m_sites;
  ofstream& m_ofsJustNegLog10PscaledAndNumOccs;
};

void SiteManager::addSite(const SiteRange& s)
{
  m_sites.push_back(s);
}

inline void SiteManager::writeLastUnreportedSite()
{
  if (!m_sites.empty())
    {
      cout << idxFromChrom(m_sites.front().chrom) << '\t'
	   << m_sites.front().begPos << '\t'
	   << m_sites.front().endPos - m_sites.front().begPos << '\t'
	   << m_sites.front().negLog10P_scaled;
#ifdef DEBUG
      cout << m_sites.front().sampled;
#endif
      cout << '\n';
      m_ofsJustNegLog10PscaledAndNumOccs << m_sites.front().negLog10P_scaled << '\t'
					 << m_sites.front().endPos - m_sites.front().begPos << '\n';
    }
}


void SiteManager::processPvalue(const long double& pval
#ifdef DEBUG
				, const bool& sampled
#endif
				)
{
  //  static const int MAX_VALUE = static_cast<int>(floor(-log10(numeric_limits<double>::min()) * CHANGE_OF_SCALE + 0.5));
  int negLog10P_scaled;
  if (pval < 0 || pval > 1)
    negLog10P_scaled = 0;
  else
    {
      //      if (pval < numeric_limits<double>::min())
      //	negLog10P_scaled = MAX_VALUE;
      //      else
      negLog10P_scaled = static_cast<int>(floor(-log10(pval) * CHANGE_OF_SCALE + 0.5));
    }
  deque<SiteRange>::iterator itCurSiteNeedingPval(m_sites.begin());
  while (itCurSiteNeedingPval != m_sites.end() && itCurSiteNeedingPval->hasPval)
    itCurSiteNeedingPval++;
  if (m_sites.end() == itCurSiteNeedingPval)
    {
      cerr << "Error:  line " << __LINE__ << ", m_sites is empty or already filled with P-values"
	   << endl << endl;
      exit(2);
    }

  itCurSiteNeedingPval->pval = pval;
  itCurSiteNeedingPval->negLog10P_scaled = negLog10P_scaled;
  itCurSiteNeedingPval->hasPval = true;
#ifdef DEBUG
  itCurSiteNeedingPval->sampled = sampled;
#endif

  if (itCurSiteNeedingPval != m_sites.begin())
    {
      deque<SiteRange>::iterator it_prev = itCurSiteNeedingPval;
      it_prev--;
      if (it_prev->chrom == itCurSiteNeedingPval->chrom && it_prev->endPos + 1 == itCurSiteNeedingPval->endPos &&
#ifdef DEBUG
	  it_prev->sampled == itCurSiteNeedingPval->sampled &&
#endif
	  it_prev->negLog10P_scaled == itCurSiteNeedingPval->negLog10P_scaled)
	itCurSiteNeedingPval->begPos = it_prev->begPos;
      else
	writeLastUnreportedSite();
      m_sites.pop_front();
    }
}

struct SiteData {
  int pos;
  int count;
  bool hasPval;
  bool sampled;
};

struct StatsForCount {
  int numOccs; // number of occurrences
  long double pmf; // probability mass function
  long double pval; // P-value, probability of observing a count this large or larger
  int MAxN; // moving average of number-of-occurrences, but not divided by N
};

class BackgroundRegionManager {
public:
  BackgroundRegionManager(const int& samplingInterval, const int& MAlength);
  void setBounds(const string* pChrom, const int posL, const int posR);
  const int& getRightEdge(void) const { return m_posR; };
  const bool& isSliding(void) const { return m_sliding; };
  void add(const SiteRange& s);
  void computePandFlush(SiteManager& sm);
  void slideAndCompute(const SiteRange& s, SiteManager& sm);

private:
  BackgroundRegionManager(void); // require use of the constructor with 2 arguments
  BackgroundRegionManager(const BackgroundRegionManager&); // ditto
  void findCutoff(void);
  void computeStats(const int& this_k);
  long double getPvalue(const unsigned int& k);
  int m_posL;
  int m_posR;
  int m_posC;
  int m_runningSum_count;
  int m_runningSum_countSquared;
  int m_numPtsInNullRegion;
  int m_runningSum_count_duringPrevComputation;
  int m_runningSum_countSquared_duringPrevComputation;
  int m_numPtsInNullRegion_duringPrevComputation;

  int m_MAlength;
  long double m_thresholdRatio;
  vector<StatsForCount> m_distn; // distribution of observed counts
  deque<SiteData> m_sitesInRegion_leftHalf; // endPos values of the sites in the region and whether P has been assigned
  deque<SiteData> m_sitesInRegion_rightHalf; // endPos values of the sites in the region and whether P has been assigned
  int m_modeXval;
  int m_modeYval;
  int m_kcutoff;
  int m_kTrendReversal;
  bool m_sliding;
  bool m_needToUpdate_kcutoff;
  set<int> m_kvalsWithMinMAxN;
  int m_minMAxN;
  int m_prev_k;

  int m_samplingInterval;
  int m_nextPosToSample;
  int m_sampledDataDistnSize;

  long double (*m_pmf)(const int&, const long double&, const vector<long double>&); // Make these member variables, not local variables,
  vector<long double> m_pmfParams; // so they can be accessed outside computeStats() for debugging.
  const string* m_pCurChrom;
};

BackgroundRegionManager::BackgroundRegionManager(const int& samplingInterval, const int& MAlength)
{
  m_posL = m_posC = m_posR = -1;
  m_runningSum_count = m_runningSum_countSquared = m_numPtsInNullRegion = 0;
  m_runningSum_count_duringPrevComputation = m_runningSum_countSquared_duringPrevComputation = m_numPtsInNullRegion_duringPrevComputation = 0;
  m_kcutoff = m_modeYval = m_modeXval = -1;
  m_MAlength = MAlength; // see explanation in findCutoff(); 5 is good when samplingInterval = 1, 15 is good when windowSize/samplingInterval ~= 250
  m_thresholdRatio = 1.33; // see explanation in findCutoff(); could instead try 1.4. 1.5 seems to be too high, 1.2 seems to be too low.
  m_sliding = false;
  m_needToUpdate_kcutoff = true;

  m_samplingInterval = samplingInterval;
  m_nextPosToSample = -1;
  m_sampledDataDistnSize = 0;

  m_minMAxN = m_kTrendReversal = -1;
  m_prev_k = -1;
  m_pmf = NULL;
  m_pCurChrom = NULL;
}

void BackgroundRegionManager::setBounds(const string* pChrom, const int posL, const int posR)
{
  if (posR <= posL)
    {
      cerr << "Error:  BackgroundRegionManager::setBounds() received posL = "
           << posL << " and posR = " << posR << "; must have posL <= posR."
           << endl
           << endl;
      exit(1);
    }
  m_pCurChrom = pChrom; // used solely for warning and error messages and debugging
  m_posL = posL;
  m_posR = posR;
  m_posC = m_posL + (m_posR - m_posL) / 2; // integer division
  m_nextPosToSample = m_posL;
}

void BackgroundRegionManager::add(const SiteRange& s)
{
  if (s.endPos > m_posR)
    {
      cerr << "Coding error:  BRM::add(), region " << *m_pCurChrom << ':' << m_posL << '-' << m_posR
           << " received out-of-bounds position " << s.endPos
           << " (line " << __LINE__ << " of the code)." << endl
           << endl;
      exit(1);
    }
  SiteData sd;
  sd.pos = s.endPos;
  sd.count = s.count;
  sd.hasPval = false;
  sd.sampled = false;

  if (1 != m_samplingInterval && s.endPos > m_nextPosToSample)
    m_nextPosToSample = s.endPos;

  // Alternatively, "sampling" could be performed strictly at multiples of m_samplingInterval.
  // If we intend to sample at, say, positions 200, 400, 600, 800, 1000, etc.,
  // and regions 390-402 and 598-620 are excluded (unmappable),
  // then the current code will sample at 200, 403, 621, 821, etc.
  // Sampling strictly at multiples of m_samplingInterval would in this case
  // sample at 200, 800, 1000, etc.
  // Use the following code if strict multiples of m_samplingInterval are desired.
  //
  // m_nextPosToSample += static_cast<int>(ceil(static_cast<double>(s.endPos - m_nextPosToSample)/static_cast<double>(m_samplingInterval))) * m_samplingInterval;
  //
  // The above is a one-step version of "while(s.endPos > m_nextPosToSample){m_nextPosToSample += m_samplingInterval;}."
  // Typically, m_nextPosToSample will only need to increase by m_samplingInterval;
  // it will need to increase by a greater multiple of m_samplingInterval when a large gap is present in the data
  // (an unmappable region or otherwise restricted/withheld region).

  if (1 == m_samplingInterval || s.endPos == m_nextPosToSample)
    {
      sd.sampled = true;
      m_nextPosToSample += m_samplingInterval;
      // Add the incoming site's count to the distribution of counts observed in this region.
      if (s.count < static_cast<int>(m_distn.size()))
        m_distn[s.count].numOccs++;
      else
        {
          StatsForCount sc;
          sc.numOccs = 0;
          sc.pmf = sc.pval = -1.;
          sc.MAxN = -1; // don't compute moving averages until/unless we're sliding, for efficiency's sake
          while (static_cast<int>(m_distn.size()) < s.count)
            {
              m_distn.push_back(sc); // create bins for unobserved interior values, e.g., count = 5 but only 0,1,2 have been observed so far
              m_sampledDataDistnSize++;
            }
          sc.numOccs = 1;
          m_distn.push_back(sc);
          m_sampledDataDistnSize++;
        }
      if ( static_cast<int>(m_distn.size()) >= m_MAlength && m_distn[s.count].MAxN != -1 && m_distn[s.count].MAxN >= m_modeYval)
        {
          m_modeXval = s.count;
          m_modeYval = m_distn[s.count].MAxN;
        }
    }

  // Add the incoming site to the appropriate deque of observed positions.
  if (s.endPos < m_posC)
    m_sitesInRegion_leftHalf.push_back(sd);
  else
    m_sitesInRegion_rightHalf.push_back(sd);
}

void BackgroundRegionManager::findCutoff()
{
  // In all three noise models employed (negative binomial (NB), binomial, and Poisson,
  // though the NB is essentially always used because the data essentially always has variance > mean),
  // the distribution monotonically decreases beyond its mode.
  // Our idea is to identify a cutoff, when one exists,
  // which essentially marks the upper bound of observed monotonic decrease.
  // (Observations beyond the cutoff are then "predominantly signal," rather than noise.)
  // We do this by starting just beyond the mode of the
  // moving averages (MAs) of occurrence counts.
  // Whenever we encounter a global minimum (among values seen so far) in the MAs,
  // we measure whether a subsequent MA is substantially greater,
  // i.e. above some threshold (e.g., an increase of 33.3% or more).
  // If such an extreme increase is observed,
  // then we declare that "global minimum so far" to be the cutoff value.
  // If we exhaust the set of observations and no such increase is detected,
  // the highest observation is returned, i.e., no points are actually "cut off."
  // If a new "global minimum so far" is observed, below the previous "global minimum so far,"
  // the previous global minimum is ignored and the search process begins anew at the new local minimum.
  //
  // If no "trend reversal" is detected, but we detect a contiguous set of empty histogram bins,
  // we set the cutoff at the value immediately preceding those empty histogram bins.

  int kcutoff_uponEntry = m_kcutoff;

  if (m_sampledDataDistnSize < m_MAlength)
    {
      // Too few distinct count values were observed to compute a moving average of length MAlength;
      // set the "cutoff" to the largest observed count, i.e., use all data, don't "cut off" any points.
      m_kcutoff = m_sampledDataDistnSize - 1;
      m_minMAxN = m_kTrendReversal = -1;
      m_modeXval = m_modeYval = -1;
      m_kvalsWithMinMAxN.clear();
      m_needToUpdate_kcutoff = false;

      if (m_kcutoff != kcutoff_uponEntry)
        {
        UpdateNullRegionStatsAndExit:
          // The boundary of the null region has changed.
          // Recompute the running sums and the number of observations in the null region.
          if (kcutoff_uponEntry <= 0)
            {
              // Compute from scratch.
              m_runningSum_count = m_runningSum_countSquared = m_numPtsInNullRegion = 0;
              for (int kk = 0; kk <= m_kcutoff; kk++)
                {
                  m_runningSum_count += m_distn[kk].numOccs * kk;
                  m_runningSum_countSquared += m_distn[kk].numOccs * kk * kk;
                  m_numPtsInNullRegion += m_distn[kk].numOccs;
                }
            }
          else
            {
              if (m_kcutoff > kcutoff_uponEntry)
                {
                  // The null region has expanded; increase the values accordingly.
                  for (int kk = kcutoff_uponEntry + 1; kk <= m_kcutoff; kk++)
                    {
                      m_runningSum_count += m_distn[kk].numOccs * kk;
                      m_runningSum_countSquared += m_distn[kk].numOccs * kk * kk;
                      m_numPtsInNullRegion += m_distn[kk].numOccs;
                    }
                }
              else
                {
                  // The null region has contracted; decrease the values accordingly.
                  for (int kk = kcutoff_uponEntry; kk > m_kcutoff; kk--)
                    {
                      m_runningSum_count -= m_distn[kk].numOccs * kk;
                      m_runningSum_countSquared -= m_distn[kk].numOccs * kk * kk;
                      m_numPtsInNullRegion -= m_distn[kk].numOccs;
                    }
                }
            }
        }

      return;
    }

  int sum(0), idxL(m_modeXval + 1), idxR(m_modeXval + 1 + m_MAlength - 1); // these idxL and idxR values will be ignored/replaced if !m_sliding
  int idxC = m_modeXval + 1 + m_MAlength / 2; // integer division; see commendt directly above re: idxL and idxR
  pair<int, int> xyCurMAxN;
  bool useGlobMin(false);

  if (!m_sliding)
    {
      // Need to compute the moving averages (MAs).
      // Technically, these aren't moving averages, they're sums,
      // because we're not dividing sums by the number of terms.
      // But these "MA x N" values function as MAs.
      // Benefits:  No unnecessary division and casting from int to double,
      // and we gain the ability to perform exact tests for equality
      // (int == int instead of double == double).
      idxL = 0;
      idxC = idxL + m_MAlength / 2;
      idxR = idxL + m_MAlength - 1;
      // Compute the initial moving average (times the number of terms).
      for (int i = idxL; i <= idxR; i++)
        sum += m_distn[i].numOccs;
      m_distn[idxC].MAxN = sum;
      if (-1 == m_modeYval || m_distn[idxC].MAxN > m_modeYval)
        {
          m_modeXval = idxC;
          m_modeYval = m_distn[idxC].MAxN;
        }
      while (idxR < m_sampledDataDistnSize - 1)
        {
          idxC++;
          idxR++;
          sum -= m_distn[idxL].numOccs;
          sum += m_distn[idxR].numOccs;
          m_distn[idxC].MAxN = sum;
          if (m_distn[idxC].MAxN > m_modeYval)
            {
              m_modeXval = idxC;
              m_modeYval = m_distn[idxC].MAxN;
            }
          idxL++;
        }
      if (idxL > m_modeXval + 1)
        {
          idxL = m_modeXval + 1;
          idxC = idxL + m_MAlength / 2;
          idxR = idxL + m_MAlength - 1;
        }
      // else the "while" loop below won't get executed
    }
  m_kvalsWithMinMAxN.clear();
  if (idxR != m_sampledDataDistnSize - 1)
    {
      m_kvalsWithMinMAxN.insert(idxC); // idxC == m_modeXval+1 + m_MAlength/2 here, whether !m_sliding or m_sliding==true
      m_minMAxN = m_distn[idxC].MAxN;
    }
  long double minMAxN = static_cast<long double>(m_minMAxN);

  idxR++;
  while (idxR < m_sampledDataDistnSize)
    {
      idxC++;
      xyCurMAxN.first = idxC;
      xyCurMAxN.second = m_distn[idxC].MAxN;
      if (static_cast<long double>(xyCurMAxN.second) > m_thresholdRatio * minMAxN)
        {
          useGlobMin = true;
          break;
        }
      else
        {
          if (0 == xyCurMAxN.second)
            {
              // We've detected a contiguous stretch of at least m_MAlength empty histogram bins.
              // Set the global minimum here, and break out of the loop.
              // m_kcutoff will be set following the exit from the loop.
              m_kvalsWithMinMAxN.clear();
              m_kvalsWithMinMAxN.insert(xyCurMAxN.first); // k == idxC
              m_minMAxN = 0;
              useGlobMin = true;
              break;
            }
        }
      if (xyCurMAxN.second <= m_minMAxN)
        {
          if (xyCurMAxN.second < m_minMAxN)
            {
              // This is a unique, new minimum.
              m_kvalsWithMinMAxN.clear();
              m_minMAxN = xyCurMAxN.second;
              minMAxN = static_cast<long double>(m_minMAxN);
            } // else it's a duplicate occurrence of the existing minimum
          m_kvalsWithMinMAxN.insert(xyCurMAxN.first); // k == idxC
        }
      idxL++;
      idxR++;
    }

  // If useGlobMin is false, then we exhausted all observed values
  // without finding an "extreme" global minimum in the moving averages.
  // Return the maximum observed count as the "cutoff."
  if (!useGlobMin)
    {
      m_kcutoff = m_sampledDataDistnSize - 1;
      m_kTrendReversal = -1;
      m_kvalsWithMinMAxN.clear();
      m_minMAxN = -1;
      m_needToUpdate_kcutoff = false;
      if (m_kcutoff != kcutoff_uponEntry)
        goto UpdateNullRegionStatsAndExit;
      return;
    }
  // Otherwise, return the "global minimum so far" as the cutoff.
  // If there are ties (multiple k values with the same MAxN values), return the highest of these k values.
  set<int>::const_iterator it = m_kvalsWithMinMAxN.end();
  int k = *(--it);
  m_kcutoff = k;

  if (m_minMAxN > 0)
    m_kTrendReversal = xyCurMAxN.first; // where we determined the monotonically decreasing trend to have reversed
  else
    m_kTrendReversal = -1;

  m_needToUpdate_kcutoff = false;
  if (m_kcutoff != kcutoff_uponEntry)
    goto UpdateNullRegionStatsAndExit;
}

void BackgroundRegionManager::computeStats(const int& this_k)
{
  if (1 == m_sampledDataDistnSize)
    {
      // All observations within this region were of count == 0.
      m_distn[0].pmf = m_distn[0].pval = 1.;
      m_runningSum_count_duringPrevComputation = m_runningSum_count;
      m_runningSum_countSquared_duringPrevComputation = m_runningSum_countSquared;
      m_numPtsInNullRegion = m_distn[0].numOccs;
      m_numPtsInNullRegion_duringPrevComputation = m_numPtsInNullRegion;
      m_prev_k = 0;
      m_pmf = &nextProbPoisson;
      m_pmfParams.clear();
      m_pmfParams.push_back(0.);
      return;
    }

  int k;
  // Make *pmf() and params private member variables,
  // so they can be accessed elsewhere during debugging.
  long double prob0; // probability of observing 0 counts by random chance
  long double m; // mean
  long double v; // variance
  long double N(static_cast<long double>(m_numPtsInNullRegion));
  m = static_cast<long double>(m_runningSum_count) / N;
  v = (static_cast<long double>(m_runningSum_countSquared) - N * m * m) / (N - 1.);
  m_pmfParams.clear();
  static bool warningAlreadyIssued(false);

  // Set up the negative binomial model.
  // Double-check that m < v; if m >= v, which is extremely unlikely,
  // use a more appropriate model (binomial or Poisson).
  // Use the definitions of variance and mean to perform an integer (exact) test for m == v.
  if (m_numPtsInNullRegion * m_runningSum_countSquared - m_runningSum_count*m_runningSum_count == m_runningSum_count*(m_numPtsInNullRegion - 1))
    {
      // Poisson, m = v
      prob0 = exp(-m);
      m_pmfParams.push_back(m);
      m_pmf = &nextProbPoisson;
      if ((0 == m_runningSum_count || 1 == m_runningSum_count) && !warningAlreadyIssued)
        {
          cerr << "Warning:  In region " << *m_pCurChrom << ':' << m_posL << '-' << m_posR
               << ", all counts used for statistics were 0, or all were 0 except one was 1.\n"
               << "This generally should not happen.  If this region is unmappable or problematic for other reasons,\n"
               << "it would almost certainly be best to filter it out of the input.\n"
               << "There may be other such regions in the input; this warning will only be issued once during this run."
               << endl;
          warningAlreadyIssued = true;
        }
    }
  else
    {
      if (m < v)
        {
          // negative binomial
          long double r = m * m / (v - m);
          m_pmfParams.push_back(m);
          m_pmfParams.push_back(r);
          prob0 = pow(r / (r + m), r);
          m_pmf = &nextProbNegativeBinomial;
        }
      else // m > v (this is very unlikely)
        {
          // binomial
          long double n(floor(m * m / (m - v) + 0.5)); // estimate n of the fit from the observed mean and variance
          long double kk_max(static_cast<long double>(m_kcutoff)); // BUGBUG if m_kcutoff has 0 observations, do we still want to use it here??!?
          if (kk_max > n + 0.1) // + 0.1 to account for roundoff error, e.g. k=16.00001 and n=15.99999
            {
              // We don't expect this to happen, but if the estimated n is smaller than the observed m_kcutoff,
              // then we need to set n = m_kcutoff to allow m_kcutoff to be generated by the binomial model.
              n = kk_max;
            }
          // Now that we've set n, there's inconsistency among n, m, and v, with respect to the binomial.
          // We choose to keep m as observed, and update the variance parameter v so that consistency is achieved.
          v = m * (1. - m / n); // Now m*m/(m-v) = the integer n. Example: m=1.8952, v=1.6747, n=16.2893-->16, v-->1.6701.
          prob0 = pow(v / m, n);
          m_pmfParams.push_back(m);
          m_pmfParams.push_back(v);
          m_pmfParams.push_back(n);
          if (v > 1.0e-8)
            m_pmf = &nextProbBinomial;
          else
            {
              m_pmf = NULL;
              return;
            }
        }
    }

  // Now compute the minimum necessary number of pmf values, from the null model (e.g. negative binomial fit).

  // BUGBUG potential speed-up:  Allow the user to specify a minimum pmf,
  // such that no pmf values below it will be computed and
  // the P-values reported will asymptote at the P-value corresponding to that cutoff.
  // E.g., the user might specify a minimum pmf of numeric_limits<double>::epsilon().
  // If a P-value is, say, 1.23456e-218, do we really want to know that,
  // or are we content to know that it's something smaller than, say, 1e-40?

  long double curPMF(prob0);
  int k_begin(1), k_end(m_sampledDataDistnSize - 1); // default (-1 == this_k):  compute for all k observed in the background window
  const int MAlenOver2(m_MAlength / 2);
  if (-1 == this_k)
    m_distn[0].pmf = curPMF;
  else
    {
      k_end = this_k;
      if (-1 == m_prev_k) // compute for 0 <= k <= this_k.
        m_distn[0].pmf = curPMF;
      else // compute for (highest k previously handled) < k <= this_k.
        {
          curPMF = m_distn[m_prev_k].pmf;
          k_begin = m_prev_k + 1;
        }
    }
  for (k = k_begin; k <= k_end; k++)
    {
      curPMF = m_pmf(k, curPMF, m_pmfParams); // note:  if pmf == binomial and k > binomial's n, 0 is returned
      if (static_cast<int>(m_distn.size()) == k)
        {
          StatsForCount sc;
          sc.numOccs = 0;
          sc.pmf = sc.pval = -1.;
          sc.MAxN = -1; // don't compute moving averages until/unless we're sliding, for efficiency's sake
          m_distn.push_back(sc); // be sure not to increment m_sampledDataDistnSize...
          if (k >= m_MAlength)
            m_distn[k - MAlenOver2].MAxN = m_distn[k - MAlenOver2 - 1].MAxN - m_distn[k - m_MAlength].numOccs + m_distn[k].numOccs;
          else
            {
              if (m_MAlength - 1 == k)
                {
                  int sum(0);
                  for (int j = 0; j <= k; j++)
                    sum += m_distn[j].numOccs;
                  m_distn[MAlenOver2].MAxN = sum;
                }
            }
        }
      m_distn[k].pmf = curPMF;
    }

  // Because pmf(k) happens to be a multiple of pmf(k-1) for every k for each model (negative binomial, binomial, Poisson),
  // we can write the P-value for observing k counts as
  //
  // pval(k) = pmf(k) * (1 + term(k+1) + term(k+2) + term(k+3) + ...),
  //
  // where the terms in the sum decrease monotonically with k.
  // We can cap the sum at, say, 5 significant digits for pval(k),
  // and then work backwards, filling in pval(k-1) = pmf(k-1) + pval(k),
  // ending with pval(0) = 1.
  k--; // reset so that k corresponds to the last pmf computed
  int j = k;
  long double sum(1.), prevTerm(1.), curTerm;
  const long double SMALL_VALUE(5.0e-7); // restrict the correctness of the final P-value to ~5 significant digits
  while ((curTerm = m_pmf(++j, prevTerm, m_pmfParams)) > SMALL_VALUE)
    {
      sum += curTerm;
      prevTerm = curTerm;
    }
  m_distn[k].pval = m_distn[k].pmf * sum;
  //  if (m_distn[k].pval < numeric_limits<double>::min())
  //    m_distn[k].pval = numeric_limits<double>::min();
  while (k > 1)
    {
      m_distn[k - 1].pval = m_distn[k - 1].pmf + m_distn[k].pval;
      //      if (m_distn[k].pval < numeric_limits<double>::min())
      //	m_distn[k].pval = numeric_limits<double>::min();
      k--;
    }
  m_distn[0].pval = 1.; // explicitly set it to 1, to avoid potential round-off error

  if (-1 != this_k)
    m_prev_k = this_k;
  else
    m_prev_k = m_sampledDataDistnSize - 1;

  m_runningSum_count_duringPrevComputation = m_runningSum_count;
  m_runningSum_countSquared_duringPrevComputation = m_runningSum_countSquared;
  m_numPtsInNullRegion_duringPrevComputation = m_numPtsInNullRegion;
}

long double BackgroundRegionManager::getPvalue(const unsigned int& k)
{
  if (k < m_distn.size()) // Yes, m_distn.size(), not m_sampledDataDistnSize
    return m_distn[k].pval;

  // k is greater than all count values in the distribution, so we need to add bin(s) for it.
  // This can only happen if the user has specified a sampling interval
  // and k happens to have not been sampled.

  if (NULL == m_pmf || m_pmfParams.empty())
    {
      cerr << "Coding error:  BRM::getPvalue(" << k << ") was called when m_pmf == NULL and/or m_pmfParams.empty() == true." << endl;
      exit(1);
    }

  StatsForCount sc;
  sc.numOccs = 0;
  sc.pmf = sc.pval = -1.;
  sc.MAxN = -1;

  const int prev_max_k(static_cast<int>(m_distn.size()) - 1), MAlenOver2(m_MAlength / 2); // Yes, m_distn.size(), not m_sampledDataDistnSize.
  int kk = prev_max_k;
  long double curPMF(m_distn[kk].pmf);
  while (k >= m_distn.size()) // Yes, k, not kk.  We're growing the vector until k fits into its highest bin.
    {
      curPMF = m_pmf(++kk, curPMF, m_pmfParams);
      sc.pmf = curPMF;
      m_distn.push_back(sc);
      if (kk >= m_MAlength)
        m_distn[kk - MAlenOver2].MAxN = m_distn[kk - MAlenOver2 - 1].MAxN - m_distn[kk - m_MAlength].numOccs + m_distn[kk].numOccs;
      else
        {
          if (m_MAlength - 1 == kk)
            {
              int sum(0);
              for (int j = 0; j <= kk; j++)
                sum += m_distn[j].numOccs;
              m_distn[MAlenOver2].MAxN = sum;
            }
        }
    }
  // Compute the P-value.  See comments in method computeStats() for further info.
  // kk == k at this point.
  long double sum(1.), prevTerm(1.), curTerm;
  const long double SMALL_VALUE(5.0e-7); // restrict the correctness of the final P-value to ~5 significant digits
  while ((curTerm = m_pmf(++kk, prevTerm, m_pmfParams)) > SMALL_VALUE)
    {
      sum += curTerm;
      prevTerm = curTerm;
    }
  m_distn[k].pval = m_distn[k].pmf * sum;
  //  if (m_distn[k].pval < numeric_limits<double>::min())
  //    m_distn[k].pval = numeric_limits<double>::min();
  // Fill in P-values for any bins that were added between prev_max_k and k.
  kk = k;
  while (kk > prev_max_k + 1)
    {
      m_distn[kk - 1].pval = m_distn[kk - 1].pmf + m_distn[kk].pval;
      //      if (m_distn[k].pval < numeric_limits<double>::min())
      //	m_distn[k].pval = numeric_limits<double>::min();
      kk--;
    }

  return m_distn[k].pval;
}

// This method gets called at the end of a chromosome, at the end of a file,
// and anytime there's a gap in the data that's wider than half the background window width.
// It computes P-values for all sites in the current background window that need them,
// passes them along, and "flushes" the distribution (deletes it, resets accompanying variables).
void BackgroundRegionManager::computePandFlush(SiteManager& sm)
{
 if (m_distn.empty())
    return;

  if (!m_sliding || m_needToUpdate_kcutoff)
    {
      // Moving averages need to be computed, m_kcutoff needs to be determined,
      // mean and variance need to be computed, and all pmfs need to be computed.
      findCutoff();
    }

  computeStats(-1); // -1 means "for all observed values of k"

  while (!m_sitesInRegion_leftHalf.empty())
    {
      if (!m_sitesInRegion_leftHalf.front().hasPval)
        {
          long double pval;
          if (m_pmf != NULL)
            pval = getPvalue(m_sitesInRegion_leftHalf.front().count);
          else
            pval = 999.;
	  sm.processPvalue(pval // pass this P-value along to the corresponding site
#ifdef DEBUG
			   , m_sitesInRegion_leftHalf.front().sampled
#endif
			   );
	}
      m_sitesInRegion_leftHalf.pop_front();
    }
  while (!m_sitesInRegion_rightHalf.empty())
    {
      if (!m_sitesInRegion_rightHalf.front().hasPval)
        {
          long double pval;
          if (m_pmf != NULL)
            pval = getPvalue(m_sitesInRegion_rightHalf.front().count);
          else
            pval = 999.;
	  sm.processPvalue(pval // pass this P-value along to the corresponding site
#ifdef DEBUG
			   , m_sitesInRegion_rightHalf.front().sampled
#endif
			   );
        }
      m_sitesInRegion_rightHalf.pop_front();
    }

  m_distn.clear();
  m_posL = m_posC = m_posR = -1;
  m_pCurChrom = NULL;
  m_runningSum_count = m_runningSum_countSquared = m_numPtsInNullRegion = 0;
  m_runningSum_count_duringPrevComputation = m_runningSum_countSquared_duringPrevComputation = m_numPtsInNullRegion_duringPrevComputation = 0;
  m_modeYval = m_modeXval = m_kcutoff = m_kTrendReversal = -1;
  m_sliding = false;
  m_needToUpdate_kcutoff = true;

  m_kvalsWithMinMAxN.clear();
  m_minMAxN = -1;

  m_prev_k = -1;
  m_pmf = NULL;
  m_pmfParams.clear();
  m_nextPosToSample = -1;
  m_sampledDataDistnSize = 0;
}

// This method gets called to start sliding the background window (m_sliding == false) or perform a slide.
// It computes a P-value for the site at the center of the background window, when present,
// and it slides the background rightward.
// For the initial slide (e.g., the initial background window on a chromosome),
// it also computes P-values for all sites in the left half of the background window.
//
// Strictly speaking, if m_sliding == false and the incoming site lands at the right boundary,
// no "slide" is performed, just an append operation and computations.
void BackgroundRegionManager::slideAndCompute(const SiteRange& s, SiteManager& sm)
{
  if (m_posR > s.endPos)
    {
      cerr << "Coding error:  line " << __LINE__ << ", slideAndCompute window = " << *m_pCurChrom << ':'
           << m_posL << '-' << m_posR << '\n'
           << "erroneously received incoming position = " << s.endPos << "." << endl
           << endl;
      exit(1);
    }

  SiteData sd;
  sd.pos = s.endPos;
  sd.count = s.count;
  sd.hasPval = false;
  sd.sampled = false;

  if (!m_sliding)
    {
      // When m_sliding == false, this function, slideAndCompute,
      // really functions as "add and compute."
      if (s.endPos == m_posR)
        {
          add(s); // append this site
          sm.addSite(s);
        }

      findCutoff(); // because !m_sliding, findCutoff() will compute all moving averages
      computeStats(-1); // -1 means compute "for all count values"

      // P-values have been computed for all counts observed in this region.
      // Assign these P-values to the counts observed in the left half of this region
      // (i.e., assign to all points to the left of the central bp of this region).
      for (deque<SiteData>::iterator it = m_sitesInRegion_leftHalf.begin();
           it != m_sitesInRegion_leftHalf.end();
           it++)
        {
          long double pval;
          if (m_pmf != NULL)
            pval = getPvalue(it->count);
          else
            pval = 999.;
	  sm.processPvalue(pval
#ifdef DEBUG
			   , it->sampled
#endif
			   );
          it->hasPval = true;
        }

      // The region is centered on a specific position.
      // If data was observed for that position (this is usually true),
      // assign a P-value for that count at that position.
      // P-values for positions to the right of this position
      // will be determined later, one position at a time,
      // as the region "slides" rightward and new positions, in turn,
      // become the central position.
      if (m_sitesInRegion_rightHalf.front().pos == m_posC)
        {
          long double pval;
          if (m_pmf != NULL)
            pval = getPvalue(m_sitesInRegion_rightHalf.front().count);
          else
            pval = 999.;
	  
	  // pass this P-value along for the corresponding site
	  sm.processPvalue(pval
#ifdef DEBUG
			   , m_sitesInRegion_rightHalf.front().sampled
#endif
			   );
          m_sitesInRegion_rightHalf.front().hasPval = true;
        }
      m_sliding = true;

      if (s.endPos == m_posR)
        return;

      // else proceed to slide this region
    } // End of if (!m_sliding).  Note that m_sliding is now true.

  int idxMin(-1), idxMax(-1);
  bool updateDistn(false);

  // If we reach here, we're about to perform a 1bp slide.
  // Either we need to perform several 1bp slide events
  // until the next 1bp slide event brings in pos == s.endPos,
  // or we're now going to slide 1bp and bring in pos == s.endPos.
  bool needToComputePMFs(false);

  if (m_needToUpdate_kcutoff)
    {
      cerr << "Coding error:  slideAndCompute(), m_sliding == true, line "
           << __LINE__ << ", expected m_needToUpdate_kcutoff = false, but it's true.\n"
           << "Region = ";
      cerr << *m_pCurChrom << ':'
           << "[" << m_posL << ',' << m_posC << ',' << m_posR
           << "], incoming pos = " << s.endPos << ", k_c = " << m_kcutoff
           << endl;
      exit(1);
    }

  while (m_posR + 1 < s.endPos)
    {
      // Pop from the left half and update if necessary.
      if (m_sitesInRegion_leftHalf.front().pos == m_posL)
        {
          int prevModeXval(m_modeXval);
          const int k = m_sitesInRegion_leftHalf.front().count;
          updateDistn = m_sitesInRegion_leftHalf.front().sampled;

          m_sitesInRegion_leftHalf.pop_front();

          if (updateDistn)
            {
              // Update the values used to compute the mean and variance of the estimated null distribution.
              if (k <= m_kcutoff)
                {
                  m_runningSum_count -= k;
                  m_runningSum_countSquared -= k * k;
                  m_numPtsInNullRegion--;
                }
              m_distn[k].numOccs--;
              // Update moving averages (technically, moving sums, not averages, because we're not dividing them by N).
              idxMin = max(k - m_MAlength / 2, m_MAlength / 2);
              idxMax = min(k + m_MAlength / 2, m_sampledDataDistnSize - 1 - m_MAlength / 2);
              for (int i = idxMin; i <= idxMax; i++)
                {
                  m_distn[i].MAxN -= 1;
                  if (i == m_modeXval) // note that m_modeXval could be -1, in which case i can't equal it
                    {
                      m_modeYval--;
                      // Check whether this subtraction reveals a new mode to the left or right of it.
                      for (int j = m_MAlength / 2; j < m_sampledDataDistnSize - m_MAlength / 2; j++)
                        {
                          if (m_distn[j].MAxN > m_modeYval)
                            {
                              m_modeXval = j;
                              m_modeYval = m_distn[j].MAxN;
                              m_needToUpdate_kcutoff = true;
                            }
                        }
                    }
                }

              if (0 == m_distn[k].numOccs && k == m_sampledDataDistnSize - 1)
                {
                  // The bin at the end of the count distribution/histogram is now empty.
                  // Delete it, and delete any empty bins immediately preceding it,
                  // so that the highest bin contains at least one observation.
                  // If m_samplingInterval != 1 and m_distn.size() > m_sampledDataDistnSize,
                  // delete all the excess bins in the upper tail for simplicity's sake.
                  // (Such bins, each with numOccs == 0, will exist when, e.g.,
                  // only observations with k <= 18 have been "sampled" for use
                  // in the distribution, but k == 23, unsampled, was nonetheless observed,
                  // and a P-value was computed for it.
                  while (!m_distn.empty() && 0 == m_distn.back().numOccs)
                    m_distn.pop_back();
                  m_sampledDataDistnSize = static_cast<int>(m_distn.size());
                  // Because we've deleted 1+ bins from the end of m_distn,
                  // 1+ moving averages at the end of m_distn are now undefined.
                  // (This will occur infrequently.)
                  // Mark them as such for bookkeeping's sake.
                  for (int i = m_sampledDataDistnSize - 1; i > m_sampledDataDistnSize - 1 - m_MAlength / 2 && i > -1; i--)
                    m_distn[i].MAxN = -1;
                  // If we deleted the bin corresponding to m_kcutoff,
                  // update m_kcutoff so that it's within range.
                  // Let findCutoff() do this, so all appropriate variables will get updated.
                  if (m_kcutoff >= m_sampledDataDistnSize)
                    m_needToUpdate_kcutoff = true;
                  // Bring idxMax back within range if necessary, now that m_distn.size() has decreased.
                  if (idxMax > m_sampledDataDistnSize - 1 - m_MAlength / 2)
                    idxMax = m_sampledDataDistnSize - 1 - m_MAlength / 2; // idxMax might now be < idxMin; ok if so
                }

              if (!m_needToUpdate_kcutoff) // then perform some additional tests and maybe set m_needToUpdate_kcutoff = true
                {
                  if (m_modeXval != prevModeXval || m_sampledDataDistnSize - 1 == m_kcutoff)
                    m_needToUpdate_kcutoff = true; // for safety's sake, at least
                  else
                    {
                      const int halfMAlength = m_MAlength / 2;
                      if (0 == m_minMAxN)
                        {
                          if (k + halfMAlength > m_modeXval + 1 && k - halfMAlength < m_kcutoff)
                            {
                              // Check whether the subtraction has created a new instance of 0 == m_minMAxN at some k < m_kcutoff.
                              idxMin = max(k - halfMAlength, halfMAlength);
                              idxMax = min(k + halfMAlength, m_sampledDataDistnSize - 1 - halfMAlength);
                              for (int i = idxMax; i >= idxMin; i--)
                                {
                                  if (0 == m_distn[i].MAxN)
                                    {
                                      m_needToUpdate_kcutoff = true;
                                      break;
                                    }
                                }
                            }
                        }
                      else
                        {
                          if (-1 == m_kTrendReversal)
                            {
                              cerr << "Coding error:  m_kTrendReversal should NOT be -1 on line " << __LINE__ << "." << endl;
                              exit(1);
                            }
                          if (k + halfMAlength > m_modeXval + 1 && k - halfMAlength <= m_kTrendReversal)
                            m_needToUpdate_kcutoff = true; // possibly overkill, but worth doing for safety's sake
                        }
                    }
                }
            } // end of "if (updateDistn)"
        } // end of "if m_sitesInRegion_leftHalf.front().pos == m_posL"

      // When we have an observation for the central position,
      // move it from the leftmost position in the right half
      // to the rightmost position in the left half.
      if (m_sitesInRegion_rightHalf.front().pos == m_posC)
        {
          m_sitesInRegion_leftHalf.push_back(m_sitesInRegion_rightHalf.front());
          m_sitesInRegion_rightHalf.pop_front();
        }
      m_posL++;
      m_posC++;
      m_posR++;

      if (m_needToUpdate_kcutoff)
        {
          findCutoff(); // sets m_needToUpdate_kcutoff = false
          needToComputePMFs = true;
        }

      // Compute/assign P-value for m_posC if necessary.
      if (m_sitesInRegion_rightHalf.front().pos == m_posC)
        {
          if (m_runningSum_countSquared_duringPrevComputation != m_runningSum_countSquared || m_runningSum_count_duringPrevComputation != m_runningSum_count || m_numPtsInNullRegion_duringPrevComputation != m_numPtsInNullRegion)
            needToComputePMFs = true; // mean and/or variance have changed; won't change if removed k > m_kcutoff
          if (needToComputePMFs)
            {
              m_prev_k = -1; // compute pmfs from k=0 through current k
              computeStats(m_sitesInRegion_rightHalf.front().count); // sets m_prev_k = m_sitesInRegion_rightHalf.front().count
              needToComputePMFs = false;
            }
          else // We've calculated pmfs for this distribution, but it's possible we haven't computed one for a k this large.
              if (m_prev_k != -1 && m_sitesInRegion_rightHalf.front().count > m_prev_k)
            {
              computeStats(m_sitesInRegion_rightHalf.front().count); // sets m_prev_k = m_sitesInRegion_rightHalf.front().count
              needToComputePMFs = false;
            }
          long double pval;
          if (m_pmf != NULL)
            pval = getPvalue(m_sitesInRegion_rightHalf.front().count);
          else
            pval = 999.;
	  // pass this P-value along for the corresponding site
	  sm.processPvalue(pval
#ifdef DEBUG
			   , m_sitesInRegion_rightHalf.front().sampled
#endif
			   );
          m_sitesInRegion_rightHalf.front().hasPval = true;
        }
    } // end of "while sliding and not bringing in any new observations because there's missing data there"

  // If we reach here,
  // we're about to slide 1bp and bring in pos == m_posR + 1.

  m_sitesInRegion_rightHalf.push_back(sd); // process its addition to m_distn below
  sm.addSite(s);
  if (1 != m_samplingInterval)
    {
      if (s.endPos > m_nextPosToSample)
        m_nextPosToSample = s.endPos;

      // Alternatively, "sampling" could be performed strictly at multiples of m_samplingInterval.
      // If we intend to sample at, say, positions 200, 400, 600, 800, 1000, etc.,
      // and regions 390-402 and 598-620 are excluded (unmappable),
      // then the current code will sample at 200, 403, 621, 821, etc.
      // Sampling strictly at multiples of m_samplingInterval would in this case
      // sample at 200, 800, 1000, etc.
      // Use the following code if strict multiples of m_samplingInterval are desired.
      //
      // m_nextPosToSample += static_cast<int>(ceil(static_cast<double>(s.endPos - m_nextPosToSample)/static_cast<double>(m_samplingInterval))) * m_samplingInterval;
      //
      // The above is a one-step version of "while(s.endPos > m_nextPosToSample){m_nextPosToSample += m_samplingInterval;}."
      // Typically, m_nextPosToSample will only need to increase by m_samplingInterval;
      // it will need to increase by a greater multiple of m_samplingInterval when a large gap is present in the data
      // (an unmappable region or otherwise restricted/withheld region).

      if (s.endPos == m_nextPosToSample)
        {
          m_sitesInRegion_rightHalf.back().sampled = true;
          m_nextPosToSample += m_samplingInterval;
        }
    }
  else
    m_sitesInRegion_rightHalf.back().sampled = true;
  // increment m_posR below

  if (m_needToUpdate_kcutoff)
    {
      cerr << "Coding error:  slideAndCompute(), m_sliding == true, line "
           << __LINE__ << ", successfully slid through missing data, expected m_needToUpdate_kcutoff = false, but it's true.\n"
           << "Region = ";
      cerr << *m_pCurChrom << ':'
           << "[" << m_posL << ',' << m_posC << ',' << m_posR
           << "], incoming pos = " << s.endPos << ", k_c = " << m_kcutoff
           << endl;
      exit(1);
    }

  // If there's an observation exiting the region whose count
  // equals that of the observation entering the region,
  // and both were sampled,
  // the distribution remains unchanged, and no calculations need to be made,
  // unless the necessary calculations were postponed during a previous execution of this method.
  if (m_sitesInRegion_leftHalf.front().pos == m_posL && m_sitesInRegion_leftHalf.front().count == s.count && m_sitesInRegion_leftHalf.front().sampled && m_sitesInRegion_rightHalf.back().sampled)
    {
      m_sitesInRegion_leftHalf.pop_front();
      // When we have an observation for the central position,
      // move it from the leftmost position in the right half
      // to the rightmost position in the left half.
      if (m_sitesInRegion_rightHalf.front().pos == m_posC)
        {
          m_sitesInRegion_leftHalf.push_back(m_sitesInRegion_rightHalf.front());
          m_sitesInRegion_rightHalf.pop_front();
        }
      m_posL++;
      m_posC++;
      m_posR++;
      // Assign P-value for m_posC if necessary.
      if (m_sitesInRegion_rightHalf.front().pos == m_posC)
        {
          if (m_runningSum_countSquared_duringPrevComputation != m_runningSum_countSquared || m_runningSum_count_duringPrevComputation != m_runningSum_count || m_numPtsInNullRegion_duringPrevComputation != m_numPtsInNullRegion)
            needToComputePMFs = true; // Should only be true if we previously pop_fronted k < m_kcutoff without push_backing
          // (due to missing data at the right edge), and additionally,
          // there was missing data at m_posC, so no pmfs were computed.
          if (needToComputePMFs) // very unlikely to be true here
            {
              m_prev_k = -1; // compute pmfs for k=0 through current k
              computeStats(m_sitesInRegion_rightHalf.front().count); // sets m_prev_k = m_sitesInRegion_rightHalf.front().count
              needToComputePMFs = false;
            }
          else if (m_prev_k != -1 && m_sitesInRegion_rightHalf.front().count > m_prev_k)
            {
              // Only pmfs for k <= m_prev_k have been computed and stored, to save time.
              // Stil need to compute the pmfs for m_prev_k < k <= this k.
              computeStats(m_sitesInRegion_rightHalf.front().count); // sets m_prev_k = m_sitesInRegion_rightHalf.front().count
            }
          long double pval;
          if (m_pmf != NULL)
            pval = getPvalue(m_sitesInRegion_rightHalf.front().count);
          else
            pval = 999.;
	  // pass this P-value along for the corresponding site
	  sm.processPvalue(pval
#ifdef DEBUG
			   , m_sitesInRegion_rightHalf.front().sampled
#endif
			   );
          m_sitesInRegion_rightHalf.front().hasPval = true;
        }

      return;
    }

  // If we reach here,
  // we're adding a site at the right boundary.
  // We might be removing one from the left boundary (probably so).
  // If an observed site slides into the center position,
  // we'll have to compute a P-value for it.
  // (Or assign one from the current distribution
  // in the unlikely event of incoming and outgoing sites
  // both having high counts in what has been declared an upper tail beyond the null.)

  int k_incoming = s.count;
  int k_outgoing = -1;
  int origModeXval(m_modeXval);
  int origDistnSize(m_sampledDataDistnSize);

  if (m_sitesInRegion_leftHalf.front().pos == m_posL)
    {
      if (m_sitesInRegion_leftHalf.front().sampled)
        k_outgoing = m_sitesInRegion_leftHalf.front().count;
      m_sitesInRegion_leftHalf.pop_front();
    }
  if (k_outgoing != -1)
    {
      // Update the values used to compute the mean and variance of the estimated null distribution.
      if (k_outgoing <= m_kcutoff)
        {
          m_runningSum_count -= k_outgoing;
          m_runningSum_countSquared -= k_outgoing * k_outgoing;
          m_numPtsInNullRegion--;
        }
      m_distn[k_outgoing].numOccs--;
      if (0 == m_distn[k_outgoing].numOccs && k_outgoing == m_sampledDataDistnSize - 1)
        {
          // The bin at the end of the count distribution/histogram is now empty.
          // Delete it, and delete any empty bins immediately preceding it,
          // so that the highest bin contains at least one observation.
          // If m_samplingInterval != 1 and m_distn.size() > m_sampledDataDistnSize,
          // delete all the excess bins in the upper tail for simplicity's sake.
          // (Such bins, each with numOccs == 0, will exist when, e.g.,
          // only observations with k <= 18 have been "sampled" for use
          // in the distribution, but k == 23, unsampled, was nonetheless observed,
          // and a P-value was computed for it.
          while (!m_distn.empty() && 0 == m_distn.back().numOccs)
            m_distn.pop_back();
          m_sampledDataDistnSize = static_cast<int>(m_distn.size());
          // Because we've deleted 1+ bins from the end of m_distn,
          // 1+ moving averages at the end of m_distn are now undefined.
          // (This will happen infrequently.)
          // Mark them as such for bookkeeping's sake.
          for (int i = m_sampledDataDistnSize - 1; i > m_sampledDataDistnSize - 1 - m_MAlength / 2 && i > -1; i--)
            m_distn[i].MAxN = -1;
          if (origDistnSize >= m_MAlength && m_sampledDataDistnSize < m_MAlength)
            {
              // There are now too few bins to compute a MAxN of length m_MAlength, so the mode is undefined.
              m_modeXval = m_modeYval = -1;
            }
          // If we deleted the bin corresponding to m_kcutoff,
          // update m_kcutoff so that it's within range.
          // Let findCutoff() do this, so all appropriate variables will get updated.
          if (m_kcutoff >= m_sampledDataDistnSize)
            m_needToUpdate_kcutoff = true;
          // Now redefine "origDistnSize" so that it's correct with respect to k_incoming, which might replace deleted bin(s) and add new ones.
          origDistnSize = m_sampledDataDistnSize;
        }

      // Update moving averages (technically, moving sums, not averages, because we're not dividing them by N).
      idxMin = max(k_outgoing - m_MAlength / 2, m_MAlength / 2);
      idxMax = min(k_outgoing + m_MAlength / 2, m_sampledDataDistnSize - 1 - m_MAlength / 2);
      for (int i = idxMin; i <= idxMax; i++)
        {
          m_distn[i].MAxN -= 1;
          if (i == m_modeXval)
            {
              m_modeYval--;
              // There's a small chance that this subtraction has moved the mode leftward or rightward.
              for (int j = m_MAlength / 2; j < m_sampledDataDistnSize - m_MAlength / 2; j++)
                {
                  if (m_distn[j].MAxN > m_modeYval || (j > m_modeXval && m_distn[j].MAxN == m_modeYval))
                    {
                      m_modeXval = j;
                      m_modeYval = m_distn[j].MAxN;
                    }
                }
            }
        }
    }
  m_posL++;

  // Note:  The incoming site was pushed onto the end of m_sitesInRegion_rightHalf above, hence no need to do it here.
  m_posR++;

  if (m_sitesInRegion_rightHalf.back().sampled)
    {
      // Update the values used to compute the mean and variance of the estimated null distribution.
      if (k_incoming <= m_kcutoff)
        {
          m_runningSum_count += k_incoming;
          m_runningSum_countSquared += k_incoming * k_incoming;
          m_numPtsInNullRegion++;
        }
      if (k_incoming < m_sampledDataDistnSize)
        m_distn[k_incoming].numOccs++; // NOTE:  Still need to update MAxN values; will do that below.
      else
        {
          // Add bins to m_distn.         NOTE:  MAxN values get updated here in this case.
          int startHere = m_sampledDataDistnSize - m_MAlength / 2;
          StatsForCount sc;
          sc.numOccs = 0;
          sc.pmf = sc.pval = -1.;
          sc.MAxN = -1;

          while (m_sampledDataDistnSize < k_incoming && m_sampledDataDistnSize < static_cast<int>(m_distn.size()))
            m_sampledDataDistnSize++;
          if (m_sampledDataDistnSize < static_cast<int>(m_distn.size()))
            {
              m_distn[m_sampledDataDistnSize].numOccs = 1;
              m_sampledDataDistnSize++;
            }
          else
            {
              while (static_cast<int>(m_distn.size()) < k_incoming)
                m_distn.push_back(sc); // create bins for unobserved interior values, e.g., count = 5 but only 0,1,2 have been observed so far
              sc.numOccs = 1;
              m_distn.push_back(sc);
              m_sampledDataDistnSize = static_cast<int>(m_distn.size());
            }
          if (startHere >= m_MAlength / 2) // then we have at least one valid MAxN value that we will now update
            {
              int sum(0);
              int idxL(startHere - m_MAlength / 2), idxC(startHere), idxR(startHere + m_MAlength / 2);
              for (int i = idxL; i <= idxR; i++)
                sum += m_distn[i].numOccs;
              m_distn[idxC].MAxN = sum;
              if (m_distn[idxC].MAxN > m_modeYval) // also true when m_modeYval == m_modeXval == -1
                {
                  m_modeXval = idxC;
                  m_modeYval = m_distn[idxC].MAxN;
                }
              idxR++;
              while (idxR < m_sampledDataDistnSize)
                {
                  sum -= m_distn[idxL++].numOccs;
                  sum += m_distn[idxR++].numOccs;
                  idxC++;
                  m_distn[idxC].MAxN = sum;
                  if (m_distn[idxC].MAxN > m_modeYval) // also true when m_modeYval == m_modeXval == -1
                    {
                      m_modeXval = idxC;
                      m_modeYval = m_distn[idxC].MAxN;
                    }
                }
            }
          else if (m_sampledDataDistnSize >= m_MAlength)
            {
              // m_distn contained very few bins before an observation of k_incoming slid into the current region,
              // too few to compute any moving averages (times N, MAxN), but k_incoming is large enough
              // that now, m_distn contains enough bins to compute at least one MAxN, maybe even several
              // (e.g., m_distn had bins for k=0,1,2, and then k_incoming=6 suddenly slid into the region).
              // So we need to fill in the rightmost MAxN value and work leftwards from there.
              int stopHere = max(startHere, m_MAlength / 2 - 1);
              startHere = m_sampledDataDistnSize - 1 - m_MAlength / 2;
              int sum(0);
              int idxL(startHere - m_MAlength / 2), idxC(startHere), idxR(startHere + m_MAlength / 2);
              for (int i = idxL; i <= idxR; i++)
                sum += m_distn[i].numOccs;
              m_distn[idxC].MAxN = sum;
              if (m_distn[idxC].MAxN > m_modeYval) // recall m_modeYval == -1 if there had been too few bins to compute a MAxN value
                {
                  m_modeXval = idxC;
                  m_modeYval = m_distn[idxC].MAxN;
                }
              idxC--;
              idxL--;
              while (idxC != stopHere)
                {
                  sum -= m_distn[idxR--].numOccs;
                  sum += m_distn[idxL--].numOccs;
                  m_distn[idxC].MAxN = sum;
                  if (m_distn[idxC].MAxN > m_modeYval) // >, not >=, because in the event of a tie, we want to choose the rightmost mode
                    {
                      m_modeXval = idxC;
                      m_modeYval = m_distn[idxC].MAxN;
                    }
                  idxC--;
                }
            }
        }

      // If one or more bins were added to the end of m_distn because k_incoming was >= m_distn.size(),
      // even if bins were deleted due to k_outgoing and replaced due to k_incoming,
      // then the MAxN values are up-to-date with respect to k_incoming.
      // Otherwise, update moving averages (technically, moving sums, not averages, because we're not dividing them by N)
      // to reflect the addition of k_incoming.
      if (m_sampledDataDistnSize <= origDistnSize) // <=, not ==, because k_outgoing could have caused shrinkage of m_distn.
        {
          idxMin = max(k_incoming - m_MAlength / 2, m_MAlength / 2);
          idxMax = min(k_incoming + m_MAlength / 2, m_sampledDataDistnSize - 1 - m_MAlength / 2);
          for (int i = idxMin; i <= idxMax; i++)
            {
              m_distn[i].MAxN += 1;
              if (i == m_modeXval)
                m_modeYval++;
              else
                {
                  if (m_distn[i].MAxN > m_modeYval)
                    {
                      m_modeXval = i;
                      m_modeYval = m_distn[i].MAxN;
                    }
                  else if (i > m_modeXval && m_distn[i].MAxN == m_modeYval)
                    {
                      // The addition has moved the mode rightward.
                      m_modeXval = i;
                      m_modeYval = m_distn[i].MAxN;
                    }
                }
            }
        }
      else
        {
          // Bins were added to the end of m_distn; if m_kcutoff encompassed all of m_distn before, it needs to encompass the new bin(s) too.
          if (origDistnSize - 1 == m_kcutoff)
            m_needToUpdate_kcutoff = true;
        }

      // At this point, all MAxN values have been updated.

      if (!m_needToUpdate_kcutoff)
        {
          // If the mode has changed, recompute everything.
          if (m_modeXval != origModeXval)
            m_needToUpdate_kcutoff = true;
          else
            {
              if (m_sampledDataDistnSize - 1 == m_kcutoff)
                m_needToUpdate_kcutoff = true; // Probably m_kcutoff won't change, but do a full check to make sure.
              else
                {
                  if (0 == m_minMAxN)
                    {
                      if (m_distn[m_kcutoff].MAxN != 0)
                        {
                          // m_kcutoff probably won't change, but it no longer has MAxN==0, so recompute.
                          m_needToUpdate_kcutoff = true;
                        }
                      else
                        {
                          int halfMAlength = m_MAlength / 2;
                          if (k_outgoing != -1 && k_outgoing + halfMAlength > m_modeXval + 1 && k_outgoing - halfMAlength < m_kcutoff)
                            {
                              // There's a tiny chance that the subtraction caused a bin left of m_kcutoff
                              // to get its MAxN value reduced to 0, thereby moving m_kcutoff leftward.
                              idxMin = max(k_outgoing - halfMAlength, halfMAlength);
                              idxMax = min(k_outgoing + halfMAlength, m_sampledDataDistnSize - 1 - halfMAlength);
                              for (int i = idxMax; i >= idxMin; i--)
                                {
                                  if (0 == m_distn[i].MAxN)
                                    {
                                      m_needToUpdate_kcutoff = true;
                                      break;
                                    }
                                }
                            }
                          // There's a tiny chance that an addition far to the left of m_kcutoff (where MAxN == 0)
                          // raised an MAxN value high enough that it is now higher than an MAxN value to its left
                          // by more than the threshold.  This has been observed: m_kcutoff = 86 where MAxN = 0,
                          // k = 47 has MAxN = 7, k = 51 has MAxN = 9 with 9 < 7*m_thresholdRatio, and then
                          // k_incoming raises 51's MAxN from 9 to 10 without touching 47's MAxN,
                          // so that 10 > 7*m_thesholdRatio and m_kcutoff needs to get set to k=47.
                          if (m_modeXval + 1 < k_incoming - halfMAlength && k_incoming + halfMAlength < m_kcutoff)
                            m_needToUpdate_kcutoff = true;
                        }
                    }
                  else
                    {
                      if (-1 == m_kTrendReversal)
                        {
                          cerr << "Coding error:  m_kTrendReversal should NOT be -1 on line " << __LINE__ << " of BRM::slideAndCompute()." << endl;
                          cerr << "k_c = " << m_kcutoff << ", m_minMAxN = " << m_minMAxN << ", m_modeXval = "
                               << m_modeXval << ", m_modeYval = " << m_modeYval
                               << ", k_out = " << k_outgoing << ", k_in = " << k_incoming
                               << ", region = ";
                          cerr << *m_pCurChrom << ':'
                               << "[" << m_posL << ',' << m_posC + 1 << ',' << m_posR << ']' << endl;
                          cerr << "m_distn = {{0," << m_distn[0].numOccs << ',' << m_distn[0].MAxN;
                          for (unsigned int q = 1; q < m_distn.size(); q++)
                            {
                              if (0 == (q + 1) % 5)
                                cerr << "},\n{" << q << ',' << m_distn[q].numOccs << ',' << m_distn[q].MAxN;
                              else
                                cerr << "}, {" << q << ',' << m_distn[q].numOccs << ',' << m_distn[q].MAxN;
                            }
                          cerr << "}}" << endl;
                          exit(1);
                        }
                      // Prior to k_incoming and, when present, k_outgoing,
                      // m_kcutoff's MAxN was a "global minimum so far" and m_kTrendReversal's MAxN
                      // was sufficiently higher to define a "trend reversal."
                      int halfMAlength = m_MAlength / 2;
                      if (k_incoming + halfMAlength > m_modeXval + 1 && k_incoming - halfMAlength < m_kTrendReversal)
                        m_needToUpdate_kcutoff = true; // possibly overkill, but worth doing for safety's sake
                      if (k_outgoing != -1 && k_outgoing + halfMAlength > m_modeXval + 1 && k_outgoing - halfMAlength <= m_kTrendReversal)
                        m_needToUpdate_kcutoff = true; // possibly overkill, but worth doing for safety's sake
                    }
                }
            }
        }
    } // end of "if (m_sitesInRegion_rightHalf.back().sampled)"

  // When we have an observation for the central position,
  // move it from the leftmost position in the right half
  // to the rightmost position in the left half.
  if (m_sitesInRegion_rightHalf.front().pos == m_posC)
    {
      m_sitesInRegion_leftHalf.push_back(m_sitesInRegion_rightHalf.front());
      m_sitesInRegion_rightHalf.pop_front();
    }
  m_posC++;

  if (m_needToUpdate_kcutoff)
    {
      findCutoff(); // sets m_needToUpdate_kcutoff = false
      needToComputePMFs = true;
    }

  // Compute/assign P-value for m_posC if necessary.
  if (m_sitesInRegion_rightHalf.front().pos == m_posC)
    {
      if (m_runningSum_countSquared_duringPrevComputation != m_runningSum_countSquared || m_runningSum_count_duringPrevComputation != m_runningSum_count || m_numPtsInNullRegion_duringPrevComputation != m_numPtsInNullRegion)
        needToComputePMFs = true; // mean and/or variance have changed; won't change if added and removed k > m_kcutoff
      if (needToComputePMFs)
        {
          m_prev_k = -1; // compute pmfs from k=0 through current k
          computeStats(m_sitesInRegion_rightHalf.front().count); // sets m_prev_k = m_sitesInRegion_rightHalf.front().count
          needToComputePMFs = false;
        }
      else // We've calculated pmfs for this distribution, but it's possible we haven't computed one for a k this large.
          if (m_prev_k != -1 && m_sitesInRegion_rightHalf.front().count > m_prev_k)
        {
          computeStats(m_sitesInRegion_rightHalf.front().count); // sets m_prev_k = m_sitesInRegion_rightHalf.front().count
          needToComputePMFs = false;
        }
      long double pval;
      if (m_pmf != NULL)
        pval = getPvalue(m_sitesInRegion_rightHalf.front().count);
      else
        pval = 999.;
      sm.processPvalue(pval
#ifdef DEBUG
		       , m_sitesInRegion_rightHalf.front().sampled
#endif
		       );
      m_sitesInRegion_rightHalf.front().hasPval = true;
    }

  else
    {
      // There's no observation corresponding to the central position in the region.
      // Therefore we don't need to assign, or compute, a P-value.
      // If m_needToUpdate_kcutoff == true, the computation will be done later.
      // If m_needToUpdate_kcutoff == false but needToComputePMFs == true,
      // we need to ensure that needToCompute will get triggered
      // the next time any P-value needs to be assigned,
      // either in a subsequent call to this method or in a call to the "compute and flush" one.
      if (!m_needToUpdate_kcutoff && needToComputePMFs)
        {
          if (m_runningSum_countSquared == m_runningSum_countSquared_duringPrevComputation)
            m_runningSum_countSquared_duringPrevComputation = -1; // this will ensure pmfs get computed next time
        }
    }
}

bool parseAndProcessInput(const int& windowSize, const int& samplingInterval, const int& MAlength,
			  ofstream& ofsPvalData);
bool parseAndProcessInput(const int& windowSize, const int& samplingInterval, const int& MAlength,
			  ofstream& ofsPvalData)
{
  const int BUFSIZE(1000);
  char buf[BUFSIZE], *p;
  long linenum(0);
  int fieldnum;
  const int halfWindowSize(windowSize / 2); // integer division
  SiteRange curSite, prevSite;
  
  BackgroundRegionManager brm(samplingInterval, MAlength);
  SiteManager sm(ofsPvalData);

  prevSite.chrom = NULL;
  prevSite.endPos = -1;
  curSite.ID = NULL;
  curSite.hasPval = false;
  curSite.pval = -1.;
#ifdef DEBUG
  curSite.sampled = false;
#endif

  long start, end;

  while (cin.getline(buf, BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      p = strtok(buf, "\t");
      curSite.chrom = intern(string(p));
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        {
        MissingField:
          cerr << "Error:  Missing required field " << fieldnum
               << " on line " << linenum << "." << endl
               << endl;
          return false;
        }
      start = atol(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      end = atol(p);
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      //curSite.ID = intern(string(p));
      fieldnum++;
      if (!(p = strtok(NULL, "\t")))
        goto MissingField;
      curSite.count = atoi(p);
      // ignore any further fields

      // When contiguous stretches of sites with identical counts are observed
      // within a line of input, they need to be processed one site at a time,
      // for statistical reasons.
      for (int siteEnd = start + 1; siteEnd <= end; siteEnd++)
        {
          curSite.endPos = siteEnd;
	  curSite.begPos = curSite.endPos - 1;

          if (curSite.chrom != prevSite.chrom || curSite.endPos > prevSite.endPos + halfWindowSize)
            {
              brm.computePandFlush(sm); // Compute P-values for all unprocessed sites in the window.
	      // Writes values to disk.
              // This method removes all count data from brm.
              brm.setBounds(curSite.chrom, curSite.endPos, curSite.endPos + windowSize - 1);
            }

          if (!brm.isSliding() && curSite.endPos < brm.getRightEdge())
            {
              brm.add(curSite);
              sm.addSite(curSite);
            }
          else
            brm.slideAndCompute(curSite, sm); // calls sm.addSite(curSite)

          prevSite = curSite;
        }
    }

  brm.computePandFlush(sm); // See explanatory comment above.
  sm.writeLastUnreportedSite();

  return true;
}

int main(int argc, char* argv[])
{

  // Option defaults
  int background_size = 50001;
  int sampling_interval = 1;
  int smoothing_parameter = 5; // recommend ca. 15 when the maximum # of sampled observations is ca. 250
  int print_help = 0;
  int print_version = 0;
  string infilename = "";
  string outfilename = "";
  string outfilenameChromNames = "";
  string outfilenamePvals = "";

  // Long-opt definitions
  static struct option long_options[] = {
    { "background_size", required_argument, 0, 'b' },
    { "sampling_interval", required_argument, 0, 'n' },
    { "smoothing_parameter", required_argument, 0, 'm' },
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "outputChromlist", required_argument, 0, 'c' },
    { "outputPvals", required_argument, 0, 'p' }, // note we had been using 'p' for num_pvals
    { "help", no_argument, &print_help, 1 },
    { "version", no_argument, &print_version, 1 },
    { 0, 0, 0, 0 }
  };

  // Parse options
  char c;
  stringstream ss; // Used for parsing doubles (allows scientific notation)
  while ((c = getopt_long(argc, argv, "b:f:m:n:p:s:i:o:c:hvV", long_options, NULL)) != -1)
    {
      switch (c)
        {
        case 'b':
          background_size = atoi(optarg);
          break;
        case 'n':
          sampling_interval = atoi(optarg);
          break;
        case 'm':
          smoothing_parameter = atoi(optarg);
          break;
        case 'p':
	  outfilenamePvals = optarg;
          break;
        case 'i':
          infilename = optarg;
          break;
        case 'o':
          outfilename = optarg;
          break;
	case 'c':
          outfilenameChromNames = optarg;
          break;
	case 'h':
          print_help = 1;
          break;
        case 'v':
        case 'V':
          print_version = 1;
          break;
        case 0:
          // long option received, do nothing
          break;
        default:
          print_help = 1;
        }
    }

  if (!print_help && !print_version && outfilenameChromNames.empty())
    {
      cerr << "Error:  No filename supplied for (temporary) output file of integer-to-chromosomeName mapping."
	   << endl
	   << endl;
      print_help = 1;
    }
  if (!print_help && !print_version && outfilenamePvals.empty())
    {
      cerr << "Error:  No filename supplied for (temporary) output file of scaled -log10(P) values."
	   << endl
	   << endl;
      print_help = 1;
    }
  
  // Print usage and exit if necessary
  if (print_help)
    {
      cerr << "Usage:  " << argv[0] << " [options] < in.cutcounts.bed > out.pvalues.bed\n"
           << "\n"
           << "Options: \n"
           << "  -b, --background_size=SIZE     The size of the background region (50001)\n"
           << "  -n, --sampling_interval=INT    How often (bp) to sample for null modeling (1)\n"
           << "  -m, --smoothing_prameter=INT   Smoothing parameter used in null modeling (5)\n"
           << "  -i, --input=FILE               A file to read input from (STDIN)\n"
           << "  -o, --output=FILE              A file to write output to (STDOUT)\n"
	   << "  -c, --outputChromlist=FILE     Output file to store chromName-to-int mapping\n"
	   << "  -p, --outputPvals=FILE         Output file to store scaled -log10(P) values and # occurrences\n"
	   << "  -v, --version                  Print the version information and exit\n"
           << "  -h, --help                     Display this helpful help\n"
           << "\n"
           << " output (sent to stdout) will be a .bed4 file with a scaled transformation of the P-value in field 4\n"
           << " input (received from stdin) requires IDs in field 4 and counts in field 5.\n"
           << endl
           << endl;
      return -1;
    }

  if (print_version)
    {
      cout << argv[0] << " version " << hotspot2_VERSION_MAJOR
           << '.' << hotspot2_VERSION_MINOR << endl;
      return 0;
    }

  ios_base::sync_with_stdio(false); // calling this static method in this way turns off checks, speeds up I/O

  if (!infilename.empty() && infilename != "-")
    {
      if (freopen(infilename.c_str(), "r", stdin) == NULL)
        {
          cerr << "Error: Couldn't open input file " << infilename << endl;
          return 1;
        }
    }
  if (!outfilename.empty() && outfilename != "-")
    {
      if (freopen(outfilename.c_str(), "w", stdout) == NULL)
        {
          cerr << "Error: Couldn't open output file " << outfilename << " for writing" << endl;
          return 1;
        }
    }

  ofstream ofsIntToChrnameMapping(outfilenameChromNames.c_str());
  if (!ofsIntToChrnameMapping)
    {
      cerr << "Error:  Unable to open file \"" << outfilenameChromNames << "\" for write."
	   << endl
	   << endl;
      return -1;
    }
  ofstream ofsPvals(outfilenamePvals.c_str());
  if (!ofsPvals)
    {
      cerr << "Error:  Unable to open file \"" << outfilenamePvals << "\" for write."
	   << endl
	   << endl;
      return -1;
    }
  
  if (!parseAndProcessInput(background_size, sampling_interval, smoothing_parameter,
			    ofsPvals))
    return -1;

  for (map<string*, int>::const_iterator it = chromAsInt.begin(); it != chromAsInt.end(); it++)
    ofsIntToChrnameMapping << it->second << '\t' << *(it->first) << '\n';
  
  return 0;
}
