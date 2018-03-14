#include <iostream>
#include <vector>
#include <list>
#include <cstdlib>
#include <cstring>
#include <string>

using namespace std;

struct Peak {
  string chrom;
  long beg;
  long end;
  string id;
  float score;
  string score_string;
};

// function for sorting Peaks on the same chromosome in sort-bed order
bool Peak_coordLT(const Peak& a, const Peak& b);
bool Peak_coordLT(const Peak& a, const Peak& b)
{
  if (a.beg == b.beg)
    {
      if (a.end == b.end)
	{
	  if (a.id == b.id)
	    return a.score < b.score;
	  else
	    return a.id < b.id;
	}
      else
	return a.end < b.end;
    }
  return a.beg < b.beg;
}

// function for sorting Peaks in descending order by score
bool Peak_scoreGT(const Peak& a, const Peak& b);
bool Peak_scoreGT(const Peak& a, const Peak& b)
{
  return a.score > b.score;
}

// Function for resolving overlaps:
// Any Peak that overlaps the strongest Peak gets removed.  ("Strongest" = highest score.)
// In the resulting list, any Peak that overlaps the then-2nd-strongest Peak gets removed,
// and so on until we've compared all Peaks.
void resolveOverlaps(list<Peak>& peaks);
void resolveOverlaps(list<Peak>& peaks)
{
  peaks.sort(Peak_scoreGT);

  for (list<Peak>::iterator it1 = peaks.begin(); it1 != peaks.end(); it1++)
    {
      list<Peak>::iterator it2(it1);
      it2++;
      while (it2 != peaks.end())
	{
	  if (it2->beg < it1->end && it2->end > it1->beg)
	    {
	      list<Peak>::iterator it3(it2);
	      it2++;
	      peaks.erase(it3);
	    }
	  else
	    it2++;
	}
    }
  // restore to genomic, sort-bed order
  peaks.sort(Peak_coordLT);
}

// function to parse/extract the data for each peak
// into a Peak data structure
bool PeakfromString(char *s, Peak& d);
bool PeakfromString(char *s, Peak& d)
{
  char *p = strtok(s, "\t");
  int fieldnum(1);
  if (!p)
    {
    MissingField:
      cerr << "Error:  Missing required field " << fieldnum
	   << " in " << s << "." << endl << endl;
      return false;
    }
  d.chrom = string(p);
  fieldnum++;
  if (!(p = strtok(NULL, "\t")))
    goto MissingField;
  d.beg = atol(p);
  fieldnum++;
  if (!(p = strtok(NULL, "\t")))
    goto MissingField;
  d.end = atol(p);
  fieldnum++;
  if (!(p = strtok(NULL, "\t")))
    goto MissingField;
  d.id = string(p);
  fieldnum++;
  if (!(p = strtok(NULL, "\t")))
    goto MissingField;
  d.score = atof(p);
  d.score_string = string(p);
  return true;
}

// This function reads overlap data from stdin (cin),
// parses it, resolves overlaps between peaks,
// and writes the non-overlapping subsets to stdout (cout).
// Peak data is expected to contain 5 fields, with a score in field 5.
// The input is expected to be the output from:
//    bedops -m <peakFile> | bedmap --echo --echo-map - <peakFile>
// Each line thus contains merged peaks (the outer bounds of overlapping peaks
// or of a peak that doesn't overlap any others), then "|",
// then a semicolon-delimited list of the peaks whose bounds are on the left of "|".
// See the comment above the resolveOverlaps function for a description
// of how the overlaps are resolved.
bool parseInputWriteOutput(void);
bool parseInputWriteOutput(void)
{
  const int BUFSIZE(500000);
  char buf[BUFSIZE], *p;
  long linenum(0);

  while (cin.getline(buf, BUFSIZE))
    {
      list<Peak> Peaks;
      vector<string> vec; // this vector is used for parsing/extracting Peak data
      linenum++;
      if (!(p = strtok(buf, "|")))
	{
	  cerr << "Error:  Failed to find |-delimited input on line " << linenum
	       << " of the input file." << endl << endl;
	  return false;
	}
      while ((p = strtok(NULL, ";")))
	vec.push_back(string(p));
      for (int i = 0; i < static_cast<int>(vec.size()); i++)
	{
	  Peak d;
	  char s[1000];
	  strcpy(s,vec[i].c_str());
	  if (!PeakfromString(s, d))
	    return false;
	  Peaks.push_back(d);
	}
      resolveOverlaps(Peaks);
      for (list<Peak>::const_iterator it = Peaks.begin();
	   it != Peaks.end(); it++)
	cout << it->chrom << '\t'
	     << it->beg << '\t'
	     << it->end << '\t'
	     << it->id << '\t'
	     << it->score_string << endl; // write the score exactly as we received it (same # of digits, etc.)
    }
  
  return true;
}

int main(int argc, const char* argv[])
{
  ios_base::sync_with_stdio(false); // calling this static method in this way turns off checks, and can speed up I/O

  if (!parseInputWriteOutput()) // see comment accompanying the function re: input and output
    return -1;

  return 0;
}

