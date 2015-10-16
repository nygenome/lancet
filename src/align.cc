#include "align.hh"

/******************************************************************
** Align.cc
**
** Gapped alignment based on the Smith-Waterman algorithm
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

using namespace std;

int MATCH = 2;
int MISMATCH = -4;

int INDEL = -2;

int GAP_OPEN = -8;
int GAP_EXTEND = -1;


int cmp(char s, char t)
{
  if (s == t) { return MATCH; }
  return MISMATCH;
}

int max(int a, int b)
{
  return (a > b) ? a : b;
}

int maxscore(int x, int y, int m, char & t)
{
  int max = m; t = '\\';

  if (x > max) { max = x; t = '<'; }
  if (y > max) { max = y; t = '^'; }

  return max;
}



typedef struct cell
{
public:
  cell(int s=0, char t='*') : score(s), tb(t) {}

  int score;
  char tb;
} cell;

cell mk_cell(int s, char t)
{
  return cell(s, t);
}

cell maxx(int a, int b)
{
  if (a > b) { return cell(a, '-'); }
               return cell(b, '<');
}

cell maxy(int a, int b)
{
  if (a > b) { return cell(a, '|'); }
               return cell(b, '^');
}

cell maxscorexy(cell scorex, cell scorey, int scorez)
{
  cell r = mk_cell(scorez, '\\');

  if (scorex.score > r.score) { r.score = scorex.score; r.tb = '<'; }
  if (scorey.score > r.score) { r.score = scorey.score; r.tb = '^'; }

  return r;
}

typedef struct aligned
{
  char s;
  char t;
  //float c;
  int score;
} aligned;


void global_align(const string & S, const string & T, 
                  string & S_aln, string & T_aln, 
                  int endfree, int V)
{
  S_aln.clear();
  T_aln.clear();

  //if (endfree) { V = 1; }

  int n = S.length();
  int m = T.length();

  if (V) { cout << "S: " << S << " " << n << endl; }
  if (V) { cout << "T: " << T << " " << m << endl; }

  std::vector< std::vector<int>  > scores   (n+2, std::vector<int>(m+2));
  std::vector< std::vector<char> > trackback(n+2, std::vector<char>(m+2));

  //int scores[n+2][m+2];
  //char trackback[n+2][m+2];

  int i, j;

  for (i = 0; i <= n; i++) { scores[i][0] = i*INDEL; trackback[i][0] = '<'; }
  for (j = 0; j <= m; j++) { scores[0][j] = j*INDEL; trackback[0][j] = '^'; }

  if (V) { cout << " "; for (int i = 1; i <= n; i++) { cout << "    " << S[i-1]; } cout << endl; }

  for (int j = 1; j <= m; j++)
  {
    if (V) { cout << T[j-1]; }

    for (int i = 1; i <= n; i++)
    {
      int leftscore  = scores[i-1][j] + INDEL;
      int upscore    = scores[i][j-1] + INDEL;
      int matchscore = scores[i-1][j-1] + cmp(S[i-1], T[j-1]);
      char tb;

      int score = maxscore(leftscore, upscore, matchscore, tb);

      scores[i][j] = score;
      trackback[i][j] = tb;

      if (V) { printf(" % 3d%c", scores[i][j], tb); }
    }

    if (V) { cout << endl; }
  }

  if (V) { cout << "Score is " << scores[n][m] << endl; }

  vector<aligned> trace;

  i = n; j = m;

  if (endfree)
  {
    int maxval = scores[0][m];
    i = 0;
    for (int q = 0; q < n; q++)
    {
      if (scores[q][m] > maxval)
      {
        i = q;
        maxval = scores[q][m];
      }
    }

    if (V) { cout << "Tracking back from: " << i << " score: " << maxval << endl; }
  }

  while (i > 0 || j > 0)
  {
   // cout << "at: " << i << "," << j << endl;
    int s = scores[i][j];
    char t = trackback[i][j];

    aligned a;
    a.score = s;

    if      (t == '*')  { break; }
    else if (t == '\\') { a.s = S[i-1]; a.t = T[j-1]; i--; j--; }
    else if (t == '<')  { a.s = S[i-1]; a.t = '-';    i--; }
    else if (t == '^')  { a.s = '-';    a.t = T[j-1]; j--; }

    trace.push_back(a);
  }

  for (int k = trace.size() - 1; k >= 0; k--)
  {
    if (V) { cout << "   " << trace[k].s; }
    S_aln.push_back(trace[k].s);
  }
  if (V) { cout << endl; }

  for (int k = trace.size() - 1; k >= 0; k--)
  {
    if (V) { cout << "   " << trace[k].t; }
    T_aln.push_back(trace[k].t);
  }
  if (V) { cout << endl; }

  if (V) 
  {
    for (int k = trace.size() - 1; k >= 0; k--)
    {
      printf(" %3d", trace[k].score);
    }
    cout << endl;

    cout << "1: " << S_aln << endl;
    cout << "2: " << T_aln << endl;
  }
}



void global_align_aff(const string & S, const string & T, 
	string & S_aln, string & T_aln,
	int endfree, int V)
{
  S_aln.clear();
  T_aln.clear();

  //if (endfree) { V = 1; }

  int n = S.length();
  int m = T.length();

  if (V) { cout << "S: " << S << " " << n << endl; }
  if (V) { cout << "T: " << T << " " << m << endl; }

  std::vector< std::vector<cell> > M (n+2, std::vector<cell>(m+2));
  std::vector< std::vector<cell> > X (n+2, std::vector<cell>(m+2));
  std::vector< std::vector<cell> > Y (n+2, std::vector<cell>(m+2));

  // base conditions
  int i, j;

  for (j = 0; j <= m; j++) { X[0][j] = mk_cell(GAP_OPEN + j*GAP_EXTEND, '^'); M[0][j] = X[0][j]; }
  for (i = 0; i <= n; i++) { Y[i][0] = mk_cell(GAP_OPEN + i*GAP_EXTEND, '<'); M[i][0] = Y[i][0]; }
  M[0][0] = mk_cell(0, '*');

  if (V) { cout << " "; for (int i = 1; i <= n; i++) { cout << "     " << S[i-1]; } cout << endl; }

  for (int j = 1; j <= m; j++)
  {
    if (V) { cout << T[j-1]; }

    for (int i = 1; i <= n; i++)
    {
      X[i][j] = maxx(X[i-1][j].score+GAP_EXTEND, M[i-1][j].score+GAP_OPEN);

      Y[i][j] = maxy(Y[i][j-1].score+GAP_EXTEND, M[i][j-1].score+GAP_OPEN);

      M[i][j] = maxscorexy(X[i][j], 
                           Y[i][j],
                           M[i-1][j-1].score+cmp(S[i-1],T[j-1]));

      if (V) { printf(" % 3d%c%c", M[i][j].score, M[i][j].tb, (i==j)?'*':' '); }
    }

    if (V) { cout << endl; }
  }

  if (V) { cout << "Score is " << M[n][m].score << endl; }

  vector<aligned> trace;

  i = n; j = m;

  if (endfree)
  {
    int maxval = M[0][m].score;
    i = 0;
    for (int q = 0; q < n; q++)
    {
      if (M[q][m].score > maxval)
      {
        i = q;
        maxval = M[q][m].score;
      }
    }

    if (V) { cout << "Tracking back from: " << i << " score: " << maxval << endl; }
  }

  bool forcey = false;
  bool forcex = false;

  while (i > 0 || j > 0)
  {
    int s = M[i][j].score;
    char t = M[i][j].tb;
    char x = X[i][j].tb;
    char y = Y[i][j].tb;
    char z = t;

    aligned a;
    a.score = s;

    if      (t == '*')  { break; }

    else if (forcex)    { a.s = S[i-1]; a.t = '-'; z=x; if (X[i][j].tb == '<') { forcex = false; } i--; }
    else if (t == '<')  { a.s = S[i-1]; a.t = '-'; z=x; if (X[i][j].tb == '-') { forcex = true;  } i--; }

    else if (forcey)    { a.s = '-'; a.t = T[j-1]; z=y; if (Y[i][j].tb == '^') { forcey = false; } j--; }
    else if (t == '^')  { a.s = '-'; a.t = T[j-1]; z=y; if (Y[i][j].tb == '|') { forcey = true;  } j--; }

    else if (t == '\\') { a.s = S[i-1]; a.t = T[j-1]; i--; j--; }
    
    else { cerr << "WTF!" << endl; exit(1); }

    if (V) { cerr << "M[" << i << "," << j << "]:\t" 
                  << s << " " << t << "," << x << (forcex ? '*' : ' ') << "," << y << (forcey ? '*' : ' ') << "," << z
                  << " | " << a.s << " " << a.t << endl; }

    trace.push_back(a);
  }

  for (int k = trace.size() - 1; k >= 0; k--)
  {
    if (V) { cout << "   " << trace[k].s; }
    S_aln.push_back(trace[k].s);
  }
  if (V) { cout << endl; }

  for (int k = trace.size() - 1; k >= 0; k--)
  {
    if (V) { cout << "   " << trace[k].t; }
    T_aln.push_back(trace[k].t);
  }
  if (V) { cout << endl; }

  if (V) 
  {
    for (int k = trace.size() - 1; k >= 0; k--)
    {
      printf(" %3d", trace[k].score);
    }
    cout << endl;

    cout << "S': " << S_aln << endl;
    cout << "T': " << T_aln << endl;

  }
}

/*
void global_cov_align_aff(const string & S, const string & T, const vector<float> & C, 
	string & S_aln, string & T_aln, string & C_aln,
	int endfree, int V)
{
  S_aln.clear();
  T_aln.clear();
  C_aln.clear();

  //if (endfree) { V = 1; }

  int n = S.length();
  int m = T.length();

  if (V) { cout << "S: " << S << " " << n << endl; }
  if (V) { cout << "T: " << T << " " << m << endl; }

  std::vector< std::vector<cell> > M (n+2, std::vector<cell>(m+2));
  std::vector< std::vector<cell> > X (n+2, std::vector<cell>(m+2));
  std::vector< std::vector<cell> > Y (n+2, std::vector<cell>(m+2));

  // base conditions
  int i, j;

  for (j = 0; j <= m; j++) { X[0][j] = mk_cell(GAP_OPEN + j*GAP_EXTEND, '^'); M[0][j] = X[0][j]; }
  for (i = 0; i <= n; i++) { Y[i][0] = mk_cell(GAP_OPEN + i*GAP_EXTEND, '<'); M[i][0] = Y[i][0]; }
  M[0][0] = mk_cell(0, '*');

  if (V) { cout << " "; for (int i = 1; i <= n; i++) { cout << "     " << S[i-1]; } cout << endl; }

  for (int j = 1; j <= m; j++)
  {
    if (V) { cout << T[j-1]; }

    for (int i = 1; i <= n; i++)
    {
      X[i][j] = maxx(X[i-1][j].score+GAP_EXTEND, M[i-1][j].score+GAP_OPEN);

      Y[i][j] = maxy(Y[i][j-1].score+GAP_EXTEND, M[i][j-1].score+GAP_OPEN);

      M[i][j] = maxscorexy(X[i][j], 
                           Y[i][j],
                           M[i-1][j-1].score+cmp(S[i-1],T[j-1]));

      if (V) { printf(" % 3d%c%c", M[i][j].score, M[i][j].tb, (i==j)?'*':' '); }
    }

    if (V) { cout << endl; }
  }

  if (V) { cout << "Score is " << M[n][m].score << endl; }

  vector<aligned> trace;

  i = n; j = m;

  if (endfree)
  {
    int maxval = M[0][m].score;
    i = 0;
    for (int q = 0; q < n; q++)
    {
      if (M[q][m].score > maxval)
      {
        i = q;
        maxval = M[q][m].score;
      }
    }

    if (V) { cout << "Tracking back from: " << i << " score: " << maxval << endl; }
  }

  bool forcey = false;
  bool forcex = false;

  while (i > 0 || j > 0)
  {
    int s = M[i][j].score;
    char t = M[i][j].tb;
    char x = X[i][j].tb;
    char y = Y[i][j].tb;
    char z = t;

    aligned a;
    a.score = s;

    if      (t == '*')  { break; }

    else if (forcex)    { a.s = S[i-1]; a.t = '-'; a.c = '-'; z=x; if (X[i][j].tb == '<') { forcex = false; } i--; }
    else if (t == '<')  { a.s = S[i-1]; a.t = '-'; a.c = '-'; z=x; if (X[i][j].tb == '-') { forcex = true;  } i--; }

    else if (forcey)    { a.s = '-'; a.t = T[j-1]; a.c = C[j-1]; z=y; if (Y[i][j].tb == '^') { forcey = false; } j--; }
    else if (t == '^')  { a.s = '-'; a.t = T[j-1]; a.c = C[j-1]; z=y; if (Y[i][j].tb == '|') { forcey = true;  } j--; }

	else if (t == '\\') { a.s = S[i-1]; a.t = T[j-1]; a.c = C[j-1]; i--; j--; }
    
    else { cerr << "WTF!" << endl; exit(1); }

    if (V) { cerr << "M[" << i << "," << j << "]:\t" 
                  << s << " " << t << "," << x << (forcex ? '*' : ' ') << "," << y << (forcey ? '*' : ' ') << "," << z
                  << " | " << a.s << " " << a.t << endl; }

    trace.push_back(a);
  }

  for (int k = trace.size() - 1; k >= 0; k--)
  {
    if (V) { cout << "   " << trace[k].s; }
    S_aln.push_back(trace[k].s);
  }
  if (V) { cout << endl; }

  for (int k = trace.size() - 1; k >= 0; k--)
  {
    if (V) { cout << "   " << trace[k].t; }
    T_aln.push_back(trace[k].t);
  }
  if (V) { cout << endl; }

  for (int k = trace.size() - 1; k >= 0; k--)
  {
    if (V) { cout << "   " << trace[k].c; }
    //C_aln.push_back(trace[k].c);
	stringstream ss;
	ss << trace[k].c;
    C_aln += ss.str();
	if(k!=0) { C_aln += " "; }
  }
  if (V) { cout << endl; }

  if (V) 
  {
    for (int k = trace.size() - 1; k >= 0; k--)
    {
      printf(" %3d", trace[k].score);
    }
    cout << endl;

    cout << "S': " << S_aln << endl;
    cout << "T': " << T_aln << endl;
    cout << "C': " << C_aln << endl;

  }
}
*/

