#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "htslib/sam.h"
//#include "htslib/faidx.h"
//#include "bam.h"
#include <deque>
#include <cctype>

using namespace std;

///////////////////////

struct Pileup {
  // this is a basic data structure used to build a queue and report pileups
  int pos; // genomic coordinate
  char ref; // reference (genomic) nucleotide
  double nreads; // read coverage (overlap)
  double nreads_fwd;
  double nreads_rev;
  int nreads_uniq;
  char strand; // strand indicator

  Pileup() : pos(0), ref('N'), nreads(0), nreads_uniq(0), strand('+') {}
  Pileup(int p) : pos(p), ref('N'), nreads(0), nreads_uniq(0), strand('+'){ }
};

void process_queue(deque<Pileup> &q, int upto_pos, bool process_all,
		   const string &ref_id) {
  if (!q.empty())
  {
    // coordinates pos in queue are 0-based
    // output 0-based half-open interval!!!
    int end_pos = q.back().pos + 1;
    if (upto_pos < end_pos) end_pos = upto_pos;
    if (process_all) end_pos = q.back().pos+1;
    //cout << "upto: " << upto_pos << endl;
    //cout << "q_front_pos: " << q.front().pos << endl;
    //cout << "q_back_pos: " << q.back().pos << endl;
    //cout << "end_pos: " << end_pos << endl;
    double cov = q.front().nreads;
    int block_start = q.front().pos;
    double cov_prev = cov;
    int cov_uniq_prev = q.front().nreads_uniq;
    int i;
    char strand = q.front().strand;
    int q_front = q.front().pos;
    // one caveat: consecutively output intervals may have same value i.e. need to be merged!
    for (i=q.front().pos; i<end_pos; ++i)
    {
      cov = q.front().nreads;
      if (cov != cov_prev)
      {
        cout << ref_id << "\t" 
         << block_start << "\t"
         << i << "\t"
         << cov_prev << "\t" << cov_uniq_prev << "\t" << strand << endl;
        block_start = i; 
      }
      cov_prev = q.front().nreads;
      cov_uniq_prev = q.front().nreads_uniq;
      q.pop_front();
    }
    if ( q_front < i)
    cout << ref_id << "\t" 
         << block_start << "\t"
         << i << "\t"
         << cov_prev << "\t" << cov_uniq_prev << "\t" << strand << endl;
 
  }

  while( (!q.empty()) &&
	 (process_all || (q.front().pos < upto_pos))) {
    // output all completed positions (to the left of upto)
    // OR all queued positions if chr changed (process_all == TRUE)
    //starting with the leftmost position (queue front)
    // we should never end up here!!!!
    cout << "HERE" << endl;
    q.pop_front();
  }
}

/////////////////////
class DNAComplementer {
  char c[256];
public:
  DNAComplementer() { 
    for(int i=0; i<256; ++i)
      c[i] = 'n';
    c['a'] = 't'; c['c'] = 'g';  c['g'] = 'c'; c['t'] = 'a';
    c['m'] = 'k'; c['r'] = 'y';  c['k'] = 'm'; c['y'] = 'r';
    c['s'] = 's'; c['w'] = 'w';
    c['v'] = 'b'; c['h'] = 'd'; c['b'] = 'v'; c['d'] = 'h';
    for(int i='a'; i<='z'; ++i)
      c[toupper(i)] = char(toupper(c[int(i)]));
  }
  char operator () (char x) const { return c[int(x)]; }
};

/////////////////////

int main(int argc, char **argv) {
  vector<string> opts;
  vector<string> args;

  for(int i=1; i < argc; ++i) {
    string a(argv[i]);
    if (a.size() > 1 && a[0] == '-')
      opts.push_back(a);
    else
      args.push_back(a);
  }
  
  bool opts_valid = true;
  bool no_ss = false;         // not strand specific
  bool pair_ends = false;    //paired end reads

  bool filter5pClipped = true; // filter reads 5p-clipped by more than 1 nt
  int maxReadLength = 44; // filter reads with length > maxReadLength
  int minReadLength = 15; // filter reads with length < minReadLength
  
  for(int i=0; i < opts.size(); ++i) {
    if (opts[i] == "--noss")
      no_ss = true;
    else	if (opts[i] == "--paired")
		pair_ends = true;
    else if (opts[i] == "--keep5pClipped")
       filter5pClipped=false;
    else
      opts_valid = false;
  }
  //cerr << "filter5pClipped=" << filter5pClipped << endl;
  //return 1;
  if(args.size() < 1 || !opts_valid) {
    cerr << "USAGE: " << argv[0] << "[options] reads.bam <minReadLength> <maxReadLength>\n";
    cerr << "\t--noss\t\tNot strand-specific (convert everything to +)\n\t--paired\t\tPaired end sequencing (filter anything that is not properly paired)\n\t--keep5pClipped\t\tdefault is to throw away 5p clipped reads. Use this option to keep 5p clipped reads.\n";
    return 1;
  }

  string bam_fn( args[0] );
  if ( args.size() > 1 )
  {
     minReadLength = atoi( args[1].c_str() );
  }
  if ( args.size() > 2 )
  { 
     maxReadLength = atoi( args[2].c_str() );
  }

  //cout << "maxReadLength=" << maxReadLength << endl;
  //return 1;
  //string fas_fn( args[1] );

  // index the fasta file by finding out where each chr starts
  //ifstream file(fas_fn.c_str());
  //if (!file.is_open()) {
  //  cerr << "Failed to open FASTA file " << fas_fn << "\n";
  //  return 1;
  //}

  //faidx_t *fai = fai_load(fas_fn.c_str());


  samFile *bam_file;
  bam_hdr_t *bam_hdr;

  if ((bam_file = sam_open(bam_fn.c_str(), "r")) == 0) {
    cerr << "Failed to open BAM file " << bam_fn << "\n";
    return 1;
  }

  bam_hdr = sam_hdr_read(bam_file);
  bam1_t *bam = bam_init1();

/*typedef struct {
    int32_t n_targets, ignore_sam_err;
    uint32_t l_text;
    uint32_t *target_len;
    int8_t *cigar_tab;
    char **target_name;
    char *text;
    void *sdict;
} bam_hdr_t;
*/
   //cout << "numChr" << "\t" << bam_hdr->n_targets << endl;
   for (int i = 0; i < bam_hdr->n_targets; ++i)
   {
        char *chrName = bam_hdr->target_name[i];
        int chrLen = bam_hdr->target_len[i];
        cout << "chrInfo"
             << "\t" << chrName
             << "\t" << chrLen
             << endl;
   }

  string curr_ref;
  string prev_ref;
  char *ref_seq = NULL;
  int ref_len(0);
  bool changed_ref = false;

  string bases(16, 'X');
  bases[1] = 'A';
  bases[2] = 'C';
  bases[4] = 'G';
  bases[8] = 'T';
  bases[15] = 'N';

  //ReverseComplementer rev_comp;

  deque<Pileup> q; // double-ended queue; (forward strand)
                   // holds pileups for positions being processed

  deque<Pileup> qr; // double-ended queue; (reverse strand)
                   // holds pileups for positions being processed


  vector<double> interval_cnt(maxReadLength+1,0);
  vector<double> interval_cnt_rev(maxReadLength+1,0);
  vector<int> interval_idx(maxReadLength+1,-1);
  vector<int> interval_idx_rev(maxReadLength+1,-1);
  int cnt_interval = 0;
  int cnt_interval_rev = 0;
  int read_pos_prev = 0;
  char read_strand_prev;
  // read BAM file
  // bam variable points at current alignment record
  while( sam_read1(bam_file,bam_hdr,bam) > 0) {
    // bam = read alignment bam1_t
    changed_ref = false;

    // skip non-unique reads
    //int num_hits = bam_aux2i( bam_aux_get(bam, "NH") );
    //if (num_hits > 1)
    //  continue;

    // skip unmapped reads
    if ((bam->core.flag & 0x4) > 0)
      continue;
	//skip if unpaired, secondary or duplicate if pair_ends is set
	if(pair_ends&&((bam->core.flag & 0x2) == 0 || bam->core.flag & 0x100 || bam->core.flag & 0x400)){
		continue;
	}	
    // check for soft-clipped reads (clipped seq is present in SEQ field but missing from reference)
    // on both ends
    bool has_indels(false);
    bool has_skip(false); // CIGAR will have N (reference skip) in it
    int nclipstart(0); // left soft-clip size 
    int nclipend(0); // right soft-clip size

    for(int i=0; i < bam->core.n_cigar; ++i) {
      if (bam_cigar_op(bam_get_cigar(bam)[i]) == BAM_CINS || 
	  bam_cigar_op(bam_get_cigar(bam)[i]) == BAM_CDEL )  {
	has_indels = true; 
	break;
      }
      
      if (bam_cigar_op(bam_get_cigar(bam)[i]) == BAM_CREF_SKIP)
      {
          has_skip = true; // N in CIGAR string, spliced alignment!!!
          break; // skip spliced alignments!!!
      }

      if (bam_cigar_op(bam_get_cigar(bam)[i]) == BAM_CSOFT_CLIP) {
	if (i == 0)
	  nclipstart = bam_cigar_oplen(bam_get_cigar(bam)[i]);
	else
	  nclipend = bam_cigar_oplen(bam_get_cigar(bam)[i]);
      }
    }

    // discard reads with non-continuous (spliced or with indels) alignments
    if (has_indels || has_skip)
      continue;

    int read_len(bam->core.l_qseq); // length of the read
    read_len -= nclipstart;
    read_len -= nclipend; 
    if ( read_len > maxReadLength || read_len < minReadLength )
       continue; // filter reads > maxReadLen or < minReadLength   


    bool rev_strand = bam_is_rev(bam);
    char read_strand = (rev_strand ? '-' : '+');
    // skip reads clipped at 5p by more than 1 nucleotide
    if (filter5pClipped && rev_strand && nclipend>1)
        continue;
    if (filter5pClipped && !rev_strand && nclipstart>1)
        continue;

    // get chr name for this bam line
    string ref( bam_hdr->target_name[bam->core.tid] );
    //cout << "ref=" << ref << endl;
    // load chr sequence 
    if (ref != curr_ref) {
      cerr << "Encountered new ref: " << ref << "; loading...";
      if (prev_ref.length()>0) {
	//free(ref_seq);
	changed_ref = true;
	prev_ref = curr_ref;
      } else {
	// initialize prev_ref to ref
	prev_ref = ref;
      }
      curr_ref = ref;
      //ref_seq = fai_fetch(fai, ref.c_str(), &ref_len);
      //for(int i=0; i < ref_len; ++i)
      //	ref_seq[i] = toupper(ref_seq[i]);
      cerr << " loaded\n";
      cerr.flush();
    }
    
    int read_pos(bam->core.pos); // 0-based leftmost coordinate (forward strand)
    int read_end = read_pos + read_len - 1;                             
    // process queue: output all positions to the left of current read
    // i.e. with pos < read_pos OR if chr changed, output all positions
    // currently in the queue

    
    if (!rev_strand || changed_ref)
      process_queue(q, read_pos, changed_ref, 
		  changed_ref ? prev_ref : curr_ref);
    if (rev_strand || changed_ref)
      process_queue(qr, read_pos, changed_ref, 
		  changed_ref ? prev_ref : curr_ref);
   
    int numHits = bam_aux2i( bam_aux_get(bam, "NH") );

    double read_weight = 1.0 / numHits;
     
    if ( !changed_ref && (read_pos < read_pos_prev) )
    {
      cerr << "ERROR: Input is unsorted: " << endl;
      cerr << curr_ref << ":" << read_pos << " < " << read_pos_prev << endl;
      return 1;
    }

    if (read_pos != read_pos_prev || changed_ref)
    {
      string int_ref = (changed_ref ? prev_ref : curr_ref);
      // output intervals starting at read_pos_prev
      for (int j=0; j < cnt_interval; ++j)
      {
         int int_end = interval_idx[j];
         double int_cnt = interval_cnt[int_end];
         interval_cnt[int_end]=0; // reset read cnt
         cout << "interval" << "\t"
              << int_ref << "\t"
              << read_pos_prev << "\t" 
              << int_end << "\t" 
              << int_cnt << "\t"
              << "+"
              << endl;
      }
      for (int j=0; j < cnt_interval_rev; ++j)
      {
         int int_end = interval_idx_rev[j];
         double int_cnt = interval_cnt_rev[int_end];
         interval_cnt_rev[int_end]=0; // reset read cnt
         cout << "interval" << "\t"
              << int_ref << "\t"
              << read_pos_prev << "\t" 
              << int_end << "\t" 
              << int_cnt << "\t"
              << "-"
              << endl;
      }
 

      // reset read count
      cnt_interval = 0;
      cnt_interval_rev = 0;
      //memset(&interval_cnt[0], 0, interval_cnt.size() * sizeof(interval_cnt[0]));
      //memset(&interval_cnt_rev[0], 0, interval_cnt_rev.size() * sizeof(interval_cnt_rev[0]));
      
    }   
    //cout << "read_len=" << read_len << endl; 
    //cout << interval_cnt[read_len] << endl;
    if (!rev_strand)
    {
      if (interval_cnt[read_len]==0)
      {
        //cout << "here" << endl;
        interval_idx[cnt_interval++]=read_len;
        interval_cnt[read_len]=0;
      }
      interval_cnt[read_len]+=read_weight;
    }
    else
    {
      if (interval_cnt_rev[read_len]==0)
      {
        //cout << "here" << endl;
        interval_idx_rev[cnt_interval_rev++]=read_len;
        interval_cnt_rev[read_len]=0;
      }
      interval_cnt_rev[read_len]+=read_weight;
    }

    bool rev_read = false;
    if(pair_ends){
    	if(bam->core.flag & 0x80 ){
		rev_read = true;		
	}
    }
    deque<Pileup>::iterator q_it;
    if (!rev_strand)
      q_it = q.begin();
    else
      q_it = qr.begin();
    //deque<Pileup>::iterator q_it( qr.begin() );
    for(int i=0; i < read_len; ++i, ++q_it) {
      int g(read_pos + i);  // genomic position

      if (!rev_strand) 
      {
        if (q_it == q.end()) {
          // grow queue if reached the end of the current queuer
	  q.push_back(Pileup(g)); // add position (read_pos+i) to the queue
	  q_it = q.end()-1; // point q to the queue end
        }
        
      }
      else
      {
        if (q_it == qr.end()) {
	  qr.push_back(Pileup(g)); // add position (read_pos+i) to the queue
	  q_it = qr.end()-1; // point q to the queue end
        }
      }

      if (q_it->pos != g) {
	  cerr << "ASSERT: read pos " << g
	     << " != queue pos " << q_it->pos << "\n";
 	  return 1;
      }

      if(rev_strand^rev_read)
        q_it->strand = '-';
      
      if (!rev_strand)
        q_it->nreads_fwd+=read_weight;
      else
        q_it->nreads_rev+=read_weight;
      q_it->nreads+=read_weight;
      if (numHits==1) q_it->nreads_uniq+=1;
      //++(q_it->nreads); // read coverage
      
    } // for each position in the read
    read_pos_prev = read_pos;  
    read_strand_prev = read_strand; 
  }

  // process all remaining positions in the queue
  process_queue(q, 0, true, curr_ref);
  process_queue(qr, 0, true, curr_ref);

  // output intervals starting at read_pos_prev
  for (int j=0; j < cnt_interval; ++j)
  {
         int int_end = interval_idx[j];
         double int_cnt = interval_cnt[int_end];
         interval_cnt[int_end]=0;
         cout << "interval" << "\t"
              << (changed_ref ? prev_ref : curr_ref) << "\t"
              << read_pos_prev << "\t" 
              << int_end << "\t" 
              << int_cnt << "\t"
              << "+"
              << endl;
  }
  for (int j=0; j < cnt_interval_rev; ++j)
  {
         int int_end = interval_idx_rev[j];
         double int_cnt = interval_cnt_rev[int_end];
         interval_cnt_rev[int_end]=0;
         cout << "interval" << "\t"
              << (changed_ref ? prev_ref : curr_ref) << "\t"
              << read_pos_prev << "\t" 
              << int_end << "\t" 
              << int_cnt << "\t"
              << "-"
              << endl;
  }


  return 0;
} 
