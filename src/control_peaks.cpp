// Control-peak engine

#include <Rcpp.h> //Communication between cpp and R
#include <algorithm> //Contains algs for interval searchign
#include <array>    //Fixed Size
#include <climits>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>      //Dynamic Array
#include <atomic>
#include <thread>
#include <exception>
#include <mutex>


//Constants
static constexpr int MIN_REGION_LENGTH    = 50;
static constexpr int SLIDING_WINDOW_SIZE  = 200;
static constexpr int SLIDING_WINDOW_STEP  = 50;

//Types
enum RegionID : int {
  REG_UTR3   = 0,
  REG_UTR5   = 1,
  REG_CDS    = 2,
  REG_EXON   = 3,
  REG_INTRON = 4,
  REG_N      = 5
};

//Stores genomic interval
struct FlatRegion {
  std::vector<int> starts;
  std::vector<int> ends;
  bool empty() const { return starts.empty(); }
};

//All info for 1 transcript
struct Transcript {
  std::array<FlatRegion, REG_N> regions; //Stores all intervals for biotype
  std::array<bool, REG_N>       has{}; //Stores if transcriopt has biotype
  std::vector<int>              splice_sites;
  int   tx_start = -1;
  int   tx_end   = -1;
  char  strand   = '*';
  std::string gene;             //gene name
};

//Transcript Name -> Transcript Object
using TranscriptMap = std::unordered_map<std::string, Transcript>;

//Gene interval
//Gene Name -> Gene Object
struct GeneRec { int start; int end; };
using GeneMap = std::unordered_map<std::string, GeneRec>;

//Input per chromosome, built by main thread from Rcpp inputs
//For thread safety
struct ChromInput {
  std::vector<int>         anno_start;
  std::vector<int>         anno_end;
  std::vector<std::string> anno_transcript;
  std::vector<std::string> anno_region;
  std::vector<std::string> anno_strand;
  std::vector<std::string> anno_gene;
  std::vector<std::string> gene_name;
  std::vector<int>         gene_start;
  std::vector<int>         gene_end;
  std::vector<int>         peak_start;
  std::vector<int>         peak_end;
  std::vector<std::string> peak_FC;
  std::vector<std::vector<std::string>> peak_transcripts;
  uint32_t seed = 0;
};

//Output per chromosome, converted to Rcpp by main thread
struct ChromOutput {
  std::vector<int>         peak_start;
  std::vector<int>         peak_end;
  std::vector<int>         control_start;
  std::vector<int>         control_end;
  std::vector<std::string> name;
  std::vector<std::string> strand;
};


// Region to name helpers
static const char* region_cstr(RegionID r) {
  switch (r) {
    case REG_UTR3:   return "UTR3";
    case REG_UTR5:   return "UTR5";
    case REG_CDS:    return "CDS";
    case REG_EXON:   return "exon";
    case REG_INTRON: return "intron";
    default:         return "?";
  }
}

static RegionID region_from_cstr(const char* s) {
  if (!s) return REG_N;
  if (std::strcmp(s, "UTR3") == 0) return REG_UTR3;
  if (std::strcmp(s, "UTR5") == 0) return REG_UTR5;
  if (std::strcmp(s, "CDS")  == 0) return REG_CDS;
  if (std::strcmp(s, "exon") == 0) return REG_EXON;
  return REG_N; // unknown / not stored
}


// Flat-vector interval helpers
// Checks if 1 interval overlaps against any other interval
static bool flat_overlaps_any_scalar(int qs, int qe,
                                     const std::vector<int>& starts,
                                     const std::vector<int>& ends) {
  const int n = (int)starts.size();
  if (n == 0) return false;

  // upper_bound gives count of elements <= qs.
  int i = (int)(std::upper_bound(starts.begin(), starts.end(), qs) - starts.begin());
  if (i >= 1 && ends[i - 1] > qs) return true;
  if (i + 1 <= n && starts[i] < qe) return true;
  return false;
}

// Check if interval is fully contained in another
static bool flat_contained_scalar(int qs, int qe,
                                  const std::vector<int>& starts,
                                  const std::vector<int>& ends) {
  const int n = (int)starts.size();
  if (n == 0) return false;
  int i = (int)(std::upper_bound(starts.begin(), starts.end(), qs) - starts.begin());
  if (i < 1) return false;
  return qe <= ends[i - 1];
}

// Batch Contained Check
static std::vector<char> flat_contained_batch(const std::vector<int>& qs,
                                              const std::vector<int>& qe,
                                              const std::vector<int>& starts,
                                              const std::vector<int>& ends) {
  std::vector<char> out(qs.size(), 0);
  const int n = (int)starts.size();
  if (n == 0) return out;
  for (size_t k{}; k < qs.size(); ++k) {
    int i = (int)(std::upper_bound(starts.begin(), starts.end(), qs[k]) - starts.begin());
    if (i >= 1 && qe[k] <= ends[i - 1]) out[k] = 1;
  }
  return out;
}

// splice_slice: take every-other element of splice_sites starting at
// zero_based (0/1).
static std::vector<int> splice_slice(const std::vector<int>& splice_sites,
                                     int zero_based) {
  std::vector<int> out;
  const int n = (int)splice_sites.size();
  if (n == 0) return out;
  int first = zero_based;
  if (first >= n) return out;
  out.reserve((n - first + 1) / 2);
  for (int i = first; i < n; i += 2) out.push_back(splice_sites[i]);
  return out;
}

// Intersect a single (a_s, a_e) with many intervals
static FlatRegion flat_intersect_one(int a_s, int a_e,
                                     const std::vector<int>& b_s,
                                     const std::vector<int>& b_e) {
  FlatRegion out;
  if (b_s.empty()) return out;
  for (size_t i{}; i < b_s.size(); ++i) {
    if (b_s[i] < a_e && b_e[i] > a_s) {
      out.starts.push_back(std::max(b_s[i], a_s));
      out.ends.push_back(std::min(b_e[i], a_e));
    }
  }
  return out;
}

// Set difference: a minus b.
static FlatRegion flat_setdiff(const std::vector<int>& a_s,
                               const std::vector<int>& a_e,
                               const std::vector<int>& b_s,
                               const std::vector<int>& b_e) {
  FlatRegion out;
  const size_t na = a_s.size();
  if (na == 0) return out;
  if (b_s.empty()) {
    out.starts = a_s;
    out.ends   = a_e;
    return out;
  }
  out.starts.reserve(na + b_s.size());
  out.ends.reserve(na + b_s.size());

  for (size_t i{}; i < na; ++i) {
    int cur = a_s[i];
    const int a_end = a_e[i];

    // ov is all B intervals that overlap with this A interval
    std::vector<size_t> ov;
    for (size_t jk = 0; jk < b_s.size(); ++jk) {
      if (b_s[jk] < a_end && b_e[jk] > a_s[i]) ov.push_back(jk);
    }

    if (ov.empty()) {
      out.starts.push_back(cur);
      out.ends.push_back(a_end);
      continue;
    }
    //Loop through each overlapping interval, keep space before interval B
    for (size_t jk : ov) {
      const int js = b_s[jk];
      const int je = b_e[jk];
      if (js > cur) {
        out.starts.push_back(cur);
        out.ends.push_back(js);
      }
      if (je > cur) cur = je;
      if (cur >= a_end) break;  //if it consumed all of A then stop
    }
    //Keep leftover space after interval B
    if (cur < a_end) {
      out.starts.push_back(cur);
      out.ends.push_back(a_end);
    }
  }
  return out;
}

// Merge a small unsorted set of intervals into disjoint sorted intervals.
// order(s, e) ascending.
//Reduce Function
static FlatRegion flat_reduce(std::vector<int> s, std::vector<int> e) {
  FlatRegion out;
  const size_t n = s.size();
  if (n == 0) return out;

  //Sort indicies by start and then by end
  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
    if (s[a] != s[b]) return s[a] < s[b];
    return e[a] < e[b];
  });

  int cur_s = s[idx[0]];
  int cur_e = e[idx[0]];
  //Loop through intervals in sorted order and merge intervals that overlap
  for (size_t k{1}; k < n; ++k) {
    const int sk = s[idx[k]];
    const int ek = e[idx[k]];
    if (sk <= cur_e) {
      if (ek > cur_e) cur_e = ek;
    } else {
      out.starts.push_back(cur_s);
      out.ends.push_back(cur_e);
      cur_s = sk;
      cur_e = ek;
    }
  }
  out.starts.push_back(cur_s);
  out.ends.push_back(cur_e);
  return out;
}


// Sliding windows
struct Windows { std::vector<int> starts; std::vector<int> ends; };
//Splits long peaks into windows
static Windows sliding_window_python(int start, int end,
                                     int window_size = SLIDING_WINDOW_SIZE,
                                     int min_len     = SLIDING_WINDOW_STEP) {
  Windows w;
  const int seq_len_total = end - start + 1;
  if (seq_len_total <= 0) return w;

  //Loop through sequence and create leftover winodws
  int last_start = INT_MIN;
  int last_end   = INT_MIN;
  int i = 0;
  while (i < seq_len_total) {
    if (i + window_size > seq_len_total) {
      if (seq_len_total - i > min_len) {
        last_start = start + i;
        last_end   = end;
      }
    } else {
      last_start = start + i;
      last_end   = start + i + window_size;
    }
    w.starts.push_back(last_start);
    w.ends.push_back(last_end);
    i += window_size;
  }
  return w;
}


// Annotation builders
//Recieves column vectors from R
static TranscriptMap build_transcript_annotations(
    const std::vector<std::string>& anno_transcript,
    const std::vector<std::string>& anno_region,
    const std::vector<int>&         anno_start,
    const std::vector<int>&         anno_end,
    const std::vector<std::string>& anno_strand,
    const std::vector<std::string>& anno_gene) {
  TranscriptMap out;
  const size_t n = anno_transcript.size();
  if (n == 0) return out;

  // Temp transcript: raw (unreduced) intervals per non-intron
  // region + the raw exon starts/ends used to derive splice sites + introns.
  struct TmpTx {
    std::array<std::vector<int>, REG_N> rs;
    std::array<std::vector<int>, REG_N> re;
    std::array<bool, REG_N> has_raw{};
    std::vector<int> exon_s;
    std::vector<int> exon_e;
    char strand = '*';
    std::string gene;
  };
  std::unordered_map<std::string, TmpTx> tmp;
  tmp.reserve(n / 4);

  // Pass 1: file order. gene gets last-write-wins to mirror Python.
  for (size_t i{}; i < n; ++i) {
    const std::string& tx_s   = anno_transcript[i];
    const std::string& reg_s  = anno_region[i];
    const std::string& gen_s  = anno_gene[i];
    const std::string& str_s  = anno_strand[i];
    const int   s      = anno_start[i];
    const int   e      = anno_end[i];

    TmpTx& t = tmp[tx_s]; //Map tmp object
    if (t.strand == '*' && !str_s.empty()) t.strand = str_s[0]; //Store strand
    t.gene = gen_s; //Store gene name

    RegionID rid = region_from_cstr(reg_s.c_str()); //Convert Region String to ID
    if (rid == REG_N) continue; // skip unknown / non-stored regions

    //Save Intervals
    t.rs[rid].push_back(s);
    t.re[rid].push_back(e);
    t.has_raw[rid] = true;

    //Save splice regions if its an Exon
    if (rid == REG_EXON) {
      t.exon_s.push_back(s);
      t.exon_e.push_back(e);
    }
  }

  // Pass 2: reduce per (tx, region), derive splice sites + introns.
  for (auto& kv : tmp) {
    Transcript& t = out[kv.first];
    t.strand = kv.second.strand;
    t.gene   = kv.second.gene;

    for (int r{}; r < REG_INTRON; ++r) {
      if (!kv.second.has_raw[r]) continue;
      t.regions[r] = flat_reduce(std::move(kv.second.rs[r]),
                                 std::move(kv.second.re[r]));
      t.has[r] = !t.regions[r].empty();
    }

    if (!kv.second.exon_s.empty()) {
      // Splice sites = sort(concat(exon_starts, exon_ends)).
      std::vector<int> sites;
      sites.reserve(kv.second.exon_s.size() * 2);
      sites.insert(sites.end(), kv.second.exon_s.begin(), kv.second.exon_s.end());
      sites.insert(sites.end(), kv.second.exon_e.begin(), kv.second.exon_e.end());
      std::sort(sites.begin(), sites.end());
      t.splice_sites = std::move(sites);
      t.tx_start = t.splice_sites.front();
      t.tx_end   = t.splice_sites.back();

      // Introns: pairs taken from sites[1..n-2]
      if (t.splice_sites.size() > 2) {
        std::vector<int> i_s, i_e;
        for (size_t k = 1; k + 1 < t.splice_sites.size(); ++k) {
          const size_t mid_idx = k - 1;
          if (mid_idx % 2 == 0) i_s.push_back(t.splice_sites[k]);
          else                  i_e.push_back(t.splice_sites[k]);
        }
        const size_t pair_n = std::min(i_s.size(), i_e.size());
        i_s.resize(pair_n);
        i_e.resize(pair_n);
        if (!i_s.empty()) {
          t.regions[REG_INTRON] = flat_reduce(std::move(i_s), std::move(i_e));
          t.has[REG_INTRON] = !t.regions[REG_INTRON].empty();
        }
      }
    }
  }

  return out;
}

//Build Gene Map
static GeneMap build_gene_map(const std::vector<std::string>& gene_name,
                              const std::vector<int>&         gene_start,
                              const std::vector<int>&         gene_end) {
  GeneMap out;
  const size_t n = gene_name.size();
  out.reserve(n);
  for (size_t i{}; i < n; ++i) {
    if (gene_name[i].empty()) continue;
    out[gene_name[i]] = GeneRec{gene_start[i], gene_end[i]};
  }
  return out;
}


// Per-peak qualification
static std::vector<RegionID> find_overlap_for_peak(int peak_start, int peak_end,
                                                   const Transcript& tx) {

  std::vector<RegionID> out;
  static const RegionID order_arr[] = { REG_UTR3, REG_UTR5, REG_CDS, REG_EXON, REG_INTRON };
  for (RegionID r : order_arr) {
    if (!tx.has[r]) continue;
    const FlatRegion& rf = tx.regions[r];
    if (rf.empty()) continue;
    if (r == REG_INTRON) {
      if (flat_overlaps_any_scalar(peak_start, peak_end, rf.starts, rf.ends)) {
        out.push_back(r);
      }
    } else {
      if (flat_contained_scalar(peak_start, peak_end, rf.starts, rf.ends)) {
        out.push_back(r);
      }
    }
  }
  return out;
}

//Peaks choosen in order CDS > UTR3 > UTR5 > exon > intron
static RegionID choose_peak_region(const std::array<bool, REG_N>& present) {
  if (present[REG_CDS])  return REG_CDS;
  if (present[REG_UTR3]) return REG_UTR3;
  if (present[REG_UTR5]) return REG_UTR5;
  if (present[REG_EXON]) return REG_EXON;
  return REG_INTRON;
}


// Control sampling
struct Control {
  bool ok = false;
  int control_start = 0;
  int control_end   = 0;
  std::string region_type;
};

struct PickFlat {
  bool ok = false;
  int start = 0;
  int end   = 0;
};

//Random integer in [0, n).
static int rng_pick(int n, std::mt19937 &rng) {
  std::uniform_int_distribution<int> dist(0, n - 1);
  return dist(rng);
}

//Picks random interval from [starts,ends)
static PickFlat sample_uniform_flat(const std::vector<int>& starts,
                                    const std::vector<int>& ends, std::mt19937 &rng) {
  PickFlat p;
  if (starts.empty()) return p;
  const int i = rng_pick((int)starts.size(), rng);
  p.start = starts[i];
  p.end   = ends[i];
  p.ok    = true;
  return p;
}

//Simple control sampling for UTR/Exon Regions
static Control sample_simple_control(int peak_start, int peak_end, int peak_length,
                                     const std::vector<int>& cl_starts,
                                     const std::vector<int>& cl_ends,
                                     const char* final_region,
                                     std::mt19937 &rng) {
  Control c;
  if (cl_starts.empty()) return c;

  std::vector<int> ks;  // candidate region starts
  std::vector<int> kn;  // n_starts per kept region
  ks.reserve(cl_starts.size());
  kn.reserve(cl_starts.size());
  //Check if peak length can fit within interval
  for (size_t i{}; i < cl_starts.size(); ++i) {
    const int n_s = cl_ends[i] - peak_length - cl_starts[i];
    if (n_s > 0) {
      ks.push_back(cl_starts[i]);
      kn.push_back(n_s);
    }
  }
  if (ks.empty()) return c;

  //Total number of start positions possible across all intervals
  long long total = 0;
  for (int v : kn) total += v;
  if (total <= 0 || total > INT_MAX) return c;


  //Choose random start position in [0, total).
  const int chosen_global = rng_pick((int)total, rng);

  long long cumulative = 0;
  int chosen_row = -1;
  long long previous_total = 0;
  //Find which interval contains randomly choosen position
  for (size_t i = 0; i < kn.size(); ++i) {
    const long long prev = cumulative;
    cumulative += kn[i];
    if (cumulative > chosen_global) {
      chosen_row = (int)i;
      previous_total = prev;
      break;
    }
  }
  if (chosen_row < 0) return c;

  const int offset = chosen_global - (int)previous_total;
  c.control_start = ks[chosen_row] + offset;
  c.control_end   = c.control_start + peak_length;
  c.region_type   = final_region;
  c.ok            = true;
  (void)peak_start; (void)peak_end;
  return c;
}

static Control get_cds_control(int peak_start, int peak_end, int peak_length,
                               const std::vector<int>& cl_starts,
                               const std::vector<int>& cl_ends,
                               const FlatRegion& rf, const Transcript& tx, std::mt19937 &rng) {
  Control c;
  // Last overlapping CDS interval (Python's "last hit wins").
  int k = -1;
  for (int i = (int)rf.starts.size() - 1; i >= 0; --i) {
    if (rf.starts[i] < peak_end && rf.ends[i] > peak_start) { k = i; break; }
  }
  if (k < 0) return c;

  const int ec_s = rf.starts[k];
  const int ec_e = rf.ends[k];
  const int dist_to_start = peak_start - ec_s;
  const int dist_to_end   = ec_e - peak_end;

  std::vector<int> ss;
  int distance;
  if (dist_to_start < dist_to_end) {
    ss = splice_slice(tx.splice_sites, 1);
    distance = dist_to_start;
  } else {
    ss = splice_slice(tx.splice_sites, 0);
    distance = dist_to_end;
  }
  if (ss.empty()) return c;

  std::vector<int> cs(ss.size()), ce(ss.size());
  for (size_t i = 0; i < ss.size(); ++i) {
    cs[i] = ss[i] + distance;
    ce[i] = cs[i] + peak_length;
  }

  auto ok = flat_contained_batch(cs, ce, cl_starts, cl_ends);
  std::vector<int> keep_s, keep_e;
  for (size_t i = 0; i < ok.size(); ++i) {
    if (ok[i]) { keep_s.push_back(cs[i]); keep_e.push_back(ce[i]); }
  }
  if (keep_s.empty()) return c;

  FlatRegion red = flat_reduce(std::move(keep_s), std::move(keep_e));
  PickFlat p = sample_uniform_flat(red.starts, red.ends, rng);
  if (!p.ok) return c;

  c.control_start = p.start;
  c.control_end   = p.end;
  c.region_type   = "CDS";
  c.ok            = true;
  return c;
}

static Control get_intron_control(int peak_start, int peak_end, int peak_length,
                                  const std::vector<int>& cl_starts,
                                  const std::vector<int>& cl_ends,
                                  const FlatRegion& rf, const Transcript& tx, std::mt19937 &rng) {
  Control c;
  int region_start = -1, region_end = -1;
  int intersect_length = 0;
  bool is_exon_overlap = false;
  bool is_5_prime_ss   = false;

  // Last overlapping intron.
  for (int i = (int)rf.starts.size() - 1; i >= 0; --i) {
    if (rf.starts[i] < peak_end && rf.ends[i] > peak_start) {
      region_start = rf.starts[i];
      region_end   = rf.ends[i];
      intersect_length = std::min(peak_end, rf.ends[i]) - std::max(peak_start, rf.starts[i]);
      break;
    }
  }

  if (intersect_length == 0) {
    // Try the last overlapping exon instead.
    if (tx.has[REG_EXON]) {
      const FlatRegion& erf = tx.regions[REG_EXON];
      for (int i = (int)erf.starts.size() - 1; i >= 0; --i) {
        if (erf.starts[i] < peak_end && erf.ends[i] > peak_start) {
          region_start = erf.starts[i];
          region_end   = erf.ends[i];
          intersect_length = std::min(peak_end, erf.ends[i]) - std::max(peak_start, erf.starts[i]);
          is_exon_overlap = true;
          break;
        }
      }
    }
  }

  if (intersect_length == 0 || region_start < 0) return c;

  const int dist_to_start = peak_start - region_start;
  const int dist_to_end   = region_end - peak_end;

  std::vector<int> ss;
  std::string final_region;
  int distance = INT_MIN;
  bool have_dist = false;
  const char p_plus = (tx.strand == '+');

  if (peak_length == intersect_length) {
    if (dist_to_start < dist_to_end) {
      if (is_exon_overlap) {
        ss = splice_slice(tx.splice_sites, 0);
        final_region = p_plus ? "3pexon" : "5pexon";
        distance = dist_to_start; have_dist = true;
      } else {
        ss = splice_slice(tx.splice_sites, 1);
        final_region = p_plus ? "5pintron" : "3pintron";
        distance = dist_to_start; have_dist = true;
      }
    } else {
      if (is_exon_overlap) {
        ss = splice_slice(tx.splice_sites, 1);
        final_region = p_plus ? "5pexon" : "3pexon";
        distance = dist_to_end; have_dist = true;
        is_5_prime_ss = true;
      } else {
        ss = splice_slice(tx.splice_sites, 2);
        final_region = p_plus ? "3pintron" : "5pintron";
        distance = dist_to_end; have_dist = true;
        is_5_prime_ss = true;
      }
    }
  } else {
    // Two independent ifs — both can fire (the second overwrites the first).
    if (dist_to_start < 0) {
      ss = splice_slice(tx.splice_sites, 1);
      final_region = p_plus ? "5pss" : "3pss";
      distance = dist_to_start; have_dist = true;
    }
    if (dist_to_end < 0) {
      ss = splice_slice(tx.splice_sites, 2);
      final_region = p_plus ? "3pss" : "5pss";
      distance = dist_to_end; have_dist = true;
      is_5_prime_ss = true;
    }
  }

  if (final_region.empty() || ss.empty() || !have_dist) return c;

  std::vector<int> cs(ss.size()), ce(ss.size());
  if (is_5_prime_ss) {
    for (size_t i = 0; i < ss.size(); ++i) {
      const int cp = ss[i] - distance;
      cs[i] = cp - peak_length;
      ce[i] = cp;
    }
  } else {
    for (size_t i = 0; i < ss.size(); ++i) {
      const int cp = ss[i] + distance;
      cs[i] = cp;
      ce[i] = cp + peak_length;
    }
  }

  auto ok = flat_contained_batch(cs, ce, cl_starts, cl_ends);
  std::vector<int> keep_s, keep_e;
  for (size_t i = 0; i < ok.size(); ++i) {
    if (ok[i]) { keep_s.push_back(cs[i]); keep_e.push_back(ce[i]); }
  }
  if (keep_s.empty()) return c;

  FlatRegion red = flat_reduce(std::move(keep_s), std::move(keep_e));
  PickFlat p = sample_uniform_flat(red.starts, red.ends, rng);
  if (!p.ok) return c;

  c.control_start = p.start;
  c.control_end   = p.end;
  c.region_type   = std::move(final_region);
  c.ok            = true;
  return c;
}

static Control get_control(int peak_start, int peak_end, RegionID region,
                           const std::vector<int>& cl_starts,
                           const std::vector<int>& cl_ends,
                           const Transcript& tx, std::mt19937 &rng) {
  Control c;
  const int peak_length = peak_end - peak_start;
  if (peak_length < MIN_REGION_LENGTH) return c;
  if (peak_length <= 0) return c;
  if (cl_starts.empty()) return c;

  static const FlatRegion empty_rf;
  const FlatRegion& rf = tx.has[region] ? tx.regions[region] : empty_rf;

  if (region == REG_UTR3 || region == REG_UTR5 || region == REG_EXON) {
    return sample_simple_control(peak_start, peak_end, peak_length,
                                 cl_starts, cl_ends, region_cstr(region), rng);
  }
  if (region == REG_CDS) {
    return get_cds_control(peak_start, peak_end, peak_length,
                           cl_starts, cl_ends, rf, tx, rng);
  }
  if (region == REG_INTRON) {
    return get_intron_control(peak_start, peak_end, peak_length,
                              cl_starts, cl_ends, rf, tx, rng);
  }
  return c;
}

// -----------------------------
// Sorted insert helper (mirrors R's append at findInterval position).
// -----------------------------

static void sorted_insert(std::vector<int>& starts, std::vector<int>& ends,
                          int cs, int ce) {
  // findInterval(cs, starts) returns count of elements <= cs. Insert after.
  const auto it = std::upper_bound(starts.begin(), starts.end(), cs);
  const size_t pos = it - starts.begin();
  starts.insert(starts.begin() + pos, cs);
  ends.insert(ends.begin() + pos, ce);
}

// -----------------------------
// Per-chromosome driver
// -----------------------------
// Pure-C++ per-chromosome kernel. Touches no R / Rcpp state — safe to call
// from worker threads concurrently as long as each call gets its own input.
static ChromOutput process_chromosome_core(const ChromInput& in) {
  ChromOutput out;
  std::mt19937 rng(in.seed);

  const int n_peaks = (int)in.peak_start.size();
  if (n_peaks == 0 || in.anno_transcript.empty()) return out;

  TranscriptMap tx_map = build_transcript_annotations(
      in.anno_transcript, in.anno_region, in.anno_start, in.anno_end,
      in.anno_strand, in.anno_gene);
  GeneMap gene_map = build_gene_map(in.gene_name, in.gene_start, in.gene_end);

  // Pass 1: build peak_region_dict + per-strand peak totals.
  struct Bucket {
    std::array<std::vector<std::string>, REG_N> tx_per_region;
    std::array<bool, REG_N> present{};
    bool any = false;
  };
  std::vector<Bucket> buckets(n_peaks);
  std::vector<int> peak_order; peak_order.reserve(n_peaks);

  std::vector<int> plus_s, plus_e, minus_s, minus_e;
  plus_s.reserve(n_peaks);  plus_e.reserve(n_peaks);
  minus_s.reserve(n_peaks); minus_e.reserve(n_peaks);

  for (int p_i = 0; p_i < n_peaks; ++p_i) {
    const int ps = in.peak_start[p_i];
    const int pe = in.peak_end[p_i];
    if (pe <= ps) continue;

    for (const std::string& tx_name : in.peak_transcripts[p_i]) {
      auto it = tx_map.find(tx_name);
      if (it == tx_map.end()) continue;
      const Transcript& tx = it->second;

      if (tx.strand == '+') { plus_s.push_back(ps);  plus_e.push_back(pe); }
      else                  { minus_s.push_back(ps); minus_e.push_back(pe); }

      auto qualifying = find_overlap_for_peak(ps, pe, tx);
      if (qualifying.empty()) continue;

      Bucket& b = buckets[p_i];
      if (!b.any) {
        b.any = true;
        peak_order.push_back(p_i);
      }
      for (RegionID r : qualifying) {
        b.present[r] = true;
        b.tx_per_region[r].push_back(tx_name);
      }
    }
  }

  if (peak_order.empty()) return out;

  FlatRegion tpp = flat_reduce(std::move(plus_s),  std::move(plus_e));
  FlatRegion tpm = flat_reduce(std::move(minus_s), std::move(minus_e));

  // Pass 2: per peak, pick a region, walk windows, sample.
  std::vector<int> tcp_s, tcp_e, tcm_s, tcm_e;
  std::unordered_set<long long> processed_windows;

  out.peak_start.reserve(peak_order.size());
  out.peak_end.reserve(peak_order.size());
  out.control_start.reserve(peak_order.size());
  out.control_end.reserve(peak_order.size());
  out.name.reserve(peak_order.size());
  out.strand.reserve(peak_order.size());

  for (int p_i : peak_order) {
    const Bucket& b = buckets[p_i];
    const RegionID peak_region = choose_peak_region(b.present);

    const int peak_s = in.peak_start[p_i];
    const int peak_e = in.peak_end[p_i];
    const std::string& fc = in.peak_FC[p_i];
    const int peak_length = peak_e - peak_s;

    const std::vector<std::string>& peak_txs = b.tx_per_region[peak_region];

    for (const std::string& tx_name : peak_txs) {
      auto it = tx_map.find(tx_name);
      if (it == tx_map.end()) continue;
      const Transcript& tx = it->second;

      const bool plus = (tx.strand == '+');
      const std::vector<int>& tp_s = plus ? tpp.starts : tpm.starts;
      const std::vector<int>& tp_e = plus ? tpp.ends   : tpm.ends;

      std::vector<int> bl_s, bl_e;
      auto gene_it = gene_map.find(tx.gene);
      if (gene_it != gene_map.end()) {
        FlatRegion bl = flat_intersect_one(gene_it->second.start,
                                           gene_it->second.end,
                                           tp_s, tp_e);
        bl_s = std::move(bl.starts);
        bl_e = std::move(bl.ends);
      } else {
        bl_s = tp_s;
        bl_e = tp_e;
      }

      std::vector<int> ri_s, ri_e;
      if (peak_region == REG_INTRON) {
        if (tx.tx_start < 0 || tx.tx_end < 0) continue;
        ri_s.push_back(tx.tx_start);
        ri_e.push_back(tx.tx_end);
      } else if (tx.has[peak_region]) {
        ri_s = tx.regions[peak_region].starts;
        ri_e = tx.regions[peak_region].ends;
      }

      FlatRegion cl = flat_setdiff(ri_s, ri_e, bl_s, bl_e);

      Windows ws;
      if (peak_length <= 400) {
        ws.starts.push_back(peak_s);
        ws.ends.push_back(peak_e);
      } else {
        ws = sliding_window_python(peak_s, peak_e);
      }

      bool peak_processed = false;

      for (size_t w_i = 0; w_i < ws.starts.size(); ++w_i) {
        const int ws_s = ws.starts[w_i];
        const int ws_e = ws.ends[w_i];
        const long long key = (long long)ws_s * 100000000LL + (long long)ws_e;
        if (processed_windows.count(key)) continue;

        Control ctrl = get_control(ws_s, ws_e, peak_region, cl.starts, cl.ends, tx, rng);
        if (!ctrl.ok) continue;
        if (ctrl.control_start == 0) continue;

        const int cs = ctrl.control_start;
        const int ce = ctrl.control_end;

        if (flat_overlaps_any_scalar(cs, ce, tp_s, tp_e)) continue;
        if (plus) {
          if (flat_overlaps_any_scalar(cs, ce, tcp_s, tcp_e)) continue;
        } else {
          if (flat_overlaps_any_scalar(cs, ce, tcm_s, tcm_e)) continue;
        }

        if (plus) sorted_insert(tcp_s, tcp_e, cs, ce);
        else      sorted_insert(tcm_s, tcm_e, cs, ce);

        processed_windows.insert(key);

        out.peak_start.push_back(ws_s);
        out.peak_end.push_back(ws_e);
        out.control_start.push_back(cs);
        out.control_end.push_back(ce);
        out.name.push_back(tx.gene + "_" + ctrl.region_type + "_" + fc);
        out.strand.push_back(std::string(1, tx.strand));

        peak_processed = true;
      }

      if (peak_processed) break;
    }
  }

  return out;
}

/**
 * Converting between RCpp and C++ core to allow for multhreading
 */

//Rcpp to plain conversion helpers
static std::vector<std::string> to_string_vec(const Rcpp::CharacterVector x){
  std::vector<std::string> out;
  const int n = x.size();
  out.reserve(n);
  for (int i{}; i<n; ++i){
    if (SEXP(x[i]) == NA_STRING) out.emplace_back();
    else out.emplace_back(Rcpp::as<std::string>(x[i]));
  }
  return out;
}

static std::vector<int> to_int_vec(const Rcpp::IntegerVector& x) {
  return std::vector<int>(x.begin(), x.end());
}

//Build C++ Chrom Input from Rcpp
static ChromInput build_chrom_input(
    const Rcpp::CharacterVector& anno_transcript,
    const Rcpp::CharacterVector& anno_region,
    const Rcpp::IntegerVector&   anno_start,
    const Rcpp::IntegerVector&   anno_end,
    const Rcpp::CharacterVector& anno_strand,
    const Rcpp::CharacterVector& anno_gene,
    const Rcpp::CharacterVector& gene_name,
    const Rcpp::IntegerVector&   gene_start,
    const Rcpp::IntegerVector&   gene_end,
    const Rcpp::IntegerVector&   peak_start,
    const Rcpp::IntegerVector&   peak_end,
    const Rcpp::CharacterVector& peak_FC,
    const Rcpp::List&            peak_transcripts,
    uint32_t                     seed) {
  ChromInput in;
  in.anno_transcript = to_string_vec(anno_transcript);
  in.anno_region     = to_string_vec(anno_region);
  in.anno_start      = to_int_vec(anno_start);
  in.anno_end        = to_int_vec(anno_end);
  in.anno_strand     = to_string_vec(anno_strand);
  in.anno_gene       = to_string_vec(anno_gene);
  in.gene_name       = to_string_vec(gene_name);
  in.gene_start      = to_int_vec(gene_start);
  in.gene_end        = to_int_vec(gene_end);
  in.peak_start      = to_int_vec(peak_start);
  in.peak_end        = to_int_vec(peak_end);
  in.peak_FC         = to_string_vec(peak_FC);

  const int n_peaks = peak_transcripts.size();
  in.peak_transcripts.resize(n_peaks);
  for (int i{}; i < n_peaks; ++i) {
    Rcpp::CharacterVector v = peak_transcripts[i];
    in.peak_transcripts[i] = to_string_vec(v);
  }
  in.seed = seed;
  return in;
}

//Build Rcpp Chrom Output from C++
static Rcpp::List wrap_chrom_output(const ChromOutput& out) {
  const int nr = (int)out.peak_start.size();
  Rcpp::IntegerVector   o_ps(nr), o_pe(nr), o_cs(nr), o_ce(nr);
  Rcpp::CharacterVector o_nm(nr), o_st(nr);
  for (int i = 0; i < nr; ++i) {
    o_ps[i] = out.peak_start[i];
    o_pe[i] = out.peak_end[i];
    o_cs[i] = out.control_start[i];
    o_ce[i] = out.control_end[i];
    o_nm[i] = out.name[i];
    o_st[i] = out.strand[i];
  }
  return Rcpp::List::create(
    Rcpp::Named("peak_start")    = o_ps,
    Rcpp::Named("peak_end")      = o_pe,
    Rcpp::Named("name")          = o_nm,
    Rcpp::Named("strand")        = o_st,
    Rcpp::Named("control_start") = o_cs,
    Rcpp::Named("control_end")   = o_ce
  );
}




// Check the per-chromosome list shape on the main thread before spawning
// workers. Names the chromosome and the bad field in the error.
static void validate_chrom_input(const Rcpp::List& el) {
  static const char* required[] = {
    "anno_transcript", "anno_region", "anno_start", "anno_end",
    "anno_strand", "anno_gene",
    "gene_name", "gene_start", "gene_end",
    "peak_start", "peak_end", "peak_FC", "peak_transcripts", "seed", "chrom"
  };
  Rcpp::CharacterVector names = el.names();
  std::unordered_set<std::string> have;
  have.reserve(names.size());
  for (int i = 0; i < names.size(); ++i) have.insert(Rcpp::as<std::string>(names[i]));

  std::string chrom = have.count("chrom")
    ? Rcpp::as<std::string>(el["chrom"]) : std::string("<unknown>");

  for (const char* nm : required) {
    if (!have.count(nm)) {
      Rcpp::stop("control peak engine: chromosome '%s' is missing field '%s'.",
                 chrom, nm);
    }
  }

  auto check_len = [&](const char* a, R_xlen_t la,
                       const char* b, R_xlen_t lb) {
    if (la != lb) {
      Rcpp::stop("control peak engine: chromosome '%s' has length mismatch "
                 "between '%s' (%lld) and '%s' (%lld).",
                 chrom, a, (long long)la, b, (long long)lb);
    }
  };

  const R_xlen_t n_anno = Rf_xlength(el["anno_transcript"]);
  check_len("anno_transcript", n_anno, "anno_region",  Rf_xlength(el["anno_region"]));
  check_len("anno_transcript", n_anno, "anno_start",   Rf_xlength(el["anno_start"]));
  check_len("anno_transcript", n_anno, "anno_end",     Rf_xlength(el["anno_end"]));
  check_len("anno_transcript", n_anno, "anno_strand",  Rf_xlength(el["anno_strand"]));
  check_len("anno_transcript", n_anno, "anno_gene",    Rf_xlength(el["anno_gene"]));

  const R_xlen_t n_gene = Rf_xlength(el["gene_name"]);
  check_len("gene_name", n_gene, "gene_start", Rf_xlength(el["gene_start"]));
  check_len("gene_name", n_gene, "gene_end",   Rf_xlength(el["gene_end"]));

  const R_xlen_t n_peak = Rf_xlength(el["peak_start"]);
  check_len("peak_start", n_peak, "peak_end",         Rf_xlength(el["peak_end"]));
  check_len("peak_start", n_peak, "peak_FC",          Rf_xlength(el["peak_FC"]));
  check_len("peak_start", n_peak, "peak_transcripts", Rf_xlength(el["peak_transcripts"]));

  if (TYPEOF(el["peak_transcripts"]) != VECSXP) {
    Rcpp::stop("control peak engine: chromosome '%s' 'peak_transcripts' must "
               "be a list, got SEXP type %d.",
               chrom, TYPEOF(el["peak_transcripts"]));
  }
}

//Multi-threaded entry point
// [[Rcpp::export]]
Rcpp::List process_chromosomes_threaded_cpp(Rcpp::List per_chrom, int max_threads) {
  const int n_chrom = per_chrom.size();
  if (n_chrom == 0) return Rcpp::List::create();

  //Main thread Rcpp to C++ conversion. After this loop `inputs` is read-only;
  //workers only read from it.
  std::vector<ChromInput> inputs;
  inputs.reserve(n_chrom);
  for (int i{}; i < n_chrom; ++i) {
    Rcpp::List el = per_chrom[i];
    validate_chrom_input(el);
    inputs.push_back(build_chrom_input(
        Rcpp::as<Rcpp::CharacterVector>(el["anno_transcript"]),
        Rcpp::as<Rcpp::CharacterVector>(el["anno_region"]),
        Rcpp::as<Rcpp::IntegerVector>(el["anno_start"]),
        Rcpp::as<Rcpp::IntegerVector>(el["anno_end"]),
        Rcpp::as<Rcpp::CharacterVector>(el["anno_strand"]),
        Rcpp::as<Rcpp::CharacterVector>(el["anno_gene"]),
        Rcpp::as<Rcpp::CharacterVector>(el["gene_name"]),
        Rcpp::as<Rcpp::IntegerVector>(el["gene_start"]),
        Rcpp::as<Rcpp::IntegerVector>(el["gene_end"]),
        Rcpp::as<Rcpp::IntegerVector>(el["peak_start"]),
        Rcpp::as<Rcpp::IntegerVector>(el["peak_end"]),
        Rcpp::as<Rcpp::CharacterVector>(el["peak_FC"]),
        Rcpp::as<Rcpp::List>(el["peak_transcripts"]),
        (uint32_t)Rcpp::as<int>(el["seed"])
    ));
  }

  //Pick thread count: min(requested, hw_concurrency, n_chrom). At least 1.
  unsigned int hw_threads = std::thread::hardware_concurrency();
  if (hw_threads == 0) hw_threads = 1;
  int n_workers = (max_threads > 0) ? max_threads : (int)hw_threads;
  if (n_workers > (int)hw_threads) n_workers = (int)hw_threads;
  if (n_workers > n_chrom)         n_workers = n_chrom;
  if (n_workers < 1)               n_workers = 1;

  //Output: pre-sized; each slot owned by exactly one worker (no lock needed).
  std::vector<ChromOutput> outputs(n_chrom);

  //Atomic counter is the entire synchronization mechanism for work claiming.
  std::atomic<int> next_idx{0};

  //Error capture for the cold path.
  std::atomic<bool> failed{false};
  std::mutex err_mu;
  std::string err_msg;

  auto worker = [&]() {
    while (true) {
      const int i = next_idx.fetch_add(1, std::memory_order_relaxed);
      if (i >= n_chrom) return;
      if (failed.load(std::memory_order_relaxed)) return;
      try {
        outputs[i] = process_chromosome_core(inputs[i]);
      } catch (const std::exception& e) {
        std::lock_guard<std::mutex> lk(err_mu);
        if (err_msg.empty()) err_msg = e.what();
        failed.store(true, std::memory_order_relaxed);
        return;
      } catch (...) {
        std::lock_guard<std::mutex> lk(err_mu);
        if (err_msg.empty()) err_msg = "unknown error in chromosome worker";
        failed.store(true, std::memory_order_relaxed);
        return;
      }
    }
  };

  //Create threads. Skip spawn if only one worker.
  if (n_workers == 1) {
    worker();
  } else {
    std::vector<std::thread> threads;
    threads.reserve(n_workers);
    for (int t = 0; t < n_workers; ++t) threads.emplace_back(worker);
    for (auto& th : threads) th.join();
  }

  if (failed.load()) {
    Rcpp::stop("control peak engine worker failed: %s", err_msg);
  }

  //Turn outputs back to Rcpp on the main thread.
  Rcpp::List result(n_chrom);
  for (int i = 0; i < n_chrom; ++i) {
    result[i] = wrap_chrom_output(outputs[i]);
  }
  return result;
}







