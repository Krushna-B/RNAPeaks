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

// Global random number generator.
//Not thread safe
static std::mt19937 g_rng;

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
    const Rcpp::CharacterVector& anno_transcript,
    const Rcpp::CharacterVector& anno_region,
    const Rcpp::IntegerVector&   anno_start,
    const Rcpp::IntegerVector&   anno_end,
    const Rcpp::CharacterVector& anno_strand,
    const Rcpp::CharacterVector& anno_gene) {
  TranscriptMap out;
  const int n = anno_transcript.size();
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
  for (int i{}; i < n; ++i) {
    const char* tx_c   = anno_transcript[i];
    const char* reg_c  = anno_region[i];
    const char* gen_c  = anno_gene[i];
    const char* str_c  = anno_strand[i];
    const int   s      = anno_start[i];
    const int   e      = anno_end[i];

    TmpTx& t = tmp[std::string(tx_c)]; //Map tmp object
    if (t.strand == '*' && str_c && str_c[0]) t.strand = str_c[0]; //Store strand
    t.gene = gen_c ? gen_c : ""; //Store gene name

    RegionID rid = region_from_cstr(reg_c); //Convert Region String to ID
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
static GeneMap build_gene_map(const Rcpp::CharacterVector& gene_name,
                              const Rcpp::IntegerVector&   gene_start,
                              const Rcpp::IntegerVector&   gene_end) {
  GeneMap out;
  const int n = gene_name.size();
  out.reserve(n);
  for (int i{}; i < n; ++i) {
    const char* nm = gene_name[i];
    if (!nm) continue;
    out[std::string(nm)] = GeneRec{gene_start[i], gene_end[i]};
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

//Random Integer Picker from 1 to n
static int rng_pick_1_based(int n) {
  std::uniform_int_distribution<int> dist(1, n);
  return dist(g_rng);
}

//Picks random interval from [starts,ends)
static PickFlat sample_uniform_flat(const std::vector<int>& starts,
                                    const std::vector<int>& ends) {
  PickFlat p;
  if (starts.empty()) return p;
  const int i1 = rng_pick_1_based((int)starts.size()); // 1..n
  const int i  = i1 - 1;
  p.start = starts[i];
  p.end   = ends[i];
  p.ok    = true;
  return p;
}

//Simple control sampling for UTR/Exon Regions
static Control sample_simple_control(int peak_start, int peak_end, int peak_length,
                                     const std::vector<int>& cl_starts,
                                     const std::vector<int>& cl_ends,
                                     const char* final_region) {
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


  //Choose random start position
  const int chosen_global = rng_pick_1_based((int)total);


  long long cumulative = 0;
  int chosen_row = -1;
  long long previous_total = 0;
  //Find which interval contains randomly choosen position
  for (size_t i = 0; i < kn.size(); ++i) {
    const long long prev = cumulative;
    cumulative += kn[i];
    if (cumulative >= chosen_global) {
      chosen_row = (int)i;
      previous_total = prev;
      break;
    }
  }
  if (chosen_row < 0) return c;

  //Convert the global choice into an actual genomic start coordinate.
  //Create control interval.
  const int offset = chosen_global - (int)previous_total - 1;
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
                               const FlatRegion& rf, const Transcript& tx) {
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
  PickFlat p = sample_uniform_flat(red.starts, red.ends);
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
                                  const FlatRegion& rf, const Transcript& tx) {
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
  PickFlat p = sample_uniform_flat(red.starts, red.ends);
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
                           const Transcript& tx) {
  Control c;
  const int peak_length = peak_end - peak_start;
  if (peak_length < MIN_REGION_LENGTH) return c;
  if (peak_length <= 0) return c;
  if (cl_starts.empty()) return c;

  static const FlatRegion empty_rf;
  const FlatRegion& rf = tx.has[region] ? tx.regions[region] : empty_rf;

  if (region == REG_UTR3 || region == REG_UTR5 || region == REG_EXON) {
    return sample_simple_control(peak_start, peak_end, peak_length,
                                 cl_starts, cl_ends, region_cstr(region));
  }
  if (region == REG_CDS) {
    return get_cds_control(peak_start, peak_end, peak_length,
                           cl_starts, cl_ends, rf, tx);
  }
  if (region == REG_INTRON) {
    return get_intron_control(peak_start, peak_end, peak_length,
                              cl_starts, cl_ends, rf, tx);
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

// [[Rcpp::export]]
Rcpp::List process_chromosome_cpp(
    Rcpp::IntegerVector   anno_line_id,
    Rcpp::CharacterVector anno_transcript,
    Rcpp::CharacterVector anno_region,
    Rcpp::IntegerVector   anno_start,
    Rcpp::IntegerVector   anno_end,
    Rcpp::CharacterVector anno_strand,
    Rcpp::CharacterVector anno_gene,
    Rcpp::CharacterVector gene_name,
    Rcpp::IntegerVector   gene_start,
    Rcpp::IntegerVector   gene_end,
    Rcpp::IntegerVector   peak_start,
    Rcpp::IntegerVector   peak_end,
    Rcpp::CharacterVector peak_FC,
    Rcpp::CharacterVector peak_range,
    Rcpp::List            peak_transcripts,
    uint32_t              seed
) {
  (void)anno_line_id; // sorting is the caller's responsibility (line_id ascending)

  g_rng.seed(seed);

  const int n_peaks = peak_start.size();

  // No peaks or no anno -> nothing to do.
  if (n_peaks == 0 || anno_transcript.size() == 0) {
    return Rcpp::List::create(
      Rcpp::Named("peak_start")    = Rcpp::IntegerVector(0),
      Rcpp::Named("peak_end")      = Rcpp::IntegerVector(0),
      Rcpp::Named("name")          = Rcpp::CharacterVector(0),
      Rcpp::Named("strand")        = Rcpp::CharacterVector(0),
      Rcpp::Named("control_start") = Rcpp::IntegerVector(0),
      Rcpp::Named("control_end")   = Rcpp::IntegerVector(0)
    );
  }

  // Build per-chromosome lookups.
  TranscriptMap tx_map = build_transcript_annotations(
      anno_transcript, anno_region, anno_start, anno_end, anno_strand, anno_gene);
  GeneMap gene_map = build_gene_map(gene_name, gene_start, gene_end);

  // Per-peak transcript-name lists (materialize once).
  std::vector<std::vector<std::string>> peak_tx(n_peaks);
  for (int i = 0; i < n_peaks; ++i) {
    Rcpp::CharacterVector v = peak_transcripts[i];
    peak_tx[i].reserve(v.size());
    for (int j = 0; j < v.size(); ++j) {
      peak_tx[i].push_back(std::string(v[j]));
    }
  }

  // -----------------------------
  // Pass 1: build peak_region_dict + per-strand peak totals.
  // -----------------------------

  // bucket[p_i][region] = list of tx names qualifying that region for peak p_i.
  struct Bucket {
    std::array<std::vector<std::string>, REG_N> tx_per_region;
    std::array<bool, REG_N> present{};
    bool any = false;
  };
  std::vector<Bucket> buckets(n_peaks);
  std::vector<int> peak_order; peak_order.reserve(n_peaks);

  std::vector<int> plus_s, plus_e, minus_s, minus_e;
  plus_s.reserve(n_peaks); plus_e.reserve(n_peaks);
  minus_s.reserve(n_peaks); minus_e.reserve(n_peaks);

  for (int p_i = 0; p_i < n_peaks; ++p_i) {
    const int ps = peak_start[p_i];
    const int pe = peak_end[p_i];
    if (ps == NA_INTEGER || pe == NA_INTEGER || pe <= ps) continue;
    const int peak_length = pe - ps;
    (void)peak_length;

    for (const std::string& tx_name : peak_tx[p_i]) {
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

  if (peak_order.empty()) {
    return Rcpp::List::create(
      Rcpp::Named("peak_start")    = Rcpp::IntegerVector(0),
      Rcpp::Named("peak_end")      = Rcpp::IntegerVector(0),
      Rcpp::Named("name")          = Rcpp::CharacterVector(0),
      Rcpp::Named("strand")        = Rcpp::CharacterVector(0),
      Rcpp::Named("control_start") = Rcpp::IntegerVector(0),
      Rcpp::Named("control_end")   = Rcpp::IntegerVector(0)
    );
  }

  // Total peaks per strand -> sorted reduced.
  FlatRegion tpp = flat_reduce(std::move(plus_s),  std::move(plus_e));
  FlatRegion tpm = flat_reduce(std::move(minus_s), std::move(minus_e));

  // -----------------------------
  // Pass 2: per peak (in input order), pick a region, walk windows, sample.
  // -----------------------------

  // Running totals of accepted controls per strand; kept sorted by start.
  std::vector<int> tcp_s, tcp_e, tcm_s, tcm_e;
  // processed_windows: shared across all peaks within this chromosome,
  // mirrors `new.env(...)` in R created once before the peak loop.
  std::unordered_set<long long> processed_windows;

  // Output buffers.
  std::vector<int>         res_pstart, res_pend, res_cstart, res_cend;
  std::vector<std::string> res_name, res_strand;
  res_pstart.reserve(peak_order.size());
  res_pend.reserve(peak_order.size());
  res_cstart.reserve(peak_order.size());
  res_cend.reserve(peak_order.size());
  res_name.reserve(peak_order.size());
  res_strand.reserve(peak_order.size());

  for (int p_i : peak_order) {
    const Bucket& b = buckets[p_i];
    const RegionID peak_region = choose_peak_region(b.present);

    const int peak_s = peak_start[p_i];
    const int peak_e = peak_end[p_i];
    const std::string fc(Rcpp::as<const char*>(peak_FC[p_i]));
    const int peak_length = peak_e - peak_s;

    const std::vector<std::string>& peak_txs = b.tx_per_region[peak_region];

    for (const std::string& tx_name : peak_txs) {
      auto it = tx_map.find(tx_name);
      if (it == tx_map.end()) continue;
      const Transcript& tx = it->second;

      const bool plus = (tx.strand == '+');
      const std::vector<int>& tp_s = plus ? tpp.starts : tpm.starts;
      const std::vector<int>& tp_e = plus ? tpp.ends   : tpm.ends;

      // Blacklist (peak-occupied space, restricted to this gene's window).
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

      // Region in-set for setdiff. For intron we collapse to a single
      // [tx_start, tx_end] interval, matching the R reference.
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

      // Windows: peak itself if short, otherwise the Python sliding window.
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

        Control ctrl = get_control(ws_s, ws_e, peak_region, cl.starts, cl.ends, tx);
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

        // Accept: insert into the per-strand running totals.
        if (plus) sorted_insert(tcp_s, tcp_e, cs, ce);
        else      sorted_insert(tcm_s, tcm_e, cs, ce);

        processed_windows.insert(key);

        // Output row.
        res_pstart.push_back(ws_s);
        res_pend.push_back(ws_e);
        res_cstart.push_back(cs);
        res_cend.push_back(ce);
        res_name.push_back(tx.gene + "_" + ctrl.region_type + "_" + fc);
        res_strand.push_back(std::string(1, tx.strand));

        peak_processed = true;
      }

      if (peak_processed) break;
    }
  }

  // -----------------------------
  // Assemble Rcpp output.
  // -----------------------------

  const int nr = (int)res_pstart.size();
  Rcpp::IntegerVector   o_ps(nr), o_pe(nr), o_cs(nr), o_ce(nr);
  Rcpp::CharacterVector o_nm(nr), o_st(nr);
  for (int i = 0; i < nr; ++i) {
    o_ps[i] = res_pstart[i];
    o_pe[i] = res_pend[i];
    o_cs[i] = res_cstart[i];
    o_ce[i] = res_cend[i];
    o_nm[i] = res_name[i];
    o_st[i] = res_strand[i];
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
