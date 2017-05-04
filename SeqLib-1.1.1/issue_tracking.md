* In `RefGenome.cpp`, the ctor is better to use initializer list rather than 
assignment in the function body, following best practices. In fact, considering the issues raise with these raw pointers, it may be better to use `std::shared_ptr<faidx_t>` and call `get()` whenever the pointer is requested.
  The proposed temporary fix is (comments are removed because the codes are obvious and 
  excessive comments may become maintenance nightmare):

      RefGenome::RefGenome(const std::string& file ): index(nullptr) {

        if (!read_access_test(file)) {
          index = nullptr;
          throw std::invalid_argument("RefGenome: file not found - " + file);
        }

        index = fai_load(file.c_str());
      }

* In `RefGenome::queryRegion`, the automatic variable `int len` was left uninitialized, while requested by `faidx_fetch_seq` defined in `htslib`.
  Proposal:

      int len{0};

* In `RefGenome:queryRegion`, the check for returned `char *` is throwing, and the data member `index` is non-null (because the check precedes this one passes).
This leads to resource leak.
  Proposal:

      if (!f){
        index = nullptr; // added line for exception safety
        throw std::invalid_argument("RefGenome::queryRegion - Could not find valid sequence");
      }

* In `SnowUtils.h`, `SnowTools::AddCommas` has an issue that if template parameter is floating point, format is off. The proposal is to change the name of the function to `SnowTools::AddCommasToInt`, to avoid applying the function to floating point types.

* In `GenomicRegion.h`, the `friend` declaration at the beginning is unnecessary because right now there isn't any private data.

* In `GenomicRegion.h`, why ctor `GenomicRegion(int32_t, int32_t, int32_t, char)` doesn't check for positivity of chromosome range, and start, end?
    In addition, why aren't the members `chr`, `pos1`, `pos2` unsigned and private?

* In `GenomicRegion.h`, ctor 
    `GenomicRegion(const std::string&, const std::string&, const std::string&, bam_hdr_t*)`, why subtract 1 from the string converted chromosome number?
    Is it better to have an `else` block than early return statement?

* In `GenomicRegion.h`, for function `GenomicRegion::random`, the bound value of 25 in the `for` loop is better to be replaced with a macro, as the tool itself is not necessarily bound to human studies.
    In addition, why an assertion that `k>0`? This immediately fails when `k` starts from 0, as the `for` loop does.

* In `GenomicRegion.h`, why is function `GenomicRegion::chrToString` marked `static`? Is it possible to have an utility function extracted and put it `SnowUtils.h`? Again, the implementation assumes human.

* In `GenomicRegion.h`, why does function `GenomicRegion::isEmpty` insist on `chr=-1`. Is it possible that `chr != 1`, but `pos1 == pos2`, so that it is actually empty?
  Is empty `GenomicRegion` always invalid (because validity check as implemented in `GenomicRegion::valid()`, which by the way should be renamed as `is_valid`, checks for positivity of chromosome)?

* In `GenomicRegion.h`, function `GenomicRegion::cetromereOverlap` is declared but not defined.

* In `GenomicRegion.h` and `GenomicRegion.cpp`, function `GenomicRegion::getOverlap` has a typo: it should be taking references instead of a copy.

* In `GenomicRegion.cpp`, the check in function `GenomicRegion::pad` is possibly wrong. If the argument provided is indeed negative, then multiplying by 2 is unnecessary and comparing with width is also unnecessary. Second, it is actually better to change the function parameter from `int32_t` to `uint32_t` because negative padding is counter-intuitive. Third, the check didn't check for the case that padding to the left/right goes out of range (actually commented out).

* In `Fractions.cpp`, the `getline` function in the `while` loop shouldn't use "\n" as the delimiter, as it will be different on other OS (Windows for sure).
  There also is no member function declared as `bool valid() const` for `FracRegion`, but is used here.
  Last, it is better to have a `stringSplit` function defined in `SnowUtils.h` for breaking lines of delimited files (eg. CSV, tab) into fields.

* In `gChain.h`, why isn't default ctor deleted? It left data members uninitialized.

* In `gChain.cpp`, the function definitions are not wrapped in `namespace SnowTools`.

* In `BamStats.h`, not ctors defined although there are data members. 
  Proposal:
    
      BamStats() : m_group_map(){}

* In `BamStats.cpp`, function `BamReadGroup::addRead` should take a `const` reference instead of reference to signal the read is not being modified. 
    Similarly for `BamStats::addRead`.

*

*

*

*

*

*

* In `run_snowman2.cpp` in the `SnowmanSV` repo, line 25 `#define MATE_LOOKUP_MIN 3` and line 146 `static int32_t mate_lookup_min = 3;` potentially conflicts.
