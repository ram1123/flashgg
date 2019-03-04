#include "flashgg/Taggers/interface/HHWWggCandidateDumper.h"

#include "flashgg/Taggers/interface/PluggableAnalyzer.h"

namespace flashgg {
  namespace fwlite {
    PLUGGABLE_ANALYZER( HHWWggCandidateDumper );
    PLUGGABLE_ANALYZER( CutBasedHHWWggCandidateDumper );
  }
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
