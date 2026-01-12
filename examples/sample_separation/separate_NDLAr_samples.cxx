#include "duneanasel/nd/Selections.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TChain.h"

#include <memory>

int main(int argc, char const *argv[]) {

  TChain pch("cafmaker/cafTree");
  pch.Add(argv[1]);
  caf::StandardRecordProxy srp(&pch, "rec");

  Long64_t ents = pch.GetEntries();
  std::cout << "Input tree has " << ents << " entries." << std::endl;

  pch.GetEntry(0);

  // open a second copy of the file to copy full records out of
  TChain ch("cafmaker/cafTree");
  ch.Add(argv[1]);

  caf::StandardRecord *SR = nullptr;
  ch.SetBranchAddress("rec", &SR);

  std::map<sel::beam::Sample, int> selections;
  for (Long64_t i = 0; i < ents; ++i) {
    pch.GetEntry(i);

    for (auto const &nd_int : srp.common.ixn.pandora) {
      selections[sel::beam::ndlar::numode::ApplySelection(nd_int)]++;
    }
  }
}
