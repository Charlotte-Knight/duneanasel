#include "duneanasel/nd/Selections.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TChain.h"
#include "TH1D.h"

#include <memory>

int main(int argc, char const *argv[]) {

  TChain pch("cafmaker/cafTree");
  pch.Add(argv[1]);
  caf::StandardRecordProxy srp(&pch, "rec");

  Long64_t ents = pch.GetEntries();
  std::cout << "Input tree has " << ents << " entries." << std::endl;

  pch.GetEntry(0);

  TFile hout("mutracklen.root", "RECREATE");
  TH1D mutracklen("mutracklen", ";Track Length [cm]; Count", 100, 0, 500);
  for (Long64_t i = 0; i < ents; ++i) {
    pch.GetEntry(i);

    for (auto const &nd_int : srp.common.ixn.pandora) {
      if (sel::beam::ndlar::numode::ApplySelection(nd_int) !=
          sel::beam::kNuMuCCLikeContained) {
        continue;
      }

      mutracklen.Fill(sel::beam::ndlar::ParticleLength(
          *sel::beam::ndlar::GetLongestParticle(nd_int)));
    }
  }

  hout.Write();
  hout.Close();
}
