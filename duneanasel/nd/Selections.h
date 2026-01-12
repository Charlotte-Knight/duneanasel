#pragma once

#include "duneanaobj/StandardRecord/SRInteraction.h"

#include "duneanasel/common/type_traits.h"

//The selection in this file was originally developed by Eva Sabater, what follows
//  is a transcription from the original CAFAna code.

namespace sel {
namespace beam {

enum Sample { kRejected = 0, kNuMuCCLikeContained, kNuECCLikeContained };

namespace ndlar {

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool InFV(T const &nd_int) {

  constexpr double NDLArXLo = -346.9;
  constexpr double NDLArXHi = 346.9;
  constexpr double NDLArYLo = -215.5;
  constexpr double NDLArYHi = 81.7;
  constexpr double NDLArZLo = 418.2;
  constexpr double NDLArZHi = 913.3;

  float x = nd_int.vtx.x;
  if ((x <= (NDLArXLo + 25)) || (x >= (NDLArXHi - 25))) {
    return false;
  }
  float y = nd_int.vtx.y;
  if ((y <= (NDLArYLo + 25)) || (y >= (NDLArYHi - 25))) {
    return false;
  }
  float z = nd_int.vtx.z;
  if ((z <= (NDLArZLo + 25)) || (z >= (NDLArZHi - 25))) {
    return false;
  }

  return true;
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool AllPartContained(T const &nd_int) {
  for (auto const &p : nd_int.part.dlp) {
    if ((p.primary == 1) && (p.contained != 1)) {
      return false;
    }
  }
  return true;
}

template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
inline double ParticleLength(T const &p) {
  return std::hypot(p.end.x - p.start.x, p.end.y - p.start.y,
                    p.end.z - p.start.z);
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline caf::SRRecoParticleProxy const *GetLongestParticle(T const &nd_int) {
  if (!nd_int.part.dlp.size()) {
    return nullptr;
  }

  double longest = -std::numeric_limits<double>::max();
  size_t longest_idx = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < nd_int.part.dlp.size(); ++i) {
    const auto &p = nd_int.part.dlp[i];
    if (p.primary != 1) {
      continue;
    }
    double const &plen = ParticleLength(p);
    if (plen > longest) {
      longest = plen;
      longest_idx = i;
    }
  }
  return &nd_int.part.dlp[longest_idx];
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool HasParticleWithReconstructedPID(T const &nd_int, int pid) {
  for (auto const &p : nd_int.part.dlp) {
    if ((p.primary == 1) && (p.pdg == pid)) {
      return false;
    }
  }
  return false;
}

namespace numode {

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool NuMuCCLikeContained(T const &nd_int) {
  if (!InFV(nd_int) || !AllPartContained(nd_int)) {
    return false;
  }
  auto plong = GetLongestParticle(nd_int);
  if (!plong || (plong->pdg != 13)) {
    return false;
  }
  return true;
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool NuECCLikeContained(T const &nd_int) {
  if (!InFV(nd_int) || !AllPartContained(nd_int) ||
      !HasParticleWithReconstructedPID(nd_int, 11)) {
    return false;
  }
  return true;
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline sel::beam::Sample ApplySelection(T const &fd_int) {

  if (!sel::beam::ndlar::InFV(fd_int)) {
    return sel::beam::kRejected;
  }

  if (sel::beam::ndlar::numode::NuMuCCLikeContained(fd_int)) {
    return sel::beam::kNuMuCCLikeContained;
  } else if (sel::beam::ndlar::numode::NuECCLikeContained(fd_int)) {
    return sel::beam::kNuECCLikeContained;
  }
  return sel::beam::kRejected;
}

} // namespace numode
} // namespace ndlar

} // namespace beam
} // namespace sel
