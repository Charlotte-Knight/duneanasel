add_executable(separate_FD_samples separate_FD_samples.cxx)
target_link_libraries(separate_FD_samples duneanasel::all)


add_executable(separate_NDLAr_samples separate_NDLAr_samples.cxx)
target_link_libraries(separate_NDLAr_samples duneanasel::all)

install(TARGETS
  separate_FD_samples
  separate_NDLAr_samples
  DESTINATION bin)
