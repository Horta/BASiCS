add_library(basics_objs OBJECT posterior.cpp phi.cpp s.cpp gig.cpp)
set_target_properties(basics_objs
  PROPERTIES
    POSITION_INDEPENDENT_CODE ON)
