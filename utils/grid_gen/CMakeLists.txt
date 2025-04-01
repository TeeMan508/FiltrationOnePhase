project(GridGen)
set(SOURCE main.cpp)

find_package(inmost REQUIRED)
link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

add_executable(GridGen ${SOURCE})

target_link_libraries(GridGen ${INMOST_LIBRARIES})

if(USE_MPI)
  message("linking GridGen with MPI")
  target_link_libraries(GridGen ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(GridGen PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)


install(TARGETS GridGen EXPORT inmost-targets RUNTIME DESTINATION bin)