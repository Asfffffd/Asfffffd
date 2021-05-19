dune_add_test(NAME exercise_fractures
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              SOURCES main.cc)

# add a symlink for each input file
add_input_file_links()
# add a symlink for the grids folder
dune_symlink_to_source_files(FILES grids)

#install sources
install(FILES
exercisefractures.cc
matrixproblem.hh
fractureproblem.hh
matrixspatialparams.hh
fracturespatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/exercises/exercise-fractures)
