# CMake generated Testfile for 
# Source directory: /Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/tests
# Build directory: /Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/build/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(bamtools_help "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/build/src/bamtools" "--help")
set_tests_properties(bamtools_help PROPERTIES  _BACKTRACE_TRIPLES "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/tests/CMakeLists.txt;1;add_test;/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/tests/CMakeLists.txt;0;")
add_test(bamtools_stats "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/build/src/bamtools" "stats" "-in" "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/tests/data/sam_spec_example.bam")
set_tests_properties(bamtools_stats PROPERTIES  _BACKTRACE_TRIPLES "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/tests/CMakeLists.txt;6;add_test;/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/tests/CMakeLists.txt;0;")
