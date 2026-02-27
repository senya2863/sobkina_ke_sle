add_test([=[progonka.bebebe]=]  /home/ksenya/sobkina_ke_sle/tridiagonal/build/DiagThreeTest [==[--gtest_filter=progonka.bebebe]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[progonka.bebebe]=]  PROPERTIES WORKING_DIRECTORY /home/ksenya/sobkina_ke_sle/tridiagonal/build SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  DiagThreeTest_TESTS progonka.bebebe)
