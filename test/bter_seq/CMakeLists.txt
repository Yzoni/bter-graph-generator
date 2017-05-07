include_directories(
        ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}
        ${gmock_SOURCE_DIR}/include ${gmock_SOURCE_DIR})

add_executable(runTests
        SetupTest.cpp)

target_link_libraries(runTests gtest gtest_main gmock gmock_main)
target_link_libraries(runTests bter_seq)
