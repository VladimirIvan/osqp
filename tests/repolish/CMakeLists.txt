get_directory_property(headers
                        DIRECTORY ${PROJECT_SOURCE_DIR}/tests
                        DEFINITION headers)
                        
set(headers ${headers}
    ${CMAKE_CURRENT_SOURCE_DIR}/test_repolish.h PARENT_SCOPE)

