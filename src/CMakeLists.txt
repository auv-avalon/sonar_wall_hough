rock_library(sonar_wall_hough
    SOURCES Dummy.cpp
    HEADERS Dummy.hpp)

rock_executable(sonar_wall_hough_bin Main.cpp
    DEPS sonar_wall_hough)

