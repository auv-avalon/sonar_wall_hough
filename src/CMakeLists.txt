rock_library(sonar_wall_hough
    SOURCES Hough.cpp Houghspace.cpp SonarPeak.cpp
    HEADERS Line.hpp Hough.hpp SonarPeak.hpp Config.hpp Houghspace.hpp
    DEPS_PKGCONFIG base-types)
