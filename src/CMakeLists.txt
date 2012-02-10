rock_library(sonar_wall_hough
    SOURCES Line.cpp Hough.cpp Houghspace.cpp FilterHelper.cpp
    HEADERS Line.hpp Hough.hpp SonarPeak.hpp Config.hpp Houghspace.hpp FilterHelper.hpp
    DEPS_PKGCONFIG base-types)
