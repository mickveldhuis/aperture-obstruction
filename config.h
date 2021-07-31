#pragma once

namespace cfg {
    // Mount + telescope geometry
    inline constexpr double LENGTH_1 = 0.814;
    inline constexpr double LENGTH_2 = 1.098;
    inline constexpr double LENGTH_3 = 0.439;
    inline constexpr double LENGTH_4 = 0.356;
    inline constexpr double LENGTH_5 = 0.200;

    inline constexpr double OFFSET_ANGLE = 45;

    inline constexpr double APERTURE_DIAMETER = 0.4;
    inline constexpr double GUIDER_DIAMETER = 0.15;
    inline constexpr double FINDER_DIAMETER = 0.05;

    inline constexpr double APERTURE_SEC_DIAMETER = 0.175;
    inline constexpr double GUIDER_SEC_DIAMETER = 0.044;

    // Dome geometry
    inline constexpr double DOME_SLIT_WIDTH = 1.84;
    inline constexpr double DOME_EXTENT = 1.66;
    inline constexpr double DOME_DIAMETER = 6.25;
    
    // Observatory location
    inline constexpr double LATITUDE = 53.24;
    inline constexpr double LONGITUDE = 6.54;
}